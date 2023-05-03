use std::{
    thread,
    sync::Arc,
    io::{self, Write},
};
use crate::{
    ext::vec::Tuples,
};
use super::{
    NamedSeq,
    wfa::{Aligner, Penalties},
    cigar::{Cigar, Operation},
};

/// Calculates all pairwise divergences between all sequences, writes alignments to PAF file,
/// and returns condensed distance matrix (see `kodama` crate).
pub fn pairwise_divergences(
    entries: &[NamedSeq],
    mut paf_writer: impl Write,
    penalties: &Penalties,
    threads: u16,
) -> io::Result<Vec<f64>> {
    if threads == 1 {
        let n = entries.len();
        let mut divergences = Vec::with_capacity(n * (n - 1) / 2);
        let aligner = Aligner::new(penalties);
        for (i, entry1) in entries.iter().enumerate() {
            for entry2 in entries[i + 1..].iter() {
                let (cigar, score) = aligner.align(entry1.seq(), entry2.seq());
                let divergence = write_paf(&mut paf_writer, entry2, entry1, score, &cigar)?;
                divergences.push(divergence);
            }
        }
        Ok(divergences)
    } else {
        divergences_multithread(entries, paf_writer, penalties, threads)
    }
}

fn divergences_multithread(
    entries: &[NamedSeq],
    mut paf_writer: impl Write,
    penalties: &Penalties,
    threads: u16,
) -> io::Result<Vec<f64>> {
    let threads = usize::from(threads);
    let seqs = Arc::new(entries.iter().map(|entry| entry.seq.clone()).collect::<Vec<Vec<u8>>>());
    let n = seqs.len();
    let pairs = Arc::new(Tuples::combinations(&(0..n).collect::<Vec<_>>(), 2));
    let n_pairs = pairs.len();
    let mut handles = Vec::with_capacity(threads);

    let mut start = 0;
    for i in 0..threads {
        if start == n_pairs {
            break;
        }
        let rem_workers = threads - i;
        // Ceiling division.
        let end = start + ((n_pairs - start) + rem_workers - 1) / rem_workers;
        // Closure with cloned data.
        {
            let pairs = Arc::clone(&pairs);
            let seqs = Arc::clone(&seqs);
            let penalties = penalties.clone();
            handles.push(thread::spawn(move || {
                assert!(start < end);
                let aligner = Aligner::new(&penalties);
                pairs.iter_range(start..end)
                    .map(|pair| aligner.align(&seqs[pair[0]], &seqs[pair[1]]))
                    .collect::<Vec<_>>()
            }));
        }
        start = end;
    }
    assert_eq!(start, n_pairs);

    let mut ix = 0;
    let mut divergences = Vec::with_capacity(n_pairs);
    for handle in handles.into_iter() {
        for (cigar, score) in handle.join().expect("Worker process failed") {
            match &pairs[ix] {
                [i, j] => divergences.push(write_paf(&mut paf_writer, &entries[*j], &entries[*i], score, &cigar)?),
                _ => unreachable!(),
            }
            ix += 1;
        }
    }
    assert_eq!(ix, n_pairs);
    Ok(divergences)
}

/// Writes the alignment to PAF file and returns sequence divergence (edit distance / total aln length),
/// where total alignment length includes gaps into both sequences.
fn write_paf(
    writer: &mut impl Write,
    query: &NamedSeq,
    refer: &NamedSeq,
    score: i32,
    cigar: &Cigar,
) -> io::Result<f64>
{
    let qlen = cigar.query_len();
    let rlen = cigar.ref_len();
    assert_eq!(query.seq().len() as u32, cigar.query_len());
    assert_eq!(refer.seq().len() as u32, cigar.ref_len());
    let qname = query.name();
    let rname = refer.name();
    write!(writer, "{qname}\t{qlen}\t0\t{qlen}\t+\t{rname}\t{rlen}\t0\t{rlen}\t")?;

    let mut nmatches = 0;
    let mut total_size = 0;
    for item in cigar.iter() {
        total_size += item.len();
        if item.operation() == Operation::Equal {
            nmatches += item.len();
        }
    }
    let edit_dist = total_size - nmatches;
    let divergernce = f64::from(edit_dist) / f64::from(total_size);
    writeln!(writer, "{nmatches}\t{total_size}\t60\t\
        NM:i:{edit_dist}\tAS:i:{score}\tdv:f:{divergernce:.7}\tcg:Z:{cigar}")?;
    Ok(divergernce)
}
