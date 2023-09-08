use std::{
    thread,
    sync::Arc,
    io::{self, Write},
};
use fx::FxHashMap;
use smallvec::SmallVec;
use super::{
    NamedSeq,
    wfa::{Aligner, Penalties},
    cigar::{Cigar, Operation},
};

const CAPACITY: usize = 4;
type KmerCache<'a> = FxHashMap<&'a [u8], SmallVec<[u32; CAPACITY]>>;

/// Creates a HashMap containing all the k-mers in the sequence.
/// A good rolling hash function should speed up the code.
fn cache_kmers(seq: &[u8], k: usize) -> KmerCache {
    let mut map = KmerCache::default();
    for (i, kmer) in seq.windows(k).enumerate() {
        map.entry(kmer)
            .or_insert_with(SmallVec::new)
            .push(i as u32);
    }
    map
}

/// Finds all matches between k-mers in two caches.
fn kmer_matches(cache1: &KmerCache, cache2: &KmerCache) -> Vec<(u32, u32)> {
    let mut matches = Vec::new();
    for (kmer1, positions1) in cache1.iter() {
        if let Some(positions2) = cache2.get(kmer1) {
            for &pos1 in positions1 {
                for &pos2 in positions2 {
                    matches.push((pos1, pos2));
                }
            }
        }
    }
    matches.sort_unstable();
    matches
}

// /// Aligns two sequences between k-mer matches.
// fn align(
//     seq1: &[u8],
//     cache1: &KmerCache,
//     seq2: &[u8],
//     cache2: &KmerCache,
//     aligner: &Aligner,
//     penalties: &Penalties,
//     k: usize,
// ) -> Result<(Cigar, i32), Error>
// {
//     let matches = kmer_matches(cache1, cache2);

// }

fn safe_align(aligner: &Aligner, entry1: &NamedSeq, entry2: &NamedSeq) -> (Cigar, i32) {
    unimplemented!()
    // match aligner.align(entry1.seq(), entry2.seq()) {
    //     Ok(cigar_and_score) => cigar_and_score,
    //     Err((ch, raw_cigar)) => {
    //         log::error!("Could not align {} and {}. Violating CIGAR character '{}' ({}) in {:?}",
    //             entry1.name(), entry2.name(), char::from(ch), ch, raw_cigar);
    //         (Cigar::new(), i32::MIN)
    //     }
    // }
}

/// Calculates all pairwise divergences between all sequences, writes alignments to PAF file,
/// and returns condensed distance matrix (see `kodama` crate).
pub fn pairwise_divergences(
    entries: &[NamedSeq],
    mut paf_writer: impl Write,
    penalties: &Penalties,
    threads: u16,
) -> io::Result<Vec<f64>> {
    if threads == 1 {
        log::debug!("        Aligning sequences in 1 thread");
        let n = entries.len();
        let mut divergences = Vec::with_capacity(n * (n - 1) / 2);
        let aligner = Aligner::new(penalties.clone());
        for (i, entry1) in entries.iter().enumerate() {
            for entry2 in entries[i + 1..].iter() {
                let (cigar, score) = safe_align(&aligner, entry1, entry2);
                let divergence = write_paf(&mut paf_writer, entry2, entry1, &cigar, score)?;
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
    log::debug!("        Aligning sequences in {} threads", threads);
    let entries = Arc::new(entries.to_vec());
    let n = entries.len() as u32;
    let pairs: Arc<Vec<(u32, u32)>> = Arc::new(
        (0..n - 1).flat_map(|i| (i + 1..n).map(move |j| (i, j))).collect());
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
            let entries = Arc::clone(&entries);
            let penalties = penalties.clone();
            handles.push(thread::spawn(move || {
                assert!(start < end);
                let aligner = Aligner::new(penalties);
                pairs[start..end].iter()
                    .map(|&(i, j)| safe_align(&aligner, &entries[i as usize], &entries[j as usize]))
                    .collect::<Vec<_>>()
            }));
        }
        start = end;
    }
    assert_eq!(start, n_pairs);

    let mut pairs_iter = pairs.iter();
    let mut divergences = Vec::with_capacity(n_pairs);
    for handle in handles.into_iter() {
        for (cigar, score) in handle.join().expect("Worker process failed") {
            let &(i, j) = pairs_iter.next().expect("Number of solutions does not match the number of tasks");
            divergences.push(write_paf(&mut paf_writer, &entries[j as usize], &entries[i as usize], &cigar, score)?);
        }
    }
    assert!(pairs_iter.next().is_none(), "Number of solutions does not match the number of tasks");
    Ok(divergences)
}

/// Writes the alignment to PAF file and returns sequence divergence (edit distance / total aln length),
/// where total alignment length includes gaps into both sequences.
fn write_paf(
    writer: &mut impl Write,
    query: &NamedSeq,
    refer: &NamedSeq,
    cigar: &Cigar,
    score: i32,
) -> io::Result<f64>
{
    let qname = query.name();
    let rname = refer.name();
    let qlen = query.seq().len() as u32;
    let rlen = refer.seq().len() as u32;
    write!(writer, "{qname}\t{qlen}\t0\t{qlen}\t+\t{rname}\t{rlen}\t0\t{rlen}\t")?;

    if cigar.is_empty() {
        writeln!(writer, "0\t0t\t0")?;
        return Ok(1.0);
    }

    assert_eq!(qlen, cigar.query_len());
    assert_eq!(rlen, cigar.ref_len());
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
