use std::{
    thread,
    sync::Arc,
    io::{self, Write},
};
use fx::FxHashMap;
use smallvec::SmallVec;
use crate::{
    seq::{
        NamedSeq,
        wfa::{Aligner, Penalties},
        cigar::{Cigar, Operation},
    },
    err::{Error, add_path},
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
fn find_kmer_matches(cache1: &KmerCache, cache2: &KmerCache) -> Vec<(u32, u32)> {
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

/// Aligns two sequences between k-mer matches.
/// First sequence: reference, sequence: query.
fn align(
    aligner: &Aligner,
    seq1: &[u8],
    seq2: &[u8],
    kmer_matches: &[(u32, u32)],
    k: usize,
) -> Result<(Cigar, i32), Error>
{
    let penalties = aligner.penalties();
    let sparse_aln = bio::alignment::sparse::lcskpp(&kmer_matches, k);
    let mut cigar = Cigar::new();
    let mut score = 0;
    let mut i1 = 0;
    let mut j1 = 0;
    let k = k as u32;
    for ix in sparse_aln.path.into_iter() {
        let (i2, j2) = kmer_matches[ix];
        if i1 > i2 {
            continue;
        }
        let jump1 = i2 - i1;
        let jump2 = j2.checked_sub(j1).unwrap();
        if jump1 > 0 && jump2 > 0 {
            score += aligner.align(&seq1[i1 as usize..i2 as usize], &seq2[j1 as usize..j2 as usize], &mut cigar)?;
        } else if jump1 > 0 {
            score -= penalties.gap_open + jump1 as i32 * penalties.gap_extend;
            cigar.push_checked(Operation::Del, jump1);
        } else if jump2 > 0 {
            score -= penalties.gap_open + jump2 as i32 * penalties.gap_extend;
            cigar.push_checked(Operation::Ins, jump2);
        }

        cigar.push_checked(Operation::Equal, k);
        i1 = i2 + k;
        j1 = j2 + k;
    }

    let n = seq1.len() as u32;
    let m = seq2.len() as u32;
    let jump1 = n.checked_sub(i1).unwrap();
    let jump2 = m.checked_sub(j1).unwrap();
    if jump1 > 0 && jump2 > 0 {
        score += aligner.align(&seq1[i1 as usize..n as usize], &seq2[j1 as usize..m as usize], &mut cigar)?;
    } else if jump1 > 0 {
        score -= penalties.gap_open + jump1 as i32 * penalties.gap_extend;
        cigar.push_checked(Operation::Del, jump1);
    } else if jump2 > 0 {
        score -= penalties.gap_open + jump2 as i32 * penalties.gap_extend;
        cigar.push_checked(Operation::Ins, jump2);
    }
    Ok((cigar, score))
}

/// Calculates all pairwise divergences between all sequences, writes alignments to PAF file,
/// and returns condensed distance matrix (see `kodama` crate).
pub fn pairwise_divergences(
    entries: &[NamedSeq],
    mut paf_writer: impl Write,
    penalties: &Penalties,
    k: usize,
    threads: u16,
) -> Result<Vec<f64>, Error> {
    let caches: Vec<_> = entries.iter().map(|entry| cache_kmers(entry.seq(), k)).collect();
    if threads == 1 {
        log::debug!("        Aligning sequences in 1 thread");
        let n = entries.len();
        let mut divergences = Vec::with_capacity(n * (n - 1) / 2);
        let aligner = Aligner::new(penalties.clone());
        for (i, (entry1, cache1)) in entries.iter().zip(&caches).enumerate() {
            for (entry2, cache2) in entries[i + 1..].iter().zip(&caches[i + 1..]) {
                let kmer_matches = find_kmer_matches(cache1, cache2);
                let (cigar, score) = align(&aligner, entry1.seq(), entry2.seq(), &kmer_matches, k)?;
                let divergence = write_paf(&mut paf_writer, entry2, entry1, &cigar, score)
                    .map_err(add_path!(!))?;
                divergences.push(divergence);
            }
        }
        Ok(divergences)
    } else {
        divergences_multithread(entries, &caches, paf_writer, penalties, k, threads)
    }
}

fn divergences_multithread(
    entries: &[NamedSeq],
    caches: &[KmerCache],
    mut paf_writer: impl Write,
    penalties: &Penalties,
    k: usize,
    threads: u16,
) -> Result<Vec<f64>, Error> {
    let threads = usize::from(threads);
    log::debug!("        Aligning sequences in {} threads", threads);
    let n = entries.len() as u32;
    let pairs: Arc<Vec<(u32, u32)>> = Arc::new(
        (0..n - 1).flat_map(|i| (i + 1..n).map(move |j| (i, j))).collect());
    let n_pairs = pairs.len();
    let mut handles = Vec::with_capacity(threads);

    let entries = Arc::new(entries.to_vec());
    // Need to precompute all kmer matches as `caches` cannot be passed between threads.
    let all_kmer_matches = Arc::new(pairs.iter()
        .map(|&(i, j)| find_kmer_matches(&caches[i as usize], &caches[j as usize])).collect::<Vec<_>>());
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
            let all_kmer_matches = Arc::clone(&all_kmer_matches);
            let penalties = penalties.clone();
            handles.push(thread::spawn(move || {
                assert!(start < end);
                let aligner = Aligner::new(penalties);
                pairs[start..end].iter().zip(all_kmer_matches[start..end].iter()).map(|(&(i, j), kmer_matches)|
                        align(&aligner, &entries[i as usize].seq(), entries[j as usize].seq(), kmer_matches, k))
                    .collect::<Result<Vec<_>, Error>>()
            }));
        }
        start = end;
    }
    assert_eq!(start, n_pairs);

    // Collect all results so that even if there is a problem in one of the sequences, all handles are closed.
    let results: Vec<_> = handles.into_iter().map(|handle| handle.join().expect("Worker process failed")).collect();
    let mut pairs_iter = pairs.iter();
    let mut divergences = Vec::with_capacity(n_pairs);
    for res in results.into_iter() {
        for (cigar, score) in res?.into_iter() {
            let &(i, j) = pairs_iter.next().expect("Number of solutions does not match the number of tasks");
            divergences.push(write_paf(&mut paf_writer, &entries[j as usize], &entries[i as usize], &cigar, score)
                .map_err(add_path!(!))?);
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
