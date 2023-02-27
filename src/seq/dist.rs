use std::{
    cmp::min,
    io::{self, Write},
};
use once_cell::unsync::OnceCell;
use nohash::{IntSet, IntMap};
use bio::alignment::sparse::lcskpp;
use super::{
    contigs::{ContigId, ContigNames},
    kmers::{KmerCounts, canonical_kmers, N_KMER},
};

struct ContigKmers<'a> {
    contig_id: ContigId,
    seq: &'a [u8],

    k: u8,
    /// Key: k-mer hash, value: number of appearances.
    kmer_multiset: IntMap<u64, u16>,
    /// Store positions of rare k-mer in static vectors, which are spilled into heap if needed.
    kmer_positions: IntMap<u64, Vec<u32>>,
    /// Total number of rare k-mers, same as `kmer_multiset.values().sum`.
    rare_kmers: u32,
    /// Maximal number of non-overlapping rare k-mers.
    nonoverl_kmers: u32,
}

impl<'a> ContigKmers<'a> {
    fn new(contig_id: ContigId, seq: &'a [u8], kmer_counts: &'a KmerCounts, threshold: u8) -> Self {
        let ref_occurances = kmer_counts.get(contig_id);
        let k: u8 = kmer_counts.k().try_into().unwrap();
        let canon_kmers = canonical_kmers(seq, k);
        assert_eq!(canon_kmers.len(), ref_occurances.len());

        // TODO: kmer_positions are enough, no need to also store kmer_multiset.
        let mut kmer_multiset = IntMap::default();
        let mut kmer_positions = IntMap::default();
        let mut rare_kmers = 0;
        let mut nonoverl_kmers = 0;
        let mut next_nonoverl_pos = 0;

        // println!("New for {}", contig_id);
        for (pos, (kmer, &ref_count)) in canon_kmers.into_iter().zip(ref_occurances).enumerate() {
            let pos = pos as u32;
            // println!("    Pos {:5}  count {}   ({})", pos, ref_count, kmer != N_KMER && ref_count <= threshold);
            if kmer != N_KMER && ref_count <= threshold {
                *kmer_multiset.entry(kmer).or_insert(0) += 1;
                kmer_positions.entry(kmer).or_insert_with(|| Vec::with_capacity(2)).push(pos);
                rare_kmers += 1;
                if pos >= next_nonoverl_pos {
                    // println!("        Non-overl");
                    next_nonoverl_pos = pos + u32::from(k);
                    nonoverl_kmers += 1;
                }
            }
        }
        // println!("In total: {} rare kmers,  {} non-overl", rare_kmers, nonoverl_kmers);
        Self {
            contig_id, seq, k,
            kmer_multiset, kmer_positions,
            rare_kmers, nonoverl_kmers,
        }
    }

    /// Compares `self` v `othr` and writes to a CSV stream `f`.
    /// First three columns: contig1, contig2, strand.
    ///
    /// Next, three groups of two columns of type `inters/size1,size2   dist`, where distance
    /// is calculated as `1 - jaccard_index`.
    /// Three groups correspond to multisets of k-mers; LCS; and non-overlapping LCS.
    fn compare<W: Write>(&self, othr: &Self, contigs: &ContigNames, f: &mut W) -> io::Result<()> {
        // println!("Compare {} & {}", contigs.get_name(self.contig_id), contigs.get_name(othr.contig_id));
        write!(f, "{}\t{}\t", contigs.get_name(self.contig_id), contigs.get_name(othr.contig_id))?;
        let mset_inters = multiset_intersection(&self.kmer_multiset, &othr.kmer_multiset);

        let fw_matches = find_kmer_matches(&self.kmer_positions, &othr.kmer_positions,
            min(self.rare_kmers, othr.rare_kmers) as usize);
        // let fw_lcs = lcskpp(&fw_matches, usize::from(self.k)).path;
        let fw_lcs = find_lcs(self.seq.len(), othr.seq.len(), &fw_matches);
        let len2 = othr.seq.len() as u32;
        let mut rv_matches: Vec<_> = fw_matches.iter().rev().map(|(pos1, pos2)| (*pos1, len2 - 1 - *pos2)).collect();
        rv_matches.sort_unstable();
        // let rv_lcs = lcskpp(&rv_matches, usize::from(self.k)).path;
        let rv_lcs = find_lcs(self.seq.len(), othr.seq.len(), &rv_matches);

        let (lcs, matches) = if fw_lcs.len() >= rv_lcs.len() {
            write!(f, "+")?;
            (fw_lcs, fw_matches)
        } else {
            write!(f, "-")?;
            (rv_lcs, rv_matches)
        };
        // let lcs_indices: std::collections::HashSet<_> = lcs.iter().cloned().collect();
        // println!("All matches:");
        // for (i, (pos1, pos2)) in matches.iter().enumerate() {
            // println!("    {}, {}  ({})", pos1, pos2, lcs_indices.contains(&i));
        // }
        // let nonoverl_lcs = nonoverl_lcs_size(&matches, &lcs, u32::from(self.k));
        let nonoverl_lcs = nonoverl_lcs_size2(&lcs, u32::from(self.k));

        write_dist(f, mset_inters, self.rare_kmers, othr.rare_kmers)?;
        write_dist(f, lcs.len() as u32, self.rare_kmers, othr.rare_kmers)?;
        write_dist(f, nonoverl_lcs, self.nonoverl_kmers, othr.nonoverl_kmers)?;
        writeln!(f)
    }
}

/// Calculates the size of the multiset intersection.
/// `counts*`: key = canonical k-mer, value = number of k-mer occurances in the contig sequence.
fn multiset_intersection(counts1: &IntMap<u64, u16>, counts2: &IntMap<u64, u16>) -> u32 {
    let mut inters = 0;
    for (key, &val1) in counts1.iter() {
        if let Some(&val2) = counts2.get(key) {
            inters += u32::from(min(val1, val2));
        }
    }
    inters
}

/// Write to f: `inters/size1,size2  jaccard_distance`.
fn write_dist<W: Write>(f: &mut W, inters: u32, size1: u32, size2: u32) -> io::Result<()> {
    assert!(inters <= min(size1, size2));
    let dist = 1.0 - f64::from(inters) / f64::from(size1 + size2 - inters);
    write!(f, "\t{}/{},{}\t{:.8}", inters, size1, size2, dist)
}

/// Finds all matching positions between two sets.
/// `positions*`: key = canonical k-mer, value =
fn find_kmer_matches(
        positions1: &IntMap<u64, Vec<u32>>,
        positions2: &IntMap<u64, Vec<u32>>,
        start_capacity: usize,
) -> Vec<(u32, u32)> {
    let mut matches = Vec::with_capacity(start_capacity);
    for (kmer, kmer_pos1) in positions1.iter() {
        if let Some(kmer_pos2) = positions2.get(kmer) {
            matches.extend(kmer_pos1.iter()
                .flat_map(move |&pos1| kmer_pos2.iter().map(move |&pos2| (pos1, pos2))));
        }
    }
    matches.sort_unstable();
    matches
}

fn find_lcs(n: usize, m: usize, matches: &[(u32, u32)]) -> Vec<(u32, u32)> {
    let n1 = n + 1;
    let matches: IntSet<u32> = matches.iter().map(|&(i, j)| j * n1 as u32 + i).collect();
    let m1 = m + 1;
    let mut a = vec![0_u32; n1 * m1];
    let mut p = vec![usize::MAX; n1 * m1];

    for ix in n1 + 1..n1 * m1 {
        if ix % n1 == 0 {
            continue;
        }
        let ix_up = ix - 1;
        let ix_left = ix - n1;
        let ix_diag = ix_left - 1;
        if a[ix_up] > a[ix_left] {
            a[ix] = a[ix_up];
            p[ix] = ix_up;
        } else {
            a[ix] = a[ix_left];
            p[ix] = ix_left;
        }

        if a[ix_diag] + 1 > a[ix] && matches.contains(&(ix as u32)) {
            a[ix] = a[ix_diag] + 1;
            p[ix] = ix_diag;
        }
    }

    let mut ix = n1 * m1 - 1;
    let mut history = Vec::with_capacity(a[ix] as usize);
    while p[ix] != usize::MAX {
        if p[ix] == ix - n1 - 1 {
            history.push(((ix / n1) as u32, (ix % n1) as u32));
        }
        ix = p[ix as usize];
    }
    history.reverse();
    history
}

/// Calculates the size of non-overlapping LCS.
fn nonoverl_lcs_size(matches: &[(u32, u32)], path: &[usize], k: u32) -> u32 {
    let mut next_pos1 = 0;
    let mut next_pos2 = 0;
    let mut count = 0;
    for &i in path.iter() {
        let (pos1, pos2) = matches[i];
        if pos1 >= next_pos1 && pos2 >= next_pos2 {
            next_pos1 = pos1 + k;
            next_pos2 = pos2 + k;
            count += 1;
        }
    }
    count
}

/// Calculates the size of non-overlapping LCS.
fn nonoverl_lcs_size2(lcs: &[(u32, u32)], k: u32) -> u32 {
    let mut next_pos1 = 0;
    let mut next_pos2 = 0;
    let mut count = 0;
    for &(pos1, pos2) in lcs.iter() {
        if pos1 >= next_pos1 && pos2 >= next_pos2 {
            next_pos1 = pos1 + k;
            next_pos2 = pos2 + k;
            count += 1;
        }
    }
    count
}

pub fn find_differences<W, I>(
    f: &mut W,
    contigs: &ContigNames,
    ref_seqs: &[Vec<u8>],
    kmer_counts: &KmerCounts,
    rare_threshold: u8,
    pairs: I,
) -> io::Result<()>
where W: Write,
      I: Iterator<Item = (ContigId, ContigId)>,
{
    writeln!(f, "# k: {}", kmer_counts.k())?;
    writeln!(f, "# threshold: {}", rare_threshold)?;
    writeln!(f, "contig1\tcontig2\tstrand\tmultisets\tmultiset_dist\tlcs\tlcs_dist\tnonoverl_lcs\tnonoverl_lcs_dist")?;
    let init_fn = |contig_id| {
        move || ContigKmers::new(contig_id, &ref_seqs[contig_id.ix()], kmer_counts, rare_threshold)
    };
    let contig_kmers: Vec<_> = (0..contigs.len()).map(|_| OnceCell::new()).collect();
    for (contig1, contig2) in pairs {
        let contig_kmers1 = contig_kmers[contig1.ix()].get_or_init(init_fn(contig1));
        let contig_kmers2 = contig_kmers[contig2.ix()].get_or_init(init_fn(contig2));
        contig_kmers1.compare(contig_kmers2, contigs, f)?;
    }
    Ok(())
}
