use std::{
    cmp::min,
    io::{self, Write},
};
use once_cell::unsync::OnceCell;
use nohash::IntMap;
use bio::alignment::sparse::lcskpp;
use super::{
    contigs::{ContigId, ContigNames},
    kmers::{KmerCounts, KmerCount, canonical_kmers, N_KMER},
};

struct ContigKmers<'a> {
    contig_id: ContigId,
    seq: &'a [u8],

    k: u8,
    /// Store positions of rare k-mer in static vectors, which are spilled into heap if needed.
    kmer_positions: IntMap<u64, Vec<u32>>,
    /// Total number of rare k-mers, same as the sum over all `kmer_positions.values().len()`.
    rare_kmers: u32,
    /// Maximal number of non-overlapping rare k-mers.
    nonoverl_kmers: u32,
}

impl<'a> ContigKmers<'a> {
    fn new(contig_id: ContigId, seq: &'a [u8], kmer_counts: &'a KmerCounts, threshold: KmerCount) -> Self {
        let ref_occurances = kmer_counts.get(contig_id);
        let k: u8 = kmer_counts.k().try_into().unwrap();
        let canon_kmers = canonical_kmers(seq, k);
        assert_eq!(canon_kmers.len(), ref_occurances.len());

        let mut kmer_positions = IntMap::default();
        let mut rare_kmers = 0;
        let mut nonoverl_kmers = 0;
        let mut next_nonoverl_pos = 0;

        for (pos, (kmer, &ref_count)) in canon_kmers.into_iter().zip(ref_occurances).enumerate() {
            let pos = pos as u32;
            if kmer != N_KMER && ref_count <= threshold {
                kmer_positions.entry(kmer).or_insert_with(|| Vec::with_capacity(2)).push(pos);
                rare_kmers += 1;
                if pos >= next_nonoverl_pos {
                    next_nonoverl_pos = pos + u32::from(k);
                    nonoverl_kmers += 1;
                }
            }
        }
        Self {
            contig_id, seq, k, kmer_positions,
            rare_kmers, nonoverl_kmers,
        }
    }

    /// Compares `self` v `othr` and writes to a CSV stream `f`.
    /// First three columns: contig1, contig2, strand.
    ///
    /// Next, four groups of two columns of type `inters/size1,size2   dist`, where distance
    /// is calculated as `1 - jaccard_index`.
    /// Four groups correspond to sets of k-mers; multisets of k-mers; LCS; and non-overlapping LCS.
    fn compare<W: Write>(&self, othr: &Self, contigs: &ContigNames, f: &mut W) -> io::Result<()> {
        log::debug!("Compare {} & {}", contigs.get_name(self.contig_id), contigs.get_name(othr.contig_id));
        write!(f, "{}\t{}\t", contigs.get_name(self.contig_id), contigs.get_name(othr.contig_id))?;
        let (set_inters, mset_inters) = multiset_intersection(&self.kmer_positions, &othr.kmer_positions);

        let fw_matches = find_kmer_matches(&self.kmer_positions, &othr.kmer_positions,
            min(self.rare_kmers, othr.rare_kmers) as usize);
        let fw_lcs = lcskpp(&fw_matches, usize::from(self.k)).path;
        let len2 = othr.seq.len() as u32;
        let mut rv_matches: Vec<_> = fw_matches.iter().rev().map(|(pos1, pos2)| (*pos1, len2 - 1 - *pos2)).collect();
        rv_matches.sort_unstable();
        let rv_lcs = lcskpp(&rv_matches, usize::from(self.k)).path;

        let lcs = if fw_lcs.len() >= rv_lcs.len() {
            write!(f, "+")?;
            improve_lcs(&fw_matches, &fw_lcs)
        } else {
            write!(f, "-")?;
            improve_lcs(&rv_matches, &rv_lcs)
        };
        let nonoverl_lcs = nonoverl_lcs_size(&lcs, u32::from(self.k));

        write_dist(f, set_inters, self.kmer_positions.len() as u32, othr.kmer_positions.len() as u32)?;
        write_dist(f, mset_inters, self.rare_kmers, othr.rare_kmers)?;
        write_dist(f, lcs.len() as u32, self.rare_kmers, othr.rare_kmers)?;
        write_dist(f, nonoverl_lcs, self.nonoverl_kmers, othr.nonoverl_kmers)?;
        writeln!(f)
    }
}

/// Calculates the size of the set and multiset intersection.
/// `counts*`: key = canonical k-mer, value = number of k-mer occurances in the contig sequence.
fn multiset_intersection(positions1: &IntMap<u64, Vec<u32>>, positions2: &IntMap<u64, Vec<u32>>) -> (u32, u32) {
    let mut set_inters = 0;
    let mut mset_inters = 0;
    for (kmer, kmer_pos1) in positions1.iter() {
        if let Some(kmer_pos2) = positions2.get(kmer) {
            set_inters += 1;
            mset_inters += min(kmer_pos1.len(), kmer_pos2.len()) as u32;
        }
    }
    (set_inters, mset_inters)
}

/// Write to f: `inters/size1,size2  jaccard_distance`.
fn write_dist<W: Write>(f: &mut W, inters: u32, size1: u32, size2: u32) -> io::Result<()> {
    assert!(inters <= min(size1, size2));
    let dist = 1.0 - f64::from(inters) / f64::from(size1 + size2 - inters);
    write!(f, "\t{}/{},{}\t{:.8}", inters, size1, size2, dist)
}

/// Finds all matching positions between two sets.
/// `positions*`: key = canonical k-mer, value = k-mer positions.
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

/// LCSk++ finds an approximate LCS solution.
/// This function improves found LCS by greedily adding matches where possible.
fn improve_lcs(matches: &[(u32, u32)], path: &[usize]) -> Vec<(u32, u32)> {
    let mut possible_pos1 = 0;
    let mut possible_pos2 = 0;
    let mut last_ix = 0;
    let mut lcs = Vec::new();
    for &ix in path {
        let (new_pos1, new_pos2) = matches[ix];
        for &(pos1, pos2) in matches[last_ix..ix].iter() {
            if pos1 >= possible_pos1 && pos2 >= possible_pos2 && pos1 < new_pos1 && pos2 < new_pos2 {
                lcs.push((pos1, pos2));
                possible_pos1 = pos1 + 1;
                possible_pos2 = pos2 + 1;
            }
        }
        last_ix = ix;
    }

    for &(pos1, pos2) in matches[last_ix..].iter() {
        if pos1 >= possible_pos1 && pos2 >= possible_pos2 {
            lcs.push((pos1, pos2));
            possible_pos1 = pos1 + 1;
            possible_pos2 = pos2 + 1;
        }
    }
    lcs
}

/// Calculates the size of non-overlapping LCS.
fn nonoverl_lcs_size(lcs: &[(u32, u32)], k: u32) -> u32 {
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
    rare_threshold: KmerCount,
    pairs: I,
) -> io::Result<()>
where W: Write,
      I: Iterator<Item = (ContigId, ContigId)>,
{
    writeln!(f, "# k: {}", kmer_counts.k())?;
    writeln!(f, "# threshold: {}", rare_threshold)?;
    write!(f, "contig1\tcontig2\tstrand\tsets\tset_dist\tmultisets\tmultiset_dist\t")?;
    writeln!(f, "lcs\tlcs_dist\tnonoverl_lcs\tnonoverl_lcs_dist")?;
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
