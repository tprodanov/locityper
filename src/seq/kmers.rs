use std::{
    io::{self, Error, ErrorKind},
    cmp::{min, max},
};
use super::contigs::{ContigId, ContigNames};

/// Stores k-mer counts for each input k-mer.
pub struct KmerCounts {
    k: u32,
    counts: Vec<Vec<u8>>,
}

impl KmerCounts {
    /// Load k-mer counts from a string.
    /// First line is "k=<number>". All consecutive lines contain just a single number.
    /// Must contain exact number of k-mer as `contigs`.
    pub fn load(full_contents: &str, contigs: &ContigNames) -> io::Result<Self> {
        assert!(!contigs.is_empty(), "Cannot load k-mer counts for empty contigs set!");
        let mut split = full_contents.split('\n');
        let first = split.next().ok_or_else(|| Error::new(ErrorKind::InvalidData,
            "Empty file with k-mer counts!"))?;
        let k = match first.split_once('=') {
            Some(("k", v)) => v.parse().map_err(|e: std::num::ParseIntError|
                Error::new(ErrorKind::InvalidData, e))?,
            _ => return Err(Error::new(ErrorKind::InvalidData,
                format!("Incorrect k-mer counts format, first line ({}) must be in format \"k=integer\"", first))),
        };

        let mut counts = Vec::with_capacity(contigs.len());
        for contig_len in contigs.lengths() {
            assert!(contig_len >= k, "Contig too short!");
            let n_kmers = contig_len - k + 1;
            let mut curr_counts = Vec::with_capacity(n_kmers as usize);
            for _ in 0..n_kmers {
                let val: u8 = split.next()
                    .ok_or_else(|| Error::new(ErrorKind::InvalidData, "Not enough k-mer counts!"))?
                    .parse().map_err(|e: std::num::ParseIntError| Error::new(ErrorKind::InvalidData, e))?;
                // We assume that each k-mer appears at least once.
                curr_counts.push(max(val, 1));
            }
            counts.push(curr_counts);
        }
        match split.next() {
            Some("") | None => Ok(Self { k, counts }),
            _ => Err(Error::new(ErrorKind::InvalidData, "Too many k-mer counts!")),
        }
    }

    /// Returns k-mer size.
    pub fn k(&self) -> u32 {
        self.k
    }

    /// Returns all k-mer counts for the corresponding contig.
    /// Length: contig len - k + 1.
    pub fn get(&self, contig_id: ContigId) -> &[u8] {
        &self.counts[contig_id.ix()]
    }
}

/// k-mer that contains Ns.
pub const N_KMER: u64 = u64::MAX;

/// Iterates over canonical k-mers in the sequence.
/// Returns N_KMER for k-mers containing N.
/// K-mer size `k` must be at most 32.
///
/// Canonical k-mer is calculated as minimal between the k-mer and its rev.complement.
pub fn canonical_kmers(seq: &[u8], k: u8) -> Vec<u64> {
    assert!(0 < k && k <= 32, "k-mer size must be within [1, 32].");
    let k_usize = usize::from(k);
    assert!(seq.len() >= k_usize, "Sequence is too short!");
    let mask: u64 = if k == 32 { -1_i64 as u64 } else { (1_u64 << 2 * k) - 1 };
    let rv_shift = 2 * k - 2;
    let mut fw_kmer: u64 = 0;
    let mut rv_kmer: u64 = 0;
    let mut reset = k_usize - 1;
    let mut kmers = Vec::with_capacity(seq.len() - k_usize + 1);

    for (i, &nt) in seq.iter().enumerate() {
        let fw_enc: u64 = match nt {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                reset = i + k_usize;
                if i + 1 >= k_usize {
                    kmers.push(N_KMER);
                }
                continue;
            },
        };
        let rv_enc = 3 - fw_enc;
        fw_kmer = (fw_kmer << 2) | fw_enc;
        rv_kmer = (rv_enc << rv_shift) | (rv_kmer >> 2);

        if i >= reset {
            kmers.push(min(fw_kmer & mask, rv_kmer));
        } else if i + 1 >= k_usize {
            kmers.push(N_KMER);
        }
    }
    debug_assert_eq!(kmers.len(), seq.len() - k_usize + 1);
    kmers
}
