use std::{
    io::{self, Error, ErrorKind},
    cmp::max,
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
    pub fn kmer_size(&self) -> u32 {
        self.k
    }

    /// Returns all k-mer counts for the corresponding contig.
    /// Length: contig len - k + 1.
    pub fn counts(&self, contig_id: ContigId) -> &[u8] {
        &self.counts[contig_id.ix()]
    }
}