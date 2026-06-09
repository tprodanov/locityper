//! Transfer alignments from one haplotype to another.

use std::{
    path::Path,
};

use crate::{
    seq::cigar::Cigar,
};

/// Haplotype-haplotype alignments.
pub struct HapHapAlns {
    /// Alignments for each contig, sorted by edit distance.
    alns: Vec<Vec<Cigar>>,
}

impl HapHapAlns {
    /// Loads alignments from a PAF file, keep <= `n_alns` best alignments per haplotype.
    pub fn from_file(filename: &Path, n_alns: usize) -> crate::Result<Self> {
        todo!()
    }
}