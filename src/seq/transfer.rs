//! Transfer alignments from one haplotype to another.

use std::{
    path::Path,
    cmp::{min, max},
};
use crate::{
    seq::{
        ContigId, ContigNames,
        cigar::SearchableCigar,
        dist::PafFile,
    },
    ext::{
        self,
    },
    algo::IntSet,
};

/// Haplotype-haplotype alignments.
pub struct HapHapAlns {
    /// Alignments from each contig. Tuples contain (second contig, total number of matches, cigar).
    alns: Vec<Vec<(ContigId, u32, SearchableCigar)>>,
}

impl HapHapAlns {
    /// Loads alignments from a PAF file, keep <= `n_alns` best alignments per haplotype.
    pub fn from_file(
        filename: &Path,
        contigs: &ContigNames,
        max_div: f64,
    ) -> crate::Result<Self>
    {
        // [TODO] Replace leave-out haplotypes.
        let mut alns = vec![Vec::new(); contigs.len()];
        let min_simil = 1.0 - max_div;
        let mut seen = IntSet::default();
        let mut file = ext::sys::open(filename).map(PafFile::new)?;
        while let Some(entry) = file.next() {
            let entry = entry?;
            let hap1: &str = entry.query_name();
            let hap2: &str = entry.target_name();
            let Some(id1) = contigs.try_get_id(hap1) else { continue };
            let Some(id2) = contigs.try_get_id(hap2) else { continue };
            if id1 == id2 { continue };
            let key = (u32::from(min(id1.get(), id2.get())) << 16) | u32::from(max(id1.get(), id2.get()));
            if !seen.insert(key) { continue };

            if !entry.full_positive_alignment()? {
                log::warn!("Alignment between {} and {} is on the reverse strand or does not fully cover both sequences",
                hap1, hap2);
                continue
            }
            let n_matches = entry.n_matches()?;
            let simil = f64::from(n_matches) / f64::from(entry.aln_len()?);
            if simil < min_simil { continue };
            let Some(cigar) = entry.cigar().transpose()? else { continue };
            let cigar = SearchableCigar::new(&cigar);
            alns[id1.ix()].push((id2, n_matches, cigar.invert()));
            alns[id2.ix()].push((id1, n_matches, cigar));
        }

        for hap_alns in &mut alns {
            hap_alns.sort_by(|a, b| b.1.cmp(&a.1));
        }
        Ok(Self { alns })
    }
}