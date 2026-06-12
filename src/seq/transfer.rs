//! Transfer alignments from one haplotype to another.

use std::{
    path::Path,
    ops::Range,
};
use crate::{
    seq::{
        contigs::{ContigId, ContigNames, ContigSet},
        cigar::SearchableCigar,
        dist::PafFile,
        aln::Alignment,
        wfa::Aligner,
    },
    model::locs::MateData,
    ext::{
        self,
        TriangleMatrix,
    },
};

/// Haplotype-haplotype alignments.
pub struct HapAlns {
    /// Alignments from each contig to other contigs.
    /// matrix[i, j] contains alignment where contig i is query, j is reference.
    aln_matrix: TriangleMatrix<Option<SearchableCigar>>,

    /// For each contig, contains list of other indices and number of matches.
    /// Inner lists are ordered by decreasing edit distance.
    best_ixs: Vec<Vec<(ContigId, u32)>>,
}

impl HapAlns {
    /// Loads alignments from a PAF file, keep <= `n_alns` best alignments per haplotype.
    pub fn load(
        filename: &Path,
        contigs: &ContigNames,
        max_div: f64,
    ) -> crate::Result<Self>
    {
        // [TODO] Replace leave-out haplotypes.
        let mut aln_matrix = TriangleMatrix::new(contigs.len(), None);
        let mut best_ixs = vec![Vec::new(); contigs.len()];
        let min_simil = 1.0 - max_div;
        let mut file = ext::sys::open(filename).map(PafFile::new)?;
        while let Some(entry) = file.next() {
            let entry = entry?;
            let hap1: &str = entry.query_name();
            let hap2: &str = entry.target_name();
            let Some(id1) = contigs.try_get_id(hap1) else { continue };
            let Some(id2) = contigs.try_get_id(hap2) else { continue };
            if id1 == id2 { continue };
            let matrix_cell = aln_matrix.get_symmetric_mut(id1.ix(), id2.ix());
            if matrix_cell.is_some() { continue };

            if !entry.full_positive_alignment()? {
                log::warn!("Alignment between {} and {} is on the reverse strand or does not fully cover both sequences",
                hap1, hap2);
                continue
            }
            let n_matches = entry.n_matches()?;
            let simil = f64::from(n_matches) / f64::from(entry.aln_len()?);
            if simil < min_simil { continue };
            let Some(cigar) = entry.cigar().transpose()? else { continue };
            *matrix_cell = Some(SearchableCigar::new(&cigar, id1 > id2));
            best_ixs[id1.ix()].push((id2, n_matches));
            best_ixs[id2.ix()].push((id1, n_matches));
        }

        best_ixs.iter_mut().for_each(|v| v.sort_by(|a, b| b.1.cmp(&a.1)));
        Ok(Self { aln_matrix, best_ixs })
    }

    pub fn transfer_alignments(
        &self,
        alns: &mut Vec<Alignment>,
        start_ix: usize,
        mate_data: &MateData,
        contig_set: &ContigSet,
        aligner: &Aligner,
    ) -> crate::Result<()>
    {
        let mut intervals = iset::IntervalSet::new();
        // log::debug!("    Initial alignments");
        for aln in &alns[start_ix..] {
            // log::debug!("        * {:15} {:?} {}", aln.interval(), aln.cigar(), aln.strand());
            intervals.insert(aln_interval(aln));
        }
        let first_aln = &alns[start_ix];
        log::debug!("    First alignment: {} {:?} {}", first_aln.interval(), first_aln.cigar(), first_aln.strand());
        let first_seq = mate_data.get_seq(first_aln.strand());
        let first_id = first_aln.contig_id();
        for &(oth_contig_id, _) in &self.best_ixs[first_id.ix()] {
            let cigar = self.aln_matrix.get_symmetric(first_id.ix(), oth_contig_id.ix())
                .as_ref().expect("Alignment between haplotypes must exist");
            let query_to_ref = first_id < oth_contig_id;
            let cursor = cigar.find_position_dynamic(first_aln.interval().start(), query_to_ref);
            log::debug!("        to {}: cursor {:?}", contig_set.contigs().get_name(oth_contig_id), cursor);
            let (_rstart, _rend, new_cigar) = cigar.transfer_alignment_dynamic(cursor, first_aln.cigar(), first_seq,
                contig_set.get_seq(oth_contig_id), aligner, query_to_ref);
            log::debug!("        -> {:?}", new_cigar);
        }
        Ok(())
    }
}

/// Convert alignment intervals into ranges of u64, useful for interval set.
#[inline(always)]
fn aln_interval(aln: &Alignment) -> Range<u64> {
    let prefix = u64::from(aln.contig_id().get()) << 32;
    (prefix | u64::from(aln.interval().start()))..(prefix | u64::from(aln.interval().end()))
}
