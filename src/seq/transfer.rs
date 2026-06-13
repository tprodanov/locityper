//! Transfer alignments from one haplotype to another.

use std::{
    path::Path,
    ops::Range,
    sync::Arc,
};
use crate::{
    seq::{
        Interval,
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

    fn transfer_alignment<const TO_REF: bool>(
        &self,
        source_id: ContigId,
        target_id: ContigId,
        read_source_aln: &Alignment,
        read_seq: &[u8],
        contig_set: &ContigSet,
        aligner: &Aligner,
    ) {
        let target_seq = contig_set.get_seq(target_id);
        let haps_cigar = self.aln_matrix.get_symmetric(source_id.ix(), target_id.ix())
            .as_ref().expect("Alignment between haplotypes must exist");
        let source_aln_start = read_source_aln.interval().start();
        let (cigar_lbound, cigar_rbound, approx_pos) = haps_cigar.find_approx_position::<TO_REF>(source_aln_start);
        log::debug!("        : Approx. start = {},  cigar bounds: {}..{}", approx_pos, cigar_lbound, cigar_rbound);
        let (start, end, new_cigar) = haps_cigar.transfer_alignment::<TO_REF>(
            source_aln_start, cigar_lbound, cigar_rbound, read_source_aln.cigar(), read_seq, target_seq, aligner);
        log::debug!("        : {:?}", new_cigar);
        let new_interval = Interval::new(Arc::clone(contig_set.contigs()), target_id, start, end);
        log::debug!("        : {}", new_interval);
    }

    pub fn transfer_alignments(
        &self,
        alns: &mut Vec<Alignment>,
        start_ix: usize,
        mate_data: &MateData,
        contig_set: &ContigSet,
        aligner: &Aligner,
    ) {
        let mut intervals = iset::IntervalSet::new();
        // log::debug!("    Initial alignments");
        for aln in &alns[start_ix..] {
            // log::debug!("        * {:15} {:?} {}", aln.interval(), aln.cigar(), aln.strand());
            intervals.insert(aln_interval(aln));
        }
        let first_aln = &alns[start_ix];
        log::debug!("    First alignment: {} {:?} {}", first_aln.interval(), first_aln.cigar(), first_aln.strand());
        let read_seq = mate_data.get_seq(first_aln.strand());
        let source_id = first_aln.contig_id();
        for &(target_id, _) in &self.best_ixs[source_id.ix()] {
            log::debug!("        Transferring to [{}] {}", target_id.get(), contig_set.contigs().get_name(target_id));
            if source_id < target_id {
                self.transfer_alignment::<true>(source_id, target_id, first_aln, read_seq, contig_set, aligner);
            } else {
                self.transfer_alignment::<false>(source_id, target_id, first_aln, read_seq, contig_set, aligner);
            }
        }
    }
}

/// Convert alignment intervals into ranges of u64, useful for interval set.
#[inline(always)]
fn aln_interval(aln: &Alignment) -> Range<u64> {
    let prefix = u64::from(aln.contig_id().get()) << 32;
    (prefix | u64::from(aln.interval().start()))..(prefix | u64::from(aln.interval().end()))
}
