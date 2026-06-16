//! Transfer alignments from one haplotype to another.

use std::{
    path::Path,
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
    bg::ErrorProfile,
    model::locs::{ReadData, PrelimAlignments},
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
    ) -> crate::Result<Self> {
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
            // if std::cmp::min(id1, id2).get() != 270 || std::cmp::max(id1, id2).get() != 356 { continue };
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

    /// Try to identify any new alignments for a given read/read pair.
    /// Returns the number of new alignments.
    pub(crate) fn transfer_alignments(
        &self,
        prelim_alignments: &mut PrelimAlignments,
        read_data: &ReadData,
        contig_set: &ContigSet,
        aligner: &Aligner,
        err_prof: &ErrorProfile,
    ) -> usize {
        let n = prelim_alignments.len();
        let mut seen = vec![false; n];

        for i in 0..n {
            if std::mem::replace(&mut seen[i], true) { continue }
            let source_aln = &prelim_alignments[i];
            // Copy various variables before `source_aln` is dropped to make `prelim_alignments` mutable.
            let source_contig_id = source_aln.contig_id();
            let source_aln_start = source_aln.interval().start();
            let source_cigar = source_aln.cigar().clone();
            let source_strand = source_aln.strand();
            let read_end = source_aln.read_end();
            let read_seq = read_data.mate_data(read_end).get_seq(source_strand);
            let passable_dist = prelim_alignments.passable_dist(read_end);

            for &(target_contig_id, _) in &self.best_ixs[source_contig_id.ix()] {
                let target_seq = contig_set.get_seq(target_contig_id);
                let haps_cigar = self.aln_matrix.get_symmetric(source_contig_id.ix(), target_contig_id.ix())
                    .as_ref().expect("Alignment between haplotypes must exist");
                let query_to_ref = source_contig_id < target_contig_id;

                let (min_cigar_ix, max_cigar_ix, approx_start) = if query_to_ref {
                    haps_cigar.find_approx_position::<true>(source_aln_start)
                } else {
                    haps_cigar.find_approx_position::<false>(source_aln_start)
                };
                if let Some(ix) = prelim_alignments.pos_collection().get(read_end, target_contig_id, approx_start) {
                    // Similar position is observed in alignment `ix`. Do not transfer it later.
                    if let Some(v) = seen.get_mut(ix as usize) {
                        *v = true;
                    }
                    continue;
                }

                let (new_start, new_cigar) = if query_to_ref {
                    haps_cigar.transfer_alignment::<true>(
                        source_aln_start, min_cigar_ix, max_cigar_ix, &source_cigar, read_seq, target_seq, aligner)
                } else {
                    haps_cigar.transfer_alignment::<false>(
                        source_aln_start, min_cigar_ix, max_cigar_ix, &source_cigar, read_seq, target_seq, aligner)
                };
                const MIN_ALN_SIZE: u32 = 50;
                let cigar_ref_len = new_cigar.ref_len();
                // Either too high edit distance or too short alignment.
                if cigar_ref_len.abs_diff(new_cigar.query_len()) > passable_dist || cigar_ref_len < MIN_ALN_SIZE {
                    continue
                }

                let interval = Interval::new(Arc::clone(contig_set.contigs()), target_contig_id,
                    new_start, new_start + cigar_ref_len);
                let new_aln = Alignment::new(interval, new_cigar, source_strand, read_end);
                prelim_alignments.push(new_aln, contig_set.contigs(), err_prof);
            }
        }
        prelim_alignments.len() - n
    }
}
