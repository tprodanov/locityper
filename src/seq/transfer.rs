//! Transfer alignments from one haplotype to another.

use std::{
    sync::Arc,
    cmp::min,
};
use crate::{
    seq::{
        Interval,
        contigs::{ContigId, ContigNames, ContigSet},
        cigar::{Cigar, CigarIndex, QueryToRef, RefToQuery},
        paf::{PafEntry},
        aln::Alignment,
        wfa::Aligner,
    },
    bg::ErrorProfile,
    model::locs::{ReadData, PrelimAlignments},
    ext::{
        TriangleMatrix,
        opt::With,
    }
};

/// Haplotype-haplotype alignments.
pub struct HapAlns {
    /// Alignments from each contig to other contigs.
    /// matrix[i, j] contains alignment where contig i is query, j is reference.
    aln_matrix: TriangleMatrix<Option<(Cigar, CigarIndex)>>,
    transfer_fails: u32,

    local_step: u32,
    /// Haplotype -> Window -> offset in `locally_similar_ixs`. Inner vector has size of (n_windows + 1).
    contig_offsets: Vec<Vec<usize>>,
    /// Slice [contig_offset[contig][window] .. contig_offset[contig][window + 1]] will contain a sorted list of
    /// other locally similar contigs.
    locally_similar_ixs: Vec<ContigId>,
}

impl HapAlns {
    pub fn new(n_contigs: usize, transfer_fails: u32) -> Self {
        Self {
            aln_matrix: TriangleMatrix::new(n_contigs, None),
            transfer_fails,
            local_step: 0,
            contig_offsets: Vec::new(),
            locally_similar_ixs: Vec::new(),
        }
    }

    pub fn add(&mut self, entry: &mut PafEntry, contigs: &ContigNames, max_div: f64) {
        let id1 = entry.query_id();
        let id2 = entry.target_id();
        let matrix_cell = self.aln_matrix.get_symmetric_mut(id1.ix(), id2.ix());
        if matrix_cell.is_some() { return };

        if !entry.full_positive_alignment() {
            log::warn!("Alignment between {} and {} is on the reverse strand or does not fully cover both sequences",
                contigs.get_name(id1), contigs.get_name(id2));
            return;
        }
        if entry.divergence().unwrap_or(f64::INFINITY) > max_div { return };
        let Some(mut cigar) = entry.take_cigar() else { return };
        if id1 > id2 {
            cigar = cigar.invert();
        }
        let cigar_index = CigarIndex::new(&cigar);
        *matrix_cell = Some((cigar, cigar_index));
    }

    /// Identify closest local haplotypes.
    pub fn create_locally_similar_lists(
        &mut self,
        contigs: &ContigNames,
        basis_contigs: &[ContigId],
        window: u32,
        step: u32,
        max_div: f64,
    ) {
        let max_edit = (f64::from(window) * max_div).floor() as u32;
        let n_contigs = contigs.len();
        let mut contig_offsets = vec![Vec::new(); n_contigs];
        let mut locally_similar_ixs = Vec::new();

        let mut buffer = Vec::new();
        let max_len = contigs.lengths().iter().max().expect("At least one contig must be present");
        // Temporary structure: window index -> [(haplotype index, edit size)]
        let mut contig_locally_similar = vec![Vec::new(); max_len.saturating_sub(window).div_ceil(step) as usize + 1];
        for &contig_id in basis_contigs {
            let contig_len = contigs.get_len(contig_id);
            let n_windows = contig_len.saturating_sub(window).div_ceil(step) as usize + 1;
            // This way we will for certain not get out of bounds.
            let curr_locally_similar = &mut contig_locally_similar[..n_windows];

            let i = contig_id.ix();
            for j in 0..n_contigs {
                if i == j { continue };
                let Some((cigar, _)) = self.aln_matrix.get_symmetric(i, j).as_ref() else { continue };
                let in_query = i < j;
                cigar.locally_similar::<(u32, u32)>(in_query, window, step, max_edit, &mut buffer);
                for &(window, edit) in &buffer {
                    curr_locally_similar[window as usize].push((j, edit));
                }
            }

            let curr_offsets = &mut contig_offsets[i];
            curr_offsets.reserve_exact(n_windows + 1);
            curr_offsets.push(locally_similar_ixs.len());
            for similar_at_window in curr_locally_similar {
                similar_at_window.sort_unstable_by_key(|tup| tup.1);
                locally_similar_ixs.extend(similar_at_window.iter().map(|tup| ContigId::new(tup.0)));
                curr_offsets.push(locally_similar_ixs.len());
                similar_at_window.clear();
            }
            assert!(curr_offsets.len() == n_windows + 1);
        }

        self.local_step = step;
        self.contig_offsets = contig_offsets;
        self.locally_similar_ixs = locally_similar_ixs;
    }

    /// Try to identify any new alignments for a given read/read pair.
    /// Returns the number of new alignments and a total number of attempted transfers.
    pub(crate) fn transfer_alignments(
        &self,
        prelim_alignments: &mut PrelimAlignments,
        read_data: &ReadData,
        contig_set: &ContigSet,
        aligner: &Aligner,
        err_prof: &ErrorProfile,
    ) -> (usize, usize) {
        let contigs = contig_set.contigs();
        let n = prelim_alignments.len();
        let mut seen = vec![false; n];
        let mut total_transfers = 0;

        for i in 0..n {
            if std::mem::replace(&mut seen[i], true) {
                continue;
            }
            let source_aln = &prelim_alignments[i];
            // Copy various variables before `source_aln` is dropped to make `prelim_alignments` mutable.
            let source_contig_id = source_aln.contig_id();
            let source_aln_start = source_aln.interval().start();
            let source_aln_end = source_aln.interval().end();
            let source_cigar = source_aln.cigar().clone();
            let source_strand = source_aln.strand();
            let read_end = source_aln.read_end();
            let read_seq = read_data.mate_data(read_end).get_seq(source_strand);
            let mut n_fails = 0..self.transfer_fails;

            let (left_clipping, right_clipping) = source_cigar.soft_clipping();

            let curr_offsets = &self.contig_offsets[source_contig_id.ix()];
            assert!(!curr_offsets.is_empty(), "Trying to reconstruct alignments not from a basis haplotype");
            let n_windows = curr_offsets.len().strict_sub(1);
            // If we are in the last window, start / step can put us beyond the size of the vector.
            let k = min(n_windows, (source_aln_start / self.local_step + 1) as usize);
            // Safety: we know that k is at most n_windows, and length = n_windows + 1.
            let start_offset = unsafe { *curr_offsets.get_unchecked(k.strict_sub(1)) };
            let end_offset = unsafe { *curr_offsets.get_unchecked(k) };
            for &target_contig_id in &self.locally_similar_ixs[start_offset..end_offset] {
                let target_seq = contig_set.get_seq(target_contig_id);
                let (haps_cigar, cigar_index) = self.aln_matrix
                    .get_symmetric(source_contig_id.ix(), target_contig_id.ix())
                    .as_ref().expect("Alignment between haplotypes must exist");
                let query_to_ref = source_contig_id < target_contig_id;

                let offset = if query_to_ref {
                    cigar_index.find_cigar_offset(source_aln_start, QueryToRef)
                } else {
                    cigar_index.find_cigar_offset(source_aln_start, RefToQuery)
                };
                if let Some(prev_aln_entry) = prelim_alignments.pos_collection()
                    .get(read_end, target_contig_id, offset.rpos())
                {
                    if let Some(seen_val) = seen.get_mut(prev_aln_entry.index() as usize) {
                        *seen_val = true;
                    }
                    if !prev_aln_entry.allowed_to_improve() {
                        continue;
                    }
                }

                let (new_start, new_cigar) = if query_to_ref {
                    haps_cigar.transfer_read_alignment(offset, &source_cigar, read_seq, target_seq,
                        aligner, QueryToRef)
                } else {
                    haps_cigar.transfer_read_alignment(offset, &source_cigar, read_seq, target_seq,
                        aligner, RefToQuery)
                };

                {
                    // How much to add to the sides of the reference sequence
                    const PADDING: u32 = 10;
                    // Inside the region corresponding to the reference alignment,
                    // allow this size of free indels on the sides.
                    const FREE_END: u32 = 10;
                    let target_len = target_seq.len() as u32;
                    let target_end = if query_to_ref {
                        cigar_index.find_cigar_offset(source_aln_end - 1, QueryToRef)
                    } else {
                        cigar_index.find_cigar_offset(source_aln_end - 1, RefToQuery)
                    }.rpos() + 1;
                    let padded_target_start = offset.rpos().saturating_sub(PADDING);
                    let padded_target_end = min(target_end + PADDING, target_len);
                    let target_subseq = &target_seq[padded_target_start as usize..padded_target_end as usize];
                    let target_free_left = FREE_END + offset.rpos() - padded_target_start;
                    let target_free_right = FREE_END + padded_target_end.saturating_sub(target_end);
                    let mut alt_cigar = Cigar::new();
                    let alt_start = padded_target_start + aligner.align_free_ends(target_subseq, read_seq, &mut alt_cigar,
                        With((target_free_left as i32, (FREE_END + left_clipping) as i32)),
                        With((target_free_right as i32, (FREE_END + right_clipping) as i32)));
                    alt_cigar.boundary_ins_to_soft();
                    log::debug!("{} R{} to {}.   Approx {}-{}", read_data.name(), read_end.ix() + 1,
                        contigs.get_name(target_contig_id), offset.rpos() + 1, target_end);
                    log::debug!("    {:4}  {}", new_start, new_cigar);
                    log::debug!("    {:4}  {}", alt_start, alt_cigar);
                    if new_cigar != alt_cigar {
                        log::warn!("    Different");
                    }
                }

                total_transfers += 1;
                let cigar_ref_len = new_cigar.ref_len();
                const MIN_ALN_SIZE: u32 = 50;
                // Either too high edit distance or too short alignment.
                if cigar_ref_len < MIN_ALN_SIZE {
                    prelim_alignments.record_bad_alignment(read_end, target_contig_id, new_start, false);
                    // Break if produced too many fails for this specific alignment.
                    if n_fails.next().is_none() { break } else { continue };
                }

                let interval = Interval::new(Arc::clone(contigs), target_contig_id,
                    new_start, new_start + cigar_ref_len);
                let mut new_aln = Alignment::new(interval, new_cigar, source_strand, read_end);
                new_aln.calculate_edit_and_prob(contigs, err_prof);
                if let Err(opt_prev_ix) = prelim_alignments.insert(new_aln, false) {
                    if let Some(j) = opt_prev_ix && let Some(seen_val) = seen.get_mut(j as usize) {
                        *seen_val = true;
                    }
                    if n_fails.next().is_none() {
                        break;
                    }
                }
            }
        }
        (prelim_alignments.len() - n, total_transfers)
    }
}
