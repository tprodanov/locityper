use std::rc::Rc;
use crate::{
    math::Ln,
    seq::contigs::{ContigId, ContigNames},
};
use super::locs::{TwoIntervals, ReadPairAlignments};

/// First window in `ContigWindows` represents an unmapped window.
pub(crate) const UNMAPPED_WINDOW: u32 = 0;
/// First contig will have windows starting from index 1.
pub(crate) const INIT_WSHIFT: u32 = 1;

/// Alignment of a read pair to a specific contig windows.
pub struct ReadWindows {
    /// Index in the list of read-pair alignments.
    ix: u32,
    /// Window for each read-end (UNMAPPED_WINDOW if unmapped).
    windows: (u32, u32),
    /// ln-probability of this location.
    ln_prob: f64,
}

impl ReadWindows {
    fn new(ix: u32, windows: (u32, u32), ln_prob: f64) -> Self {
        Self { ix, windows, ln_prob }
    }
}

/// Stores the contigs and windows corresponding to the windows.
pub struct ContigWindows {
    /// Window size.
    window: u32,

    /// Ignore left- and right-most `boundary_size` bp.
    boundary: u32,

    /// All contig names.
    contigs: Rc<ContigNames>,

    /// Windows correspond to these contigs.
    contig_ids: Vec<ContigId>,

    /// ln(contig_ids.len()).
    ln_ploidy: f64,

    /// Starts and ends within each contig, after removing boundary regions.
    bounded_starts: Vec<u32>,
    bounded_ends: Vec<u32>,

    /// Start index for each contig id (length: n-contigs + 1).
    /// Contig with index `i` will have windows between `wshifts[i]..wshifts[i + 1]`.
    wshifts: Vec<u32>,
}

impl ContigWindows {
    pub fn new(window: u32, boundary: u32, contig_ids: Vec<ContigId>, contigs: Rc<ContigNames>) -> Self {
        let n = contig_ids.len();
        let mut bounded_starts = Vec::with_capacity(n);
        let mut bounded_ends = Vec::with_capacity(n);

        let mut wshifts = Vec::with_capacity(n + 1);
        let mut curr_wshift = INIT_WSHIFT;
        wshifts.push(curr_wshift);

        for &contig_id in contig_ids.iter() {
            let contig_len = contigs.get_len(contig_id);
            let n_windows = contig_len.saturating_sub(2 * boundary) / window;
            assert!(n_windows > 0, "Contig {} is too short (len = {})", contigs.get_name(contig_id), contig_len);

            let start = (contig_len - window * n_windows) / 2;
            let end = start + window * n_windows;
            assert!(start >= boundary && end + boundary <= contig_len,
                "Incorrect calculations! len {}, window {}, boundary {} -> n_windows {}, start {}, end {}",
                contig_len, window, boundary, n_windows, start, end);
            bounded_starts.push(start);
            bounded_ends.push(end);
            curr_wshift += n_windows;
            wshifts.push(curr_wshift);
        }
        let ln_ploidy = (contig_ids.len() as f64).ln();
        Self { window, boundary, contigs, contig_ids, ln_ploidy, bounded_starts, bounded_ends, wshifts }
    }

    /// Given all pair-alignments for a single read pair,
    /// finds all alignments corresponding to `self`, and appends them to `out_alns`.
    /// Returns the number of added alignments.
    ///
    /// Output pair-alignments with ln-probability worse than `best_prob - prob_diff` are discarded.
    /// Remaining alignments have random order, probabilities are not normalized.
    ///
    /// Any alignments that were in the vector before, stay as they are and in the same order.
    pub fn read_alignments(&self, pair_alns: &ReadPairAlignments,
        out_alns: &mut Vec<ReadWindows>,
        prob_diff: f64,
    ) -> usize {
        let start_len = out_alns.len();
        // Probability of being unmapped to any of the contigs.
        let mut unmapped_prob = pair_alns.unmapped_prob() + self.ln_ploidy;
        // Current threshold, is updated during the for-loop.
        let mut thresh_prob = unmapped_prob - prob_diff;
        for (i, &contig_id) in self.contig_ids.iter().enumerate() {
            let wshift = self.wshifts[i];
            let start = self.bounded_starts[i];
            let end = self.bounded_ends[i];

            let (start_ix, palns) = pair_alns.contig_alns(contig_id);
            for (aln_ix, paln) in (start_ix..).zip(palns) {
                let windows = match paln.intervals() {
                    TwoIntervals::Both(interval1, interval2) => {
                        let middle1 = interval1.middle();
                        let middle2 = interval2.middle();
                        if middle1 < start || middle1 >= end || middle2 < start || middle2 >= end {
                            unmapped_prob = Ln::add(unmapped_prob, paln.ln_prob());
                            continue;
                        } else {
                            (wshift + (middle1 - start) / self.window, wshift + (middle2 - start) / self.window)
                        }
                    },
                    TwoIntervals::First(interval1) => {
                        let middle1 = interval1.middle();
                        if middle1 < start || middle1 >= end {
                            unmapped_prob = Ln::add(unmapped_prob, paln.ln_prob());
                            continue;
                        } else {
                            (wshift + (middle1 - start) / self.window, UNMAPPED_WINDOW)
                        }
                    },
                    TwoIntervals::Second(interval2) => {
                        let middle2 = interval2.middle();
                        if middle2 < start || middle2 >= end {
                            unmapped_prob = Ln::add(unmapped_prob, paln.ln_prob());
                            continue;
                        } else {
                            (UNMAPPED_WINDOW, wshift + (middle2 - start) / self.window)
                        }
                    },
                };
                let ln_prob = paln.ln_prob();
                if ln_prob >= thresh_prob {
                    thresh_prob = thresh_prob.max(ln_prob - prob_diff);
                    out_alns.push(ReadWindows::new(aln_ix as u32, windows, ln_prob));
                }
            }
        }

        // As threshold was updated during the for-loop, some alignments in the beginning may need to removed.
        let mut i = start_len;
        while i < out_alns.len() {
            if out_alns[i].ln_prob < thresh_prob {
                out_alns.swap_remove(i);
            } else {
                i += 1;
            }
        }
        if unmapped_prob >= thresh_prob {
            out_alns.push(ReadWindows::new(u32::MAX, (UNMAPPED_WINDOW, UNMAPPED_WINDOW), unmapped_prob));
        }
        out_alns.len() - start_len
    }
}
