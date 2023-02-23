use std::{
    rc::Rc,
    fmt::Write,
};
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
    #[allow(dead_code)]
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

    /// Returns two windows, to which the read pair is aligned.
    /// If unmapped, one of the values if UNMAPPED_WINDOW.
    pub fn windows(&self) -> (u32, u32) {
        self.windows
    }

    /// Returns ln-probability of the alignment.
    pub fn ln_prob(&self) -> f64 {
        self.ln_prob
    }
}

/// Stores the contigs and windows corresponding to the windows.
pub struct ContigWindows {
    /// All contig names.
    contigs: Rc<ContigNames>,
    /// Windows correspond to these contigs.
    ids: Vec<ContigId>,
    /// Copy number of each contig.
    cns: Vec<u8>,
    /// `ln(sum(cns))`.
    ln_ploidy: f64,

    /// Window size.
    window: u32,
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

        let mut dedup_ids = Vec::with_capacity(n);
        let mut cns = Vec::with_capacity(n);
        for &contig_id in contig_ids.iter() {
            if let Some(i) = dedup_ids.iter().position(|&id| id == contig_id) {
                cns[i] += 1;
                continue;
            }
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
            dedup_ids.push(contig_id);
            cns.push(1);
        }
        Self {
            ids: dedup_ids,
            ln_ploidy: (n as f64).ln(),
            contigs, cns, window, bounded_starts, bounded_ends, wshifts,
        }
    }

    pub fn contig_names(&self) -> &ContigNames {
        &self.contigs
    }

    /// Total number of contigs. Can be less than ploidy, as contigs can be repeated.
    pub fn n_contigs(&self) -> usize {
        self.ids.len()
    }

    /// Total number of windows in the contig group.
    pub fn n_windows(&self) -> u32 {
        self.wshifts[self.wshifts.len() - 1]
    }

    pub fn window_size(&self) -> u32 {
        self.window
    }

    pub fn ids(&self) -> impl Iterator<Item = ContigId> + '_ {
        self.ids.iter().cloned()
    }

    /// Returns string with all contig names through a comma.
    pub fn ids_str(&self) -> String {
        let mut s = String::new();
        let mut first = true;
        for (id, cn) in self.contigs_cns() {
            for _ in 0..cn {
                if first {
                    first = false;
                } else {
                    write!(s, ",").unwrap();
                }
                write!(s, "{}", self.contigs.get_name(id)).unwrap();
            }
        }
        s
    }

    /// Returns iterator over pairs `(contig_id, contig_cn)`.
    pub fn contigs_cns(&self) -> impl Iterator<Item = (ContigId, u8)> + '_ {
        self.ids().zip(self.cns.iter().cloned())
    }

    /// Returns the number of windows corresponding to `i`-th contig.
    pub(crate) fn get_n_windows(&self, i: usize) -> u32 {
        self.wshifts[i + 1] - self.wshifts[i]
    }

    /// Returns window shift for the `i`-th contig.
    pub(crate) fn get_wshift(&self, i: usize) -> u32 {
        self.wshifts[i]
    }

    /// Given all pair-alignments for a single read pair,
    /// finds all alignments corresponding to `self`, and appends them to `out_alns`.
    /// Returns the number of added alignments.
    ///
    /// Output pair-alignments with ln-probability worse than `best_prob - prob_diff` are discarded.
    /// Remaining alignments have random order, probabilities are not normalized.
    ///
    /// Any alignments that were in the vector before, stay as they are and in the same order.
    pub fn read_windows(&self,
        pair_alns: &ReadPairAlignments,
        out_alns: &mut Vec<ReadWindows>,
        prob_diff: f64,
    ) -> usize {
        let start_len = out_alns.len();
        // Probability of being unmapped to any of the contigs.
        let mut unmapped_prob = pair_alns.unmapped_prob() + self.ln_ploidy;
        // Current threshold, is updated during the for-loop.
        let mut thresh_prob = unmapped_prob - prob_diff;
        for (i, &contig_id) in self.ids.iter().enumerate() {
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
                            None
                        } else {
                            Some((wshift + (middle1 - start) / self.window, wshift + (middle2 - start) / self.window))
                        }
                    },
                    TwoIntervals::First(interval1) => {
                        let middle1 = interval1.middle();
                        if middle1 < start || middle1 >= end {
                            None
                        } else {
                            Some((wshift + (middle1 - start) / self.window, UNMAPPED_WINDOW))
                        }
                    },
                    TwoIntervals::Second(interval2) => {
                        let middle2 = interval2.middle();
                        if middle2 < start || middle2 >= end {
                            None
                        } else {
                            Some((UNMAPPED_WINDOW, wshift + (middle2 - start) / self.window))
                        }
                    },
                };

                let ln_prob = paln.ln_prob();
                if let Some(windows) = windows {
                    if ln_prob >= thresh_prob {
                        thresh_prob = thresh_prob.max(ln_prob - prob_diff);
                        out_alns.push(ReadWindows::new(aln_ix as u32, windows, ln_prob));
                    }
                } else {
                    unmapped_prob = Ln::add(unmapped_prob, paln.ln_prob());
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
