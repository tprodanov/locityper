use std::{
    cmp::min,
    fmt::Write,
};
use crate::{
    math::Ln,
    ext::vec::F64Ext,
    bg::depth::ReadDepth,
    seq::{
        self, ContigId, ContigNames, ContigSet,
        kmers::KmerCounts,
    },
};
use super::{
    locs::{TwoIntervals, ReadPairAlignments},
    dp_cache::{CachedDepthDistrs, DistrBox},
};

/// First window represents unmapped reads.
pub(crate) const UNMAPPED_WINDOW: u32 = 0;
/// Second window represents reads alignments within the boundary.
pub(crate) const BOUNDARY_WINDOW: u32 = 1;
/// Regular windows have indices starting with 2.
pub(crate) const REG_WINDOW_SHIFT: u32 = 2;

/// Alignment of a read pair to a specific contig windows.
pub struct ReadWindows {
    /// Index in the list of read-pair alignments.
    ix: u32,
    /// Window for each read-end.
    windows: (u32, u32),
    /// ln-probability of this location.
    ln_prob: f64,
}

impl ReadWindows {
    fn new(ix: u32, windows: (u32, u32), ln_prob: f64) -> Self {
        Self { ix, windows, ln_prob }
    }

    /// Returns two windows, to which the read pair is aligned.
    pub fn windows(&self) -> (u32, u32) {
        self.windows
    }

    /// Returns ln-probability of the alignment.
    pub fn ln_prob(&self) -> f64 {
        self.ln_prob
    }

    /// Index of the read-pair alignment across all alignments for the read pair.
    pub(crate) fn ix(&self) -> u32 {
        self.ix
    }
}

/// Calculates weight of a positive value x.
/// If x <= c1, returns 1.
/// If x = c2, returns 0.5 (returns ln 0.5).
/// For values between c1 and infinity, the weight is distributed according to a sech function.
/// ```
/// f(x) = { 1                              if x in [0, c1].
///         /      x - c1
///        { sech ------- * ln(2 + sqrt3)   if x > c1.
///         \     c2 - c1
/// ```
fn sech_weight(x: f64, c1: f64, c2: f64) -> f64 {
    if x <= c1 {
        return 1.0;
    }
    // ln(2 + sqrt(3)).
    const LN2_SQRT3: f64 = 1.31695789692481;
    let t = (x - c1) / (c2 - c1) * LN2_SQRT3;
    let expt = t.exp();
    2.0 * expt / (expt * expt + 1.0)
}

/// Set of windows for one contig.
#[derive(Clone)]
pub struct ContigWindows {
    contig_id: ContigId,
    window: u32,
    /// Start of the windows within contig (ignoring boundary).
    start: u32,
    /// End of the windows within contig (ignoring boundary). `end - start` divides window size.
    end: u32,
    /// GC-content for each window within contig.
    window_gcs: Vec<u8>,
    /// k-mer-based window weights.
    window_weights: Vec<f64>,
}

impl ContigWindows {
    pub fn new(
        contig_id: ContigId,
        contigs: &ContigNames,
        seq: &[u8],
        kmer_counts: &KmerCounts,
        depth: &ReadDepth,
        params: &super::Params,
    ) -> Self
    {
        let contig_len = seq.len() as u32;
        let window = depth.window_size();
        let window_padding = depth.window_padding();
        assert!(contig_len > window + 2 * params.boundary_size,
            "Contig {} is too short (len = {})", contigs.get_name(contig_id), contig_len);
        debug_assert_eq!(contig_len, contigs.get_len(contig_id));
        let n_windows = (f64::from(contig_len - 2 * params.boundary_size) / f64::from(window))
            .floor() as u32;
        let sum_len = n_windows * window;
        let start = (contig_len - sum_len) / 2;
        let end = start + sum_len;
        debug_assert!(start >= window_padding && end + window_padding <= contig_len);

        let k = kmer_counts.k();
        let halfk = k / 2;
        let contig_kmer_counts = kmer_counts.get(contig_id);

        let mut window_gcs = Vec::with_capacity(n_windows as usize);
        let mut window_weights = Vec::with_capacity(n_windows as usize);
        for j in 0..n_windows {
            let padded_start = start + window * j - window_padding;
            let padded_end = start + window * (j + 1) + window_padding;
            window_gcs.push(seq::gc_content(&seq[padded_start as usize..padded_end as usize]).round() as u8);

            let mean_kmer_freq = F64Ext::mean(&contig_kmer_counts[padded_start.saturating_sub(halfk) as usize
                ..min(padded_end - halfk, contig_len - k + 1) as usize]);
            window_weights.push(sech_weight(mean_kmer_freq, params.rare_kmer, params.semicommon_kmer));
        }
        Self { contig_id, window, start, end, window_gcs, window_weights }
    }

    /// Creates a set of contig windows for each contig.
    pub fn new_all(set: &ContigSet, depth: &ReadDepth, params: &super::Params) -> Vec<Self> {
        let contigs = set.contigs();
        contigs.ids().zip(set.seqs())
            .map(|(id, seq)| Self::new(id, contigs, seq, set.kmer_counts(), depth, params))
            .collect()
    }

    pub fn n_windows(&self) -> u32 {
        self.window_gcs.len() as u32
    }
}

/// Stores the contigs and windows corresponding to the windows.
pub struct MultiContigWindows {
    by_contig: Vec<ContigWindows>,
    cns: Vec<u8>,
    /// `ln(sum(cns))`.
    ln_ploidy: f64,
    /// Start index for each contig id (length: n-contigs + 1).
    /// Contig with index `i` will have windows between `wshifts[i]..wshifts[i + 1]`.
    wshifts: Vec<u32>,
    window: u32,
}

impl MultiContigWindows {
    pub fn new(contig_ids: &[ContigId], contig_windows: &[ContigWindows]) -> Self {
        let n = contig_ids.len();
        let mut by_contig = Vec::<ContigWindows>::with_capacity(n);
        let mut cns = Vec::<u8>::with_capacity(n);
        let mut wshifts = Vec::with_capacity(n + 1);
        let mut curr_wshift = REG_WINDOW_SHIFT;
        wshifts.push(curr_wshift);

        for &id in contig_ids.iter() {
            if let Some(i) = by_contig.iter().position(|el| el.contig_id == id) {
                cns[i] = cns[i].checked_add(1).unwrap();
                continue;
            }
            let curr_contig = contig_windows[id.ix()].clone();
            curr_wshift += curr_contig.n_windows();
            wshifts.push(curr_wshift);
            by_contig.push(curr_contig);
            cns.push(1);
        }
        Self {
            ln_ploidy: (n as f64).ln(),
            window: by_contig[0].window,
            by_contig, cns, wshifts,
        }
    }

    /// Total number of contigs. Can be less than ploidy, as contigs can be repeated.
    pub fn n_contigs(&self) -> usize {
        self.by_contig.len()
    }

    /// Total number of windows in the contig group.
    pub fn total_windows(&self) -> u32 {
        *self.wshifts.last().unwrap()
    }

    pub fn window_size(&self) -> u32 {
        self.window
    }

    pub fn ids(&self) -> impl Iterator<Item = ContigId> + '_ {
        self.by_contig.iter().map(|el| el.contig_id)
    }

    /// Returns string with all contig names through a comma.
    pub fn ids_str(&self, contigs: &ContigNames) -> String {
        let mut s = String::new();
        for (id, cn) in self.contigs_cns() {
            for _ in 0..cn {
                if !s.is_empty() {
                    write!(s, ",").unwrap();
                }
                write!(s, "{}", contigs.get_name(id)).unwrap();
            }
        }
        s
    }

    /// Returns iterator over pairs `(contig_id, contig_cn)`.
    pub fn contigs_cns(&self) -> impl Iterator<Item = (ContigId, u8)> + '_ {
        self.ids().zip(self.cns.iter().copied())
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
        for (i, contig_id) in self.ids().enumerate() {
            let wshift = self.wshifts[i];
            let start = self.by_contig[i].start;
            let end = self.by_contig[i].end;

            let (start_ix, palns) = pair_alns.contig_alns(contig_id);
            for (aln_ix, paln) in (start_ix..).zip(palns) {
                let w1 = match paln.intervals() {
                    TwoIntervals::Both(interval1, _) | TwoIntervals::First(interval1) => {
                        let middle1 = interval1.middle();
                        if start <= middle1 && middle1 < end {
                            wshift + (middle1 - start) / self.window
                        } else {
                            BOUNDARY_WINDOW
                        }
                    }
                    TwoIntervals::Second(_) => UNMAPPED_WINDOW,
                };
                let w2 = match paln.intervals() {
                    TwoIntervals::Both(_, interval2) | TwoIntervals::Second(interval2) => {
                        let middle2 = interval2.middle();
                        if start <= middle2 && middle2 < end {
                            wshift + (middle2 - start) / self.window
                        } else {
                            BOUNDARY_WINDOW
                        }
                    }
                    TwoIntervals::First(_) => UNMAPPED_WINDOW,
                };

                let ln_prob = paln.ln_prob();
                if w1 >= REG_WINDOW_SHIFT || w2 >= REG_WINDOW_SHIFT {
                    if ln_prob >= thresh_prob {
                        thresh_prob = thresh_prob.max(ln_prob - prob_diff);
                        out_alns.push(ReadWindows::new(aln_ix as u32, (w1, w2), ln_prob));
                    }
                } else {
                    // If the read is out of bounds, we will reduce the penalty to be unmapped.
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

    pub(crate) fn get_distributions(&self, cached_distrs: &CachedDepthDistrs) -> Vec<DistrBox> {
        const _: () = assert!(UNMAPPED_WINDOW == 0 && BOUNDARY_WINDOW == 1 && REG_WINDOW_SHIFT == 2,
            "Constants were changed!");
        let mut distrs: Vec<DistrBox> = Vec::with_capacity(self.total_windows() as usize);
        distrs.push(Box::new(cached_distrs.unmapped_distr()));
        distrs.push(Box::new(cached_distrs.boundary_distr()));

        for (curr_contig, &contig_cn) in self.by_contig.iter().zip(&self.cns) {
            for (&gc_content, &weight) in curr_contig.window_gcs.iter().zip(&curr_contig.window_weights) {
                distrs.push(cached_distrs.get_distribution(gc_content, contig_cn, weight));
            }
        }
        distrs
    }

    /// Returns all window weights for i-th contig.
    pub(crate) fn get_weights(&self, i: usize) -> &[f64] {
        &self.by_contig[i].window_weights
    }

    /// Returns all GC-contents for i-th contig.
    pub(crate) fn get_gc_contents(&self, i: usize) -> &[u8] {
        &self.by_contig[i].window_gcs
    }
}
