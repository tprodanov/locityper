use std::{
    cmp::min,
    io::Write,
};
use rand::Rng;
use crate::{
    seq::{
        self, Interval,
        contigs::{ContigId, ContigNames, ContigSet, Genotype},
        kmers::KmerCounts,
    },
    math::distr::DistrBox,
    bg::ReadDepth,
    err::{Error, add_path},
};
use super::{
    locs::{GrouppedAlignments, PairAlignment},
    dp_cache::CachedDepthDistrs,
};

/// Structure that, for predefined region boundaries and window size,
/// returns windows for each alignment boundaries.
/// This function is used both in read depth calculation and during the locus genotyping.
#[derive(Clone)]
pub struct WindowGetter {
    start: u32,
    end: u32,
    window: u32,
    // Half window: size that needs to be covered in `get_significant` - 1.
    halfw: u32,
}

impl WindowGetter {
    pub fn new(start: u32, end: u32, window: u32) -> Self {
        Self {
            start, end, window,
            halfw: window.checked_sub(1).unwrap() / 2,
        }
    }

    /// Returns range of windows (start_ix..end_ix),
    /// covered by the alignment (at least by 1 bp).
    pub fn covered_any(&self, aln_start: u32, aln_end: u32) -> (u32, u32) {
        (
            min(aln_start, self.end).saturating_sub(self.start) / self.window,
            (min(aln_end, self.end) + self.window - 1).saturating_sub(self.start) / self.window,
        )
    }

    /// Returns range of windows (start_ix..end_ix),
    /// middle of which are covered by the alignment.
    pub fn covered_middle(&self, aln_start: u32, aln_end: u32) -> (u32, u32) {
        (
            (min(aln_start, self.end) + self.halfw).saturating_sub(self.start) / self.window,
            (min(aln_end, self.end) + self.halfw).saturating_sub(self.start) / self.window,
        )
    }

    /// Returns window that is covered by the alignment middle.
    pub fn middle_window(&self, aln_middle: u32) -> Option<u32> {
        if self.start <= aln_middle && aln_middle < self.end {
            Some((aln_middle - self.start) / self.window)
        } else {
            None
        }
    }
}

/// First window represents unmapped reads (or windows on the boundary).
pub(crate) const UNMAPPED_WINDOW: u32 = 0;
/// Second window represents reads alignments within the boundary.
pub(crate) const BOUNDARY_WINDOW: u32 = 1;
/// Regular windows have indices starting with 2.
pub(crate) const REG_WINDOW_SHIFT: u32 = 2;

/// Alignment of a read pair to a genotype.
pub struct ReadGtAlns {
    /// Index in the list of read-pair alignments.
    aln_ix: u32,
    /// ln-probability of this location.
    ln_prob: f64,
    /// Index of the contig (within the genotype).
    contig_ix: Option<u8>,
    /// Middle of the first read alignment.
    middle1: Option<u32>,
    /// Middle of the second read alignment.
    middle2: Option<u32>,
    /// Windows, to which the read is aligned to.
    /// Can change based on the random tweaking.
    windows: [u32; 2],
}

impl ReadGtAlns {
    fn new(aln_ix: u32, contig_ix: u8, paln: &PairAlignment) -> Self {
        Self {
            aln_ix,
            contig_ix: Some(contig_ix),
            ln_prob: paln.ln_prob(),
            middle1: paln.intervals().first().map(Interval::middle),
            middle2: paln.intervals().second().map(Interval::middle),
            windows: [UNMAPPED_WINDOW; 2],
        }
    }

    fn both_unmapped(ln_prob: f64) -> Self {
        Self {
            ln_prob,
            aln_ix: u32::MAX,
            contig_ix: None,
            middle1: None,
            middle2: None,
            windows: [UNMAPPED_WINDOW; 2],
        }
    }

    pub fn define_windows_determ(&mut self, mcontigs: &GenotypeWindows) {
        if let Some(i) = self.contig_ix {
            let contig = &mcontigs.by_contig[i as usize];
            let shift = mcontigs.wshifts[i as usize];
            self.windows = [
                contig.get_shifted_window(shift, self.middle1),
                contig.get_shifted_window(shift, self.middle2)
            ];
        }
    }

    pub fn define_windows_random(&mut self, mcontigs: &GenotypeWindows, tweak: u32, rng: &mut impl Rng) {
        if let Some(i) = self.contig_ix {
            let contig = &mcontigs.by_contig[i as usize];
            let shift = mcontigs.wshifts[i as usize];
            let r = rng.next_u64();
            let tweak1 = (r >> 32) as u32 % (2 * tweak + 1);
            let tweak2 = r as u32 % (2 * tweak + 1);

            self.windows = [
                contig.get_shifted_window(shift, self.middle1.map(|middle| middle + tweak1)),
                contig.get_shifted_window(shift, self.middle2.map(|middle| middle + tweak2)),
            ];
        }
    }

    /// Returns range of windows, to which the first and the second read ends are aligned.
    pub fn windows(&self) -> [u32; 2] {
        self.windows
    }

    /// Returns ln-probability of the alignment.
    pub fn ln_prob(&self) -> f64 {
        self.ln_prob
    }

    /// Index of the read-pair alignment across all alignments for the read pair.
    pub(crate) fn aln_ix(&self) -> u32 {
        self.aln_ix
    }
}

/// Calculate weight function [0, 1] -> [0, 1].
struct WeightCalculator {
    coeff: f64,
    power: f64,
}

impl WeightCalculator {
    /// Creates new calculator based on two parameters:
    /// - breakpoint in (0, 1): weight(breakpoint) = 1/2.
    /// - power > 0: determines how fast do the value change near the breakpoint.
    ///   The bigger the power, the faster the change.
    fn new(breakpoint: f64, power: f64) -> Self {
        assert!(0.0 < breakpoint && breakpoint < 1.0 && power > 0.0);
        Self {
            coeff: (1.0 - breakpoint).powf(power) / breakpoint.powf(power),
            power,
        }
    }

    /// Calculates weight according to the function $t^2 / (t^2 + 1)$
    /// where $t = C x^p / (1-x)^p$.
    fn get(&self, x: f64) -> f64 {
        if x == 1.0 {
            1.0
        } else {
            let t = self.coeff * x.powf(self.power) / (1.0 - x).powf(self.power);
            let t2 = t * t;
            t2 / (t2 + 1.0)
        }
    }
}

/// Set of windows for one contig.
#[derive(Clone)]
pub struct ContigWindows {
    /// Structure, that stores start, end and window.
    window_getter: WindowGetter,
    /// GC-content for each window within contig.
    window_gcs: Vec<u8>,
    /// k-mer-based window weights.
    window_weights: Vec<f64>,
    /// Is the window used (based on the weight?).
    use_window: Vec<bool>,
    /// Read depth distributions for each window.
    depth_distrs: Option<Vec<DistrBox>>,
}

impl ContigWindows {
    fn new(
        contig_id: ContigId,
        contigs: &ContigNames,
        seq: &[u8],
        kmer_counts: &KmerCounts,
        depth: &ReadDepth,
        weight_calc: &WeightCalculator,
        params: &super::Params,
        dbg_writer: &mut impl Write,
    ) -> Result<Self, Error>
    {
        let contig_len = seq.len() as u32;
        let contig_name = contigs.get_name(contig_id);
        let window = depth.window_size();
        let neighb_size = depth.neighb_size();
        if contig_len < window + 2 * params.boundary_size {
            return Err(Error::RuntimeError(format!("Contig {} is too short (len = {})",
                contigs.get_name(contig_id), contig_len)));
        }
        debug_assert_eq!(contig_len, contigs.get_len(contig_id));
        let n_windows = (contig_len - 2 * params.boundary_size) / window;
        let sum_len = n_windows * window;
        let reg_start = (contig_len - sum_len) / 2;

        let k = kmer_counts.k();
        let halfk = k / 2;
        let contig_kmer_counts = kmer_counts.get(contig_id);

        let mut window_gcs = Vec::with_capacity(n_windows as usize);
        let mut window_weights = Vec::with_capacity(n_windows as usize);
        let left_padding = (neighb_size - window) / 2;
        let right_padding = neighb_size - window - left_padding;
        for i in 0..n_windows {
            let start = reg_start + i * window;
            let end = start + window;
            let padded_start = start.saturating_sub(left_padding);
            let padded_end = (end + right_padding).min(contig_len);
            let gc = seq::gc_content(&seq[padded_start as usize..padded_end as usize]).round() as u8;
            window_gcs.push(gc);

            // What is the biggest quantile that still produces 1?
            let kmer_start = padded_start.saturating_sub(halfk) as usize;
            let kmer_end = min(padded_end - halfk, contig_len - k + 1) as usize;
            let inv_cdf1 = contig_kmer_counts[kmer_start..kmer_end].iter().filter(|&&x| x <= 1).count() as f64
                / (kmer_end - kmer_start) as f64;
            let weight = weight_calc.get(inv_cdf1);
            window_weights.push(weight);
            writeln!(dbg_writer, "{}\t{}\t{}\t{}\t{:.3}\t{:.5}", contig_name, start, end, gc, inv_cdf1, weight)
                .map_err(add_path!(!))?;
        }
        Ok(Self {
            window_getter: WindowGetter::new(reg_start, reg_start + sum_len, window),
            depth_distrs: None,
            use_window: window_weights.iter().map(|&weight| weight >= params.min_weight).collect(),
            window_gcs, window_weights,
        })
    }

    /// Creates a set of contig windows for each contig.
    pub fn new_all(
        set: &ContigSet,
        depth: &ReadDepth,
        params: &super::Params,
        mut dbg_writer: impl Write,
    ) -> Result<Vec<Self>, Error> {
        let contigs = set.contigs();
        let seqs = set.seqs();
        let weight_calc = WeightCalculator::new(params.weight_breakpoint, params.weight_power);
        writeln!(dbg_writer, "#contig\tstart\tend\tGC\tfrac_unique\tweight").map_err(add_path!(!))?;
        contigs.ids().zip(seqs)
            .map(|(id, seq)| Self::new(id, contigs, seq, set.kmer_counts(), depth,
                &weight_calc, params, &mut dbg_writer))
            .collect()
    }

    /// Defines read depth distributions for all contig windows.
    pub fn define_all_distributions(
        contig_windows: &mut [Self],
        cached_distrs: &CachedDepthDistrs<'_>,
    ) {
        for curr_windows in contig_windows {
            curr_windows.depth_distrs = Some(
                curr_windows.window_gcs.iter().zip(&curr_windows.window_weights)
                    .map(|(&gc, &w)| cached_distrs.get_distribution(gc, w))
                    .collect()
            );
        }
    }

    pub fn n_windows(&self) -> u32 {
        self.window_gcs.len() as u32
    }

    pub fn window_size(&self) -> u32 {
        self.window_getter.window
    }

    pub fn depth_distrs(&self) -> &[DistrBox] {
        &self.depth_distrs.as_ref().expect("Read depth distributions have not been initialized")
    }

    /// Returns all window weights.
    pub fn weights(&self) -> &[f64] {
        &self.window_weights
    }

    /// Is the window used?
    pub fn use_window(&self) -> &[bool] {
        &self.use_window
    }

    /// Returns window weight based on the read middle.
    /// If read is out of bounds, returns 1.0.
    pub fn get_window_weight(&self, middle: u32) -> f64 {
        self.window_getter.middle_window(middle).map(|w| self.window_weights[w as usize]).unwrap_or(1.0)
    }

    /// Returns window range within the contig based on the read alignment range.
    fn get_shifted_window(&self, shift: u32, middle: Option<u32>) -> u32 {
        match middle {
            Some(middle) => self.window_getter.middle_window(middle).map(|w| w + shift).unwrap_or(BOUNDARY_WINDOW),
            None => UNMAPPED_WINDOW,
        }
    }
}

/// Stores the contigs and windows corresponding to the windows.
pub struct GenotypeWindows {
    genotype: Genotype,
    by_contig: Vec<ContigWindows>,
    /// Start index for each contig id (length: n-contigs + 1).
    /// Contig with index `i` will have windows between `wshifts[i]..wshifts[i + 1]`.
    wshifts: Vec<u32>,
    /// Window size.
    window: u32,
}

impl GenotypeWindows {
    pub fn new(genotype: Genotype, contig_windows: &[ContigWindows]) -> Self {
        let n = genotype.ploidy();
        let mut by_contig = Vec::<ContigWindows>::with_capacity(n);
        let mut wshifts = Vec::with_capacity(n + 1);
        let mut curr_wshift = REG_WINDOW_SHIFT;
        wshifts.push(curr_wshift);

        for &id in genotype.ids() {
            let curr_contig = contig_windows[id.ix()].clone();
            curr_wshift += curr_contig.n_windows();
            wshifts.push(curr_wshift);
            by_contig.push(curr_contig);
        }
        assert!(by_contig.len() < 256, "Multi-contig collection cannot contain more than 256 entries");
        Self {
            window: by_contig[0].window_size(),
            genotype, by_contig, wshifts,
        }
    }

    pub fn genotype(&self) -> &Genotype {
        &self.genotype
    }

    pub fn window_size(&self) -> u32 {
        self.window
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
    /// Remaining alignments are sorted from best likelihood to worst.
    ///
    /// Any alignments that were in the vector before, stay as they are and in the same order.
    pub fn extend_read_gt_alns(&self,
        groupped_alns: &GrouppedAlignments,
        out_alns: &mut Vec<ReadGtAlns>,
        prob_diff: f64,
    ) -> usize {
        let start_len = out_alns.len();
        // Probability of being unmapped to any of the contigs.
        let unmapped_prob = groupped_alns.unmapped_prob();
        // Current threshold, is updated during the for-loop.
        let mut thresh_prob = unmapped_prob - prob_diff;
        for (i, &contig_id) in self.genotype.ids().iter().enumerate() {
            let contig_ix = u8::try_from(i).unwrap();
            let (start_ix, alns) = groupped_alns.contig_alns(contig_id);
            if !alns.is_empty() {
                // First alignment should have highest probability.
                thresh_prob = thresh_prob.max(alns[0].ln_prob() - prob_diff);
                for (aln_ix, aln) in (start_ix..).zip(alns) {
                    if aln.ln_prob() >= thresh_prob {
                        out_alns.push(ReadGtAlns::new(aln_ix as u32, contig_ix, aln));
                    } else {
                        break;
                    }
                }
            }
        }

        if unmapped_prob >= thresh_prob {
            out_alns.push(ReadGtAlns::both_unmapped(unmapped_prob));
        }
        let slice = &mut out_alns[start_len..];
        // Decreasing sort by ln-probability.
        slice.sort_unstable_by(|a, b| b.ln_prob.total_cmp(&a.ln_prob));
        let keep_alns = slice.partition_point(|aln| aln.ln_prob >= thresh_prob);
        out_alns.truncate(start_len + keep_alns);
        keep_alns
    }

    pub fn total_windows(&self) -> u32 {
        *self.wshifts.last().unwrap()
    }

    pub fn contig_windows(&self) -> &[ContigWindows] {
        &self.by_contig
    }
}
