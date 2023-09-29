use std::{
    cmp::min,
    io::Write,
    sync::Arc,
};
use rand::Rng;
use crate::{
    seq::{
        self,
        contigs::{ContigId, ContigNames, ContigSet, Genotype},
        kmers::KmerCounts,
    },
    bg::ReadDepth,
    err::{Error, add_path, validate_param},
    ext::vec::IterExt,
};
use super::{
    locs::{GrouppedAlignments, PairAlignment},
};

/// Structure that, for predefined region boundaries and window size,
/// returns windows for each alignment boundaries.
/// This function is used both in read depth calculation and during the locus genotyping.
#[derive(Clone)]
pub struct WindowGetter {
    start: u32,
    end: u32,
    window: u32,
    // /// Half window: size that needs to be covered in `get_significant` - 1.
    // halfw: u32,
}

impl WindowGetter {
    pub fn new(start: u32, end: u32, window: u32) -> Self {
        Self {
            start, end, window,
            // halfw: window.checked_sub(1).unwrap() / 2,
        }
    }

    /// Returns boundaries of the i-th window.
    pub fn ith_window(&self, i: u32) -> (u32, u32) {
        let start = self.start + i * self.window;
        (start, start + self.window)
    }

    // /// Returns range of windows (start_ix..end_ix),
    // /// covered by the alignment (at least by 1 bp).
    // pub fn covered_any(&self, aln_start: u32, aln_end: u32) -> (u32, u32) {
    //     (
    //         min(aln_start, self.end).saturating_sub(self.start) / self.window,
    //         (min(aln_end, self.end) + self.window - 1).saturating_sub(self.start) / self.window,
    //     )
    // }

    // /// Returns range of windows (start_ix..end_ix),
    // /// middle of which are covered by the alignment.
    // pub fn covered_middle(&self, aln_start: u32, aln_end: u32) -> (u32, u32) {
    //     (
    //         (min(aln_start, self.end) + self.halfw).saturating_sub(self.start) / self.window,
    //         (min(aln_end, self.end) + self.halfw).saturating_sub(self.start) / self.window,
    //     )
    // }

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
    /// Contig index within the genome and the parent pair alignment,
    /// Right: ln-probability, if both mates are considered unmapped.
    parent: Option<(u8, PairAlignment)>,
    /// Ln-probability (copied from the parent, if it is defined).
    ln_prob: f64,
    /// Windows, to which the read is aligned to.
    /// Can change based on the random tweaking.
    windows: [u32; 2],
}

impl ReadGtAlns {
    fn new(contig_ix: u8, paln: PairAlignment) -> Self {
        Self {
            ln_prob: paln.ln_prob(),
            parent: Some((contig_ix, paln)),
            windows: [UNMAPPED_WINDOW; 2],
        }
    }

    fn both_unmapped(ln_prob: f64) -> Self {
        Self {
            ln_prob,
            parent: None,
            windows: [UNMAPPED_WINDOW; 2],
        }
    }

    pub fn parent(&self) -> Option<&PairAlignment> {
        self.parent.as_ref().map(|(_, paln)| paln)
    }

    pub fn define_windows_determ(&mut self, mcontigs: &GenotypeWindows) {
        if let Some((contig_ix, paln)) = &self.parent {
            let contig = &mcontigs.infos[usize::from(*contig_ix)];
            let shift = mcontigs.wshifts[usize::from(*contig_ix)];
            self.windows = [
                contig.get_shifted_window_ix(shift, paln.middle1()),
                contig.get_shifted_window_ix(shift, paln.middle2())
            ];
        }
    }

    pub fn define_windows_random(&mut self, mcontigs: &GenotypeWindows, tweak: u32, rng: &mut impl Rng) {
        if let Some((contig_ix, paln)) = &self.parent {
            let contig = &mcontigs.infos[usize::from(*contig_ix)];
            let shift = mcontigs.wshifts[usize::from(*contig_ix)];
            let r = rng.next_u64();
            let tweak1 = (r >> 32) as u32 % (2 * tweak + 1);
            let tweak2 = r as u32 % (2 * tweak + 1);

            self.windows = [
                contig.get_shifted_window_ix(shift, paln.middle1().map(|middle| middle + tweak1)),
                contig.get_shifted_window_ix(shift, paln.middle2().map(|middle| middle + tweak2)),
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
}

/// Calculate weight function [0, 1] -> [0, 1].
#[derive(Clone, Debug)]
pub struct WeightCalculator {
    /// Value in (0, 1): weight(breakpoint) = 1/2.
    breakpoint: f64,
    /// Value > 0: determines how fast do the value change near the breakpoint.
    power: f64,
    /// Constant factor in the calculation.
    coeff: f64,
    /// power * 2.
    double_power: f64,
}

impl WeightCalculator {
    pub fn new(breakpoint: f64, power: f64) -> Result<Self, Error> {
        validate_param!(0.0 < breakpoint && breakpoint < 1.0,
            "Weight breakpoint ({}) must be within (0, 1)", breakpoint);
        validate_param!(power > 0.5 && power <= 50.0, "Weight power ({}) should be within [0.5, 50].", power);
        let double_power = power * 2.0;
        Ok(Self {
            coeff: ((1.0 - breakpoint) / breakpoint).powf(double_power),
            breakpoint, power, double_power,
        })
    }

    /// See supplementary methods.
    pub fn get(&self, x: f64) -> f64 {
        if x == 1.0 {
            1.0
        } else {
            let t2 = self.coeff * (x / (1.0 - x)).powf(self.double_power);
            t2 / (t2 + 1.0)
        }
    }

    /// Weight breakpoint.
    pub fn breakpoint(&self) -> f64 {
        self.breakpoint
    }

    /// Weight power.
    pub fn power(&self) -> f64 {
        self.power
    }
}

/// Set of windows for one contig.
#[derive(Clone)]
pub struct ContigInfo {
    /// Length of the contig.
    contig_len: u32,
    /// Window padding, used to find GC-content and fraction of unique k-mers.
    left_padding: u32,
    right_padding: u32,
    /// K-mer size.
    kmer_size: u32,

    /// Structure, that stores start, end and window.
    window_getter: WindowGetter,
    weight_calc: Option<WeightCalculator>,

    /// Cumulateve number of G/C across the contig. Size = len(contig) + 1.
    cumul_gc: Vec<u32>,
    /// Cumulative number of unique k-mers across the contig. Size = len(contig) + 1.
    cumul_uniq_kmers: Vec<u32>,
    /// Linguistic complexity for each moving window.
    complexities: Vec<f64>,

    /// Default window weights (without boundary tweaking).
    default_weights: Vec<f64>,
}

impl ContigInfo {
    fn new(
        contig_id: ContigId,
        contigs: &ContigNames,
        seq: &[u8],
        kmer_counts: &KmerCounts,
        depth: &ReadDepth,
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
        let reg_end = reg_start + sum_len;
        let left_padding = (neighb_size - window) / 2;
        let right_padding = neighb_size - window - left_padding;

        let cumul_uniq_kmers = IterExt::cumul_sums(kmer_counts.get(contig_id.ix()).iter().map(|&count| count <= 1));
        let cumul_gc = IterExt::cumul_sums(seq.iter().map(|&ch| ch == b'C' || ch == b'G'));
        let mut res = Self {
            contig_len, left_padding, right_padding, cumul_gc, cumul_uniq_kmers,
            kmer_size: kmer_counts.k(),
            window_getter: WindowGetter::new(reg_start, reg_end, window),
            weight_calc: params.weight_calc.clone(),
            default_weights: Vec::with_capacity(n_windows as usize),
            complexities: seq::compl::linguistic_complexity_123(seq, window as usize),
        };

        for i in 0..n_windows {
            let start = reg_start + i * window;
            let end = start + window;
            let (gc, uniq_frac, weight) = res.window_characteristics(start, end);
            res.default_weights.push(weight);
            writeln!(dbg_writer, "{}\t{}\t{}\t{}\t{:.3}\t{:.5}\t{:.5}",
                contig_name, start, end, gc, uniq_frac, res.complexities[start as usize], weight)
                .map_err(add_path!(!))?;
        }
        Ok(res)
    }

    /// Creates a set of contig windows for each contig.
    pub fn new_all(
        set: &ContigSet,
        depth: &ReadDepth,
        params: &super::Params,
        mut dbg_writer: impl Write,
    ) -> Result<Vec<Arc<ContigInfo>>, Error>
    {
        let contigs = set.contigs();
        let seqs = set.seqs();
        writeln!(dbg_writer, "#contig\tstart\tend\tGC\tfrac_unique\tcomplexity\tweight").map_err(add_path!(!))?;
        let mut all_contig_infos = Vec::with_capacity(seqs.len());
        for (id, seq) in contigs.ids().zip(seqs) {
            all_contig_infos.push(Arc::new(
                Self::new(id, contigs, seq, set.kmer_counts(), depth, params, &mut dbg_writer)?));
        }
        Ok(all_contig_infos)
    }

    /// Returns GC-content, kmer fraction and window weight for `start..end`.
    pub fn window_characteristics(&self, start: u32, end: u32) -> (u8, f64, f64) {
        let padded_start = start.saturating_sub(self.left_padding);
        let padded_end = min(self.contig_len, end + self.right_padding);
        let gc = 100.0 * f64::from(self.cumul_gc[padded_end as usize] - self.cumul_gc[padded_start as usize])
            / f64::from(padded_end - padded_start);

        let kmer_start = padded_start.saturating_sub(self.kmer_size / 2);
        let kmer_end = min(self.contig_len - self.kmer_size + 1, padded_end - self.kmer_size / 2);
        let uniq_frac = f64::from(self.cumul_uniq_kmers[kmer_end as usize] - self.cumul_uniq_kmers[kmer_start as usize])
            / f64::from(kmer_end - kmer_start);
        let weight = self.weight_calc.as_ref().map(|calc| calc.get(uniq_frac)).unwrap_or(1.0);
        (gc.round() as u8, uniq_frac, weight)
    }

    pub fn n_windows(&self) -> u32 {
        self.default_weights.len() as u32
    }

    pub fn window_size(&self) -> u32 {
        self.window_getter.window
    }

    /// Returns window weight based on the read middle.
    /// If read is out of bounds, returns 1.0.
    pub fn default_window_weight(&self, middle: u32) -> f64 {
        self.window_getter.middle_window(middle).map(|w| self.default_weights[w as usize]).unwrap_or(1.0)
    }

    /// Default window weights (can change due to random tweaking).
    pub fn default_weights(&self) -> &[f64] {
        &self.default_weights
    }

    /// Returns window range within the contig based on the read alignment range.
    fn get_shifted_window_ix(&self, shift: u32, middle: Option<u32>) -> u32 {
        match middle {
            Some(middle) => self.window_getter.middle_window(middle).map(|w| w + shift).unwrap_or(BOUNDARY_WINDOW),
            None => UNMAPPED_WINDOW,
        }
    }

    /// Generates window boundaries with given tweak size.
    pub fn generate_windows<'a>(&'a self, tweak: u32, rng: &'a mut impl Rng) -> impl Iterator<Item = (u32, u32)> + 'a {
        (0..self.n_windows()).map(move |i| {
            let (start, end) = self.window_getter.ith_window(i);
            let left_tweak = tweak.min(start) as i32;
            let right_tweak = tweak.min(self.contig_len.checked_sub(end).unwrap()) as i32;
            let r = rng.gen_range(-left_tweak..=right_tweak);
            (start.checked_add_signed(r).unwrap(), end.checked_add_signed(r).unwrap())
        })
    }
}

/// Stores the contigs and windows corresponding to the windows.
pub struct GenotypeWindows {
    genotype: Genotype,
    infos: Vec<Arc<ContigInfo>>,
    /// Start index for each contig id (length: n-contigs + 1).
    /// Contig with index `i` will have windows between `wshifts[i]..wshifts[i + 1]`.
    wshifts: Vec<u32>,
    /// Window size.
    window: u32,
}

impl GenotypeWindows {
    pub fn new(genotype: Genotype, all_infos: &[Arc<ContigInfo>]) -> Self {
        let n = genotype.ploidy();
        let mut infos = Vec::with_capacity(n);
        let mut wshifts = Vec::with_capacity(n + 1);
        let mut curr_wshift = REG_WINDOW_SHIFT;
        wshifts.push(curr_wshift);

        for &id in genotype.ids() {
            let curr_contig = Arc::clone(&all_infos[id.ix()]);
            curr_wshift += curr_contig.n_windows();
            wshifts.push(curr_wshift);
            infos.push(curr_contig);
        }
        assert!(infos.len() < 256, "Multi-contig collection cannot contain more than 256 entries");
        Self {
            window: infos[0].window_size(),
            genotype, infos, wshifts,
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
            let alns = groupped_alns.contig_aln_pairs(contig_id);
            if !alns.is_empty() {
                // First alignment should have highest probability.
                thresh_prob = thresh_prob.max(alns[0].ln_prob() - prob_diff);
                for aln in alns.iter() {
                    if aln.ln_prob() >= thresh_prob {
                        out_alns.push(ReadGtAlns::new(contig_ix, aln.clone()));
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

    pub fn contig_infos(&self) -> &[Arc<ContigInfo>] {
        &self.infos
    }
}
