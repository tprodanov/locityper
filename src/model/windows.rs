use std::{
    cmp::min,
    path::Path,
    io::{Write, BufRead},
    sync::Arc,
    ops::Index,
};
use rand::Rng;
use crate::{
    seq::{
        self,
        Interval,
        contigs::{ContigId, ContigNames, ContigSet, Genotype},
        counts::KmerCounts,
    },
    bg::ReadDepth,
    err::{error, add_path, validate_param},
    ext::{self, vec::IterExt},
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

    pub fn len(&self) -> u32 {
        self.end - self.start
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
    /// Contig index within the genotype and the parent pair alignment,
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

    /// Contig index within the
    pub fn parent(&self) -> &Option<(u8, PairAlignment)> {
        &self.parent
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
    pub fn new(breakpoint: f64, power: f64) -> crate::Result<Self> {
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

/// Explicitely set weights for one contig.
#[derive(Default, Clone)]
struct ExplicitWeights {
    /// Each weight is converted to integer. This is done for more accurate cumulative sums.
    /// Length of the vector = length of the contig.
    /// Each value = (weight, cumulative weight without the last one).
    weights: Vec<(u64, u64)>,
    sum: u64,
}

impl ExplicitWeights {
    // 2 ^ 32
    const SCALE: f64 = (1u64 << 32) as f64;

    fn push(&mut self, val: f64) {
        let i = (val * Self::SCALE) as u64;
        self.weights.push((i, self.sum));
        self.sum = self.sum.checked_add(i).unwrap();
    }

    fn extend_by(&mut self, n: usize, val: f64) {
        let i = (val * Self::SCALE) as u64;
        for _ in 0..n {
            self.weights.push((i, self.sum));
            self.sum = self.sum.checked_add(i).unwrap();
        }
    }

    /// Needs to be called at the end.
    fn finish(&mut self) {
        self.weights.push((0, self.sum));
    }

    #[inline]
    fn at(&self, i: usize) -> f64 {
        self.weights[i].0 as f64 / Self::SCALE
    }

    #[inline]
    fn sum(&self, i: usize, j: usize) -> f64 {
        (self.weights[j].1 - self.weights[i].1) as f64 / Self::SCALE
    }

    #[inline]
    fn mean(&self, i: usize, j: usize) -> f64 {
        ((self.weights[j].1 - self.weights[i].1) / (j - i) as u64) as f64 / Self::SCALE
    }

    #[inline(always)]
    fn is_empty(&self) -> bool {
        self.weights.is_empty()
    }

    /// After calling `finish` this number will be 1 bigger than the number of pushed elements.
    #[inline(always)]
    fn len(&self) -> usize {
        self.weights.len()
    }
}

/// Load explicit weights from a BED file (contig start end value).
/// Multipliers must be defined for all contigs.
///
/// For each haplotype, returns actual weights for each point, as well as moving averages.
fn load_explicit_weights(
    filename: &Path,
    contigs: &Arc<ContigNames>,
) -> crate::Result<Vec<Option<ExplicitWeights>>>
{
    let reader = ext::sys::open(&filename)?;
    log::debug!("    Loading explicit region weights from {}", ext::fmt::path(&filename));
    // Explicit weights for each basepair.
    let mut weights = vec![Some(ExplicitWeights::default()); contigs.len()];
    let mut ignored_lines = 0;
    for line in reader.lines() {
        let line = line.map_err(add_path!(filename))?;
        let mut split = line.split_whitespace();
        let interval = match Interval::parse_bed(&mut split, contigs) {
            Ok(interval) => interval,
            Err(crate::Error::ParsingError(s)) if s.starts_with("Unknown contig") => {
                ignored_lines += 1;
                continue;
            }
            Err(e) => return Err(e),
        };

        let val: f64 = split.next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| error!(ParsingError,
                "Failed to parse explicit weights ({}, line `{}`): fourth numeric column required",
                ext::fmt::path(&filename), line.replace('\t', " ")))?;
        if val < 0.0 || val > 1.0 {
            return Err(
                error!(ParsingError, "Failed to parse explicit weights ({}, line `{}`): value must be in [0, 1]",
                ext::fmt::path(&filename), line.replace('\t', " ")));
        }

        let curr_weights = weights[interval.contig_id().ix()].as_mut().unwrap();
        if curr_weights.len() != interval.start() as usize {
            return Err(error!(ParsingError,
                "Failed to parse explicit weights ({}): haplotype {} not fully covered",
                ext::fmt::path(&filename), interval.contig_name()));
        }
        curr_weights.extend_by(interval.len() as usize, val);
    }

    if ignored_lines > 0 {
        log::debug!("        Skipped {} lines in {} (unused haplotypes)",
            ignored_lines, ext::fmt::path(&filename));
    }

    for (curr_weights, (name, &length)) in weights.iter_mut().zip(contigs.names().iter().zip(contigs.lengths())) {
        let curr_weights = curr_weights.as_mut().unwrap();
        if curr_weights.is_empty() {
            return Err(error!(ParsingError,
                "Failed to parse explicit weights ({}): haplotype {} missing", ext::fmt::path(filename), name));
        } else if curr_weights.len() != length as usize {
            return Err(error!(ParsingError,
                "Failed to parse explicit weights ({}): haplotype {} not fully covered/has different length",
                ext::fmt::path(&filename), name));
        }
        curr_weights.finish();
    }
    Ok(weights)
}

/// Characteristics of a window.
pub struct WindowCharacteristics {
    pub gc_content: u8,
    pub uniq_kmer_frac: f64,
    pub ling_compl: f64,
    /// User-provided weight.
    pub explicit_weight: f64,
    /// Total weight.
    pub weight: f64,
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

    /// Window weight calculators.
    kmers_weight_calc: Option<WeightCalculator>,
    compl_weight_calc: Option<WeightCalculator>,

    /// Cumulateve number of G/C across the contig. Size = len(contig) + 1.
    cumul_gc: Vec<u32>,
    /// Cumulative number of unique k-mers across the contig. Size = len(contig) + 1.
    cumul_uniq_kmers: Vec<u32>,
    /// Linguistic complexity for each moving window.
    mov_complexities: Vec<f64>,
    /// Explicitely set weights.
    explicit_weights: Option<ExplicitWeights>,

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
        explicit_weights: Option<ExplicitWeights>,
        params: &super::Params,
        dbg_writer: &mut impl Write,
    ) -> crate::Result<Self>
    {
        let contig_len = seq.len() as u32;
        let contig_name = contigs.get_name(contig_id);
        let window = depth.window_size();
        let neighb_size = depth.neighb_size();
        if contig_len < window + 2 * params.boundary_size {
            return Err(error!(RuntimeError, "Contig {} is too short (len = {})",
                contigs.get_name(contig_id), contig_len));
        }
        debug_assert_eq!(contig_len, contigs.get_len(contig_id));
        let n_windows = (contig_len - 2 * params.boundary_size) / window;
        let sum_len = n_windows * window;
        let reg_start = (contig_len - sum_len) / 2;
        let reg_end = reg_start + sum_len;
        let left_padding = (neighb_size - window) / 2;
        let right_padding = neighb_size - window - left_padding;

        // Unique k-mers should have count = 0 because `kmer_counts` only stores off-target matches.
        let cumul_uniq_kmers = IterExt::cumul_sums(kmer_counts.get(contig_id.ix()).iter().map(|&count| count == 0));
        let cumul_gc = IterExt::cumul_sums(seq.iter().map(|&ch| ch == b'C' || ch == b'G'));
        let mut res = Self {
            contig_len, left_padding, right_padding, cumul_gc, cumul_uniq_kmers,
            kmer_size: kmer_counts.k(),
            window_getter: WindowGetter::new(reg_start, reg_end, window),
            kmers_weight_calc: params.kmers_weight_calc.clone(),
            compl_weight_calc: params.compl_weight_calc.clone(),
            default_weights: Vec::with_capacity(n_windows as usize),
            mov_complexities: seq::compl::linguistic_complexity_123(seq, neighb_size as usize),
            explicit_weights,
        };

        for i in 0..n_windows {
            let (start, end) = res.window_getter.ith_window(i);
            let chars = res.window_characteristics(start, end);
            res.default_weights.push(chars.weight);
            writeln!(dbg_writer, "{}\t{}\t{}\t{}\t{:.3}\t{:.5}\t{:.5}\t{:.5}",
                contig_name, start, end, chars.gc_content, chars.uniq_kmer_frac, chars.ling_compl,
                chars.explicit_weight, chars.weight).map_err(add_path!(!))?;
        }
        Ok(res)
    }

    /// Returns GC-content, kmer fraction, window complexity and window weight for `start..end`.
    pub fn window_characteristics(&self, start: u32, end: u32) -> WindowCharacteristics {
        let padded_start = start.saturating_sub(self.left_padding);
        let padded_end = min(self.contig_len, end + self.right_padding);
        let gc_content = (100.0 * f64::from(self.cumul_gc[padded_end as usize] - self.cumul_gc[padded_start as usize])
            / f64::from(padded_end - padded_start)).round() as u8;
        let ling_compl = self.mov_complexities[padded_start as usize];
        let explicit_weight = self.explicit_weights.as_ref()
            .map(|weights| weights.mean(start as usize, end as usize)).unwrap_or(1.0);

        let kmer_start = padded_start.saturating_sub(self.kmer_size / 2);
        let kmer_end = min(self.contig_len - self.kmer_size + 1, padded_end - self.kmer_size / 2);
        let uniq_kmer_frac =
            f64::from(self.cumul_uniq_kmers[kmer_end as usize] - self.cumul_uniq_kmers[kmer_start as usize])
            / f64::from(kmer_end - kmer_start);
        let weight = self.kmers_weight_calc.as_ref().map(|calc| calc.get(uniq_kmer_frac)).unwrap_or(1.0)
            * self.compl_weight_calc.as_ref().map(|calc| calc.get(ling_compl)).unwrap_or(1.0)
            * explicit_weight;
        WindowCharacteristics { gc_content, uniq_kmer_frac, ling_compl, explicit_weight, weight }
    }

    #[inline(always)]
    pub fn n_windows(&self) -> u32 {
        self.default_weights.len() as u32
    }

    #[inline(always)]
    pub fn window_size(&self) -> u32 {
        self.window_getter.window
    }

    /// Default window weights (can change due to random tweaking).
    #[inline(always)]
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

    /// Returns default window boundaries (without tweaking).
    pub fn default_windows<'a>(&'a self) -> impl Iterator<Item = (u32, u32)> + 'a {
        (0..self.n_windows()).map(move |i| self.window_getter.ith_window(i))
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

    #[inline(always)]
    pub fn window_getter(&self) -> &WindowGetter {
        &self.window_getter
    }

    /// Returns the weight for a read end by taking maximum value over the middle of the read,
    /// as well as half window shifts to both sides.
    pub fn read_end_weight(&self, middle: Option<u32>) -> f64 {
        let Some(i) = middle else { return 0.0 };
        let weights = self.explicit_weights.as_ref().expect("Explicit weights must be defined");
        let i = i as usize;
        let n = weights.len();
        let u = self.window_size() as usize / 2;
        weights.at(i)
            .max(weights.at(i.saturating_sub(u)))
            .max(weights.at(min(i + u, n - 1)))
    }
}

/// Contig informations for all haplotypes.
pub struct ContigInfos {
    infos: Vec<Arc<ContigInfo>>,
    has_explicit_weights: bool,
}

impl ContigInfos {
    /// Creates a set of contig windows for each contig.
    pub fn new(
        set: &ContigSet,
        weights_filename: Option<&Path>,
        depth: &ReadDepth,
        params: &super::Params,
        mut dbg_writer: impl Write,
    ) -> crate::Result<Self>
    {
        let contigs = set.contigs();
        let explicit_weights = match weights_filename {
            Some(fname) => load_explicit_weights(fname, contigs)?,
            None => vec![None; contigs.len()],
        };
        let has_explicit_weights = explicit_weights[0].is_some();

        let seqs = set.seqs();
        writeln!(dbg_writer, "#contig\tstart\tend\tGC\tfrac_unique\tcomplexity\texplicit_weight\tweight")
            .map_err(add_path!(!))?;
        let mut infos = Vec::with_capacity(seqs.len());
        for (id, (seq, curr_weights)) in contigs.ids().zip(seqs.iter().zip(explicit_weights.into_iter())) {
            infos.push(
                ContigInfo::new(id, contigs, seq, set.kmer_counts(), depth, curr_weights, params, &mut dbg_writer)
                .map(Arc::new)?);
        }
        Ok(Self { infos, has_explicit_weights })
    }

    /// If explicit weights are defined, returns average read pair weight.
    pub fn explicit_read_weight(&self, pair_alns: &[PairAlignment]) -> f64 {
        if !self.has_explicit_weights {
            return 1.0;
        }
        let mut s = 0.0;
        for pair in pair_alns {
            let info = &self.infos[pair.contig_id().ix()];
            s += f64::max(info.read_end_weight(pair.middle1()), info.read_end_weight(pair.middle2()));
        }
        s / pair_alns.len() as f64
    }
}

impl Index<usize> for ContigInfos {
    type Output = ContigInfo;

    fn index(&self, i: usize) -> &Self::Output {
        self.infos.index(i)
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
    pub fn new(genotype: Genotype, all_infos: &ContigInfos) -> Self {
        let n = genotype.ploidy();
        let mut infos = Vec::with_capacity(n);
        let mut wshifts = Vec::with_capacity(n + 1);
        let mut curr_wshift = REG_WINDOW_SHIFT;
        wshifts.push(curr_wshift);

        for &id in genotype.ids() {
            let curr_contig = Arc::clone(&all_infos.infos[id.ix()]);
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
            let alns = groupped_alns.contig_alns(contig_id);
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
