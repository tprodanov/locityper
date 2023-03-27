use std::{
    cmp::min,
    rc::Rc,
    io,
};
use once_cell::unsync::OnceCell;
#[cfg(feature = "stochastic")]
use rand::Rng;
use crate::{
    algo::vec_ext::F64Ext,
    seq::{
        seq,
        kmers::KmerCounts,
    },
    math::{
        Ln,
        distr::{DiscretePmf, WithQuantile, Mixure, NBinom, Uniform, LinearCache},
    },
    bg::{self,
        depth::GC_BINS,
    },
};
use super::{
    locs::{AllPairAlignments, TwoIntervals},
    windows::{UNMAPPED_WINDOW, INIT_WSHIFT, ReadWindows, ContigWindows},
};

/// Fake proxy distribution, that has 1.0 probability for all values.
#[derive(Clone, Copy, Debug)]
pub struct AlwaysOneDistr;

impl DiscretePmf for AlwaysOneDistr {
    fn ln_pmf(&self, _: u32) -> f64 { 0.0 }
}

/// Discrete distribution, repeated `count` times.
/// PMF is redefined in the following way:
/// `RepeatedDistr[count].pmf(k) ~= count * inner.pmf(k / count)`.
///
/// This is done as simply modifying NBinom `n` parameter gives too high probability to repeated haplotypes.
/// `2 * NBinom::new(1000, 0.99).ln_pmf(10) << NBinom::new(2000, 0.99).ln_pmf(20)`.
///
/// Not really a distribution, sum probability will be < 1.
#[derive(Clone, Debug)]
struct RepeatedDistr<T> {
    inner: T,
    count: u8,
}

impl<T: DiscretePmf + 'static> RepeatedDistr<T> {
    /// Returns a box to either `inner` (if count is 1), or to `RepeatedDistr` if count > 1.
    fn new_box(inner: T, count: u8) -> DistrBox {
        assert_ne!(count, 0, "Count cannot be 0.");
        if count == 1 {
            Box::new(inner)
        } else {
            Box::new(Self { inner, count })
        }
    }
}

impl<T: DiscretePmf> DiscretePmf for RepeatedDistr<T> {
    fn ln_pmf(&self, k: u32) -> f64 {
        let count = u32::from(self.count);
        let div = k / count;
        let rem = k % count;
        f64::from(rem) * self.inner.ln_pmf(div + 1) + f64::from(count - rem) * self.inner.ln_pmf(div)
    }
}

/// Windows with many common k-mers have a mixure of Neg.Binomial and Uniform distributions.
/// Size of the Uniform distribution is calculated as 2 x 0.99 quantile of the Neg.Binomial distribution.
const UNIFSIZE_QUANTILE: f64 = 0.99;
const UNIFSIZE_MULT: f64 = 2.0;

/// Store read depth probabilities for values between 0 and 255 for each GC content.
const CACHE_SIZE: usize = 256;

type DistrBox = Box<dyn DiscretePmf>;

/// Store cached depth distbrutions.
pub struct CachedDepthDistrs<'a> {
    /// Background read depth distribution.
    bg_depth: &'a bg::depth::ReadDepth,
    /// Multiplication coefficient: assume that read depth is `mul_coef` * bg read depth.
    mul_coef: f64,

    /// Cached read depth distributions in windows with few common k-mers (one for each GC-content).
    cached: [OnceCell<Rc<LinearCache<NBinom>>>; GC_BINS],
    /// Uniform distribution size (one for each GC-content).
    unif_size: [OnceCell<u32>; GC_BINS],
}

impl<'a> CachedDepthDistrs<'a> {
    /// Create a set of cached depth distributions.
    /// Assume that there are `mul_coef` as much reads, as in the background distribution.
    ///
    /// As background distribution counts only first read mates,
    /// we assume that `mul_coef` should be either `1.0` for unpaired reads, or `2.0` for paired reads.
    pub fn new(bg_depth: &'a bg::depth::ReadDepth, mul_coef: f64) -> Self {
        const NBINOM_CELL: OnceCell<Rc<LinearCache<NBinom>>> = OnceCell::new();
        const U32_CELL: OnceCell<u32> = OnceCell::new();
        Self {
            bg_depth, mul_coef,
            cached: [NBINOM_CELL; GC_BINS],
            unif_size: [U32_CELL; GC_BINS],
        }
    }

    /// Returns a pointer to unmapped distribution (`AlwaysOneDistr`).
    pub fn unmapped_distr(&self) -> AlwaysOneDistr {
        AlwaysOneDistr
    }

    /// Returns read depth distribution in regular windows at GC-content and contig CN.
    pub fn regular_distr(&self, gc_content: usize) -> &Rc<LinearCache<NBinom>> {
        self.cached[gc_content].get_or_init(||
            Rc::new(self.bg_depth
                .depth_distribution(gc_content)
                .mul(self.mul_coef)
                .cached(CACHE_SIZE)))
    }

    /// Returns depth bound for the given GC-content and contig CN (see `DEPTH_BOUND_QUANTILE`).
    pub fn uniform_size(&self, gc_content: usize) -> u32 {
        *self.unif_size[gc_content].get_or_init(||
            (UNIFSIZE_MULT * self.regular_distr(gc_content).quantile(UNIFSIZE_QUANTILE)) as u32)
    }

    /// Returns a box to either `RepeatedDistr`, `NBinom`, `Uniform`, or `Mixure<NBinom, Uniform>`,
    /// depending on the CN and `nbinom_weight`.
    pub fn get_distribution(&self, gc_content: usize, cn: u8, nbinom_weight: f64) -> DistrBox {
        if nbinom_weight < 0.00001 {
            RepeatedDistr::new_box(Uniform::new(0, self.uniform_size(gc_content)), cn)
        } else if nbinom_weight > 0.99999 {
            RepeatedDistr::new_box(Rc::clone(&self.regular_distr(gc_content)), cn)
        } else {
            let mixure = Mixure::new(Rc::clone(&self.regular_distr(gc_content)), nbinom_weight,
                Uniform::new(0, self.uniform_size(gc_content)));
            RepeatedDistr::new_box(mixure, cn)
        }
    }
}

/// Calculates weight of a positive value x.
/// If x <= c1, returns 1.
/// If x = c2, returns 0.5 (returns ln 0.5).
/// For values between c1 and infinity, the weight is distributed according to a sech function.
/// ```
/// f(x) = { 1                              if x in [0, c1].
///                x - c1
///        { sech ------- * ln(2 + sqrt3)   if x > c1.
///               c2 - c1
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

/// For each window in the contig group, identifies appropriate read depth distribution.
fn identify_depth_distributions(
    windows: &ContigWindows,
    cached_distrs: &CachedDepthDistrs<'_>,
    ref_seqs: &[Vec<u8>],
    kmer_counts: &KmerCounts,
    params: &Params,
) -> Vec<DistrBox>
{
    let k = kmer_counts.k();
    assert!(k % 2 == 1, "k-mer ({}) size must be odd!", k);
    let halfk = k / 2;
    let mut distrs: Vec<DistrBox> = Vec::with_capacity(windows.n_windows() as usize);
    debug_assert!(UNMAPPED_WINDOW == 0 && INIT_WSHIFT == 1, "Constants were changed!");
    distrs.push(Box::new(cached_distrs.unmapped_distr()));

    // Percentage of unique k-mers: good, adequate, bad.
    let mut window_counts: (u16, u16, u16) = (0, 0, 0);
    let window_size = cached_distrs.bg_depth.window_size();
    assert_eq!(window_size, windows.window_size());
    let gc_padding = cached_distrs.bg_depth.gc_padding();
    for (i, (contig_id, contig_cn)) in windows.contigs_cns().enumerate() {
        let curr_kmer_counts = &kmer_counts.get(contig_id);
        let n_windows = windows.get_n_windows(i);
        let contig_len = windows.contig_names().get_len(contig_id);
        let ref_seq = &ref_seqs[contig_id.ix()];
        assert_eq!(contig_len as usize, ref_seq.len(), "Contig length and reference length do not match!");

        for j in 0..n_windows {
            let start = window_size * j;
            let end = start + window_size;
            let mean_kmer_freq = if min(contig_len, end) - start >= k {
                let start_ix = start.saturating_sub(halfk) as usize;
                let end_ix = min(end - halfk, contig_len - k + 1) as usize;
                F64Ext::mean(&curr_kmer_counts[start_ix..end_ix])
            } else { 0.0 };
            let weight = sech_weight(mean_kmer_freq, params.rare_kmer, params.semicommon_kmer);
            // log::debug!("{}:{}  ({}-{})   {:.4}  -> {:.4}", contig_id, j, start, end, mean_kmer_freq, weight);
            if mean_kmer_freq <= params.rare_kmer {
                window_counts.0 += 1;
            } else if mean_kmer_freq <= params.semicommon_kmer {
                window_counts.1 += 1;
            } else {
                window_counts.2 += 1;
            }
            let gc_content = seq::gc_content(
                &ref_seq[start.saturating_sub(gc_padding) as usize..min(end + gc_padding, contig_len) as usize])
                .round() as usize;
            distrs.push(cached_distrs.get_distribution(gc_content, contig_cn, weight));
        }
    }
    log::debug!("    There are {} windows.   k-mer uniqueness: {} good, {} adequate, {} bad",
        windows.n_windows() - 1, window_counts.0, window_counts.1, window_counts.2);
    debug_assert_eq!(windows.n_windows() as usize, distrs.len());
    distrs
}

/// Read depth model parameters.
#[derive(Clone, Debug)]
pub struct Params {
    /// Boundary size: ignore left- and right-most `boundary_size` bp.
    pub boundary_size: u32,
    /// For each read pair, all alignments less probable than `best_prob - prob_diff` are discarded.
    pub prob_diff: f64,

    /// Average k-mer frequency is calculated for a window in question.
    /// If the value is less-or-equal than `rare_kmer`, the window received a weight = 1.
    /// If the value equals to `semicommon_kmer`, weight would be 0.5.
    pub rare_kmer: f64,
    pub semicommon_kmer: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            boundary_size: 1000,
            prob_diff: Ln::from_log10(8.0),
            rare_kmer: 3.0,
            semicommon_kmer: 5.0,
        }
    }
}

/// Read assignment to a vector of contigs and the corresponding likelihoods.
pub struct ReadAssignment {
    /// Contigs subset, to which reads the reads are assigned,
    /// plus the number of windows (and shifts) at each window.
    contig_windows: ContigWindows,

    /// Current read depth at each window (concatenated across different contigs).
    /// Length: `contig_windows.n_windows()`.
    depth: Vec<u32>,

    /// Read depth distribution at each window (length: `contig_windows.n_windows()`).
    depth_distrs: Vec<DistrBox>,

    /// A vector of possible read alignments to the vector of contigs (length: n_reads).
    read_windows: Vec<ReadWindows>,

    /// Store the start of the read-pair alignments in the `read_windows` vector.
    /// Length: `n_reads + 1`, first value = `0`, last value = `read_windows.len()`.
    ///
    /// For read-pair `i`, possible read locations are `read_windows[read_ixs[i]..read_ixs[i + 1]]`.
    read_ixs: Vec<usize>,

    /// Read pair indices that have > 1 possible read location (length <= n_reads).
    non_trivial_reads: Vec<usize>,

    /// Store current read assignment (length: n_reads).
    read_assgn: Vec<u16>,

    /// Total ln-probability of the current read assignments (read depth probabilities + alignment probabilities).
    likelihood: f64,
}

impl Params {
    fn check(&self) {
        assert!(self.prob_diff >= 0.0, "Probability difference ({:.4}) cannot be negative!", self.prob_diff);
        assert!(self.rare_kmer < self.semicommon_kmer, "k-mer frequency thresholds ({:.4}, {:.4}) are non-increasing",
            self.rare_kmer, self.semicommon_kmer);
    }
}

impl ReadAssignment {
    /// Creates an instance that stores read assignments to given contigs.
    /// Read assignment itself is not stored, call `init_assignments()` to start.
    pub fn new(
        contig_windows: ContigWindows,
        all_alns: &AllPairAlignments,
        cached_distrs: &CachedDepthDistrs<'_>,
        all_ref_seqs: &[Vec<u8>],
        kmer_counts: &KmerCounts,
        params: &Params,
    ) -> Self {
        params.check();
        let depth_distrs = identify_depth_distributions(&contig_windows, cached_distrs, all_ref_seqs, kmer_counts,
            params);

        let mut ix = 0;
        let mut read_windows = Vec::new();
        let mut read_ixs = vec![ix];
        let mut non_trivial_reads = Vec::new();
        for (rp, paired_alns) in all_alns.iter().enumerate() {
            let nw = contig_windows.read_windows(paired_alns, &mut read_windows, params.prob_diff);
            // rp - read pair index,
            // nw - number of possible read windows.
            debug_assert!(nw > 0, "Read pair {} has zero possible alignment locations", rp);
            assert!(nw <= usize::from(u16::MAX), "Read pair {} has too many alignment locations ({})", rp, nw);
            ix += nw;
            read_ixs.push(ix);
            if nw > 1 {
                non_trivial_reads.push(rp);
            }
        }

        Self {
            depth: vec![0; contig_windows.n_windows() as usize],
            read_assgn: vec![0; read_ixs.len() - 1],
            likelihood: f64::NAN,
            contig_windows, depth_distrs, read_windows, read_ixs, non_trivial_reads,
        }
    }

    /// Initialize read assignments and return total likelihood.
    /// Must provide function, that provides initial assignment for all read pairs with at least two possible locations.
    /// Signature: `select_init(read_windows) -> initial_index`.
    pub fn init_assignments<F>(&mut self, mut select_init: F) -> f64
    where F: FnMut(&[ReadWindows]) -> usize,
    {
        self.depth.fill(0);
        for (rp, assgn_mut) in self.read_assgn.iter_mut().enumerate() {
            let start_ix = self.read_ixs[rp];
            let end_ix = self.read_ixs[rp + 1];
            let assgn;
            if start_ix + 1 == end_ix {
                assgn = 0;
            } else {
                assgn = select_init(&self.read_windows[start_ix..end_ix]);
                assert!(start_ix + assgn < end_ix,
                    "Read pair #{}: impossible read assignment: {} ({} locations)", rp, assgn, end_ix - start_ix);
            }
            let (w1, w2) = self.read_windows[start_ix + assgn].windows();
            self.depth[w1 as usize] += 1;
            self.depth[w2 as usize] += 1;
            *assgn_mut = assgn as u16;
        }
        self.recalc_likelihood();
        self.likelihood
    }

    /// Sets current read assignments with the new ones.
    /// This triggers complete recalculation of the model likelihood, which is then returned.
    pub fn set_assignments(&mut self, new_assgn: &[u16]) -> f64 {
        self.read_assgn.clone_from_slice(new_assgn);
        self.depth.fill(0);
        for (rp, assgn) in self.read_assgn.iter().enumerate() {
            let assgn = *assgn as usize;
            let start_ix = self.read_ixs[rp];
            let end_ix = self.read_ixs[rp + 1];
            assert!(start_ix + assgn < end_ix,
                "Read pair #{}: impossible read assignment: {} ({} locations)", rp, assgn, end_ix - start_ix);
            let (w1, w2) = self.read_windows[start_ix + assgn].windows();
            self.depth[w1 as usize] += 1;
            self.depth[w2 as usize] += 1;
        }
        self.recalc_likelihood();
        self.likelihood
    }

    /// Returns the total likelihood of the read assignment.
    pub fn likelihood(&self) -> f64 {
        self.likelihood
    }

    /// Returns slice with indices of all read pairs with at least two alignment locations.
    pub fn non_trivial_reads(&self) -> &[usize] {
        &self.non_trivial_reads
    }

    /// Returns the total number of reads.
    pub fn total_reads(&self) -> usize {
        self.read_assgn.len()
    }

    /// Returns the slice of read assignments.
    pub fn read_assignments(&self) -> &[u16] {
        &self.read_assgn
    }

    /// Returns the current read-pair alignment given the read pair with index `rp`.
    pub fn current_read_aln(&self, rp: usize) -> &ReadWindows {
        &self.read_windows[self.read_ixs[rp] + self.read_assgn[rp] as usize]
    }

    /// Returns the total number of possible alignments for the read pair `rp`.
    pub fn count_possible_alns(&self, rp: usize) -> usize {
        self.read_ixs[rp + 1] - self.read_ixs[rp]
    }

    /// Returns possible read locations for the read pair `rp`.
    pub fn possible_read_alns(&self, rp: usize) -> &[ReadWindows] {
        &self.read_windows[self.read_ixs[rp]..self.read_ixs[rp + 1]]
    }

    /// Iterates over possible read reassignments, and adds likelihood differences to the vector `buffer`.
    /// This way, buffer is extended by `n` floats,
    /// where `n` is the total number of possible locations for the read pair.
    ///
    /// NAN is added for the current read assignment.
    pub fn possible_reassignments(&self, rp: usize, buffer: &mut Vec<f64>) {
        let start_ix = self.read_ixs[rp];
        let end_ix = self.read_ixs[rp + 1];
        assert!(end_ix - start_ix >= 2,
            "Read pair #{} has {} possible locations! Impossible or or useless.", rp, end_ix - start_ix);

        let curr_ix = start_ix + self.read_assgn[rp] as usize;
        debug_assert!(curr_ix < end_ix, "Read pair #{}: impossible current location!", rp);
        let curr_paln = &self.read_windows[curr_ix];
        let (w1, w2) = curr_paln.windows();

        for i in start_ix..end_ix {
            if i == curr_ix {
                buffer.push(f64::NAN);
                continue;
            }
            let new_paln = &self.read_windows[i];
            let (w3, w4) = new_paln.windows();
            let diff = new_paln.ln_prob() - curr_paln.ln_prob() + self.depth_lik_diff(w1, w2, w3, w4);
            buffer.push(diff);
        }
    }

    /// Calculates the probability of a random reassignment from a random read pair `rp` to a random location.
    /// Returns the read pair index, new assignment and the improvement in likelihood.
    /// Does not actually update any assignments.
    #[cfg(feature = "stochastic")]
    pub fn random_reassignment<R: Rng>(&self, rng: &mut R) -> (usize, u16, f64) {
        let rp = self.non_trivial_reads[rng.gen_range(0..self.non_trivial_reads.len())];
        let start_ix = self.read_ixs[rp];
        let end_ix = self.read_ixs[rp + 1];
        let total_assgns = end_ix - start_ix;
        let curr_assgn = self.read_assgn[rp] as usize;
        let new_assgn = if total_assgns < 2 {
            panic!("Read pair #{} less than 2 possible assignments.", rp)
        } else if total_assgns == 2 {
            1 - curr_assgn
        } else {
            let i = rng.gen_range(1..total_assgns);
            if i <= curr_assgn { i - 1 } else { i }
        };
        debug_assert_ne!(new_assgn, curr_assgn);

        let old_paln = &self.read_windows[start_ix + curr_assgn];
        let new_paln = &self.read_windows[start_ix + new_assgn];
        let (w1, w2) = old_paln.windows();
        let (w3, w4) = new_paln.windows();
        let diff = new_paln.ln_prob() - old_paln.ln_prob() + self.depth_lik_diff(w1, w2, w3, w4);
        (rp, new_assgn as u16, diff)
    }

    /// Reassigns read pair `rp` to a new location and returns the difference in likelihood.
    pub fn reassign(&mut self, rp: usize, new_assignment: u16) -> f64 {
        debug_assert_ne!(self.read_assgn[rp], new_assignment,
            "Read pair #{}: cannot reassign read to the same location", rp);
        let start_ix = self.read_ixs[rp];
        let end_ix = self.read_ixs[rp + 1];
        let old_ix = start_ix + self.read_assgn[rp] as usize;
        let new_ix = start_ix + new_assignment as usize;
        assert!(new_ix < end_ix, "Read pair #{}: cannot assign read to location {} (maximum {} locations)",
            rp, new_assignment, end_ix - start_ix);

        let old_paln = &self.read_windows[old_ix];
        let new_paln = &self.read_windows[new_ix];
        let (w1, w2) = old_paln.windows();
        let (w3, w4) = new_paln.windows();
        let diff = new_paln.ln_prob() - old_paln.ln_prob() + self.depth_lik_diff(w1, w2, w3, w4);
        self.likelihood += diff;
        self.read_assgn[rp] = new_assignment;
        self.depth[w1 as usize] -= 1;
        self.depth[w2 as usize] -= 1;
        self.depth[w3 as usize] += 1;
        self.depth[w4 as usize] += 1;
        diff
    }

    /// Assuming that read depth in `window` will change by `depth_change`,
    /// calculates the difference between the new and the old ln-probabilities.
    /// Positive value means that the likelihood will improve.
    /// Does not actually update the read depth.
    fn atomic_depth_lik_diff(&self, window: u32, depth_change: i32) -> f64 {
        if depth_change == 0 {
            0.0
        } else {
            let i = window as usize;
            let old_depth = self.depth[i];
            let new_depth = old_depth.checked_add_signed(depth_change).expect("Read depth became negative");
            self.depth_distrs[i].ln_pmf(new_depth) - self.depth_distrs[i].ln_pmf(old_depth)
        }
    }

    /// Calculate likelihood difference on four windows (some of them can be equal to others).
    /// Depth at w1 and w2 is decreased by one, depth at w3 and w4 is increased by one.
    fn depth_lik_diff(&self, w1: u32, w2: u32, w3: u32, w4: u32) -> f64 {
        // Variables c1..c4 store change in read depth. If cX == 0, it means that wX is the same as some other window.
        let mut c1 = -1;

        let mut c2 = if w2 == w1 {
            c1 -= 1; 0
        } else { -1 };

        let mut c3 = if w3 == w1 {
            c1 += 1; 0
        } else if w3 == w2 {
            c2 += 1; 0
        } else { 1 };

        let c4 = if w4 == w1 {
            c1 += 1; 0
        } else if w4 == w2 {
            c2 += 1; 0
        } else if w4 == w3 {
            c3 += 1; 0
        } else { 1 };

        debug_assert_eq!(c1 + c2 + c3 + c4, 0);
        self.atomic_depth_lik_diff(w1, c1) + self.atomic_depth_lik_diff(w2, c2)
            + self.atomic_depth_lik_diff(w3, c3) + self.atomic_depth_lik_diff(w4, c4)
    }

    /// Recalculates total model likelihood, and returns separately
    /// - ln-probability of read depth across all windows,
    /// - ln-probability of of all read alignments.
    pub fn recalc_likelihood(&mut self) -> (f64, f64) {
        let depth_lik = self.depth_distrs.iter()
            .zip(&self.depth)
            .map(|(distr, &depth)| distr.ln_pmf(depth))
            .sum();
        let aln_lik = self.read_ixs.iter().zip(&self.read_assgn)
            .map(|(&start_ix, &assgn)| self.read_windows[start_ix + assgn as usize].ln_prob())
            .sum();
        self.likelihood = depth_lik + aln_lik;
        (depth_lik, aln_lik)
    }

    /// Returns all information about windows.
    pub fn contig_windows(&self) -> &ContigWindows {
        &self.contig_windows
    }

    /// Returns read depth distribution for the window.
    pub fn depth_distr(&self, window: usize) -> &DistrBox {
        &self.depth_distrs[window]
    }

    /// Write read depth to a CSV file in the following format (tab-separated):
    /// General lines: `prefix  contig  window  depth  depth_lik`.
    /// Last line:     `prefix  NA      read_lik  depth_lik  sum_lik`.
    pub fn write_depth<W: io::Write>(&self, prefix: &str, f: &mut W) -> io::Result<()> {
        let unmapped_reads = self.depth[UNMAPPED_WINDOW as usize];
        let unmapped_prob = self.depth_distrs[UNMAPPED_WINDOW as usize].ln_pmf(unmapped_reads);
        let mut sum_depth_lik = unmapped_prob;
        writeln!(f, "{}\tunmapped\tNA\t{}\t{:.3}", prefix, unmapped_reads, unmapped_prob)?;

        let contigs = self.contig_windows.contig_names();
        for (i, contig_id) in self.contig_windows.ids().enumerate() {
            let curr_prefix = format!("{}\t{}\t", prefix, contigs.get_name(contig_id));
            let wshift = self.contig_windows.get_wshift(i) as usize;
            for w in wshift..self.contig_windows.get_wshift(i + 1) as usize {
                let depth = self.depth[w];
                let ln_prob = self.depth_distrs[w].ln_pmf(depth);
                writeln!(f, "{}{}\t{}\t{:.3}", curr_prefix, w - wshift, depth, ln_prob)?;
                sum_depth_lik += ln_prob;
            }
        }
        let lik = self.likelihood;
        writeln!(f, "{}\tNA\t{:.8}\t{:.8}\t{:.8}", prefix, lik - sum_depth_lik, sum_depth_lik, lik)?;
        Ok(())
    }

    /// Write reads and their assignments to a CSV file in the following format (tab-separated):
    /// `prefix  read_hash  aln1  aln2  w1  w2  prob  selected`
    pub fn write_reads<W: io::Write>(&self, prefix: &str, f: &mut W, all_alns: &AllPairAlignments) -> io::Result<()> {
        assert_eq!(all_alns.len() + 1, self.read_ixs.len());
        for (rp, paired_alns) in all_alns.iter().enumerate() {
            let hash = paired_alns.name_hash();
            let start_ix = self.read_ixs[rp];
            let end_ix = self.read_ixs[rp + 1];
            let curr_ix = start_ix + self.read_assgn[rp] as usize;

            for i in start_ix..end_ix {
                write!(f, "{}\t{:X}\t", prefix, hash)?;

                let curr_windows = &self.read_windows[i];
                let (w1, w2) = curr_windows.windows();
                if w1 == UNMAPPED_WINDOW && w2 == UNMAPPED_WINDOW {
                    write!(f, "*\t*\t")?;
                } else {
                    match paired_alns.ith_aln(curr_windows.ix() as usize).intervals() {
                        TwoIntervals::Both(aln1, aln2) => write!(f, "{}\t{}\t", aln1, aln2),
                        TwoIntervals::First(aln1) => write!(f, "{}\t*\t", aln1),
                        TwoIntervals::Second(aln2) => write!(f, "*\t{}\t", aln2),
                    }?;
                }
                writeln!(f, "{}\t{}\t{:.3}\t{}", w1, w2, curr_windows.ln_prob(), (i == curr_ix) as u8)?;
            }
        }
        Ok(())
    }
}
