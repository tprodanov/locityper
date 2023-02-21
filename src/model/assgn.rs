use std::{
    cmp::{min, max},
    io::{self, Error, ErrorKind},
    rc::Rc,
    cell::RefCell,
};
use intmap::{IntMap, Entry};
#[cfg(feature = "stochastic")]
use rand::Rng;
use crate::{
    algo::vec_ext::F64Ext,
    seq::{
        seq,
        contigs::{ContigNames, ContigId},
        interv::Interval,
    },
    math::{
        Ln,
        distr::*,
    },
    bg::{self,
        depth::GC_BINS,
    },
    model::locs::{TwoIntervals, PairAlignment, ContigsGroup, AllPairAlignments},
};

/// Stores k-mer counts for each input k-mer.
pub struct KmerCounts {
    k: u32,
    counts: Vec<Vec<u16>>,
}

impl KmerCounts {
    /// Load k-mer counts from a string.
    /// First line is "k=<number>". All consecutive lines contain just a single number.
    /// Must contain exact number of k-mer as `contigs`.
    pub fn load(s: &str, contigs: &ContigNames) -> io::Result<Self> {
        assert!(!contigs.is_empty(), "Cannot load k-mer counts for empty contigs set!");
        let mut split = s.split('\n');
        let first = split.next().ok_or_else(|| Error::new(ErrorKind::InvalidData,
            "Empty file with k-mer counts!"))?;
        let k = match first.split_once('=') {
            Some(("k", v)) => v.parse().map_err(|e: std::num::ParseIntError|
                Error::new(ErrorKind::InvalidData, e))?,
            _ => return Err(Error::new(ErrorKind::InvalidData,
                format!("Incorrect k-mer counts format, first line ({}) must be in format \"k=integer\"", first))),
        };

        let mut counts = Vec::with_capacity(contigs.len());
        for contig_len in contigs.lengths() {
            assert!(contig_len >= k, "Contig too short!");
            let n_kmers = contig_len - k + 1;
            let mut curr_counts = Vec::with_capacity(n_kmers as usize);
            for _ in 0..n_kmers {
                let val: u16 = split.next()
                    .ok_or_else(|| Error::new(ErrorKind::InvalidData, "Not enough k-mer counts!"))?
                    .parse().map_err(|e: std::num::ParseIntError| Error::new(ErrorKind::InvalidData, e))?;
                // We assume that each k-mer appears at least once.
                curr_counts.push(max(val, 1));
            }
            counts.push(curr_counts);
        }
        match split.next() {
            Some("") | None => Ok(Self { k, counts }),
            _ => Err(Error::new(ErrorKind::InvalidData, "Too many k-mer counts!")),
        }
    }
}

/// Fake proxy distribution, that has 1.0 probability for all values.
pub struct AlwaysOneDistr;

impl DiscretePmf for AlwaysOneDistr {
    fn ln_pmf(&self, _: u32) -> f64 { 0.0 }
}

impl DiscreteCdf for AlwaysOneDistr {
    fn cdf(&self, _: u32) -> f64 {
        unreachable!("CDF is called on a proxy distribution!")
    }

    fn sf(&self, _: u32) -> f64 {
        unreachable!("SF is called on a proxy distribution!")
    }
}

/// All depth values over 2x 0.99 quantile of the Neg.Binomial distribution will be binned together, in one bin.
const DEPTH_BOUND_QUANTILE: f64 = 0.99;
const DEPTH_BOUND_MULT: f64 = 2.0;

/// Store read depth probabilities for values between 0 and 127 for each GC content.
const CACHE_SIZE: usize = 256;

type DistrPtr = Rc<dyn DiscretePmfCdf>;
type DistrBox = Box<dyn DiscretePmfCdf>;

/// Store cached depth distbrutions.
pub struct CachedDepthDistrs<'a> {
    /// Background read depth distribution.
    bg_depth: &'a bg::depth::ReadDepth,
    /// Multiplication coefficient: assume that read depth is `mul_coef` * bg read depth.
    mul_coef: f64,

    /// Distribution for proxy window with unmapped reads.
    unmapped: Rc<AlwaysOneDistr>,
    /// Read depth distributions in regular windows. Key `= GC_BINS * cn + gc_content`.
    regular: RefCell<IntMap<Rc<LinearCache<NBinom>>>>,
    /// Read depth distributions in windows close to the region boundary. Key `= GC_BINS * cn + gc_content`.
    boundary: RefCell<IntMap<Rc<LinearCache<Chimeric<Uniform, NBinom>>>>>,
    /// Cached read depth bounds (see `DEPTH_BOUND_QUANTILE`). Key `= GC_BINS * cn + gc_content`.
    depth_bounds: RefCell<IntMap<u32>>,
}

impl<'a> CachedDepthDistrs<'a> {
    /// Create a set of cached depth distributions.
    /// Assume that there are `mul_coef` as much reads, as in the background distribution.
    ///
    /// As background distribution counts only first read mates,
    /// we assume that `mul_coef` should be either `1.0` for unpaired reads, or `2.0` for paired reads.
    pub fn new(bg_depth: &'a bg::depth::ReadDepth, mul_coef: f64) -> Self {
        Self {
            bg_depth, mul_coef,
            unmapped: Rc::new(AlwaysOneDistr),
            regular: RefCell::new(IntMap::new()),
            boundary: RefCell::new(IntMap::new()),
            depth_bounds: RefCell::new(IntMap::new()),
        }
    }

    /// Returns a pointer to unmapped distribution (`AlwaysOneDistr`).
    pub fn unmapped_distr(&self) -> Rc<AlwaysOneDistr> {
        Rc::clone(&self.unmapped)
    }

    /// Returns read depth distribution in regular windows at GC-content and contig CN.
    pub fn regular_distr(&self, gc_content: usize, contig_cn: usize) -> Rc<LinearCache<NBinom>> {
        let key = (contig_cn * GC_BINS + gc_content) as u64;
        match self.regular.borrow_mut().entry(key) {
            Entry::Occupied(entry) => Rc::clone(&entry.get()),
            Entry::Vacant(entry) => {
                let distr = self.bg_depth.depth_distribution(gc_content)
                    .mul(self.mul_coef * contig_cn as f64)
                    .cached(CACHE_SIZE);
                Rc::clone(&entry.insert(Rc::new(distr)))
            },
        }
    }

    /// Returns read depth distribution in boundary windows at GC-content and contig CN.
    /// Boundary distribution is created from regular Neg.Binomial on the right and uniform distribution on the left.
    pub fn boundary_distr(&self, gc_content: usize, contig_cn: usize) -> Rc<LinearCache<Chimeric<Uniform, NBinom>>> {
        let key = (contig_cn * GC_BINS + gc_content) as u64;
        match self.boundary.borrow_mut().entry(key) {
            Entry::Occupied(entry) => Rc::clone(&entry.get()),
            Entry::Vacant(entry) => {
                let bg_distr = self.bg_depth.depth_distribution(gc_content).mul(self.mul_coef * contig_cn as f64);
                let partition = bg_distr.mean().ceil() as u32;
                let half_uniform = Chimeric::new_smooth(partition, Uniform::new(0, partition), bg_distr);
                let cached_half_uniform = LinearCache::new(half_uniform, CACHE_SIZE);
                Rc::clone(&entry.insert(Rc::new(cached_half_uniform)))
            },
        }
    }

    /// Returns depth bound for the given GC-content and contig CN (see `DEPTH_BOUND_QUANTILE`).
    pub fn depth_bound(&self, gc_content: usize, contig_cn: usize) -> u32 {
        let key = (contig_cn * GC_BINS + gc_content) as u64;
        match self.depth_bounds.borrow_mut().entry(key) {
            Entry::Occupied(entry) => *entry.get(),
            Entry::Vacant(entry) => {
                let bg_distr = self.bg_depth.depth_distribution(gc_content).mul(self.mul_coef * contig_cn as f64);
                let d = (DEPTH_BOUND_MULT * bg_distr.quantile(DEPTH_BOUND_QUANTILE)).ceil() as u32;
                *entry.insert(d)
            }
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

    // Calculating ln(sech(t))
    // // ln(2).
    // const LN2: f64 = 0.69314718055995;
    // LN2 + t - (2.0 * t).exp().ln_1p()
}

/// Stores the contigs and windows corresponding to the windows.
pub struct ContigWindows {
    /// Window size.
    window_size: u32,
    /// Relevant contigs and their copy numbers.
    contigs_group: ContigsGroup,
    /// All contig names.
    contig_names: Rc<ContigNames>,
    /// Start index for each contig id (length: n-contigs + 1).
    /// Contig with index `i` will have windows between `wshifts[i]..wshifts[i + 1]`.
    wshifts: Vec<u32>,
}

/// First window in `ContigWindows` represents an unmapped window.
pub(crate) const UNMAPPED_WINDOW: u32 = 0;
/// First contig will have windows starting from index 1.
pub(crate) const INIT_WSHIFT: u32 = 1;

impl ContigWindows {
    /// Creates new `ContigWindows`.
    /// Functionally, this constructor only counts the number of windows per each contig.
    pub fn new(window_size: u32, contigs_group: ContigsGroup) -> Self {
        let mut wshifts = Vec::with_capacity(contigs_group.len() + 1);
        let mut curr_windows = INIT_WSHIFT;
        wshifts.push(curr_windows);
        let contig_names = Rc::clone(contigs_group.contigs());
        for &id in contigs_group.ids() {
            // Ceiling division.
            curr_windows += (contig_names.get_len(id) + window_size - 1) / window_size;
            wshifts.push(curr_windows);
        }
        Self { window_size, contigs_group, contig_names, wshifts }
    }

    /// Returns the number of contigs.
    pub fn len(&self) -> usize {
        self.contigs_group.len()
    }

    /// Returns the total number of windows.
    pub fn total_windows(&self) -> u32 {
        self.wshifts[self.len()]
    }

    /// Returns the window-shift for `contig` (windows for contig are stored starting with this shift).
    /// Works linearly from the number of contigs.
    pub fn window_shift(&self, contig: ContigId) -> u32 {
        if let Some(i) = self.contigs_group.ids().iter().position(|&id| id == contig) {
            self.wshifts[i]
        } else {
            panic!("Contig {} is not in the contig set", self.contig_names.get_name(contig))
        }
    }

    /// Returns the window corresponding to the middle of the interval.
    pub fn get_window_ix(&self, interval: &Interval) -> u32 {
        self.window_shift(interval.contig_id()) + interval.middle() / self.window_size
    }

    /// Returns a pair of window indices of the pair alignment.
    /// If one (or both) of the mates are unmapped, returns `UNMAPPED_WINDOW`.
    pub fn get_pair_window_ixs(&self, paln: &PairAlignment) -> (u32, u32) {
        match paln.intervals() {
            TwoIntervals::Both(aln1, aln2) => {
                debug_assert_eq!(aln1.contig_id(), aln2.contig_id(), "Read mates are mapped to diff. contigs!");
                let shift = self.window_shift(aln1.contig_id());
                (shift + aln1.middle() / self.window_size, shift + aln2.middle() / self.window_size)
            },
            TwoIntervals::First(aln1) => (self.get_window_ix(aln1), UNMAPPED_WINDOW),
            TwoIntervals::Second(aln2) => (UNMAPPED_WINDOW, self.get_window_ix(aln2)),
            TwoIntervals::None => (UNMAPPED_WINDOW, UNMAPPED_WINDOW),
        }
    }

    /// For each window in the contig group, identifies appropriate read depth distribution.
    fn identify_depth_distributions(
        &mut self,
        ref_seqs: &[Vec<u8>],
        cached_distrs: &CachedDepthDistrs<'_>,
        kmer_counts: &KmerCounts,
        params: &Params,
    ) -> Vec<DistrBox>
    {
        assert!(kmer_counts.k % 2 == 1, "k-mer ({}) size must be odd!", kmer_counts.k);
        let halfk = kmer_counts.k / 2;
        let mut distrs: Vec<DistrBox> = Vec::with_capacity(self.total_windows() as usize);
        debug_assert!(UNMAPPED_WINDOW == 0 && INIT_WSHIFT == 1, "Constants were changed!");
        distrs.push(Box::new(cached_distrs.unmapped_distr()));

        // Normal windows, boundary windows.
        let mut window_counts1 = [0_u16, 0];
        // Stratifying uniqueness: good, adequate, bad.
        let mut window_counts2 = [0_u16, 0, 0];
        let window_size = cached_distrs.bg_depth.window_size();
        let gc_padding = cached_distrs.bg_depth.gc_padding();
        for (i, (&contig_id, &contig_cn)) in
                self.contigs_group.ids().iter().zip(self.contigs_group.multiplicities()).enumerate() {
            let curr_kmer_counts = &kmer_counts.counts[contig_id.ix()];
            let n_windows = self.wshifts[i + 1] - self.wshifts[i];
            let contig_len = self.contig_names.get_len(contig_id);
            let ref_seq = &ref_seqs[contig_id.ix()];
            assert_eq!(contig_len as usize, ref_seq.len(), "Contig length and reference length do not match!");

            for j in 0..n_windows {
                let start = window_size * j;
                let end = start + window_size;
                let mean_kmer_freq = if end - start >= kmer_counts.k {
                    let start_ix = start.saturating_sub(halfk) as usize;
                    let end_ix = min(end - halfk, contig_len - kmer_counts.k + 1) as usize;
                    F64Ext::mean(&curr_kmer_counts[start_ix..end_ix])
                } else { 0.0 };
                let weight = sech_weight(mean_kmer_freq, params.rare_kmer, params.semicommon_kmer);
                // log::debug!("{}:{}  ({}-{})   {:.4}  -> {:.4}", contig_id, j, start, end, mean_kmer_freq, weight);
                if mean_kmer_freq <= params.rare_kmer {
                    window_counts2[0] += 1;
                } else if mean_kmer_freq <= params.semicommon_kmer {
                    window_counts2[1] += 1;
                } else {
                    window_counts2[2] += 1;
                }

                let gc_content = seq::gc_content(
                    &ref_seq[start.saturating_sub(gc_padding) as usize..min(end + gc_padding, contig_len) as usize])
                    .round() as usize;
                let distr: DistrPtr = if end <= params.boundary_size || start + params.boundary_size >= contig_len {
                    window_counts1[1] += 1;
                    cached_distrs.boundary_distr(gc_content, usize::from(contig_cn))
                } else {
                    window_counts1[0] += 1;
                    cached_distrs.regular_distr(gc_content, usize::from(contig_cn))
                };
                let depth_bound = cached_distrs.depth_bound(gc_content, usize::from(contig_cn));

                if weight == 0.0 {
                    distrs.push(Box::new(Uniform::new(0, depth_bound)));
                } else if weight == 1.0 {
                    distrs.push(Box::new(distr));
                } else {
                    distrs.push(Box::new(Mixure::new(distr, weight, Uniform::new(0, depth_bound))));
                }
            }
        }
        log::debug!("    There are {} windows: {} regular, {} boundary;   uniqueness: {} good, {} adequate, {} bad",
            self.total_windows() - 1, window_counts1[0], window_counts1[1],
                window_counts2[0], window_counts2[1], window_counts2[2]);
        debug_assert_eq!(self.total_windows() as usize, distrs.len());
        distrs
    }
}

/// Read depth model parameters.
#[derive(Clone, Debug)]
pub struct Params {
    /// Boundary size: in the left-most and right-most `boundary_size` bp, use boundary read depth distributions,
    /// instead of regular ones.
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
    /// Length: `contig_windows.total_windows()`.
    depth: Vec<u32>,

    /// Read depth distribution at each window (length: `contig_windows.total_windows()`).
    depth_distrs: Vec<DistrBox>,

    /// A vector of possible read alignments to the vector of contigs (length: n_reads).
    read_locs: Vec<PairAlignment>,

    /// Store the start of the read-pair alignments in the `possible_read_locs` vector.
    /// Length: `n_reads + 1`, first value = `0`, last value = `possible_read_locs.len()`.
    ///
    /// For read-pair `i`, possible read locations are `possible_read_locs[read_ixs[i]..read_ixs[i + 1]]`.
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
        contigs_group: ContigsGroup,
        all_ref_seqs: &[Vec<u8>],
        kmer_counts: &KmerCounts,
        cached_distrs: &CachedDepthDistrs<'_>,
        all_alns: &AllPairAlignments,
        params: &Params,
    ) -> Self {
        params.check();
        let mut contig_windows = ContigWindows::new(cached_distrs.bg_depth.window_size(), contigs_group);
        let depth_distrs = contig_windows.identify_depth_distributions(
            all_ref_seqs, cached_distrs, kmer_counts, params);

        let mut ix = 0;
        let mut read_locs = Vec::new();
        let mut read_ixs = vec![ix];
        let mut non_trivial_reads = Vec::new();
        for (rp, paired_alns) in all_alns.iter().enumerate() {
            let new_alns = paired_alns.multi_contig_alns(&mut read_locs, &contig_windows.contigs_group,
                params.prob_diff);
            debug_assert!(new_alns > 0, "Read pair {} has zero possible alignment locations", rp);
            assert!(new_alns <= u16::MAX as usize, "Read pair {} has too many alignment locations ({})", rp, new_alns);
            ix += new_alns;
            read_ixs.push(ix);
            if new_alns > 1 {
                non_trivial_reads.push(rp);
            }
        }

        let n_windows = contig_windows.total_windows() as usize;
        let n_reads = read_ixs.len() - 1;
        Self {
            contig_windows,
            depth: vec![0; n_windows],
            depth_distrs, read_locs, read_ixs, non_trivial_reads,
            read_assgn: vec![0; n_reads],
            likelihood: f64::NAN,
        }
    }

    /// Initialize read assignments and return total likelihood.
    /// Must provide function, that provides initial assignment for all read pairs with at least two possible locations.
    /// Signature: `select_init(single_read_locations) -> initial_index`.
    pub fn init_assignments<F>(&mut self, mut select_init: F) -> f64
    where F: FnMut(&[PairAlignment]) -> usize,
    {
        self.depth.fill(0);
        for (rp, assgn_mut) in self.read_assgn.iter_mut().enumerate() {
            let start_ix = self.read_ixs[rp];
            let end_ix = self.read_ixs[rp + 1];
            let assgn;
            if start_ix + 1 == end_ix {
                assgn = 0;
            } else {
                assgn = select_init(&self.read_locs[start_ix..end_ix]);
                assert!(start_ix + assgn < end_ix,
                    "Read pair #{}: impossible read assignment: {} ({} locations)", rp, assgn, end_ix - start_ix);
            }
            let (w1, w2) = self.contig_windows.get_pair_window_ixs(&self.read_locs[start_ix + assgn]);
            self.depth[w1 as usize] += 1;
            self.depth[w2 as usize] += 1;
            *assgn_mut = assgn.try_into().unwrap();
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
            let (w1, w2) = self.contig_windows.get_pair_window_ixs(&self.read_locs[start_ix + assgn]);
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
    pub fn current_read_aln(&self, rp: usize) -> &PairAlignment {
        &self.read_locs[self.read_ixs[rp] + self.read_assgn[rp] as usize]
    }

    /// Returns the total number of possible alignments for the read pair `rp`.
    pub fn count_possible_alns(&self, rp: usize) -> usize {
        self.read_ixs[rp + 1] - self.read_ixs[rp]
    }

    /// Returns possible read locations for the read pair `rp`.
    pub fn possible_read_alns(&self, rp: usize) -> &[PairAlignment] {
        &self.read_locs[self.read_ixs[rp]..self.read_ixs[rp + 1]]
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
        let curr_paln = &self.read_locs[curr_ix];
        let (w1, w2) = self.contig_windows.get_pair_window_ixs(curr_paln);

        for i in start_ix..end_ix {
            if i == curr_ix {
                buffer.push(f64::NAN);
                continue;
            }
            let new_paln = &self.read_locs[i];
            let (w3, w4) = self.contig_windows.get_pair_window_ixs(new_paln);
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

        let old_paln = &self.read_locs[start_ix + curr_assgn];
        let new_paln = &self.read_locs[start_ix + new_assgn];
        let (w1, w2) = self.contig_windows.get_pair_window_ixs(old_paln);
        let (w3, w4) = self.contig_windows.get_pair_window_ixs(new_paln);
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

        let old_paln = &self.read_locs[old_ix];
        let new_paln = &self.read_locs[new_ix];
        let (w1, w2) = self.contig_windows.get_pair_window_ixs(old_paln);
        let (w3, w4) = self.contig_windows.get_pair_window_ixs(new_paln);
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
            .map(|(&start_ix, &assgn)| self.read_locs[start_ix + assgn as usize].ln_prob())
            .sum();
        self.likelihood = depth_lik + aln_lik;
        (depth_lik, aln_lik)
    }

    /// Returns the contigs group, to which the reads are assigned.
    pub fn contigs_group(&self) -> &ContigsGroup {
        &self.contig_windows.contigs_group
    }

    /// Returns all information about windows.
    pub fn contig_windows(&self) -> &ContigWindows {
        &self.contig_windows
    }

    /// Returns read depth distribution for the window.
    pub fn depth_distr(&self, window: usize) -> &DistrBox {
        &self.depth_distrs[window]
    }

    /// Write CSV in the following format (separated with `\t`):
    /// First line: `prefix  NA      NA     likelihood`.
    /// Next lines: `prefix  window  depth  depth_lik`.
    pub fn write_csv<W: io::Write>(&self, prefix: &str, f: &mut W) -> io::Result<()> {
        writeln!(f, "{}\tNA\tNA\t{:.8}", prefix, self.likelihood())?;
        for (w, (&depth, depth_distr)) in self.depth.iter().zip(&self.depth_distrs).enumerate() {
            writeln!(f, "{}\t{}\t{}\t{:.5}", prefix, w, depth, depth_distr.ln_pmf(depth))?;
        }
        Ok(())
    }
}
