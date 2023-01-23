use std::rc::Rc;
use std::cmp::min;
use once_cell::unsync::OnceCell;
use statrs::distribution::Discrete;
use crate::{
    seq::{
        seq,
        contigs::{ContigNames, ContigId},
        interv::Interval,
    },
    algo::{
        nbinom::{NBinom, UniformNBinom, CachedDistr},
    },
    bg::{self,
        depth::GC_BINS,
    },
    reconstr::locs::{TwoIntervals, PairAlignment, SeveralContigs, AllPairAlignments},
};

/// Fake proxy distribution, that has 1.0 probability for all values.
struct AlwaysOneDistr;

impl Discrete<u32, f64> for AlwaysOneDistr {
    fn pmf(&self, _: u32) -> f64 { 1.0 }
    fn ln_pmf(&self, _: u32) -> f64 { 0.0 }
}

/// Store read depth probabilities for values between 0 and 127 for each GC content.
const CACHE_SIZE: usize = 256;

/// Store cached depth distbrutions.
pub struct CachedDepthDistrs<'a> {
    /// Background read depth distribution.
    bg_depth: &'a bg::depth::ReadDepth,
    /// Multiplication coefficient: assume that read depth is `mul_coef` * bg read depth.
    mul_coef: f64,
    /// Regular read depth distributions (`* mul_coef`).
    regular: Vec<OnceCell<CachedDistr<NBinom>>>,
    /// Read depth distributions at the boundary windows.
    boundary: Vec<OnceCell<CachedDistr<UniformNBinom>>>,
    /// Distribution for proxy window with unmapped reads.
    unmapped: AlwaysOneDistr,
}

impl<'a> CachedDepthDistrs<'a> {
    /// Create a set of cached depth distributions.
    /// Assume that there are `mul_coef` as much reads, as in the background distrbution.
    ///
    /// As background distribution counts only first read mates,
    /// we assume that `mul_coef` should be either `1.0` for unpaired reads, or `2.0` for paired reads.
    pub fn new(bg_depth: &'a bg::depth::ReadDepth, mul_coef: f64) -> Self {
        Self {
            bg_depth, mul_coef,
            regular: vec![OnceCell::new(); GC_BINS],
            boundary: vec![OnceCell::new(); GC_BINS],
            unmapped: AlwaysOneDistr,
        }
    }

    /// Regular read depth distribution at GC-content (`bg_depth * mul_coef`).
    pub fn regular_distr(&self, gc_content: usize) -> &CachedDistr<NBinom> {
        self.regular[gc_content].get_or_init(||
            self.bg_depth.depth_distribution(gc_content).mul(self.mul_coef).cached(CACHE_SIZE))
    }

    /// Boundary read depth distribution at GC-content (`bg_depth * mul_coef`).
    ///
    /// Designated for boundary windows, this distribution is uniform at values < mean NBinom value.
    pub fn boundary_distr(&self, gc_content: usize) -> &CachedDistr<UniformNBinom> {
        self.boundary[gc_content].get_or_init(|| {
            let distr = self.bg_depth.depth_distribution(gc_content).mul(self.mul_coef);
            let half_uniform = UniformNBinom::new(distr.mean().ceil() as u32, distr, None);
            CachedDistr::new(half_uniform, CACHE_SIZE)
        })
    }
}

/// Stores the contigs and windows corresponding to the windows.
pub struct ContigWindows {
    /// Window size.
    window_size: u32,
    /// All contig names.
    contigs: Rc<ContigNames>,
    /// Subset of contigs, to which the reads are assigned (length: n_contigs).
    ids: Vec<ContigId>,
    /// Start index for each contig id (length: n-contigs + 1).
    /// Contig with index `i` will have windows between `wshifts[i]..wshifts[i + 1]`.
    wshifts: Vec<u32>,
}

/// First window in `ContigWindows` represents an unmapped window.
const UNMAPPED_WINDOW: u32 = 0;
/// First contig will have windows starting from index 1.
const INIT_WSHIFT: u32 = 1;

impl ContigWindows {
    /// Creates new `ContigWindows`.
    /// Functionally, this constructor only counts the number of windows per each contig.
    pub fn new(window_size: u32, contigs: Rc<ContigNames>, ids: Vec<ContigId>) -> Self {
        let mut wshifts = Vec::with_capacity(ids.len() + 1);
        let mut curr_windows = INIT_WSHIFT;
        wshifts.push(curr_windows);
        for &id in ids.iter() {
            // Ceiling division.
            curr_windows += (contigs.length(id) + window_size - 1) / window_size;
            wshifts.push(curr_windows);
        }
        Self { window_size, contigs, ids, wshifts }
    }

    /// Returns the number of contigs.
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// Returns the total number of windows.
    pub fn total_windows(&self) -> u32 {
        self.wshifts[self.len()]
    }

    /// Returns the window-shift for `contig` (windows for contig are stored starting with this shift).
    /// Works linearly from the number of contigs.
    pub fn window_shift(&self, contig: ContigId) -> u32 {
        if let Some(i) = self.ids.iter().position(|&id| id == contig) {
            self.wshifts[i]
        } else {
            panic!("Contig {} is not in the contig set", self.contigs.name(contig))
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

    /// Returns read depth distributions for all windows in a set of contigs.
    fn identify_depth_distributions<'a>(
        &self,
        ref_seqs: &[Vec<u8>],
        cached_distrs: &'a CachedDepthDistrs<'a>,
        boundary_size: u32,
    ) -> DistrsVec<'a> {
        let mut distrs: DistrsVec<'a> = Vec::with_capacity(self.total_windows() as usize);
        debug_assert!(UNMAPPED_WINDOW == 0 && INIT_WSHIFT == 1, "Constants were changed!");
        distrs.push(Box::new(&cached_distrs.unmapped));

        let window_size = cached_distrs.bg_depth.window_size();
        let gc_padding = cached_distrs.bg_depth.gc_padding();
        for (i, &contig_id) in self.ids.iter().enumerate() {
            let n_windows = self.wshifts[i + 1] - self.wshifts[i];
            let contig_len = self.contigs.length(contig_id);
            let ref_seq = &ref_seqs[contig_id.ix()];
            assert_eq!(contig_len as usize, ref_seq.len(), "Contig length and reference length do not match!");

            for i in 0..n_windows {
                let start = window_size * i;
                let end = start + window_size;
                let gc_content = seq::gc_content(
                    &ref_seq[start.saturating_sub(gc_padding) as usize..min(end + gc_padding, contig_len) as usize])
                    .round() as usize;
                if end <= boundary_size || start + boundary_size >= contig_len {
                    distrs.push(Box::new(cached_distrs.boundary_distr(gc_content)));
                } else {
                    distrs.push(Box::new(cached_distrs.regular_distr(gc_content)));
                }
            }
        }
        distrs
    }
}

type DistrsVec<'a> = Vec<Box<&'a dyn Discrete<u32, f64>>>;

/// Read assignment to a vector of contigs and the corresponding likelihoods.
pub struct ReadAssignment<'a> {
    /// Contigs subset, to which reads the reads are assigned,
    /// plus the number of windows (and shifts) at each window.
    contig_windows: ContigWindows,

    /// Current read depth at each window (concatenated across different contigs).
    /// Length: `contig_windows.total_windows()`.
    depth: Vec<u32>,

    /// Read depth distribution at each window (length: `contig_windows.total_windows()`).
    depth_distrs: DistrsVec<'a>,

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

impl<'a> ReadAssignment<'a> {
    /// Creates an instance that stores read assignments to given contigs.
    /// Read assignment itself is not stored, call `init_assignments()` to start.
    ///
    /// Boundary size: in the left-most and right-most `boundary_size` bp, use boundary read depth distributions,
    /// instead of regular ones.
    ///
    /// For each read pair, all alignments less probable than `best_prob - prob_diff` are discarded.
    pub fn new<F>(
        contigs: &SeveralContigs,
        contig_names: Rc<ContigNames>,
        all_ref_seqs: &[Vec<u8>],
        cached_distrs: &'a CachedDepthDistrs<'a>,
        all_alns: &AllPairAlignments,
        boundary_size: u32,
        prob_diff: f64,
    ) -> Self {
        let contig_windows = ContigWindows::new(cached_distrs.bg_depth.window_size(), contig_names,
            contigs.contigs().to_vec());
        let depth_distrs = contig_windows.identify_depth_distributions(all_ref_seqs, cached_distrs, boundary_size);

        let mut ix = 0;
        let mut read_locs = Vec::new();
        let mut read_ixs = vec![ix];
        let mut non_trivial_reads = Vec::new();
        for (rp, paired_alns) in all_alns.iter().enumerate() {
            let new_alns = paired_alns.multi_contig_alns(&mut read_locs, &contigs, prob_diff);
            debug_assert!(new_alns > 0, "Read pair {} has zero possible alignment location", rp);
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
            depth_distrs,
            read_locs,
            read_ixs,
            non_trivial_reads,
            read_assgn: vec![0; n_reads],
            likelihood: f64::NAN,
        }
    }

    /// Initialize read assignments.
    /// Must provide function, that provides initial assignment for all read pairs with at least two possible locations.
    /// Signature: `select_init(single_read_locations) -> initial_index as u16`
    pub fn init_assignments<F>(&mut self, mut select_init: F)
    where F: FnMut(&[PairAlignment]) -> u16,
    {
        for (rp, assgn_mut) in self.read_assgn.iter_mut().enumerate() {
            let start_ix = self.read_ixs[rp];
            let end_ix = self.read_ixs[rp + 1];
            let assgn;
            if start_ix + 1 == end_ix {
                assgn = 0;
            } else {
                assgn = select_init(&self.read_locs[start_ix..end_ix]);
                assert!(start_ix + (assgn as usize) < end_ix,
                    "Read pair #{}: impossible read assignment: {} ({} locations)", rp, assgn, end_ix - start_ix);
            }
            let (w1, w2) = self.contig_windows.get_pair_window_ixs(&self.read_locs[start_ix + assgn as usize]);
            self.depth[w1 as usize] += 1;
            self.depth[w2 as usize] += 1;
            *assgn_mut = assgn;
        }
        self.recalc_likelihood();
    }

    /// Returns the total likelihood of the read assignment.
    pub fn likelihood(&self) -> f64 {
        self.likelihood
    }

    /// Returns slice with indices of all read pairs with at least two alignment locations.
    pub fn non_trivial_reads(&self) -> &[usize] {
        &self.non_trivial_reads
    }

    /// Returns the current read-pair alignment given the read pair with index `rp`.
    pub fn current_pair_alignment(&self, rp: usize) -> &PairAlignment {
        &self.read_locs[self.read_ixs[rp] + self.read_assgn[rp] as usize]
    }

    /// Returns the total number of possible alignments for the read pair `rp`.
    pub fn count_possible_alns(&self, rp: usize) -> usize {
        self.read_ixs[rp + 1] - self.read_ixs[rp]
    }

    /// Iterates over possible read reassignments, and adds likelihood differences to the vector `buffer`.
    /// This way, buffer is extended by `n` floats,
    /// where `n` is the total number of possible locations for the read pair.
    ///
    /// NAN is added for the current read assignment.
    pub fn possible_reassignments(&self, rp: usize, buffer: &mut Vec<f64>) {
        let start_ix = self.read_ixs[rp];
        let end_ix = self.read_ixs[rp + 1];
        assert!(start_ix + 1 < end_ix,
            "Read pair #{} has {} possible locations! Impossible or or useless.", rp, end_ix - start_ix);

        let curr_ix = start_ix + self.read_assgn[rp] as usize;
        debug_assert!(curr_ix < end_ix, "Read pair #{}: impossible current location!", rp);
        let curr_paln = &self.read_locs[curr_ix];
        let curr_windows = self.contig_windows.get_pair_window_ixs(curr_paln);
        // Total likelihood difference will be:
        // P(new alignment) - P(old alignment) + P(increase depth in new windows) + P(decrease depth in old windows).
        let base_diff = -curr_paln.ln_prob() + self.depth_lik_diff2(curr_windows, -1);

        for i in start_ix..end_ix {
            if i == curr_ix {
                buffer.push(f64::NAN);
                continue;
            }
            let new_paln = &self.read_locs[i];
            let new_windows = self.contig_windows.get_pair_window_ixs(new_paln);
            let diff = base_diff + new_paln.ln_prob() + self.depth_lik_diff2(new_windows, 1);
            buffer.push(diff);
        }
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
        let old_windows = self.contig_windows.get_pair_window_ixs(old_paln);
        let new_windows = self.contig_windows.get_pair_window_ixs(new_paln);
        let diff = new_paln.ln_prob() - old_paln.ln_prob()
            + self.depth_lik_diff2(old_windows, -1) + self.depth_lik_diff2(new_windows, 1);
        self.likelihood = diff;
        self.read_assgn[rp] = new_assignment;
        diff
    }

    /// Assuming that read depth in `window` will change by `depth_change`,
    /// calculates the difference between the new and the old ln-probabilities.
    /// Positive value means that the likelihood will improve.
    /// Does not actually update the read depth.
    ///
    /// Windows must not be `UNMAPPED_WINDOW`.
    fn depth_lik_diff(&self, window: u32, depth_change: i32) -> f64 {
        let i = window as usize;
        let old_depth = self.depth[i];
        let new_depth = old_depth.checked_add_signed(depth_change).expect("Read depth became negative!");
        self.depth_distrs[i].ln_pmf(new_depth) - self.depth_distrs[i].ln_pmf(old_depth)
    }

    /// Wraps `depth_lik_diff` function, and runs on a pair windows.
    fn depth_lik_diff2(&self, (w1, w2): (u32, u32), depth_change: i32) -> f64 {
        if w1 == w2 {
            self.depth_lik_diff(w1, 2 * depth_change)
        } else {
            self.depth_lik_diff(w1, depth_change) + self.depth_lik_diff(w2, depth_change)
        }
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
}
