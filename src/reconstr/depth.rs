use std::rc::Rc;
use once_cell::unsync::OnceCell;
use statrs::distribution::Discrete;
use crate::{
    seq::{
        contigs::{ContigNames, ContigId},
        interv::Interval,
    },
    algo::{
        nbinom::{NBinom, UniformNBinom, CachedDistr},
    },
    bg::{self,
        depth::GC_BINS,
    },
    reconstr::locs::{TwoIntervals, PairAlignment},
};

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

pub const UNMAPPED_WINDOW: u32 = u32::MAX;

impl ContigWindows {
    /// Creates new `ContigWindows`.
    /// Functionally, this constructor only counts the number of windows per each contig.
    pub fn new(window_size: u32, contigs: Rc<ContigNames>, ids: Vec<ContigId>) -> Self {
        let mut wshifts = Vec::with_capacity(ids.len() + 1);
        wshifts.push(0);
        let mut curr_windows = 0;
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
        for (&id, &shift) in self.ids.iter().zip(&self.wshifts) {
            if id == contig {
                return shift;
            }
        }
        panic!("Contig {} is not in the contig set", self.contigs.name(contig));
    }

    /// Returns the window corresponding to the middle of the interval.
    pub fn get_window(&self, interval: &Interval) -> u32 {
        self.window_shift(interval.contig_id()) + interval.middle() / self.window_size
    }

    /// Returns a pair of window indices of the pair alignment.
    /// If one (or both) of the mates are unmapped, returns `UNMAPPED_WINDOW`.
    pub fn get_windows(&self, paln: &PairAlignment) -> (u32, u32) {
        match paln.intervals() {
            TwoIntervals::Both(aln1, aln2) => {
                debug_assert_eq!(aln1.contig_id(), aln2.contig_id(), "Read mates are mapped to diff. contigs!");
                let shift = self.window_shift(aln1.contig_id());
                (shift + aln1.middle() / self.window_size, shift + aln2.middle() / self.window_size)
            },
            TwoIntervals::First(aln1) => (self.get_window(aln1), UNMAPPED_WINDOW),
            TwoIntervals::Second(aln2) => (UNMAPPED_WINDOW, self.get_window(aln2)),
            TwoIntervals::None => (UNMAPPED_WINDOW, UNMAPPED_WINDOW),
        }
    }
}

/// Read assignment to a vector of contigs and the corresponding likelihoods.
pub struct ReadAssignment<'a> {
    /// Contigs subset, to which reads the reads are assigned,
    /// plus the number of windows (and shifts) at each window.
    contig_windows: ContigWindows,

    /// Current read depth at each window (concatenated across different contigs).
    /// Length: `contig_windows.total_windows()`.
    depth: Vec<u32>,

    /// Read depth distribution at each window (length: `contig_windows.total_windows()`).
    depth_distrs: Vec<Box<&'a dyn Discrete<u32, f64>>>,

    /// A vector of possible read alignments to the vector of contigs (length: n_reads).
    read_locs: Vec<PairAlignment>,

    /// Store the start of the read-pair alignments in the `possible_read_locs` vector.
    /// Length: `n_reads + 1`, first value = `0`, last value = `possible_read_locs.len()`.
    ///
    /// For read-pair `i`, possible read locations are `possible_read_locs[read_ixs[i]..read_ixs[i + 1]]`.
    read_ixs: Vec<usize>,

    /// Store current read assignment (length: n_reads).
    read_assgn: Vec<u16>,

    /// Total ln-probability of the current read assignments (read depth probabilities + alignment probabilities).
    likelihood: f64,
}

impl<'a> ReadAssignment<'a> {
    /// Returns the total likelihood of the read assignment.
    pub fn likelihood(&self) -> f64 {
        self.likelihood
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
        let curr_windows = self.contig_windows.get_windows(curr_paln);
        // Total likelihood difference will be:
        // P(new alignment) - P(old alignment) + P(increase depth in new windows) + P(decrease depth in old windows).
        let base_diff = -curr_paln.ln_prob() + self.depth_lik_diff2(curr_windows, -1);

        for i in start_ix..end_ix {
            if i == curr_ix {
                buffer.push(f64::NAN);
                continue;
            }
            let new_paln = &self.read_locs[i];
            let new_windows = self.contig_windows.get_windows(new_paln);
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
        let old_windows = self.contig_windows.get_windows(old_paln);
        let new_windows = self.contig_windows.get_windows(new_paln);
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
    fn depth_lik_diff2(&self, windows: (u32, u32), depth_change: i32) -> f64 {
        match windows {
            (UNMAPPED_WINDOW, UNMAPPED_WINDOW) => 0.0,
            (w, UNMAPPED_WINDOW) | (UNMAPPED_WINDOW, w) => self.depth_lik_diff(w, depth_change),
            (w1, w2) => if w1 == w2 {
                self.depth_lik_diff(w1, 2 * depth_change)
            } else {
                self.depth_lik_diff(w1, depth_change) + self.depth_lik_diff(w2, depth_change)
            },
        }
    }

    /// Re-calculates the total likehood (ln-probability) of the read depth across all windows.
    fn recalc_depth_likelihood(&self) -> f64 {
        self.depth_distrs.iter()
            .zip(&self.depth)
            .map(|(distr, &depth)| distr.ln_pmf(depth))
            .sum()
    }

    /// Re-calculates the total likehood (ln-probability) of all read alignments.
    fn recalc_alns_likelihood(&self) -> f64 {
        self.read_ixs.iter().zip(&self.read_assgn)
            .map(|(&start_ix, &assgn)| self.read_locs[start_ix + assgn as usize].ln_prob())
            .sum()
    }
}
