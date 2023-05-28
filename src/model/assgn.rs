use std::{
    io,
};
use rand::Rng;
use crate::{
    math::Ln,
};
use super::{
    locs::AllPairAlignments,
    windows::{UNMAPPED_WINDOW, BOUNDARY_WINDOW, ReadWindows, MultiContigWindows},
    dp_cache::{CachedDepthDistrs, DistrBox},
};

/// Count how many times each alignment was selected for a read.
#[derive(Default)]
pub struct SelectedCounter {
    buffer: Vec<u16>,
    total_count: u16,
}

impl SelectedCounter {
    /// Resets counts and updates buffer length to the new read assignment.
    pub fn reset(&mut self, assgn: &ReadAssignment) {
        self.buffer.clear();
        self.buffer.resize(assgn.total_possible_alns(), 0);
        self.total_count = 0;
    }

    /// Increments selected counts.
    pub fn update(&mut self, assgn: &ReadAssignment) {
        self.total_count += 1;
        for (&read_ix, &read_assgn) in assgn.read_ixs.iter().zip(&assgn.read_assgn) {
            self.buffer[read_ix + usize::from(read_assgn)] += 1;
        }
    }

    /// Returns iterator over fractions (how many times each location was selected / number of iterations).
    pub fn fractions(&self) -> impl Iterator<Item = f64> + '_ {
        let mult = 1.0 / f64::from(self.total_count);
        self.buffer.iter().map(move |&count| mult * f64::from(count))
    }
}

/// Read assignment to a vector of contigs and the corresponding likelihoods.
pub struct ReadAssignment {
    /// Contigs subset, to which reads the reads are assigned,
    /// plus the number of windows (and shifts) at each window.
    contig_windows: MultiContigWindows,

    /// Current read depth at each window (concatenated across different contigs).
    /// Length: `contig_windows.total_windows()`.
    depth: Vec<u32>,

    /// Read depth distribution at each window (length: `contig_windows.total_windows()`).
    depth_distrs: Vec<DistrBox>,

    /// Read depth contribution, relative to read alignment likelihoods.
    depth_contrib: f64,

    /// A vector of possible read alignments to the vector of contigs.
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

impl ReadAssignment {
    /// Creates an instance that stores read assignments to given contigs.
    /// Read assignment itself is not stored, call `init_assignments()` to start.
    pub fn new(
        contig_windows: MultiContigWindows,
        all_alns: &AllPairAlignments,
        cached_distrs: &CachedDepthDistrs,
        params: &super::Params,
    ) -> Self
    {
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
            depth: vec![0; contig_windows.total_windows() as usize],
            read_assgn: vec![0; read_ixs.len() - 1],
            likelihood: f64::NEG_INFINITY,
            depth_distrs: contig_windows.get_distributions(cached_distrs),
            depth_contrib: params.depth_contrib,
            contig_windows, read_windows, read_ixs, non_trivial_reads,
        }
    }

    /// Define read windows by randomly moving read middle by at most `tweak` bp to either side.
    pub fn define_read_windows(&mut self, tweak: u32, rng: &mut impl Rng) {
        if tweak == 0 {
            self.read_windows.iter_mut().for_each(|rw| rw.define_windows_determ(&self.contig_windows));
        } else {
            self.read_windows.iter_mut().for_each(|rw| rw.define_windows_random(&self.contig_windows, tweak, rng));
        }
        self.likelihood = f64::NEG_INFINITY;
    }

    /// Try to initialize read assignments and return total likelihood.
    /// Same as `init_assignments`, but the `select_init` function may fail.
    pub fn try_init_assignments<F, E>(&mut self, mut select_init: F) -> Result<f64, E>
    where F: FnMut(&[ReadWindows]) -> Result<usize, E>,
    {
        self.depth.fill(0);
        for (rp, assgn_mut) in self.read_assgn.iter_mut().enumerate() {
            let start_ix = self.read_ixs[rp];
            let end_ix = self.read_ixs[rp + 1];
            let assgn;
            if start_ix + 1 == end_ix {
                assgn = 0;
            } else {
                assgn = select_init(&self.read_windows[start_ix..end_ix])?;
                assert!(start_ix + assgn < end_ix,
                    "Read pair #{}: impossible read assignment: {} ({} locations)", rp, assgn, end_ix - start_ix);
            }
            let (w1, w2) = self.read_windows[start_ix + assgn].windows();
            self.depth[w1 as usize] += 1;
            self.depth[w2 as usize] += 1;
            *assgn_mut = assgn as u16;
        }
        self.recalc_likelihood();
        Ok(self.likelihood)
    }

    /// Initialize read assignments and return total likelihood.
    /// Must provide function, that provides initial assignment for all read pairs with at least two possible locations.
    /// Signature: `select_init(read_windows) -> initial_index`.
    pub fn init_assignments<F>(&mut self, mut select_init: F) -> f64
    where F: FnMut(&[ReadWindows]) -> usize,
    {
        // unwrap as select_init never returns Err.
        self.try_init_assignments::<_, ()>(|windows| Ok(select_init(windows))).unwrap()
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

    /// Consumes self and returns the vector with read assignments.
    pub fn take_read_assignments(self) -> Vec<u16> {
        self.read_assgn
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

    /// Returns the total possible number of read alignments accross all read pairs.
    pub fn total_possible_alns(&self) -> usize {
        self.read_windows.len()
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
    /// Does not account for read depth contribution.
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
        self.depth_contrib * (self.atomic_depth_lik_diff(w1, c1) + self.atomic_depth_lik_diff(w2, c2)
            + self.atomic_depth_lik_diff(w3, c3) + self.atomic_depth_lik_diff(w4, c4))
    }

    /// Recalculates total model likelihood, and returns separately
    /// - ln-probability of read depth across all windows,
    /// - ln-probability of of all read alignments.
    pub fn recalc_likelihood(&mut self) -> (f64, f64) {
        let depth_lik = self.depth_distrs.iter()
            .zip(&self.depth)
            .map(|(distr, &depth)| distr.ln_pmf(depth))
            .sum::<f64>()
            * self.depth_contrib;
        let aln_lik = self.read_ixs.iter().zip(&self.read_assgn)
            .map(|(&start_ix, &assgn)| self.read_windows[start_ix + assgn as usize].ln_prob())
            .sum();
        self.likelihood = depth_lik + aln_lik;
        (depth_lik, aln_lik)
    }

    /// Returns all information about windows.
    pub fn contig_windows(&self) -> &MultiContigWindows {
        &self.contig_windows
    }

    /// Read depth contribution, relative to read alignment contribution.
    pub fn depth_contrib(&self) -> f64 {
        self.depth_contrib
    }

    /// Returns read depth distribution for the window.
    /// WARN: Need to account for `self.depth_contrib()`.
    pub fn depth_distr(&self, window: usize) -> &DistrBox {
        &self.depth_distrs[window]
    }

    pub(crate) const DEPTH_CSV_HEADER: &'static str = "contig\twindow\tdepth\tlik";

    /// Write read depth to a CSV file in the following format (tab-separated):
    /// General lines:  `prefix  contig(1|2)  window  depth     depth_lik`.
    /// Last line:      `prefix  summary      key=value key=value ...`. (key-value pairs separated by space).
    pub fn write_depth(&self, f: &mut impl io::Write, prefix: &str) -> io::Result<()> {
        // log10-likelihood of read depth.
        let mut sum_depth_lik = 0.0;
        for i in 0..self.contig_windows.n_contigs() {
            let wshift = self.contig_windows.get_wshift(i) as usize;
            let wshift_end = self.contig_windows.get_wshift(i + 1) as usize;
            for w in wshift..wshift_end {
                let depth = self.depth[w];
                let log10_prob = Ln::to_log10(self.depth_distrs[w].ln_pmf(depth));
                writeln!(f, "{}\t{}\t{}\t{}\t{:.3}", prefix, i + 1, w - wshift + 1, depth, log10_prob)?;
                sum_depth_lik += log10_prob;
            }
        }

        sum_depth_lik *= self.depth_contrib;
        let total_lik = Ln::to_log10(self.likelihood);
        writeln!(f, "{}\tsummary\treads={} unmapped={} boundary={} aln_lik={:.5} depth_lik={:.5} lik={:.5}", prefix,
            self.read_assgn.len(), self.depth[UNMAPPED_WINDOW as usize], self.depth[BOUNDARY_WINDOW as usize],
            total_lik - sum_depth_lik, sum_depth_lik, total_lik)
    }
}
