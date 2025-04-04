use std::{
    io,
};
use rand::Rng;
use crate::{
    math::Ln,
    seq::contigs::Genotype,
};
use super::{
    locs::AllAlignments,
    windows::{UNMAPPED_WINDOW, BOUNDARY_WINDOW, REG_WINDOW_SHIFT, ReadGtAlns, ContigInfos, GenotypeWindows},
    distr_cache::{DistrCache, WindowDistr},
};

/// All read alignments to a specific genotype.
pub struct GenotypeAlignments {
    /// Contigs subset, to which reads the reads are assigned,
    /// plus the number of windows (and shifts) at each window.
    gt_windows: GenotypeWindows,
    /// Read depth distribution at each window (length: `gt_windows.total_windows()`).
    depth_distrs: Vec<WindowDistr>,

    /// Read depth and read alignments contributions (sum to 2).
    depth_contrib: f64,
    aln_contrib: f64,

    /// A vector of possible read alignments.
    alns: Vec<ReadGtAlns>,
    /// Store the start of the read-pair alignments in the `alns` vector.
    /// Length: `n_reads + 1`, first value = `0`, last value = `alns.len()`.
    ///
    /// For read-pair `i`, possible read locations are `alns[read_ixs[i]..read_ixs[i + 1]]`.
    read_ixs: Vec<usize>,
    /// Indices of reads that have > 1 possible read location (length <= n_reads).
    non_trivial_reads: Vec<usize>,
}

impl GenotypeAlignments {
    /// Creates an instance that stores read assignments to a given genotype.
    /// Read assignment itself is not stored, call `init_assignments()` to start.
    pub fn new(
        genotype: Genotype,
        all_contig_infos: &ContigInfos,
        all_alns: &AllAlignments,
        params: &super::Params,
    ) -> Self
    {
        let gt_windows = GenotypeWindows::new(genotype, all_contig_infos);
        let mut ix = 0;
        let mut alns = Vec::new();
        let mut read_ixs = vec![ix];
        let mut non_trivial_reads = Vec::new();
        for (rp, paired_alns) in all_alns.reads().iter().enumerate() {
            let nw = gt_windows.extend_read_gt_alns(paired_alns, &mut alns, params.prob_diff);
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
        // 1 because `read_ixs` includes is one more than the number of reads.
        // This warning should not appear.
        if read_ixs.len() <= 1 {
            log::warn!("Zero reads can be used for {}", gt_windows.genotype());
        }

        const _: () = assert!(UNMAPPED_WINDOW == 0 && BOUNDARY_WINDOW == 1 && REG_WINDOW_SHIFT == 2,
            "Constants have changed!");
        let total_windows = gt_windows.total_windows() as usize;
        let mut depth_distrs = Vec::with_capacity(total_windows);
        for _ in 0..REG_WINDOW_SHIFT {
            depth_distrs.push(WindowDistr::TRIVIAL);
        }

        Self {
            aln_contrib: 1.0 - params.lik_skew,
            depth_contrib: 1.0 + params.lik_skew,
            non_trivial_reads, gt_windows, alns, depth_distrs, read_ixs,
        }
    }

    /// Returns slice with indices of all read pairs with at least two alignment locations.
    pub(crate) fn non_trivial_reads(&self) -> &[usize] {
        &self.non_trivial_reads
    }

    /// Returns true if there are no non-trivial reads, therefore there exists only one possible read assignment.
    pub fn trivial(&self) -> bool {
        self.non_trivial_reads.is_empty()
    }

    /// Returns the total number of reads.
    pub fn total_reads(&self) -> usize {
        self.read_ixs.len() - 1
    }

    /// Returns possible read locations for the read pair `rp`.
    pub fn possible_read_alns(&self, rp: usize) -> &[ReadGtAlns] {
        &self.alns[self.read_ixs[rp]..self.read_ixs[rp + 1]]
    }

    /// Returns the total possible number of read alignments accross all read pairs.
    pub fn total_possible_alns(&self) -> usize {
        self.alns.len()
    }

    /// Returns all information about windows.
    pub fn gt_windows(&self) -> &GenotypeWindows {
        &self.gt_windows
    }

    /// Returns the total number of windows across the genotype.
    pub fn total_windows(&self) -> usize {
        self.depth_distrs.len()
    }

    /// Read depth and read alignment contributions to the total likelihood.
    pub fn contributions(&self) -> (f64, f64) {
        (self.depth_contrib, self.aln_contrib)
    }

    /// Randomly tweak read centers and window weights, and define read depth distributions later.
    pub fn apply_tweak(
        &mut self,
        rng: &mut impl Rng,
        distr_cache: &DistrCache,
        mean_sum_weight: f64,
        params: &super::Params,
    ) {
        let tweak = params.tweak.unwrap();
        if tweak == 0 {
            self.alns.iter_mut().for_each(|rw| rw.define_windows_determ(&self.gt_windows));
        } else {
            self.alns.iter_mut().for_each(|rw| rw.define_windows_random(&self.gt_windows, tweak, rng));
        }

        self.depth_distrs.truncate(REG_WINDOW_SHIFT as usize);
        for contig_info in self.gt_windows.contig_infos() {
            for (wstart, wend) in contig_info.generate_windows(tweak, rng) {
                let chars = contig_info.window_characteristics(wstart, wend);
                if chars.weight < params.min_weight {
                    self.depth_distrs.push(WindowDistr::TRIVIAL);
                } else {
                    self.depth_distrs.push(distr_cache.get_distribution(chars.gc_content, chars.weight));
                }
            }
        }
        if params.depth_norm_power != 0.0 {
            let sum_weight = self.depth_distrs.iter().map(WindowDistr::weight).sum::<f64>();
            let norm_fct = (mean_sum_weight / sum_weight).powf(params.depth_norm_power).clamp(0.95, 1.0526);
            self.depth_contrib = (params.lik_skew + 1.0) * norm_fct;
        }
    }

    /// Returns read depth distribution for the window.
    /// WARN: Need to account for `self.depth_contrib()`.
    pub fn depth_distr(&self, window: usize) -> &WindowDistr {
        &self.depth_distrs[window]
    }

    /// Returns maximum achievable alignment likelihood (without read depth component).
    pub fn max_aln_lik(&self) -> f64 {
        // This works as read alignments are sorted by likelihood for each read pair.
        self.read_ixs[..self.read_ixs.len() - 1].iter().map(|&i| self.alns[i].ln_prob()).sum::<f64>()
    }

    /// Returns vector that stores how many times each read assignment was selected.
    pub fn create_counts(&self) -> Vec<u16> {
        vec![0; self.total_possible_alns()]
    }
}

/// Read assignment to a vector of contigs and the corresponding likelihoods.
pub struct ReadAssignment<'a> {
    /// Genotype and read alignments to it.
    parent: &'a GenotypeAlignments,
    /// Store current read assignment (length: n_reads).
    read_assgn: Vec<u16>,
    /// Current read depth at each window (concatenated across different contigs).
    /// Length: `gt_windows.total_windows()`.
    depth: Vec<u32>,

    /// Sum alignment ln-probability.
    aln_lik: f64,
    /// Sum read depth ln-probability (without read depth contribution factor).
    depth_lik: f64,
}

impl<'a> ReadAssignment<'a> {
    /// Creates an instance that stores read assignments to given genotype.
    ///
    /// Read assignments are initialized with `select_init` function, which provides initial location index
    /// for all non-trivial reads/read pairs.
    pub fn new<F>(parent: &'a GenotypeAlignments, mut select_init: F) -> Self
    where F: FnMut(&[ReadGtAlns]) -> usize,
    {
        Self::try_new::<_, ()>(parent, |locs| Ok(select_init(locs))).unwrap()
    }

    /// Same as `new`, but initialization function may produce errors.
    pub fn try_new<F, E>(parent: &'a GenotypeAlignments, mut select_init: F) -> Result<Self, E>
    where F: FnMut(&[ReadGtAlns]) -> Result<usize, E>,
    {
        let mut depth = vec![0; parent.total_windows()];
        // TODO: Replace with `array_windows`, once stable.
        let mut i = 0;
        let read_assgn = parent.read_ixs[1..].iter().map(|&j| {
            let m = j - i;
            let mut assgn = 0;
            if m > 1 {
                assgn = select_init(&parent.alns[i..j])?;
                assert!(assgn < m, "Impossible read assignment: {} ({} locations)", assgn, m);
            }
            for w in parent.alns[i + assgn].windows().into_iter() {
                depth[w as usize] += 1;
            }
            i = j;
            Ok(assgn as u16)
        }).collect::<Result<Vec<_>, E>>()?;

        let mut assgn = Self {
            parent, read_assgn, depth,
            aln_lik: f64::NEG_INFINITY,
            depth_lik: f64::NEG_INFINITY,
        };
        assgn.recalc_likelihood();
        Ok(assgn)
    }

    // /// Returns the genotype.
    // pub fn genotype(&self) -> &Genotype {
    //     self.parent.genotype()
    // }

    /// Returns the total likelihood of the read assignment.
    /// Includes contribution factors.
    pub fn likelihood(&self) -> f64 {
        self.parent.depth_contrib * self.depth_lik + self.parent.aln_contrib * self.aln_lik
    }

    /// Assuming that read depth in `window` will change by `depth_change`,
    /// calculates the difference between the new and the old ln-probabilities.
    /// Positive value means that the likelihood will improve.
    /// Does not actually update the read depth.
    /// Does not account for read depth contribution.
    fn atomic_depth_lik_diff(&self, window: u32, depth_change: i32) -> f64 {
        let w = window as usize;
        if depth_change == 0 {
            0.0
        } else {
            let distr = &self.parent.depth_distrs[w];
            let old_depth = self.depth[w];
            let new_depth = old_depth.checked_add_signed(depth_change).expect("Read depth became negative");
            distr.ln_prob(new_depth) - distr.ln_prob(old_depth)
        }
    }

    /// Calculate likelihood difference on four windows (some of them can be equal to others).
    /// Depth at w1 and w2 is decreased by one, depth at w3 and w4 is increased by one.
    /// Does not account for read depth contribution.
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

    /// Finds the best improvement, which can be achieved by moving one read pair.
    pub fn best_read_improvement(&self, read_pair: usize) -> (ReassignmentTarget, f64) {
        let start_ix = self.parent.read_ixs[read_pair];
        let end_ix = self.parent.read_ixs[read_pair + 1];
        assert!(start_ix + 1 < end_ix, "Read has only one possible assignment");

        let old_assgn = self.read_assgn[read_pair];
        let old_ix = start_ix + usize::from(old_assgn);
        let old_aln = &self.parent.alns[old_ix];
        let [w1, w2] = old_aln.windows();

        let mut best_i = 0;
        let mut best_improv = f64::NEG_INFINITY;
        let rel_contrib = self.parent.depth_contrib / self.parent.aln_contrib;
        for (i, aln) in self.parent.alns[start_ix..end_ix].iter().enumerate() {
            if i != usize::from(old_assgn) {
                let [w3, w4] = aln.windows();
                let improv = aln.ln_prob() + rel_contrib * self.depth_lik_diff(w1, w2, w3, w4);
                if improv > best_improv {
                    best_improv = improv;
                    best_i = i;
                }
            }
        }
        let improv = self.parent.aln_contrib * (best_improv - old_aln.ln_prob());
        (ReassignmentTarget {
            read_pair,
            new_assgn: best_i as u16,
            old_ix,
            new_ix: start_ix + best_i,
        }, improv)
    }

    /// Calculates improvement, achieved by reassigning read pair `rp` to a new location.
    /// Does not actually performs reassignment.
    pub fn calculate_improvement(&self, target: &ReassignmentTarget) -> f64 {
        let old_aln = &self.parent.alns[target.old_ix];
        let new_aln = &self.parent.alns[target.new_ix];
        let [w1, w2] = old_aln.windows();
        let [w3, w4] = new_aln.windows();
        self.parent.depth_contrib * self.depth_lik_diff(w1, w2, w3, w4)
            + self.parent.aln_contrib * (new_aln.ln_prob() - old_aln.ln_prob())
    }

    /// Reassigns read pair to a new location.
    pub fn reassign(&mut self, target: &ReassignmentTarget) {
        let old_aln = &self.parent.alns[target.old_ix];
        let new_aln = &self.parent.alns[target.new_ix];
        let [w1, w2] = old_aln.windows();
        let [w3, w4] = new_aln.windows();
        self.depth_lik += self.depth_lik_diff(w1, w2, w3, w4);
        self.aln_lik += new_aln.ln_prob() - old_aln.ln_prob();
        self.depth[w3 as usize] += 1;
        self.depth[w4 as usize] += 1;
        self.depth[w1 as usize] -= 1;
        self.depth[w2 as usize] -= 1;
        self.read_assgn[target.read_pair] = target.new_assgn;
    }

    /// Recalculates total model likelihood.
    pub fn recalc_likelihood(&mut self) {
        self.depth_lik = self.parent.depth_distrs.iter()
            .zip(&self.depth)
            .map(|(distr, &depth)| distr.ln_prob(depth))
            .sum::<f64>();
        self.aln_lik = self.parent.read_ixs.iter().zip(&self.read_assgn)
            .map(|(&start_ix, &assgn)| self.parent.alns[start_ix + assgn as usize].ln_prob())
            .sum();
    }

    pub(crate) const DEPTH_CSV_HEADER: &'static str = "contig\twindow\tweight\tdepth\tlik";

    /// Write read depth to a CSV file.
    pub fn write_depth(&self, f: &mut impl io::Write, prefix: &str) -> io::Result<()> {
        let gt_windows = &self.parent.gt_windows;
        for i in 0..gt_windows.genotype().ploidy() {
            let wshift = gt_windows.get_wshift(i) as usize;
            let wshift_end = gt_windows.get_wshift(i + 1) as usize;
            for w in wshift..wshift_end {
                let depth = self.depth[w];
                let log10_prob = Ln::to_log10(self.parent.depth_distrs[w].ln_prob(depth));
                writeln!(f, "{}{}\t{}\t{:.4}\t{}\t{:.3}", prefix, i + 1, w - wshift + 1,
                    self.parent.depth_distrs[w].weight(), depth, log10_prob)?;
            }
        }
        Ok(())
    }

    pub fn update_counts(&self, counts: &mut [u16]) {
        for (&read_ix, &read_assgn) in self.parent.read_ixs.iter().zip(&self.read_assgn) {
            counts[read_ix + usize::from(read_assgn)] += 1;
        }
    }

    // /// Calculates the fraction of windows with non zero read depth.
    // fn good_depth_fraction(&self) -> f64 {
    //     // log(0.01).
    //     const PROB_THRESH: f64 = -4.605170185988091;
    //     let mut sum_weight = 0.0;
    //     let mut good_depth = 0.0;
    //     for (distr, &depth) in self.parent.depth_distrs.iter().zip(&self.depth) {
    //         if let Some(inner_distr) = distr.inner() {
    //             sum_weight += distr.weight();
    //             good_depth += if inner_distr.ln_pmf(depth) >= PROB_THRESH { distr.weight() } else { 0.0 };
    //         }
    //     }
    //     if sum_weight != 0.0 { good_depth / sum_weight } else { f64::NAN }
    // }

    // /// Fraction of reads where both read ends are mapped and have good edit distance.
    // fn good_alns_fraction(&self, all_alns: &AllAlignments, is_paired_end: bool) -> f64 {
    //     let is_single_end = !is_paired_end;
    //     let mut sum_weight = 0.0;
    //     let mut good_alns = 0.0;
    //     for (groupped_alns, (start_ix, &assgn)) in all_alns.reads().iter().zip(
    //             self.parent.read_ixs.iter().zip(&self.read_assgn)) {
    //         sum_weight += groupped_alns.weight();
    //         if let Some(sel_aln) = self.parent.alns[start_ix + usize::from(assgn)].parent() {
    //             // Both read ends are mapped and have good edit distance.
    //             let good_aln = sel_aln.ix1().map(|i| groupped_alns.ith_aln(i).has_good_dist()).unwrap_or(false)
    //                 && sel_aln.ix2().map(|i| groupped_alns.ith_aln(i).has_good_dist()).unwrap_or(is_single_end);
    //                 good_alns += if good_aln { groupped_alns.weight() } else { 0.0 };
    //         }
    //     }
    //     if sum_weight != 0.0 { good_alns / sum_weight } else { f64::NAN }
    // }

    pub(crate) const SUMMARY_HEADER: &'static str = "total_reads\tunmapped\tout_of_bounds\t\
        aln_lik\tdepth_lik\tlik";

    pub fn summarize(
        &self,
        writer: &mut impl io::Write,
        prefix: &str,
    ) -> io::Result<()>
    {
        writeln!(writer, "{}{}\t{}\t{}\t{:.7}\t{:.7}\t{:.7}",
            prefix, self.read_assgn.len(), self.depth[UNMAPPED_WINDOW as usize], self.depth[BOUNDARY_WINDOW as usize],
            Ln::to_log10(self.aln_lik), Ln::to_log10(self.depth_lik), Ln::to_log10(self.likelihood()))
    }
}

/// Reassignment target: one read pair, its old location and its new location.
#[derive(Debug)]
pub struct ReassignmentTarget {
    read_pair: usize,
    new_assgn: u16,
    old_ix: usize,
    new_ix: usize,
}

impl ReassignmentTarget {
    pub fn new(assgn: &ReadAssignment, read_pair: usize, new_assgn: u16) -> Self {
        let old_assgn = assgn.read_assgn[read_pair];
        assert_ne!(old_assgn, new_assgn, "Read pair #{}: cannot change to the same assignment", read_pair);

        let start_ix = assgn.parent.read_ixs[read_pair];
        let end_ix = assgn.parent.read_ixs[read_pair + 1];
        let old_ix = start_ix + usize::from(old_assgn);
        let new_ix = start_ix + usize::from(new_assgn);
        assert!(new_ix < end_ix, "Invalid read assignment for read pair #{}", read_pair);
        Self { read_pair, new_assgn, old_ix, new_ix }
    }

    /// Creates a random reassignment target: selects one random read pair, and its possible random new location.
    pub fn random(assgn: &ReadAssignment, rng: &mut impl Rng) -> Self {
        let read_pair = assgn.parent.non_trivial_reads[rng.random_range(0..assgn.parent.non_trivial_reads.len())];
        let start_ix = assgn.parent.read_ixs[read_pair];
        let end_ix = assgn.parent.read_ixs[read_pair + 1];
        let total_assgns = end_ix - start_ix;
        debug_assert!(total_assgns > 1, "Read pair #{} has less than 2 possible assignments", read_pair);

        let old_assgn = assgn.read_assgn[read_pair];
        let new_assgn = if total_assgns == 2 {
            1 - old_assgn
        } else {
            let i = rng.random_range(1..total_assgns as u16);
            if i <= old_assgn { i - 1 } else { i }
        };

        Self {
            read_pair, new_assgn,
            old_ix: start_ix + usize::from(old_assgn),
            new_ix: start_ix + usize::from(new_assgn),
        }
    }

    /// Returns read pair index and new assignment index
    pub fn get(&self) -> (usize, u16) {
        (self.read_pair, self.new_assgn)
    }
}
