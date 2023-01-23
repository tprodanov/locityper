use std::cmp::min;
use rand::{
    SeedableRng,
    rngs::SmallRng,
    seq::SliceRandom,
};
use crate::{
    algo::vec_ext::*,
    model::assgn::ReadAssignment,
    solvers::{self, Solver},
};

pub struct GreedySolver<'a> {
    assignments: ReadAssignment<'a>,
    rng: SmallRng,
    is_finished: bool,
    /// Number of read-pairs, examined per iteration.
    sample_size: usize,
    buffer: Vec<f64>,
}

impl<'a> GreedySolver<'a> {
    pub fn new(assignments: ReadAssignment<'a>, seed: u64, sample_size: usize) -> Self {
        Self {
            sample_size: min(sample_size, assignments.non_trivial_reads().len()),
            assignments,
            rng: SmallRng::seed_from_u64(seed),
            is_finished: false,
            buffer: Vec::with_capacity(16),
        }
    }
}

impl<'a> Solver<'a> for GreedySolver<'a> {
    /// Initialize read alignments. Returns current likelihood.
    fn initialize(&mut self) -> f64 {
        self.assignments.init_assignments(solvers::init_best)
    }

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> f64 {
        if self.is_finished {
            log::warn!("GreedySolver is finished, but `step` is called one more time.")
        }
        let mut best_rp = 0;
        let mut best_assgn = 0;
        let mut best_improv = 0.0;
        for &rp in self.assignments.non_trivial_reads().choose_multiple(&mut self.rng, self.sample_size) {
            self.buffer.clear();
            self.assignments.possible_reassignments(rp, &mut self.buffer);
            let (assgn, improv) = self.buffer.argmax();
            if improv > best_improv {
                best_rp = rp;
                best_assgn = assgn as u16;
                best_improv = improv;
            }
        }
        if best_improv > 0.0 {
            let improv2 = self.assignments.reassign(best_rp, best_assgn);
            debug_assert_eq!(improv2, best_improv, "Unexpected likelihood improvement.");
        } else {
            self.is_finished = true;
        }
        best_improv
    }

    /// Returns true if the solver is finished.
    fn is_finished(&self) -> bool {
        self.is_finished
    }

    /// Return the current read assignments.
    fn current_assignments<'b>(&'b self) -> &'b ReadAssignment<'a> where 'a: 'b {
        &self.assignments
    }

    /// Finish solving, consume the solver and return the read assignments.
    fn finish(self) -> ReadAssignment<'a> {
        self.assignments
    }
}