use std::{fmt, cmp::min};
use rand::{
    Rng,
    SeedableRng,
    rngs::SmallRng,
    seq::SliceRandom,
};
use crate::{
    Error,
    algo::vec_ext::F64Ext,
    model::assgn::ReadAssignment,
};
use super::Solver;

/// Builder, that constructs `GreedySolver`.
#[derive(Clone)]
pub struct GreedyBuilder {
    sample_size: usize,
    plato_iters: u32,
}

impl Default for GreedyBuilder {
    /// Creates default GreedyBuilder: seed is not set, sample size is 100, and there can be up to 5 plato iterations.
    fn default() -> Self {
        Self {
            sample_size: 100,
            plato_iters: 5,
        }
    }
}

impl GreedyBuilder {
    pub fn set_sample_size(&mut self, sample_size: usize) -> &mut Self {
        self.sample_size = sample_size;
        self
    }

    pub fn set_plato_iters(&mut self, plato_iters: u32) -> &mut Self {
        self.plato_iters = plato_iters;
        self
    }

    /// Builds the solver.
    pub fn build(&self, assignments: ReadAssignment) -> GreedySolver {
        GreedySolver {
            sample_size: min(self.sample_size, assignments.non_trivial_reads().len()),
            rng: SmallRng::seed_from_u64(0),
            assignments,
            is_finished: false,
            buffer: Vec::with_capacity(16),
            plato_iters: self.plato_iters,
            curr_plato: 0,
        }
    }
}

/// In addition to seed, has two parameters: `sample_size` and `plato_iters`.
///
/// In one step, the solver examines `sample_size` read pairs, and selects the best read pair to switch location.
/// If no improvement was made for `plato_iters` iterations, the solver stops.
pub struct GreedySolver {
    assignments: ReadAssignment,
    rng: SmallRng,
    is_finished: bool,
    buffer: Vec<f64>,
    /// Number of read-pairs, examined per iteration.
    sample_size: usize,
    /// Number of iteration without improvement, after which the solver stops.
    plato_iters: u32,
    /// Current number of iterations without improvement.
    curr_plato: u32,
}

impl GreedySolver {
    /// Creates GreedySolver with default parameters (seed `GreedyBuilder::default`).
    pub fn new(assignments: ReadAssignment) -> Self {
        GreedyBuilder::default().build(assignments)
    }

    /// Creates GreedyBuilder.
    pub fn builder() -> GreedyBuilder {
        GreedyBuilder::default()
    }
}

impl Solver for GreedySolver {
    fn is_seedable() -> bool { true }

    fn set_seed(&mut self, seed: u64) -> Result<(), Error> {
        self.rng = SmallRng::seed_from_u64(seed);
        Ok(())
    }

    fn reset(&mut self) -> Result<(), Error> {
        self.is_finished = false;
        // self.assignments.init_assignments(solvers::init_best)
        // self.assignments.init_assignments(|_| 0)
        self.assignments.init_assignments(|alns| self.rng.gen_range(0..alns.len()));
        Ok(())
    }

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> Result<(), Error> {
        if self.is_finished {
            log::warn!("GreedySolver is finished, but `step` is called one more time.")
        }
        let mut best_rp = 0;
        let mut best_assgn = 0;
        let mut best_improv = 0.0;
        for &rp in self.assignments.non_trivial_reads().choose_multiple(&mut self.rng, self.sample_size) {
            self.buffer.clear();
            self.assignments.possible_reassignments(rp, &mut self.buffer);
            let (assgn, improv) = F64Ext::argmax(&self.buffer);
            if improv > best_improv {
                best_rp = rp;
                best_assgn = assgn as u16;
                best_improv = improv;
            }
        }
        if best_improv > 0.0 {
            self.curr_plato = 0;
            self.assignments.reassign(best_rp, best_assgn);
        } else {
            self.curr_plato += 1;
            self.is_finished = self.curr_plato > self.plato_iters;
        }
        Ok(())
    }

    /// Returns true if the solver is finished.
    fn is_finished(&self) -> bool {
        self.is_finished
    }

    /// Return the current read assignments.
    fn current_assignments(&self) -> &ReadAssignment {
        &self.assignments
    }

    fn recalculate_likelihood(&mut self) {
        self.assignments.recalc_likelihood();
    }

    fn take(self) -> ReadAssignment {
        self.assignments
    }
}

impl fmt::Display for GreedySolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Greedy({},{})", self.sample_size, self.plato_iters)
    }
}
