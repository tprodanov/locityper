use std::cmp::min;
use rand::{
    Rng,
    SeedableRng,
    rngs::SmallRng,
    seq::SliceRandom,
};
use crate::{
    algo::vec_ext::F64Ext,
    model::assgn::ReadAssignment,
    solvers::{SolverBuilder, Solver},
};

/// Builder, that constructs `GreedySolver`.
#[derive(Clone)]
pub struct GreedyBuilder {
    seed: Option<u64>,
    sample_size: usize,
    plato_iters: u32,
}

impl Default for GreedyBuilder {
    /// Creates default GreedyBuilder: seed is not set, sample size is 100, and there can be up to 5 plato iterations.
    fn default() -> Self {
        Self {
            seed: None,
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
}

impl SolverBuilder for GreedyBuilder {
    type S<'a> = GreedySolver<'a>;

    /// Sets seed.
    fn set_seed(&mut self, seed: u64) -> &mut Self {
        self.seed = Some(seed);
        self
    }

    /// Builds the solver.
    fn build<'a>(&self, assignments: ReadAssignment<'a>) -> Self::S<'a> {
        GreedySolver {
            sample_size: min(self.sample_size, assignments.non_trivial_reads().len()),
            rng: SmallRng::seed_from_u64(self.seed.expect("GreedySolver: seed is not set")),
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
pub struct GreedySolver<'a> {
    assignments: ReadAssignment<'a>,
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

impl<'a> GreedySolver<'a> {
    /// Creates GreedySolver with default parameters (seed `GreedyBuilder::default`).
    pub fn new(assignments: ReadAssignment<'a>, seed: u64) -> Self {
        Self::builder().set_seed(seed).build(assignments)
    }

    /// Creates GreedyBuilder.
    pub fn builder() -> GreedyBuilder {
        GreedyBuilder::default()
    }
}

impl<'a> Solver<'a> for GreedySolver<'a> {
    /// Returns false, the method is not determenistic.
    fn is_determenistic() -> bool { false }

    /// Initialize read alignments. Returns current likelihood.
    fn initialize(&mut self) -> f64 {
        // self.assignments.init_assignments(solvers::init_best)
        // self.assignments.init_assignments(|_| 0)
        self.assignments.init_assignments(|alns| self.rng.gen_range(0..alns.len()))
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