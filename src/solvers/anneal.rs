use std::cmp::max;
use rand::{
    Rng,
    SeedableRng,
    rngs::SmallRng,
};
use crate::{
    model::assgn::ReadAssignment,
    solvers::{self, SolverBuilder, Solver},
};

/// Builder, that constructs `SimulatedAnnealing`.
#[derive(Clone)]
pub struct AnnealingBuilder {
    seed: Option<u64>,
    /// Overall, annealing performs exactly `steps` iterations.
    steps: u32,
    /// On each iteration, it randomly checks `max_tries` elements, until one fits.
    max_tries: u32,
    /// Initialize temperature constant in such way, that initially
    /// an average negative improvement would have `init_prob` chance to pass.
    init_prob: f64,
}

impl Default for AnnealingBuilder {
    /// Creates default AnnealingBuilder: seed is not set, steps 10000, and max_tries (per step) is 50.
    fn default() -> Self {
        Self {
            seed: None,
            steps: 10000,
            max_tries: 50,
            init_prob: 0.1,
        }
    }
}

impl AnnealingBuilder {
    pub fn set_steps(&mut self, steps: u32) -> &mut Self {
        self.steps = steps;
        self
    }

    pub fn set_max_tries(&mut self, max_tries: u32) -> &mut Self {
        self.max_tries = max_tries;
        self
    }
}

impl SolverBuilder for AnnealingBuilder {
    type S<'a> = SimulatedAnnealing<'a>;

    /// Sets seed.
    fn set_seed(&mut self, seed: u64) -> &mut Self {
        self.seed = Some(seed);
        self
    }

    /// Builds the solver.
    fn build<'a>(&self, assignments: ReadAssignment<'a>) -> Self::S<'a> {
        SimulatedAnnealing {
            assignments,
            rng: SmallRng::seed_from_u64(self.seed.expect("GreedySolver: seed is not set")),
            curr_step: 0,
            steps: self.steps,
            max_tries: self.max_tries,
            init_prob: self.init_prob,
            coef: f64::NAN,
        }
    }
}

/// In addition to seed, has two parameters: `sample_size` and `plato_iters`.
///
/// In one step, the solver examines `sample_size` read pairs, and selects the best read pair to switch location.
/// If no improvement was made for `plato_iters` iterations, the solver stops.
pub struct SimulatedAnnealing<'a> {
    assignments: ReadAssignment<'a>,
    rng: SmallRng,
    curr_step: u32,
    /// Total number of iterations.
    steps: u32,
    /// Maximal number of neighbours, visited per step.
    max_tries: u32,
    init_prob: f64,
    coef: f64,
}

impl<'a> SimulatedAnnealing<'a> {
    /// Creates GreedySolver with default parameters (seed `AnnealingBuilder::default`).
    pub fn new(assignments: ReadAssignment<'a>, seed: u64) -> Self {
        Self::builder().set_seed(seed).build(assignments)
    }

    /// Creates AnnealingBuilder.
    pub fn builder() -> AnnealingBuilder {
        AnnealingBuilder::default()
    }
}

impl<'a> Solver<'a> for SimulatedAnnealing<'a> {
    /// Returns false, the method is not determenistic.
    fn is_determenistic() -> bool { false }

    /// Initialize read alignments. Returns current likelihood.
    fn initialize(&mut self) -> f64 {
        let lik = self.assignments.init_assignments(solvers::init_best);
        let mut neg_sum = 0.0;
        let mut neg_count = 0;
        for _ in 0..max(100, self.max_tries) {
            let diff = self.assignments.random_reassignment(&mut self.rng).2;
            if diff < 0.0 {
                neg_sum += diff;
                neg_count += 1;
            }
        }
        self.coef = if neg_count == 0 {
            1.0
        } else {
            self.init_prob.ln() * neg_count as f64 / neg_sum
        };
        lik
    }

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> f64 {
        if self.curr_step >= self.steps {
            log::warn!("SimulatedAnnealing is finished, but `step` is called one more time.")
        }
        let temp = 1.0 - (self.curr_step as f64 + 1.0) / self.steps as f64;
        self.curr_step += 1;
        for _ in 0..self.max_tries {
            let (rp, new_assign, improv) = self.assignments.random_reassignment(&mut self.rng);
            if improv > 0.0 || (temp > 0.0 && self.rng.gen_range(0.0..=1.0) <= (self.coef * improv / temp).exp()) {
                return self.assignments.reassign(rp, new_assign);
            }
        }
        0.0
    }

    /// Returns true if the solver is finished.
    fn is_finished(&self) -> bool {
        self.curr_step >= self.steps
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