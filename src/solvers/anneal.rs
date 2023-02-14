use std::cmp::max;
use rand::{
    Rng,
    SeedableRng,
    rngs::SmallRng,
};
use crate::model::assgn::ReadAssignment;
use super::Solver;

/// Builder, that constructs `SimulatedAnnealing`.
#[derive(Clone)]
pub struct AnnealingBuilder {
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

    /// Builds the solver.
    pub fn build(&self, assignments: ReadAssignment) -> SimulatedAnnealing {
        SimulatedAnnealing {
            assignments,
            rng: SmallRng::seed_from_u64(0),
            curr_step: 0,
            steps: self.steps,
            max_tries: self.max_tries,
            init_prob: self.init_prob,
            coeff: f64::NAN,
        }
    }
}

/// In addition to seed, has two parameters: `sample_size` and `plato_iters`.
///
/// In one step, the solver examines `sample_size` read pairs, and selects the best read pair to switch location.
/// If no improvement was made for `plato_iters` iterations, the solver stops.
pub struct SimulatedAnnealing {
    assignments: ReadAssignment,
    rng: SmallRng,
    curr_step: u32,
    /// Total number of iterations.
    steps: u32,
    /// Maximal number of neighbours, visited per step.
    max_tries: u32,
    init_prob: f64,
    coeff: f64,
}

impl SimulatedAnnealing {
    /// Creates GreedySolver with default parameters (seed `AnnealingBuilder::default`).
    pub fn default(assignments: ReadAssignment) -> Self {
        AnnealingBuilder::default().build(assignments)
    }

    /// Creates AnnealingBuilder.
    pub fn builder() -> AnnealingBuilder {
        AnnealingBuilder::default()
    }
}

impl Solver for SimulatedAnnealing {
    fn is_seedable() -> bool { true }

    fn set_seed(&mut self, seed: u64) {
        self.rng = SmallRng::seed_from_u64(seed);
    }

    fn initialize(&mut self) {
        self.assignments.init_assignments(super::init_best);
        let mut neg_sum = 0.0;
        let mut neg_count = 0;
        for _ in 0..max(100, self.max_tries) {
            let diff = self.assignments.random_reassignment(&mut self.rng).2;
            if diff < 0.0 {
                neg_sum += diff;
                neg_count += 1;
            }
        }
        self.coeff = if neg_count == 0 {
            1.0
        } else {
            self.init_prob.ln() * neg_count as f64 / neg_sum
        };
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
            if improv > 0.0 || (temp > 0.0 && self.rng.gen_range(0.0..=1.0) <= (self.coeff * improv / temp).exp()) {
                return self.assignments.reassign(rp, new_assign);
            }
        }
        0.0
    }

    fn is_finished(&self) -> bool {
        self.curr_step >= self.steps
    }

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