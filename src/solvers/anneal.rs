use std::fmt;
use rand::{
    Rng, SeedableRng,
    rngs::SmallRng,
};
use crate::{
    Error,
    model::assgn::ReadAssignment,
};
use super::Solver;

/// Builder, that constructs `SimulatedAnnealing`.
#[derive(Clone)]
pub struct AnnealingBuilder {
    /// Temperature starts at 1, and decreases by `cooling_temp` every step.
    cooling_temp: f64,
    /// Initialize temperature constant in such way, that initially
    /// an average negative improvement would have `init_prob` chance to pass.
    init_prob: f64,
    /// Solver stops if there were no improvements during the last `plato_iters`.
    plato_iters: u32,
}

impl Default for AnnealingBuilder {
    /// Creates default AnnealingBuilder.
    fn default() -> Self {
        Self {
            cooling_temp: 1e-5,
            init_prob: 0.1,
            plato_iters: 20000,
        }
    }
}

impl AnnealingBuilder {
    pub fn set_cooling_temperature(&mut self, cooling_temp: f64) -> &mut Self {
        assert!(cooling_temp > 0.0 && cooling_temp < 1.0,
            "Cooling temperature ({}) must be within (0, 1).", cooling_temp);
        self.cooling_temp = cooling_temp;
        self
    }

    pub fn set_init_probability(&mut self, init_prob: f64) -> &mut Self {
        assert!(init_prob > 0.0 && init_prob < 1.0, "Initial probability ({}) must be within (0, 1).", init_prob);
        self.init_prob = init_prob;
        self
    }

    pub fn set_plato_iters(&mut self, plato: u32) -> &mut Self {
        self.plato_iters = plato;
        self
    }

    /// Builds the solver.
    pub fn build(&self, assignments: ReadAssignment) -> SimulatedAnnealing {
        SimulatedAnnealing {
            assignments,
            rng: SmallRng::seed_from_u64(0),
            builder: self.clone(),
            curr_temp: f64::NAN,
            coeff: f64::NAN,
            curr_plato: 0,
        }
    }
}

/// Simulated annealing solver.
/// On each steps, examines at most `max_tries` neighbours, and selects one of them
/// with a probability adjusted according to the current step. addition to seed, has two parameters: `sample_size` and `plato_iters`.
///
/// In one step, the solver examines `sample_size` read pairs, and selects the best read pair to switch location.
/// If no improvement was made for `plato_iters` iterations, the solver stops.
pub struct SimulatedAnnealing {
    assignments: ReadAssignment,
    rng: SmallRng,
    /// Builder, with which the solver was built.
    builder: AnnealingBuilder,

    /// Current temperature.
    curr_temp: f64,
    /// Coefficient, by which the temperature is multiplied.
    coeff: f64,
    /// Current number of plato iterations.
    curr_plato: u32,
}

impl SimulatedAnnealing {
    /// Creates GreedySolver with default parameters (seed `AnnealingBuilder::default`).
    pub fn new(assignments: ReadAssignment) -> Self {
        AnnealingBuilder::default().build(assignments)
    }

    /// Creates AnnealingBuilder.
    pub fn builder() -> AnnealingBuilder {
        AnnealingBuilder::default()
    }
}

impl Solver for SimulatedAnnealing {
    fn is_seedable() -> bool { true }

    fn set_seed(&mut self, seed: u64) -> Result<(), Error> {
        self.rng = SmallRng::seed_from_u64(seed);
        Ok(())
    }

    fn reset(&mut self) -> Result<(), Error> {
        self.curr_temp = 1.0;
        self.curr_plato = 0;
        self.assignments.init_assignments(super::init_best);

        let mut neg_sum = 0.0;
        let mut neg_count: u32 = 0;
        const INIT_ITERS: u32 = 100;
        for _ in 0..INIT_ITERS {
            let diff = self.assignments.random_reassignment(&mut self.rng).2;
            if diff < 0.0 {
                neg_sum += diff;
                neg_count += 1;
            }
        }
        self.coeff = if neg_count == 0 {
            1.0
        } else {
            self.builder.init_prob.ln() * f64::from(neg_count) / neg_sum
        };
        Ok(())
    }

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> Result<(), Error> {
        let (rp, new_assign, improv) = self.assignments.random_reassignment(&mut self.rng);
        if improv > 0.0 || (self.curr_temp > 0.0
                && self.rng.gen_range(0.0..=1.0) <= (self.coeff * improv / self.curr_temp).exp()) {
            self.assignments.reassign(rp, new_assign);
            self.curr_plato = 0;
        } else {
            self.curr_plato += 1;
        }
        self.curr_temp -= self.builder.cooling_temp;
        Ok(())
    }

    fn is_finished(&self) -> bool {
        self.curr_plato > self.builder.plato_iters
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

impl fmt::Display for SimulatedAnnealing {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SimAnneal({:e},{:e},{})", self.builder.cooling_temp, self.builder.init_prob,
            self.builder.plato_iters)
    }
}