use std::{
    fmt,
    cmp::{min, max},
};
use rand::{
    Rng,
    seq::SliceRandom,
};
use crate::{
    err::{error, validate_param},
    ext::rand::XoshiroRng,
    model::assgn::{GenotypeAlignments, ReadAssignment, ReassignmentTarget},
};
use super::Solver;

/// Returns maximum absolute likelihood difference, observed over `count` random reassignments.
fn max_abs_random(assignments: &ReadAssignment, rng: &mut XoshiroRng, count: usize) -> f64 {
    (0..count).map(|_| assignments.calculate_improvement(&ReassignmentTarget::random(&assignments, rng)))
        .fold(0.0, |acc, x| acc.max(x.abs()))
}

const INIT_ITER: usize = 100;

/// Consider likelihood differences between these values random noise caused by float arithmetics.
fn minimum_allowed_diff(max_abs_diff: f64) -> f64 {
    (1e-10 * max_abs_diff).max(1e-14)
}

/// Assigns reads in a greedy way.
///
/// In one step, the solver examines `sample_size` read pairs, and selects the best read pair to switch location.
/// If no improvement was made for `plato_size` iterations, the solver stops.
#[derive(Clone)]
pub struct GreedySolver {
    /// If true, start with best read alignments, otherwise, assign initial alignments randomly.
    best_start: bool,
    /// Number of read-pairs, examined per iteration.
    sample_size: usize,
    /// Number of iteration without improvement, after which the solver stops.
    plato_size: usize,
}

impl Default for GreedySolver {
    fn default() -> Self {
        Self {
            best_start: true,
            sample_size: 10,
            plato_size: 100,
        }
    }
}

impl GreedySolver {
    pub fn set_best_start(&mut self, best_start: bool) {
        self.best_start = best_start;
    }

    pub fn set_sample_size(&mut self, sample_size: usize) -> crate::Result<()> {
        validate_param!(sample_size != 0, "Greedy solver: sample size must be positive");
        self.sample_size = sample_size;
        Ok(())
    }

    pub fn set_plato_size(&mut self, plato_size: usize) {
        self.plato_size = plato_size;
    }
}

impl Solver for GreedySolver {
    /// Stochastic greedy algorithm to find the best read assignment.
    fn solve_nontrivial<'a>(
        &self,
        gt_alns: &'a GenotypeAlignments,
        rng: &mut XoshiroRng,
    ) -> crate::Result<ReadAssignment<'a>>
    {
        let non_trivial_ixs = gt_alns.non_trivial_reads();
        let sample_size = min(self.sample_size, non_trivial_ixs.len());
        let mut assignments = if self.best_start {
            // 0 - index of the best alignment.
            ReadAssignment::new(gt_alns, |_| 0)
        } else {
            ReadAssignment::new(gt_alns, |alns| rng.gen_range(0..alns.len()))
        };
        let min_diff = minimum_allowed_diff(max_abs_random(&assignments, rng, INIT_ITER));

        let mut curr_plato = 0;
        let max_iter = max(100_000, self.plato_size * 100);
        for _ in 0..max_iter {
            let mut best_target = None;
            let mut best_improv = min_diff;
            for &rp in non_trivial_ixs.choose_multiple(rng, sample_size) {
                let (target, improv) = assignments.best_read_improvement(rp);
                if improv > best_improv {
                    best_target = Some(target);
                    best_improv = improv;
                }
            }
            if let Some(target) = best_target {
                curr_plato = 0;
                assignments.reassign(&target);
            } else {
                curr_plato += 1;
                if curr_plato > self.plato_size {
                    break;
                }
            }
        }
        Ok(assignments)
    }
}

impl super::SetParams for GreedySolver {
    /// Sets solver parameters.
    fn set_param(&mut self, key: &str, val: &str) -> crate::Result<()> {
        match &key.to_lowercase() as &str {
            "best" | "beststart" | "best_start" => self.set_best_start(val.parse()
                .map_err(|_| error!(InvalidInput, "Cannot parse '{}={}'", key, val))?),
            "sample" => self.set_sample_size(val.parse()
                .map_err(|_| error!(InvalidInput, "Cannot parse '{}={}'", key, val))?)?,
            "plato" => self.set_plato_size(val.parse()
                .map_err(|_| error!(InvalidInput, "Cannot parse '{}={}'", key, val))?),
            _ => log::error!("Greedy solver: unknown parameter {:?}", key),
        }
        Ok(())
    }
}

impl fmt::Display for GreedySolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Greedy({} start, {} reads/iter, plato {})",
            if self.best_start { "best" } else { "random" }, self.sample_size, self.plato_size)
    }
}

/// Simulated annealing solver.
///
/// Randomly selects direction based on the current temperature, and stops once there are no improvement `plato_size`.
#[derive(Clone)]
pub struct SimAnneal {
    /// Initialize temperature constant in such way, that initially
    /// an -abs(max diff) would have `init_prob` chance to pass.
    init_prob: f64,
    /// How many annealing steps should there be?
    anneal_steps: usize,
    /// After annealing, continue improving likelihood until encounter `plato_size` iterations without improvement.
    plato_size: usize,
}

impl Default for SimAnneal {
    fn default() -> Self {
        Self {
            init_prob: 0.5,
            anneal_steps: 20000,
            plato_size: 10000,
        }
    }
}

impl SimAnneal {
    pub fn set_init_prob(&mut self, init_prob: f64) -> crate::Result<()> {
        validate_param!(init_prob > 0.0 && init_prob <= 1.0,
            "Initial probability ({}) must be within (0, 1]", init_prob);
        self.init_prob = init_prob;
        Ok(())
    }

    pub fn set_anneal_steps(&mut self, anneal_steps: usize) -> crate::Result<()> {
        validate_param!(anneal_steps > 0, "Number of annealing steps ({}) must be positive", anneal_steps);
        self.anneal_steps = anneal_steps;
        Ok(())
    }

    pub fn set_plato_size(&mut self, plato_size: usize) {
        self.plato_size = plato_size;
    }
}

impl Solver for SimAnneal {
    /// Run simulated annealing to find the best read assignment.
    fn solve_nontrivial<'a>(
        &self,
        gt_alns: &'a GenotypeAlignments,
        rng: &mut XoshiroRng,
    ) -> crate::Result<ReadAssignment<'a>>
    {
        // 0 - index of the best alignment.
        // let mut assignments = ReadAssignment::new(gt_alns, |_| 0);
        let mut assignments = ReadAssignment::new(gt_alns, |alns| rng.gen_range(0..alns.len()));
        let max_abs = max_abs_random(&assignments, rng, INIT_ITER);
        let min_diff = minimum_allowed_diff(max_abs);
        // Solution to equation   `init_prob = exp(-max_abs / start_temp)`.
        let start_temp = (-max_abs / self.init_prob.ln()).max(1e-5);
        let temp_step = start_temp / self.anneal_steps as f64;
        let mut curr_plato = 0;

        for i in (1..=self.anneal_steps).rev() {
            let target = ReassignmentTarget::random(&assignments, rng);
            let diff = assignments.calculate_improvement(&target) - min_diff;
            if diff >= 0.0 || rng.gen::<f64>() <= (diff / (temp_step * i as f64)).exp() {
                assignments.reassign(&target);
                curr_plato = 0;
            } else {
                curr_plato += 1;
                if curr_plato >= self.plato_size {
                    break;
                }
            }
        }

        let max_iter = max(100_000, self.plato_size * 100);
        for _ in 0..max_iter {
            if curr_plato >= self.plato_size {
                break;
            }
            let target = ReassignmentTarget::random(&assignments, rng);
            let diff = assignments.calculate_improvement(&target);
            if diff > min_diff {
                assignments.reassign(&target);
                curr_plato = 0;
            } else {
                curr_plato += 1;
            }
        }
        Ok(assignments)
    }
}

impl super::SetParams for SimAnneal {
    /// Sets solver parameters.
    fn set_param(&mut self, key: &str, val: &str) -> crate::Result<()> {
        match &key.to_lowercase() as &str {
            "steps" | "annealsteps" | "anneal_steps" => self.set_anneal_steps(val.parse()
                .map_err(|_| error!(InvalidInput, "Cannot parse '{}={}'", key, val))?)?,
            "init" | "init_prob" | "initprob" => self.set_init_prob(val.parse()
                .map_err(|_| error!(InvalidInput, "Cannot parse '{}={}'", key, val))?)?,
            "plato" => self.set_plato_size(val.parse()
                .map_err(|_| error!(InvalidInput, "Cannot parse '{}={}'", key, val))?),
            _ => log::error!("Sim.Annealing: unknown parameter {:?}", key),
        }
        Ok(())
    }
}

impl fmt::Display for SimAnneal {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SimAnneal(init.prob {}, steps {}, plato {})", self.init_prob, self.anneal_steps, self.plato_size)
    }
}
