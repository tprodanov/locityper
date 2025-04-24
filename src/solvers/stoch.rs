use std::{
    fmt,
    cmp::{min, max},
};
use rand::{
    Rng,
    seq::IndexedRandom,
};
use crate::{
    ext::{
        rand::XoshiroRng,
        fmt::PrettyUsize,
    },
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
pub struct Greedy {
    /// If true, start with best read alignments, otherwise, assign initial alignments randomly.
    best_start: bool,
    /// Number of read-pairs, examined per iteration.
    sample_size: usize,
    /// Number of iteration without improvement, after which the solver stops.
    plato_size: usize,
}

impl Default for Greedy {
    fn default() -> Self {
        Self {
            best_start: true,
            sample_size: 10,
            plato_size: 100,
        }
    }
}

impl Greedy {
    pub fn set_start(&mut self, s: &str) -> Result<(), super::ParamErr> {
        self.best_start = match s {
            "b" | "best" => true,
            "r" | "rand" | "random" => false,
            _ => return Err(super::ParamErr::Invalid(format!("Invalid start value {}", s))),
        };
        Ok(())
    }

    pub fn set_sample_size(&mut self, sample_size: usize) -> Result<(), super::ParamErr> {
        if sample_size == 0 {
            Err(super::ParamErr::Invalid("Sample size must be positive".to_string()))
        } else {
            self.sample_size = sample_size;
            Ok(())
        }
    }

    pub fn set_plato_size(&mut self, plato_size: usize) {
        self.plato_size = plato_size;
    }
}

impl Solver for Greedy {
    /// Stochastic greedy algorithm to find the best read assignment.
    fn solve_nontrivial<'a>(
        &self,
        gt_alns: &'a GenotypeAlignments,
        rng: &mut XoshiroRng,
    ) -> crate::Result<ReadAssignment<'a>>
    {
        let non_trivial_reads = gt_alns.non_trivial_reads();
        let sample_size = min(self.sample_size, non_trivial_reads.len());
        let mut assignments = if self.best_start {
            // 0 - index of the best alignment.
            ReadAssignment::new(gt_alns, |_| 0)
        } else {
            ReadAssignment::new(gt_alns, |alns| rng.random_range(0..alns.len()))
        };
        let min_diff = minimum_allowed_diff(max_abs_random(&assignments, rng, INIT_ITER));

        let mut curr_plato = 0;
        let max_iter = max(100_000, self.plato_size * 100);
        for _ in 0..max_iter {
            let mut best_target = None;
            let mut best_improv = min_diff;
            for &rp in non_trivial_reads.choose_multiple(rng, sample_size) {
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

    fn describe_params(&self) -> String {
        format!("x0={},s={},p={}", if self.best_start { "best" } else { "random" },
            PrettyUsize(self.sample_size), PrettyUsize(self.plato_size))
    }
}

impl super::SetParams for Greedy {
    /// Sets solver parameters.
    fn set_param(&mut self, key: &str, val: &str) -> Result<(), super::ParamErr> {
        match &key.to_lowercase() as &str {
            "x0" | "start" => self.set_start(val)?,
            "s" | "sample" => self.set_sample_size(val.parse::<PrettyUsize>()?.0)?,
            "p" | "plato" => self.set_plato_size(val.parse::<PrettyUsize>()?.0),
            _ => return Err(super::ParamErr::Unknown),
        }
        Ok(())
    }
}

impl fmt::Display for Greedy {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Stoch.Greedy")
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
    pub fn set_init_prob(&mut self, init_prob: f64) -> Result<(), super::ParamErr> {
        if init_prob > 0.0 && init_prob <= 1.0 {
            self.init_prob = init_prob;
            Ok(())
        } else {
            Err(super::ParamErr::Invalid(format!("Initial probability ({}) must be within (0, 1]", init_prob)))
        }
    }

    pub fn set_anneal_steps(&mut self, anneal_steps: usize) -> Result<(), super::ParamErr> {
        if anneal_steps > 0 {
            self.anneal_steps = anneal_steps;
            Ok(())
        } else {
            Err(super::ParamErr::Invalid(format!("Number of annealing steps ({}) must be positive", anneal_steps)))
        }
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
        let mut assignments = ReadAssignment::new(gt_alns, |alns| rng.random_range(0..alns.len()));
        let max_abs = max_abs_random(&assignments, rng, INIT_ITER);
        let min_diff = minimum_allowed_diff(max_abs);
        // Solution to equation   `init_prob = exp(-max_abs / start_temp)`.
        let start_temp = (-max_abs / self.init_prob.ln()).max(1e-5);
        let temp_step = start_temp / self.anneal_steps as f64;
        let mut curr_plato = 0;

        for i in (1..=self.anneal_steps).rev() {
            let target = ReassignmentTarget::random(&assignments, rng);
            let diff = assignments.calculate_improvement(&target) - min_diff;
            if diff >= 0.0 || rng.random::<f64>() <= (diff / (temp_step * i as f64)).exp() {
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

    fn describe_params(&self) -> String {
        format!("n={},p={},P={:.3}", self.init_prob, self.anneal_steps, self.plato_size)
    }
}

impl super::SetParams for SimAnneal {
    /// Sets solver parameters.
    fn set_param(&mut self, key: &str, val: &str) -> Result<(), super::ParamErr> {
        match &key.to_lowercase() as &str {
            "n" | "steps" => self.set_anneal_steps(val.parse::<PrettyUsize>()?.0)?,
            "p" | "plato" => self.set_plato_size(val.parse::<PrettyUsize>()?.0),
            "P" | "prob" | "init-prob" => self.set_init_prob(val.parse()?)?,
            _ => return Err(super::ParamErr::Unknown),
        }
        Ok(())
    }
}

impl fmt::Display for SimAnneal {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Sim.Anneal")
    }
}
