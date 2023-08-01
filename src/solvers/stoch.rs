use std::{
    fmt,
    cmp::min,
};
use rand::{
    Rng,
    seq::SliceRandom,
};
use crate::{
    err::{Error, validate_param},
    ext::rand::XoshiroRng,
    model::assgn::{GenotypeAlignments, ReadAssignment, ReassignmentTarget},
};
use super::Solver;

/// Assigns reads in a greedy way.
///
/// In one step, the solver examines `sample_size` read pairs, and selects the best read pair to switch location.
/// If no improvement was made for `plato_size` iterations, the solver stops.
#[derive(Clone)]
pub struct GreedySolver {
    /// Number of read-pairs, examined per iteration.
    sample_size: usize,
    /// Number of iteration without improvement, after which the solver stops.
    plato_size: usize,
}

impl Default for GreedySolver {
    fn default() -> Self {
        Self {
            sample_size: 100,
            plato_size: 5,
        }
    }
}

impl GreedySolver {
    pub fn set_sample_size(&mut self, sample_size: usize) -> Result<(), Error> {
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
    ) -> Result<ReadAssignment<'a>, Error>
    {
        // 0 - index of the best alignment.
        let non_trivial_ixs = gt_alns.non_trivial_reads();
        let sample_size = min(self.sample_size, non_trivial_ixs.len());
        let mut assignments = ReadAssignment::new(gt_alns, |_| 0);
        let mut curr_plato = 0;

        while curr_plato <= self.plato_size {
            let mut best_target = None;
            let mut best_improv = 0.0;
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
            }
        }
        Ok(assignments)
    }
}

impl super::SetParams for GreedySolver {
    /// Sets solver parameters.
    fn set_param(&mut self, key: &str, val: &str) -> Result<(), Error> {
        match &key.to_lowercase() as &str {
            "sample" => self.set_sample_size(val.parse()
                .map_err(|_| Error::InvalidInput(format!("Cannot parse '{}={}'", key, val)))?)?,
            "plato" => self.set_plato_size(val.parse()
                .map_err(|_| Error::InvalidInput(format!("Cannot parse '{}={}'", key, val)))?),
            _ => log::error!("Greedy solver: unknown parameter {:?}", key),
        }
        Ok(())
    }
}

impl fmt::Display for GreedySolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Greedy({} reads/iter, plato {})", self.sample_size, self.plato_size)
    }
}

/// Simulated annealing solver.
///
/// Randomly selects direction based on the current temperature, and stops once there are no improvement `plato_size`.
#[derive(Clone)]
pub struct SimAnneal {
    /// Temperature starts at 1, and decreases by `cooling_temp` every step.
    cooling_temp: f64,
    /// Initialize temperature constant in such way, that initially
    /// an average negative improvement would have `init_prob` chance to pass.
    init_prob: f64,
    /// Solver stops if there were no improvements during the last `plato_size` iterations.
    plato_size: usize,
}

impl Default for SimAnneal {
    fn default() -> Self {
        Self {
            cooling_temp: 5e-5,
            init_prob: 0.1,
            plato_size: 10000,
        }
    }
}

impl SimAnneal {
    pub fn set_cooling_temp(&mut self, cooling_temp: f64) -> Result<(), Error> {
        validate_param!(cooling_temp > 0.0 && cooling_temp < 1.0,
            "Cooling temperature ({}) must be within (0, 1)", cooling_temp);
        self.cooling_temp = cooling_temp;
        Ok(())
    }

    pub fn set_init_prob(&mut self, init_prob: f64) -> Result<(), Error> {
        validate_param!(init_prob > 0.0 && init_prob < 1.0,
            "Initial probability ({}) must be within (0, 1)", init_prob);
        self.init_prob = init_prob;
        Ok(())
    }

    pub fn set_plato_size(&mut self, plato_size: usize) {
        self.plato_size = plato_size;
    }

    /// Finds temperature coefficient by checking 100 random reassignments and their probabilities.
    fn find_temperature_coeff(&self,
        assignments: &ReadAssignment,
        rng: &mut XoshiroRng
    ) -> f64 {
        let mut neg_sum = 0.0;
        let mut neg_count: u32 = 0;
        const INIT_ITERS: u32 = 100;
        for _ in 0..INIT_ITERS {
            let diff = assignments.calculate_improvement(&ReassignmentTarget::random(&assignments, rng));
            if diff < 0.0 {
                neg_sum += diff;
                neg_count += 1;
            }
        }
        if neg_count == 0 {
            1.0
        } else {
            self.init_prob.ln() * f64::from(neg_count) / neg_sum
        }
    }

    fn accept_change(diff: f64, curr_temp: f64, rng: &mut impl Rng) -> bool {
        diff > 0.0 || (diff > 0.0 && rng.gen::<f64>() <= (diff / curr_temp).exp())
    }
}

impl Solver for SimAnneal {
    /// Run simulated annealing to find the best read assignment.
    fn solve_nontrivial<'a>(
        &self,
        gt_alns: &'a GenotypeAlignments,
        rng: &mut XoshiroRng,
    ) -> Result<ReadAssignment<'a>, Error>
    {
        // 0 - index of the best alignment.
        let mut assignments = ReadAssignment::new(gt_alns, |_| 0);
        let coeff = self.find_temperature_coeff(&assignments, rng);

        let mut curr_temp = 1.0 / coeff;
        let cooling_temp = self.cooling_temp / coeff;
        let mut curr_plato = 0;
        while curr_plato <= self.plato_size {
            let target = ReassignmentTarget::random(&assignments, rng);
            if Self::accept_change(assignments.calculate_improvement(&target), curr_temp, rng) {
                assignments.reassign(&target);
                curr_plato = 0;
            } else {
                curr_plato += 1;
            }
            curr_temp -= cooling_temp;
        }
        Ok(assignments)
    }
}

impl super::SetParams for SimAnneal {
    /// Sets solver parameters.
    fn set_param(&mut self, key: &str, val: &str) -> Result<(), Error> {
        match &key.to_lowercase() as &str {
            "cooling" | "coolingtemp" | "cooling_temp" => self.set_cooling_temp(val.parse()
                .map_err(|_| Error::InvalidInput(format!("Cannot parse '{}={}'", key, val)))?)?,
            "init_prob" | "initprob" => self.set_init_prob(val.parse()
                .map_err(|_| Error::InvalidInput(format!("Cannot parse '{}={}'", key, val)))?)?,
            "plato" => self.set_plato_size(val.parse()
                .map_err(|_| Error::InvalidInput(format!("Cannot parse '{}={}'", key, val)))?),
            _ => log::error!("Sim.Annealing: unknown parameter {:?}", key),
        }
        Ok(())
    }
}

impl fmt::Display for SimAnneal {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SimAnneal(cool.temp {}, init.prob {}, plato {})",
            self.cooling_temp, self.init_prob, self.plato_size)
    }
}
