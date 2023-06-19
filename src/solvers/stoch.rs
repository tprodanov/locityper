use std::fmt;
use rand::Rng;
use crate::{
    Error,
    ext::{
        vec::IterExt,
        rand::XoshiroRng,
    },
    model::{
        windows::ReadWindows,
        assgn::{ReadAssignment, ReassignmentTarget},
    },
    bg::ser::json_get,
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
    pub fn set_sample_size(&mut self, sample_size: usize) -> &mut Self {
        assert_ne!(sample_size, 0, "Number of read-pairs (sample_size) cannot be 0.");
        self.sample_size = sample_size;
        self
    }

    pub fn set_plato_size(&mut self, plato_size: usize) -> &mut Self {
        self.plato_size = plato_size;
        self
    }
}

impl Solver for GreedySolver {
    /// Single greedy iteration to find the best read assignment.
    fn solve(&self, assignments: &mut ReadAssignment, rng: &mut XoshiroRng) -> Result<f64, Error> {
        assignments.init_assignments(|alns| rng.gen_range(0..alns.len()));
        let mut target = ReassignmentTarget::new();
        let mut curr_plato = 0;

        while curr_plato <= self.plato_size {
            let mut best_rp = 0;
            let mut best_assgn = 0;
            let mut best_improv = 0.0;
            for _ in 0..self.sample_size {
                target.set_random(assignments, rng);
                let improv = assignments.calculate_improvement(&target);
                if improv > best_improv {
                    (best_rp, best_assgn) = target.get();
                    best_improv = improv;
                }
            }
            if best_improv > 0.0 {
                curr_plato = 0;
                target.set(best_rp, best_assgn, assignments);
                assignments.reassign(&target);
            } else {
                curr_plato += 1;
            }
        }
        Ok(assignments.likelihood())
    }
}

impl super::SetParams for GreedySolver {
    /// Sets solver parameters.
    fn set_params(&mut self, obj: &json::JsonValue) -> Result<(), Error> {
        json_get!(obj -> sample_size? (as_usize), plato_size? (as_usize));
        if let Some(val) = sample_size {
            self.set_sample_size(val);
        }
        if let Some(val) = plato_size {
            self.set_plato_size(val);
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
    pub fn set_cooling_temp(&mut self, cooling_temp: f64) -> &mut Self {
        assert!(cooling_temp > 0.0 && cooling_temp < 1.0,
            "Cooling temperature ({}) must be within (0, 1).", cooling_temp);
        self.cooling_temp = cooling_temp;
        self
    }

    pub fn set_init_prob(&mut self, init_prob: f64) -> &mut Self {
        assert!(init_prob > 0.0 && init_prob < 1.0, "Initial probability ({}) must be within (0, 1).", init_prob);
        self.init_prob = init_prob;
        self
    }

    pub fn set_plato_size(&mut self, plato_size: usize) -> &mut Self {
        self.plato_size = plato_size;
        self
    }

    /// Finds temperature coefficient by checking 100 random reassignments and their probabilities.
    fn find_temperature_coeff(&self,
        target: &mut ReassignmentTarget,
        assignments: &ReadAssignment,
        rng: &mut XoshiroRng
    ) -> f64 {
        let mut neg_sum = 0.0;
        let mut neg_count: u32 = 0;
        const INIT_ITERS: u32 = 100;
        for _ in 0..INIT_ITERS {
            target.set_random(assignments, rng);
            let diff = assignments.calculate_improvement(target);
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
    /// Run simulated annealing once to find the best read assignment.
    fn solve(&self, assignments: &mut ReadAssignment, rng: &mut XoshiroRng) -> Result<f64, Error> {
        assignments.init_assignments(
            |possible_alns| IterExt::argmax(possible_alns.iter().map(ReadWindows::ln_prob)).0);
        let mut target = ReassignmentTarget::new();
        let coeff = self.find_temperature_coeff(&mut target, assignments, rng);

        let mut curr_temp = 1.0 / coeff;
        let cooling_temp = self.cooling_temp / coeff;
        let mut curr_plato = 0;
        while curr_plato <= self.plato_size {
            target.set_random(assignments, rng);
            if Self::accept_change(assignments.calculate_improvement(&target), curr_temp, rng) {
                assignments.reassign(&target);
                curr_plato = 0;
            } else {
                curr_plato += 1;
            }
            curr_temp -= cooling_temp;
        }
        Ok(assignments.likelihood())
    }
}

impl super::SetParams for SimAnneal {
    /// Sets solver parameters.
    fn set_params(&mut self, obj: &json::JsonValue) -> Result<(), Error> {
        json_get!(obj -> cooling_temp? (as_f64), init_prob? (as_f64), plato_size? (as_usize));
        if let Some(val) = cooling_temp {
            self.set_cooling_temp(val);
        }
        if let Some(val) = init_prob {
            self.set_init_prob(val);
        }
        if let Some(val) = plato_size {
            self.set_plato_size(val);
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
