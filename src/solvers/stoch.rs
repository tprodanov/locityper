use std::fmt;
use rand::{Rng, seq::SliceRandom};
use crate::{
    Error,
    ext::{
        vec::{F64Ext, IterExt},
        rand::XoshiroRng,
    },
    model::{
        windows::ReadWindows,
        assgn::ReadAssignment,
    },
    bg::ser::json_get,
};
use super::MultiTrySolver;

/// Assigns reads in a greedy way.
///
/// In one step, the solver examines `sample_size` read pairs, and selects the best read pair to switch location.
/// If no improvement was made for `plato_size` iterations, the solver stops.
#[derive(Clone)]
pub struct GreedySolver {
    /// Number of tries the solver makes to assign reads anew.
    tries: u16,
    /// Number of read-pairs, examined per iteration.
    sample_size: usize,
    /// Number of iteration without improvement, after which the solver stops.
    plato_size: usize,
}

impl Default for GreedySolver {
    fn default() -> Self {
        Self {
            tries: 3,
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

impl MultiTrySolver for GreedySolver {
    fn tries(&self) -> u16 {
        self.tries
    }

    fn set_tries(&mut self, tries: u16) -> &mut Self {
        assert_ne!(tries, 0, "Number of tries cannot be 0.");
        self.tries = tries;
        self
    }

    /// Single greedy iteration to find the best read assignment.
    fn solve_once(&self, assignments: &mut ReadAssignment, rng: &mut XoshiroRng) -> Result<(), Error> {
        let mut buffer = Vec::new();
        assignments.init_assignments(|alns| rng.gen_range(0..alns.len()));
        let mut curr_plato = 0;
        while curr_plato <= self.plato_size {
            let mut best_rp = 0;
            let mut best_assgn = 0;
            let mut best_improv = 0.0;
            for &rp in assignments.non_trivial_reads().choose_multiple(rng, self.sample_size) {
                buffer.clear();
                assignments.possible_reassignments(rp, &mut buffer);
                let (assgn, improv) = F64Ext::argmax(&buffer);
                if improv > best_improv {
                    best_rp = rp;
                    best_assgn = assgn as u16;
                    best_improv = improv;
                }
            }
            if best_improv > 0.0 {
                curr_plato = 0;
                assignments.reassign(best_rp, best_assgn);
            } else {
                curr_plato += 1;
            }
        }
        Ok(())
    }
}

impl super::SetParams for GreedySolver {
    /// Sets solver parameters.
    fn set_params(&mut self, obj: &json::JsonValue) -> Result<(), Error> {
        if obj.has_key("tries") {
            json_get!(obj -> tries (as_u16));
            self.set_tries(tries);
        }
        if obj.has_key("sample_size") {
            json_get!(obj -> sample_size (as_usize));
            self.set_sample_size(sample_size);
        }
        if obj.has_key("plato_size") {
            json_get!(obj -> plato_size (as_usize));
            self.set_plato_size(plato_size);
        }
        Ok(())
    }
}

impl fmt::Display for GreedySolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Greedy({} tries, {} reads/iter, plato {})", self.tries, self.sample_size, self.plato_size)
    }
}

/// Simulated annealing solver.
///
/// Randomly selects direction based on the current temperature, and stops once there are no improvement `plato_size`.
#[derive(Clone)]
pub struct SimAnneal {
    /// Number of tries the solver makes to assign reads anew.
    tries: u16,
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
            tries: 3,
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
    fn find_temperature_coeff(&self, assignments: &ReadAssignment, rng: &mut XoshiroRng) -> f64 {
        let mut neg_sum = 0.0;
        let mut neg_count: u32 = 0;
        const INIT_ITERS: u32 = 100;
        for _ in 0..INIT_ITERS {
            let diff = assignments.random_reassignment(rng).2;
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
}

impl MultiTrySolver for SimAnneal {
    fn tries(&self) -> u16 {
        self.tries
    }

    fn set_tries(&mut self, tries: u16) -> &mut Self {
        assert_ne!(tries, 0, "Number of tries cannot be 0.");
        self.tries = tries;
        self
    }

    /// Run simulated annealing once to find the best read assignment.
    fn solve_once(&self, assignments: &mut ReadAssignment, rng: &mut XoshiroRng) -> Result<(), Error> {
        assignments.init_assignments(
            |possible_alns| IterExt::argmax(possible_alns.iter().map(ReadWindows::ln_prob)).0);
        let coeff = self.find_temperature_coeff(assignments, rng);

        let mut curr_temp = 1.0 / coeff;
        let cooling_temp = self.cooling_temp / coeff;
        let mut curr_plato = 0;
        while curr_plato <= self.plato_size {
            let (rp, new_assign, improv) = assignments.random_reassignment(rng);
            if improv > 0.0 || (curr_temp > 0.0 && rng.gen::<f64>() <= (improv / curr_temp).exp()) {
                assignments.reassign(rp, new_assign);
                curr_plato = 0;
            } else {
                curr_plato += 1;
            }
            curr_temp -= cooling_temp;
        }
        Ok(())
    }
}

impl super::SetParams for SimAnneal {
    /// Sets solver parameters.
    fn set_params(&mut self, obj: &json::JsonValue) -> Result<(), Error> {
        // TODO: Set params only if they exist.
        json_get!(obj -> tries (as_u16), cooling_temp (as_f64), init_prob (as_f64), plato_size (as_usize));
        self.set_tries(tries)
            .set_cooling_temp(cooling_temp)
            .set_init_prob(init_prob)
            .set_plato_size(plato_size);
        Ok(())
    }
}

impl fmt::Display for SimAnneal {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SimAnneal({} tries, cool.temp {}, init.prob {}, plato {})",
            self.tries, self.cooling_temp, self.init_prob, self.plato_size)
    }
}
