use std::{fmt, cmp::min};
use rand::{
    Rng,
    SeedableRng,
    rngs::SmallRng,
    seq::SliceRandom,
};
use crate::{
    Error,
    ext::vec::F64Ext,
    model::assgn::ReadAssignment,
    bg::ser::{JsonSer, json_get},
};
use super::Solver;

/// Assigns reads in a greedy way.
///
/// In one step, the solver examines `sample_size` read pairs, and selects the best read pair to switch location.
/// If no improvement was made for `plato_iters` iterations, the solver stops.
#[derive(Clone)]
pub struct GreedySolver {
    /// Number of tries the solver makes to assign reads anew.
    tries: usize,
    /// Number of read-pairs, examined per iteration.
    sample_size: usize,
    /// Number of iteration without improvement, after which the solver stops.
    plato_iters: usize,
}

impl Default for GreedySolver {
    fn default() -> Self {
        Self {
            tries: 3,
            sample_size: 100,
            plato_iters: 5,
        }
    }
}

impl GreedySolver {
    pub fn set_tries(&mut self, tries: usize) -> &mut Self {
        assert_ne!(tries, 0, "Number of tries cannot be 0.");
        self.tries = tries;
        self
    }

    pub fn set_sample_size(&mut self, sample_size: usize) -> &mut Self {
        assert_ne!(sample_size, 0, "Number of read-pairs (sample_size) cannot be 0.");
        self.sample_size = sample_size;
        self
    }

    pub fn set_plato_iters(&mut self, plato_iters: usize) -> &mut Self {
        self.plato_iters = plato_iters;
        self
    }

    /// Single greedy iteration to find the best read assignment.
    fn solve_once(
        &self,
        assignments: &mut ReadAssignment,
        rng: &mut impl Rng,
        buffer: &mut Vec<f64>
    ) -> Result<(), Error>
    {
        assignments.init_assignments(|alns| rng.gen_range(0..alns.len()));
        let mut curr_plato = 0;
        while curr_plato <= self.plato_iters {
            let mut best_rp = 0;
            let mut best_assgn = 0;
            let mut best_improv = 0.0;
            for &rp in assignments.non_trivial_reads().choose_multiple(rng, self.sample_size) {
                buffer.clear();
                assignments.possible_reassignments(rp, buffer);
                let (assgn, improv) = F64Ext::argmax(buffer);
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

impl Solver for GreedySolver {
    fn solve(&self, assignments: &mut ReadAssignment, rng: &mut impl Rng) -> Result<(), Error> {
        let mut last_lik = f64::NEG_INFINITY;
        let mut best_lik = f64::NEG_INFINITY;
        let mut best_assgns = assignments.read_assignments().to_vec();
        let mut buffer = Vec::new();

        for _ in 0..self.tries {
            self.solve_once(assignments, rng, &mut buffer).expect("None is impossible");
            last_lik = assignments.likelihood();
            if last_lik > best_lik {
                best_assgns.copy_from_slice(assignments.read_assignments());
                best_lik = last_lik;
            }
        }
        if best_lik > last_lik {
            assignments.set_assignments(&best_assgns);
        }
        Ok(())
    }
}

impl fmt::Display for GreedySolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Greedy({} tries, {} reads/iter, plato {})", self.tries, self.sample_size, self.plato_iters)
    }
}

impl JsonSer for GreedySolver {
    fn save(&self) -> json::JsonValue {
        json::object!{
            tries: self.tries,
            sample_size: self.sample_size,
            plato_iters: self.plato_iters,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> tries (as_usize), sample_size (as_usize), plato_iters (as_usize));
        let mut res = Self::default();
        res.set_tries(tries).set_sample_size(sample_size).set_plato_iters(plato_iters);
        Ok(res)
    }
}
