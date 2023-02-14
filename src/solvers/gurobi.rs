use grb::*;
use std::result::Result;
use crate::model::assgn::ReadAssignment;
use super::Solver;

pub struct GurobiSolver {
    model: Model,
}

impl Solver for GurobiSolver {
    type Error = grb::Error;

    fn is_seedable() -> bool { true }

    fn set_seed(&mut self, seed: u64) -> Result<(), Self::Error> {
        self.model.set_param(parameter::IntParam::Seed, seed as i32)
    }

    /// Resets and initializes anew read assignments.
    fn initialize(&mut self) -> Result<(), Self::Error> {
        unimplemented!()
    }

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> Result<f64, Self::Error> {
        unimplemented!()
    }

    /// Returns true if the solver is finished.
    fn is_finished(&self) -> bool {
        unimplemented!()
    }

    /// Return the current read assignments.
    fn current_assignments(&self) -> &ReadAssignment {
        unimplemented!()
    }

    /// Recalculate likelihood and check if it matches the stored one.
    fn recalculate_likelihood(&mut self) {
        unimplemented!()
    }

    /// Consumes solver and returns the read assignments.
    fn take(self) -> ReadAssignment {
        unimplemented!()
    }
}