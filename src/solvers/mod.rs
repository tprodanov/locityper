use std::fmt::Display;
use crate::{
    model::assgn::ReadAssignment,
};

pub mod scheme;
pub mod stoch;
#[cfg(feature = "gurobi")]
pub mod gurobi;
#[cfg(feature = "highs")]
pub mod highs;

pub use self::stoch::{GreedySolver, SimAnneal};
#[cfg(feature = "gurobi")]
pub use self::gurobi::GurobiSolver;
#[cfg(feature = "highs")]
pub use self::highs::HighsSolver;
use crate::Error;

/// In order for `Solver` to be object-safe, rng type should be known in advance.
/// For that reason we use `SmallRng`, and not `impl rand::Rng`.
type SolverRng = rand::rngs::SmallRng;

pub trait SetParams {
    /// Sets solver parameters.
    fn set_params(&mut self, obj: &json::JsonValue) -> Result<(), Error>;
}

/// General trait for all solvers.
pub trait Solver : Send + SetParams + CloneSolver + Display {
    /// Distribute reads between several haplotypes in a best way.
    /// Returns likelihood of the assignment.
    fn solve(&self, assignments: &mut ReadAssignment, rng: &mut SolverRng) -> Result<f64, Error>;
}

/// Solver that tries several times and selects the best result.
pub trait MultiTrySolver {
    /// Returns the number of tries.
    fn tries(&self) -> u16;

    /// Set the number of tries.
    fn set_tries(&mut self, tries: u16) -> &mut Self;

    /// Try once.
    fn solve_once(&self, assignments: &mut ReadAssignment, rng: &mut SolverRng) -> Result<(), Error>;
}

impl<T: 'static + MultiTrySolver + SetParams + Send + Clone + Display> Solver for T {
    fn solve(&self, assignments: &mut ReadAssignment, rng: &mut SolverRng) -> Result<f64, Error> {
        let mut last_lik = f64::NEG_INFINITY;
        let mut best_lik = assignments.likelihood();
        let mut best_assgns = assignments.read_assignments().to_vec();

        for _ in 0..self.tries() {
            self.solve_once(assignments, rng)?;
            last_lik = assignments.likelihood();
            if last_lik > best_lik {
                best_assgns.copy_from_slice(assignments.read_assignments());
                best_lik = last_lik;
            }
        }
        if best_lik > last_lik {
            assignments.set_assignments(&best_assgns);
        }
        Ok(best_lik)
    }
}

/// We cannot ask `Solver` to inherit `Clone`, because then the trait is not object-safe.
pub trait CloneSolver {
    fn clone_box(&self) -> Box<dyn Solver>;
}

impl<T: 'static + Solver + Clone> CloneSolver for T {
    fn clone_box(&self) -> Box<dyn Solver> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Solver> {
    fn clone(&self) -> Box<dyn Solver> {
        self.clone_box()
    }
}