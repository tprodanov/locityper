use std::fmt::Display;
use crate::{
    model::assgn::ReadAssignment,
    ext::rand::XoshiroRng,
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

pub trait SetParams {
    /// Sets solver parameters.
    fn set_params(&mut self, obj: &json::JsonValue) -> Result<(), Error>;
}

/// General trait for all solvers.
pub trait Solver : Send + Sync + SetParams + CloneSolver + Display {
    /// Same as `solve`, but it is known that there is at least one read with more than one potential location.
    fn solve_nontrivial(&self, assignments: &mut ReadAssignment, rng: &mut XoshiroRng) -> Result<(), Error>;

    /// Distribute reads between several haplotypes in a best way.
    /// Returns likelihood of the assignment.
    ///
    /// In order for `Solver` to be object-safe, rng type should be known in advance.
    fn solve(&self, assignments: &mut ReadAssignment, rng: &mut XoshiroRng) -> Result<f64, Error> {
        if assignments.trivial() {
            assignments.init_assignments(|_| 0);
        } else {
            self.solve_nontrivial(assignments, rng)?;
        }
        Ok(assignments.likelihood())
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
