use crate::{
    model::assgn::ReadAssignment,
    bg::ser::JsonSer,
};

#[cfg(feature = "stochastic")]
pub mod stoch;
#[cfg(feature = "gurobi")]
pub mod gurobi;
#[cfg(feature = "highs")]
pub mod highs;

#[cfg(feature = "stochastic")]
pub use self::stoch::{GreedySolver, SimAnneal};
#[cfg(feature = "gurobi")]
pub use self::gurobi::GurobiSolver;
#[cfg(feature = "highs")]
pub use self::highs::HighsSolver;
use crate::Error;

/// General trait for all solvers.
pub trait Solver : Default + Send {
    /// Distribute reads between several haplotypes in a best way.
    fn solve(&self, assignments: &mut ReadAssignment, rng: &mut impl rand::Rng) -> Result<(), Error>;

    /// Sets solver parameters.
    fn set_params(&mut self, obj: &json::JsonValue) -> Result<(), Error>;
}

/// Solver that tries several times and selects the best result.
pub trait MultiTrySolver {
    /// Returns the number of tries.
    fn tries(&self) -> u16;

    /// Set the number of tries.
    fn set_tries(&mut self, tries: u16) -> &mut Self;

    /// Try once.
    fn solve_once(&self, assignments: &mut ReadAssignment, rng: &mut impl rand::Rng) -> Result<(), Error>;

    /// Sets solver parameters.
    fn set_params(&mut self, obj: &json::JsonValue) -> Result<(), Error>;
}

impl<T: MultiTrySolver + Default + Send + JsonSer> Solver for T {
    fn solve(&self, assignments: &mut ReadAssignment, rng: &mut impl rand::Rng) -> Result<(), Error> {
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
        Ok(())
    }

    fn set_params(&mut self, obj: &json::JsonValue) -> Result<(), Error> {
        self.set_params(obj)
    }
}
