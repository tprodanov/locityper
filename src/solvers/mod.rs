use std::{
    io,
    time::Instant,
    any::TypeId,
    fmt::Display,
};
use crate::{
    ext::vec::IterExt,
    model::{
        windows::ReadWindows,
        assgn::ReadAssignment,
        locs::AllPairAlignments,
    },
    bg::ser::JsonSer,
};

#[cfg(feature = "stochastic")]
pub mod stoch;
#[cfg(feature = "gurobi")]
pub mod gurobi;
#[cfg(feature = "highs")]
pub mod highs;

// #[cfg(feature = "stochastic")]
// pub use self::{
//     greedy::GreedySolver,
//     anneal::SimulatedAnnealing,
// };
// #[cfg(feature = "gurobi")]
// pub use self::gurobi::GurobiSolver;
// #[cfg(feature = "highs")]
// pub use self::highs::HighsSolver;
use crate::Error;

/// General trait for all solvers.
pub trait Solver : Default + Send + JsonSer {
    /// Distribute reads between several haplotypes in a best way.
    fn solve(&self, assignments: &mut ReadAssignment, rng: &mut impl rand::Rng) -> Result<(), Error>;
}

/// Solver that tries several times and selects the best result.
pub trait MultiTrySolver {
    /// Returns the number of tries.
    fn tries(&self) -> usize;

    /// Set the number of tries.
    fn set_tries(&mut self, tries: usize) -> &mut Self;

    /// Try once.
    fn solve_once(&self, assignments: &mut ReadAssignment, rng: &mut impl rand::Rng) -> Result<(), Error>;
}

impl<T: MultiTrySolver + Default + Send + JsonSer> Solver for T {
    fn solve(&self, assignments: &mut ReadAssignment, rng: &mut impl rand::Rng) -> Result<(), Error> {
        let mut last_lik = f64::NEG_INFINITY;
        let mut best_lik = f64::NEG_INFINITY;
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
}
