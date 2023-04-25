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
pub mod greedy;
#[cfg(feature = "stochastic")]
pub mod anneal;
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
