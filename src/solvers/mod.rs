use std::{
    fmt,
};
use crate::{
    model::assgn::{GenotypeAlignments, ReadAssignment},
    ext::rand::XoshiroRng,
};

pub mod solve;
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

pub enum ParamErr {
    Unknown,
    Parse,
    Invalid(String),
}

impl From<std::num::ParseFloatError> for ParamErr {
    fn from(_: std::num::ParseFloatError) -> Self { Self::Parse }
}

impl From<std::num::ParseIntError> for ParamErr {
    fn from(_: std::num::ParseIntError) -> Self { Self::Parse }
}

impl From<std::str::ParseBoolError> for ParamErr {
    fn from(_: std::str::ParseBoolError) -> Self { Self::Parse }
}

pub trait SetParams {
    fn set_param(&mut self, key: &str, val: &str) -> Result<(), ParamErr>;
}

/// General trait for all solvers.
pub trait Solver: Send + Sync + SetParams + CloneSolver + fmt::Display {
    /// Distribute reads between several haplotypes in the best way,
    /// when at least one read pair has several possible locations.
    fn solve_nontrivial<'a>(
        &self,
        gt_alns: &'a GenotypeAlignments,
        rng: &mut XoshiroRng,
    ) -> crate::Result<ReadAssignment<'a>>;

    /// Distribute reads between several haplotypes in the best way.
    ///
    /// In order for `Solver` to be object-safe, rng type should be known in advance.
    fn solve<'a>(
        &self,
        gt_alns: &'a GenotypeAlignments,
        rng: &mut XoshiroRng,
    ) -> crate::Result<ReadAssignment<'a>>
    {
        if gt_alns.trivial() {
            Ok(ReadAssignment::new(gt_alns, |_| unreachable!("There are no non-trivial reads")))
        } else {
            self.solve_nontrivial(gt_alns, rng)
        }
    }

    fn describe_params(&self) -> String { String::new() }
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
