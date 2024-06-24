use std::fmt::Display;
use crate::{
    model::assgn::{GenotypeAlignments, ReadAssignment},
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
use crate::err::error;

pub trait SetParams {
    /// Sets solver parameters.
    /// Each element: string `key=value`.
    fn set_params(&mut self, params: &[String]) -> crate::Result<()> {
        for param in params {
            let (key, value) = param.split_once('=')
                .ok_or_else(|| error!(InvalidInput, "Cannot parse parameter {:?} (must contain =)", param))?;
            self.set_param(key.trim(), value.trim())?;
        }
        Ok(())
    }

    fn set_param(&mut self, key: &str, val: &str) -> crate::Result<()>;
}

/// General trait for all solvers.
pub trait Solver: Send + Sync + SetParams + CloneSolver + Display {
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
