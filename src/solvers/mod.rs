use std::{
    io,
    fmt::Debug,
};
use crate::{
    algo::vec_ext::IterExt,
    model::{
        locs::PairAlignment,
        assgn::ReadAssignment,
    },
};

pub mod dbg;
#[cfg(feature = "stochastic")]
pub mod greedy;
#[cfg(feature = "stochastic")]
pub mod anneal;
#[cfg(feature = "gurobi")]
pub mod gurobi;
#[cfg(feature = "highs")]
pub mod highs;

pub use crate::solvers::dbg::{DbgWrite, NoDbg, DbgWriter};
use crate::solvers::dbg::Iteration;
#[cfg(feature = "stochastic")]
pub use crate::solvers::{
    greedy::GreedySolver,
    anneal::SimulatedAnnealing,
};
#[cfg(feature = "gurobi")]
pub use crate::solvers::gurobi::GurobiSolver;
#[cfg(feature = "highs")]
pub use crate::solvers::highs::HighsSolver;

/// Trait that distributes the reads between their possible alignments
pub trait Solver {
    type Error: Debug;

    /// Returns true if the solver can take seed.
    fn is_seedable() -> bool;

    /// Sets seed.
    /// Can panic if the seed does not fit the model, or if the solver is deterministic.
    fn set_seed(&mut self, seed: u64) -> Result<(), Self::Error>;

    /// Resets and initializes anew read assignments.
    fn initialize(&mut self) -> Result<(), Self::Error>;

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> Result<f64, Self::Error>;

    /// Returns true if the solver is finished.
    fn is_finished(&self) -> bool;

    /// Return the current read assignments.
    fn current_assignments(&self) -> &ReadAssignment;

    /// Recalculate likelihood and check if it matches the stored one.
    fn recalculate_likelihood(&mut self);

    /// Consumes solver and returns the read assignments.
    fn take(self) -> ReadAssignment;
}

/// Initialize reads based on their best alignment.
fn init_best(possible_alns: &[PairAlignment]) -> usize {
    IterExt::argmax(possible_alns.iter().map(PairAlignment::ln_prob)).0
}

#[derive(Debug)]
pub enum Error<E: Debug> {
    IoErr(io::Error),
    SolverErr(E),
}

impl<E: Debug> From<io::Error> for Error<E> {
    fn from(err: io::Error) -> Self {
        Self::IoErr(err)
    }
}

/// Distribute read assignment in at most `max_iters` iterations.
pub fn solve<S, I, W>(
    mut solver: S,
    seeds: I,
    max_iters: u32,
    dbg_writer: &mut W
) -> Result<ReadAssignment, Error<S::Error>>
where S: Solver,
      I: Iterator<Item = u64>,
      W: DbgWrite,
{
    let mut best_lik = f64::NEG_INFINITY;
    let mut last_lik = f64::NEG_INFINITY;
    let mut best_assgns: Vec<u16> = solver.current_assignments().read_assignments().to_vec();

    let mut outer = 0;
    for seed in seeds {
        if S::is_seedable() {
            solver.set_seed(seed).map_err(Error::SolverErr)?;
        } else if outer > 0 {
            break;
        }

        outer += 1;
        for inner in 0..=max_iters {
            if inner == 0 {
                solver.initialize().map_err(Error::SolverErr)?;
            } else {
                solver.step().map_err(Error::SolverErr)?;
            }
            last_lik = solver.current_assignments().likelihood();
            if last_lik > best_lik {
                best_lik = last_lik;
                best_assgns.clone_from_slice(solver.current_assignments().read_assignments());
            }

            if inner == 0 {
                dbg_writer.write(solver.current_assignments(), Iteration::Init(outer))?;
            } else if solver.is_finished() {
                dbg_writer.write(solver.current_assignments(), Iteration::Last(outer, inner))?;
                break;
            } else {
                dbg_writer.write(solver.current_assignments(), Iteration::Step(outer, inner))?;
            }
        }
        solver.recalculate_likelihood();
        let new_lik = solver.current_assignments().likelihood();
        let divergence = (last_lik - new_lik).abs();
        assert!(divergence < 1e-2, "Likelihood estimates diverged too much {} and {}", last_lik, new_lik);
        if divergence > 1e-6 {
            log::error!("Likelihood estimates diverged by {} ({} and {})", divergence, last_lik, new_lik);
        }
    }

    let mut assgns = solver.take();
    if last_lik < best_lik {
        assgns.set_assignments(&best_assgns);
    }
    dbg_writer.write(&assgns, Iteration::Best)?;
    Ok(assgns)
}
