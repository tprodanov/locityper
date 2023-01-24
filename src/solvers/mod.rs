use std::io;
use crate::{
    algo::vec_ext::IterExt,
    model::{
        locs::PairAlignment,
        assgn::ReadAssignment,
    },
};

#[cfg(feature = "stochastic")]
pub mod greedy;
pub mod dbg;

#[cfg(feature = "stochastic")]
pub use crate::solvers::greedy::GreedySolver;
pub use crate::solvers::dbg::{DbgWrite, NoDbg, DbgWriter};
use crate::solvers::dbg::Iteration;

pub trait SolverBuilder {
    type S<'a>: Solver<'a>;

    /// Sets seed.
    /// Can panic if the seed does not fit the model, or if the solver is deterministic.
    fn set_seed(&mut self, seed: u64) -> &mut Self;

    /// Builds the solver.
    fn build<'a>(&self, assignments: ReadAssignment<'a>) -> Self::S<'a>;
}

/// Trait that distributes the reads between their possible alignments
pub trait Solver<'a> {
    /// Static function that returns true, if running the method multiple times will produce the same results
    /// irrespective of the seed.
    fn is_determenistic() -> bool;

    /// Initialize read alignments. Returns current likelihood.
    fn initialize(&mut self) -> f64;

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> f64;

    /// Returns true if the solver is finished.
    fn is_finished(&self) -> bool;

    /// Return the current read assignments.
    // Comment: Lifetime 'a outlives 'b.
    fn current_assignments<'b>(&'b self) -> &'b ReadAssignment<'a> where 'a: 'b;

    /// Finish solving, consume the solver and return the read assignments.
    fn finish(self) -> ReadAssignment<'a>;
}

/// Initialize reads based on their best alignment.
fn init_best(possible_alns: &[PairAlignment]) -> usize {
    IterExt::argmax(possible_alns.iter().map(PairAlignment::ln_prob)).0
}

/// Distribute read assignment in at most `max_iters` iterations.
pub fn solve<'a, B, I, W>(
    mut assgns: ReadAssignment<'a>,
    mut solver_builder: B,
    seeds: I,
    max_iters: u32,
    dbg_writer: &mut W
) -> io::Result<ReadAssignment<'a>>
where B: SolverBuilder,
      I: Iterator<Item = u64>,
      W: DbgWrite,
{
    let mut best_lik = f64::NEG_INFINITY;
    let mut last_lik = best_lik;
    let mut best_assgns: Vec<u16> = assgns.read_assignments().to_vec();

    let mut outer = 0;
    for seed in seeds {
        outer += 1;
        let mut solver = if B::S::<'a>::is_determenistic() {
            solver_builder.build(assgns)
        } else {
            solver_builder.set_seed(seed).build(assgns)
        };

        last_lik = solver.initialize();
        dbg_writer.write(solver.current_assignments(), Iteration::Init(outer))?;
        if last_lik > best_lik {
            best_lik = last_lik;
            best_assgns.clone_from_slice(solver.current_assignments().read_assignments());
        }

        for inner in 1..=max_iters {
            solver.step();
            last_lik = solver.current_assignments().likelihood();
            if last_lik > best_lik {
                best_lik = last_lik;
                best_assgns.clone_from_slice(solver.current_assignments().read_assignments());
            }

            if solver.is_finished() {
                dbg_writer.write(solver.current_assignments(), Iteration::Last(outer, inner))?;
                break;
            } else {
                dbg_writer.write(solver.current_assignments(), Iteration::Step(outer, inner))?;
            }
        }
        assgns = solver.finish();
        if B::S::<'a>::is_determenistic() {
            break;
        }
    }

    if last_lik < best_lik {
        assgns.set_assignments(&best_assgns);
    }
    dbg_writer.write(&assgns, Iteration::Best)?;
    Ok(assgns)
}
