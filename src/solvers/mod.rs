use std::{
    io,
    time::Instant,
    any::TypeId,
    fmt::{Display, Debug},
};
use crate::{
    algo::vec_ext::IterExt,
    model::{
        locs::PairAlignment,
        assgn::ReadAssignment,
    },
};

#[cfg(feature = "stochastic")]
pub mod greedy;
#[cfg(feature = "stochastic")]
pub mod anneal;
#[cfg(feature = "gurobi")]
pub mod gurobi;
#[cfg(feature = "highs")]
pub mod highs;

#[cfg(feature = "stochastic")]
pub use self::{
    greedy::GreedySolver,
    anneal::SimulatedAnnealing,
};
#[cfg(feature = "gurobi")]
pub use self::gurobi::GurobiSolver;
#[cfg(feature = "highs")]
pub use self::highs::HighsSolver;

/// Trait that distributes the reads between their possible alignments
pub trait Solver: Display + Display {
    type Error: Debug;

    /// Returns true if the solver can take seed.
    fn is_seedable() -> bool;

    /// Sets seed.
    /// Can panic if the seed does not fit the model, or if the solver is deterministic.
    fn set_seed(&mut self, seed: u64) -> Result<(), Self::Error>;

    /// Resets the solver.
    fn reset(&mut self) -> Result<(), Self::Error>;

    /// Perform one iteration.
    fn step(&mut self) -> Result<(), Self::Error>;

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
///
/// dbg_writer writes intermediate likelihood values and runtime for each seed.
/// Format: `assignment-tag  solver  runtime  seed  likelihood`.
pub fn solve<S, F, I, W, U>(
    assgns: ReadAssignment,
    build: F,
    seeds: I,
    dbg_writer: &mut W,
    assgn_writer: &mut U,
) -> Result<ReadAssignment, Error<S::Error>>
where S: Solver,
      F: Fn(ReadAssignment) -> S,
      I: Iterator<Item = u64>,
      W: io::Write,
      U: io::Write + 'static,
{
    let mut best_lik = f64::NEG_INFINITY;
    let mut last_lik = f64::NEG_INFINITY;
    let mut best_assgns: Vec<u16> = assgns.read_assignments().to_vec();

    let mut last_instant = Instant::now();
    let mut update_timer = || {
        let new_instant = Instant::now();
        let duration = new_instant.duration_since(last_instant);
        last_instant = new_instant;
        duration
    };

    let contigs_group_str = assgns.contigs_group().to_string();
    let mut solver = build(assgns);
    let solver_str = solver.to_string();
    log::debug!("Solving {} with {}", contigs_group_str, solver_str);
    let dur = update_timer();
    let prefix = format!("{}\t{}", contigs_group_str, solver_str);
    writeln!(dbg_writer, "{}\t{}.{:06}\tstart\tNA", prefix, dur.as_secs(), dur.subsec_micros())?;
    let mut outer = 0;
    for mut seed in seeds {
        if S::is_seedable() {
            solver.set_seed(seed).map_err(Error::SolverErr)?;
        } else if outer == 0 {
            seed = 0;
        } else {
            break;
        }

        outer += 1;
        for inner in 0.. {
            if inner == 0 {
                solver.reset().map_err(Error::SolverErr)?;
            } else {
                solver.step().map_err(Error::SolverErr)?;
            }
            last_lik = solver.current_assignments().likelihood();
            if last_lik > best_lik {
                best_lik = last_lik;
                best_assgns.clone_from_slice(solver.current_assignments().read_assignments());
            }
            if solver.is_finished() {
                break;
            }
        }
        let dur = update_timer();
        writeln!(dbg_writer, "{}\t{}.{:06}\t{:X}\t{:.8}", prefix, dur.as_secs(), dur.subsec_micros(), seed, last_lik)?;
    }

    let mut assgns = solver.take();
    if last_lik < best_lik {
        assgns.set_assignments(&best_assgns);
    }
    log::debug!("Solved  {} with {}.  ln-likelihood: {:.3}", contigs_group_str, solver_str, best_lik);
    if TypeId::of::<U>() != TypeId::of::<io::Sink>() {
        assgns.write_csv(&prefix, assgn_writer)?;
    }
    Ok(assgns)
}
