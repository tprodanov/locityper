use std::{
    io,
    time::Instant,
    any::TypeId,
    fmt::Display,
};
use crate::{
    algo::vec_ext::IterExt,
    model::{
        windows::ReadWindows,
        assgn::ReadAssignment,
        locs::AllPairAlignments,
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
use crate::Error;

/// Trait that distributes the reads between their possible alignments
pub trait Solver: Display + Display {
    /// Returns true if the solver can take seed.
    fn is_seedable() -> bool;

    /// Sets seed.
    /// Returns error if the seed does not fit the model, or if the solver is deterministic.
    fn set_seed(&mut self, seed: u64) -> Result<(), Error>;

    /// Resets the solver.
    fn reset(&mut self) -> Result<(), Error>;

    /// Perform one iteration.
    fn step(&mut self) -> Result<(), Error>;

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
fn init_best(possible_alns: &[ReadWindows]) -> usize {
    IterExt::argmax(possible_alns.iter().map(ReadWindows::ln_prob)).0
}

/// Distribute read assignment in at most `max_iters` iterations.
///
/// dbg_writer writes intermediate likelihood values and runtime for each seed.
/// Format: `assignment-tag  solver  runtime  seed  likelihood`.
pub fn solve<S, F, I, W, U, V>(
    assgns: ReadAssignment,
    build: F,
    seeds: I,
    all_alns: &AllPairAlignments,
    dbg_writer: &mut W,
    depth_writer: &mut U,
    reads_writer: &mut V,
) -> Result<ReadAssignment, Error>
where S: Solver,
      F: Fn(ReadAssignment) -> S,
      I: Iterator<Item = u64>,
      W: io::Write,
      U: io::Write + 'static,
      V: io::Write + 'static,
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

    let contigs_str = assgns.contig_windows().ids_str();
    let mut solver = build(assgns);
    let solver_str = solver.to_string();
    log::debug!("    Solving {} with {}.", contigs_str, solver_str);
    let dur = update_timer();
    let prefix = format!("{}\t{}", contigs_str, solver_str);
    writeln!(dbg_writer, "{}\t{}.{:06}\tstart\tNA", prefix, dur.as_secs(), dur.subsec_micros())?;
    let mut outer = 0;
    for mut seed in seeds {
        if S::is_seedable() {
            solver.set_seed(seed)?;
        } else if outer == 0 {
            seed = 0;
        } else {
            break;
        }

        outer += 1;
        for inner in 0.. {
            if inner == 0 {
                solver.reset()?;
            } else {
                solver.step()?;
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
    log::debug!("    Solved  {} with {}.  ln-likelihood: {:.3}", contigs_str, solver_str, best_lik);
    if TypeId::of::<U>() != TypeId::of::<io::Sink>() {
        assgns.write_depth(&prefix, depth_writer)?;
    }
    if TypeId::of::<V>() != TypeId::of::<io::Sink>() {
        assgns.write_reads(&prefix, reads_writer, all_alns)?;
    }
    Ok(assgns)
}
