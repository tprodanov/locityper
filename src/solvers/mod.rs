use crate::{
    algo::vec_ext::IterExt,
    model::{
        locs::PairAlignment,
        assgn::ReadAssignment,
    },
};

pub mod greedy;

/// Trait that distributes the reads between their possible alignments
pub trait Solver<'a> {
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

/// Initialize reads based on the first alignment in the list (not necessarily best, partially random).
fn init_first(_: &[PairAlignment]) -> usize { 0 }

/// Initialize reads based on their best alignment.
fn init_best(possible_alns: &[PairAlignment]) -> usize {
    IterExt::argmax(possible_alns.iter().map(PairAlignment::ln_prob)).0
}

/// Distribute read assignment in at most `max_iter` iterations.
pub fn solve<'a, 'b, S>(mut solver: S, max_iter: usize) -> ReadAssignment<'a>
// 'a - lifetime of ReadAssignments, 'b - lifetime of Solver, 'a outlives 'b.
where 'a: 'b,
    S: Solver<'a> + 'b,
{
    let mut best_lik = solver.initialize();
    let mut last_lik = best_lik;
    let mut best_assgns: Vec<u16> = solver.current_assignments().read_assignments().to_vec();
    for _ in 0..max_iter {
        solver.step();
        last_lik = solver.current_assignments().likelihood();
        if last_lik > best_lik {
            best_lik = last_lik;
            best_assgns.clone_from_slice(solver.current_assignments().read_assignments());
        }
        if solver.is_finished() {
            break;
        }
    }
    let mut model = solver.finish();
    if last_lik < best_lik {
        model.set_assignments(&best_assgns);
    }
    model
}
