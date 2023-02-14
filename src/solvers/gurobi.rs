use std::{
    ops::Deref,
    cmp::min,
    result::Result,
};
use grb::*;
use grb::expr::GurobiSum;
use crate::model::assgn::{ReadAssignment, UNMAPPED_WINDOW, INIT_WSHIFT, split_depth};
use super::Solver;

pub struct GurobiSolver {
    model: Model,
    assignments: ReadAssignment,
    is_finished: bool,
}

impl GurobiSolver {
    pub fn new(assgns: ReadAssignment, depth_step: u32) -> Result<Self, grb::Error> {
        let mut model = Model::new(&format!("{}", assgns.contigs_group()))?;
        let contig_windows = assgns.contig_windows();
        let total_windows = contig_windows.total_windows() as usize;

        let mut by_window_vars: Vec<Vec<Var>> = vec![Vec::new(); total_windows];
        let mut window_trivial_depth = vec![0; total_windows];
        // let mut assignment_vars = Vec::new();
        // let mut assignment_coefs = Vec::new();
        // let mut depth_vars = Vec::new();
        // let mut depth_coefs = Vec::new();

        let mut objective = expr::LinExpr::new();
        let mut assignment_vars = Vec::new();
        for rp in 0..assgns.total_reads() {
            let read_locs = assgns.read_locs(rp);
            if read_locs.len() == 1 {
                let loc = &read_locs[0];
                let (w1, w2) = contig_windows.get_pair_window_ixs(loc);
                window_trivial_depth[w1 as usize] += 1;
                window_trivial_depth[w2 as usize] += 1;
                objective.add_constant(loc.ln_prob());
                continue;
            }

            assignment_vars.clear();
            for (j, loc) in read_locs.iter().enumerate() {
                // TODO: Remove names later.
                let var = add_binvar!(model, name: &format!("R{}_{}", rp, j))?;
                assignment_vars.push(var);
                objective.add_term(loc.ln_prob(), var);
                let (w1, w2) = contig_windows.get_pair_window_ixs(loc);
                if w1 != UNMAPPED_WINDOW {
                    by_window_vars[w1 as usize].push(var);
                }
                if w2 != UNMAPPED_WINDOW {
                    by_window_vars[w2 as usize].push(var);
                }
            }
            model.add_constr(&format!("R{}", rp), c!( (&assignment_vars).grb_sum() == 1 ))?;
        }

        let mut depth_vars = Vec::new();
        for w in INIT_WSHIFT as usize..total_windows {
            let depth_distr = assgns.depth_distr(w);
            let depth_bound = assgns.depth_bound(w);
            let bin_probs = split_depth(depth_distr.deref(), depth_step, depth_bound);
            let nbins = bin_probs.len();

            depth_vars.clear();
            let wsum: Expr = window_trivial_depth[w] + (&by_window_vars[w]).grb_sum();
            let mut constr1 = wsum.clone().into_linexpr()?;
            let mut constr2 = wsum.into_linexpr()?;
            for (bin, &prob) in bin_probs.iter().enumerate() {
                let min_depth;
                let max_depth;
                let var = if bin == nbins - 1 {
                    min_depth = depth_bound + 1;
                    max_depth = u32::MAX;
                    add_binvar!(model, name: &format!("D{}_{}..", w, min_depth))?
                } else {
                    min_depth = bin as u32 * depth_step;
                    max_depth = min((bin as u32 + 1) * depth_step - 1, depth_bound);
                    add_binvar!(model, name: &format!("D{}_{}..{}", w, min_depth, max_depth + 1))?
                };
                depth_vars.push(var);
                objective.add_term(prob, var);
                constr1.add_term(-f64::from(min_depth), var);
                constr2.add_term(-f64::from(max_depth), var);
            }
            model.add_constr(&format!("D{}_base", w), c!( (&depth_vars).grb_sum() == 1 ))?;
            model.add_constr(&format!("D{}_min", w), c!( constr1 >= 0 ))?;
            model.add_constr(&format!("D{}_max", w), c!( constr2 <= 0 ))?;
        }
        model.set_objective(objective, grb::ModelSense::Maximize)?;
        model.write("model.lp")?;

        unimplemented!()
        // Self {
        //     model,
        //     assgns,
        //     is_finished: false,
        // }
    }
}

impl Solver for GurobiSolver {
    type Error = grb::Error;

    fn is_seedable() -> bool { true }

    fn set_seed(&mut self, seed: u64) -> Result<(), Self::Error> {
        self.model.set_param(parameter::IntParam::Seed, seed as i32)
    }

    /// Resets and initializes anew read assignments.
    fn initialize(&mut self) -> Result<(), Self::Error> {
        self.is_finished = false;
        self.assignments.init_assignments(|_| 0);
        self.model.reset()
    }

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> Result<f64, Self::Error> {
        unimplemented!()
    }

    /// Returns true if the solver is finished.
    fn is_finished(&self) -> bool {
        self.is_finished
    }

    /// Return the current read assignments.
    fn current_assignments(&self) -> &ReadAssignment {
        &self.assignments
    }

    /// Recalculate likelihood and check if it matches the stored one.
    fn recalculate_likelihood(&mut self) {
        self.assignments.recalc_likelihood();
    }

    /// Consumes solver and returns the read assignments.
    fn take(self) -> ReadAssignment {
        self.assignments
    }
}