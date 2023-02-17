use std::result::Result;
use grb::*;
use grb::expr::GurobiSum;
use crate::model::assgn::{ReadAssignment, UNMAPPED_WINDOW, INIT_WSHIFT};
use super::Solver;

pub struct GurobiSolver {
    model: Model,
    assignments: ReadAssignment,
    is_finished: bool,
    assignment_vars: Vec<Var>,
}

impl GurobiSolver {
    pub fn new(mut assignments: ReadAssignment) -> Result<Self, grb::Error> {
        assignments.init_assignments(|_| 0);
        let mut model = Model::new(&format!("{}", assignments.contigs_group()))?;
        model.set_param(parameter::IntParam::Threads, 0)?;
        model.set_param(parameter::IntParam::LogToConsole, 0)?;

        let contig_windows = assignments.contig_windows();
        let total_windows = contig_windows.total_windows() as usize;
        let mut by_window_vars: Vec<Vec<Var>> = vec![Vec::new(); total_windows];
        let mut window_trivial_depth = vec![0; total_windows];

        let mut objective = expr::LinExpr::new();
        let mut assignment_vars = Vec::new();
        for rp in 0..assignments.total_reads() {
            let read_locs = assignments.read_locs(rp);
            if read_locs.len() == 1 {
                let loc = &read_locs[0];
                let (w1, w2) = contig_windows.get_pair_window_ixs(loc);
                window_trivial_depth[w1 as usize] += 1;
                window_trivial_depth[w2 as usize] += 1;
                objective.add_constant(loc.ln_prob());
                continue;
            }

            let prev_len = assignment_vars.len();
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
            model.add_constr(&format!("R{}", rp), c!( (&assignment_vars[prev_len..]).grb_sum() == 1 ))?;
        }

        let mut depth_vars = Vec::new();
        for w in INIT_WSHIFT as usize..total_windows {
            let depth_distr = assignments.depth_distr(w);
            let potential_min_depth = window_trivial_depth[w] as u32;
            let potential_max_depth = potential_min_depth + by_window_vars[w].len() as u32;
            if potential_min_depth == potential_max_depth {
                objective.add_constant(depth_distr.ln_pmf(potential_min_depth));
                continue;
            }

            depth_vars.clear();
            let mut constr = (&by_window_vars[w]).grb_sum().into_linexpr()?;
            for depth in potential_min_depth..=potential_max_depth {
                let var = add_binvar!(model, name: &format!("D{}_{}", w, depth))?;
                depth_vars.push(var);
                let prob = depth_distr.ln_pmf(depth);
                if depth > potential_min_depth {
                    constr.add_term(-f64::from(depth - potential_min_depth), var);
                }
                objective.add_term(prob, var);
            }
            model.add_constr(&format!("D{}_base", w), c!( (&depth_vars).grb_sum() == 1 ))?;
            model.add_constr(&format!("D{}_eq", w), c!( constr == 0 ))?;
        }
        model.set_objective(objective, grb::ModelSense::Maximize)?;
        model.write("model.lp")?;

        Ok(Self {
            model, assignments, assignment_vars,
            is_finished: false,
        })
    }

    /// Query read assignments from the ILP solution, and set them in the `self.assignments`.
    fn set_assignments(&mut self) -> Result<(), grb::Error> {
        let vals = self.model.get_obj_attr_batch(attr::X, self.assignment_vars.iter().cloned())?;
        let mut i = 0;
        self.assignments.init_assignments(|locs| {
            let j = i + locs.len();
            let new_assgn = vals[i..j].iter().position(|&v| v == 1.0).expect("Read has no assignment!");
            i = j;
            new_assgn
        });
        assert_eq!(i, vals.len(), "Numbers of total read assignments do not match");
        Ok(())
    }
}

impl Solver for GurobiSolver {
    type Error = grb::Error;

    fn is_seedable() -> bool { true }

    fn set_seed(&mut self, seed: u64) -> Result<(), Self::Error> {
        self.model.set_param(parameter::IntParam::Seed, (seed % 0x7fffffff) as i32)
    }

    /// Resets and initializes anew read assignments.
    fn initialize(&mut self) -> Result<(), Self::Error> {
        self.is_finished = false;
        self.model.reset()
    }

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> Result<f64, Self::Error> {
        let old_lik = self.assignments.likelihood();
        self.model.optimize()?;
        let ilp_lik = self.model.get_attr(attr::ObjVal)?;
        self.set_assignments()?;
        let new_lik = self.assignments.likelihood();
        let status = self.model.status()?;
        if status != Status::Optimal {
            log::error!("Gurobi achieved non-optimal status {:?}", status);
        }
        if (ilp_lik - new_lik).abs() > 1e-5 {
            log::error!("Gurobi likehood differs from the model likelihood: {} and {}", ilp_lik, new_lik);
        }
        self.is_finished = true;
        Ok(new_lik - old_lik)
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