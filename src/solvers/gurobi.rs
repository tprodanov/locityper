use std::{fmt, result::Result};
use grb::{
    *,
    expr::{LinExpr, GurobiSum},
};
use crate::{
    Error,
    model::{
        windows::UNMAPPED_WINDOW,
        assgn::ReadAssignment,
    },
};
use super::Solver;

pub struct GurobiSolver {
    model: Model,
    assignments: ReadAssignment,
    is_finished: bool,
    assignment_vars: Vec<Var>,
    seed: Option<i32>,
}

impl GurobiSolver {
    pub fn new(mut assignments: ReadAssignment) -> Result<Self, Error> {
        let contig_windows = assignments.contig_windows();
        let env = Env::new("")?;
        let mut model = Model::with_env(&contig_windows.ids_str(), &env)?;
        model.set_param(parameter::IntParam::OutputFlag, 0)?;
        model.set_param(parameter::IntParam::Threads, 1)?;

        let total_windows = contig_windows.n_windows() as usize;
        let mut window_depth = vec![(0_u32, 0_u32); total_windows];
        let mut window_depth_constrs = vec![LinExpr::new(); total_windows];
        let mut objective = LinExpr::new();
        let mut assignment_vars = Vec::new();
        for rp in 0..assignments.total_reads() {
            let read_alns = assignments.possible_read_alns(rp);
            if read_alns.len() == 1 {
                let loc = &read_alns[0];
                let (w1, w2) = loc.windows();
                window_depth[w1 as usize].0 += 1;
                window_depth[w2 as usize].0 += 1;
                objective.add_constant(loc.ln_prob());
                continue;
            }

            let prev_len = assignment_vars.len();
            for (j, loc) in read_alns.iter().enumerate() {
                let var = add_binvar!(model, name: &format!("R{:x}_{}", rp, j))?;
                assignment_vars.push(var);
                objective.add_term(loc.ln_prob(), var);
                let (w1, w2) = loc.windows();
                let inc = if w1 == w2 { 2 } else { 1 };
                for &w in &[w1, w2] {
                    if w == UNMAPPED_WINDOW {
                        window_depth[w as usize].0 += inc;
                    } else {
                        window_depth[w as usize].1 += inc;
                        window_depth_constrs[w as usize].add_term(f64::from(inc), var);
                    }
                    if inc == 2 {
                        break;
                    }
                }
            }
            model.add_constr(&format!("R{:x}", rp), c!( (&assignment_vars[prev_len..]).grb_sum() == 1 ))?;
        }

        let mut depth_vars = Vec::new();
        for (w, mut depth_constr0) in window_depth_constrs.into_iter().enumerate() {
            let depth_distr = assignments.depth_distr(w);
            let (trivial_reads, non_trivial_reads) = window_depth[w];
            if non_trivial_reads == 0 {
                objective.add_constant(depth_distr.ln_pmf(trivial_reads));
                continue;
            }

            depth_vars.clear();
            for depth_inc in 0..=non_trivial_reads {
                let depth = trivial_reads + depth_inc;
                let var = add_binvar!(model, name: &format!("D{:x}_{}", w, depth))?;
                depth_vars.push(var);
                if depth_inc > 0 {
                    depth_constr0.add_term(-f64::from(depth_inc), var);
                }
                objective.add_term(depth_distr.ln_pmf(depth), var);
            }
            model.add_constr(&format!("D{:x}_base", w), c!( (&depth_vars).grb_sum() == 1 ))?;
            model.add_constr(&format!("D{:x}_eq", w), c!( depth_constr0 == 0 ))?;
        }
        model.set_objective(objective, grb::ModelSense::Maximize)?;
        // model.write("model.lp")?;

        assignments.init_assignments(|_| 0);
        Ok(Self {
            model, assignments, assignment_vars,
            is_finished: false,
            seed: None,
        })
    }

    /// Query read assignments from the ILP solution, and set them in the `self.assignments`.
    fn set_assignments(&mut self) -> Result<(), Error> {
        let vals = self.model.get_obj_attr_batch(attr::X, self.assignment_vars.iter().cloned())?;
        let mut i = 0;
        self.assignments.try_init_assignments::<_, Error>(|locs| {
            let j = i + locs.len();
            let new_assgn = vals[i..j].iter().position(|&v| v >= 0.9999 && v <= 1.0001)
                .ok_or_else(|| Error::solver("Gurobi", "Cannot identify read assignment (output is not 0/1)"))?;
            i = j;
            Ok(new_assgn)
        })?;
        assert_eq!(i, vals.len(), "Numbers of total read assignments do not match");
        Ok(())
    }
}

impl Solver for GurobiSolver {
    fn is_seedable() -> bool { true }

    fn set_seed(&mut self, seed: u64) -> Result<(), Error> {
        self.seed = Some((seed % 0x7fffffff) as i32);
        Ok(())
    }

    /// Resets the solver.
    fn reset(&mut self) -> Result<(), Error> {
        self.is_finished = false;
        self.model.reset()?;
        self.model.set_param(parameter::IntParam::OutputFlag, 0)?;
        self.model.set_param(parameter::IntParam::Threads, 1)?;
        if let Some(seed) = self.seed.take() {
            self.model.set_param(parameter::IntParam::Seed, seed)?;
        }
        Ok(())
    }

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> Result<(), Error> {
        self.model.optimize()?;
        self.set_assignments()?;
        let status = self.model.status()?;
        if status != Status::Optimal {
            log::error!("Gurobi achieved non-optimal status {:?}", status);
        }
        let ilp_lik = self.model.get_attr(attr::ObjVal)?;
        let assgn_lik = self.assignments.likelihood();
        if (ilp_lik - assgn_lik).abs() > 1e-5 {
            log::error!("Gurobi likehood differs from the model likelihood: {} and {}", ilp_lik, assgn_lik);
        }
        self.is_finished = true;
        Ok(())
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

impl fmt::Display for GurobiSolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Gurobi")
    }
}