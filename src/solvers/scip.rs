use std::fmt;
use crate::{
    err::Error,
    model::assgn::{GenotypeAlignments, ReadAssignment},
    ext::rand::XoshiroRng,
};

/// ILP Solver using HiGHS library.
#[derive(Clone)]
pub struct ScipSolver;

impl Default for ScipSolver {
    fn default() -> Self {
        Self
    }
}

impl ScipSolver {
    /// Defines SCIP model and returns the model and vector of assignment variables.
    fn define_model(
        &self,
        gt_alns: &GenotypeAlignments
    ) -> (russcip::model::Model<russcip::model::ProblemCreated>, Vec<russcip::variable::Variable>)
    {
        use russcip::prelude::*;
        let total_windows = gt_alns.total_windows();
        // Number of trivial and non-trivial reads mapped to a window.
        let mut window_depth = vec![(0_u32, 0_u32); total_windows];
        let mut window_depth_constrs: Vec<Vec<(f64, Expr)>> = vec![Vec::new(); total_windows];

        let (depth_contrib, aln_contrib) = gt_alns.contributions();
        // let mut const_term = 0.0;
        let mut model = Model::default()
            .hide_output()
            .maximize();
        let mut assgn_constr = Vec::new();
        let mut assgn_vars = Vec::new();
        for rp in 0..gt_alns.total_reads() {
            let read_alns = gt_alns.possible_read_alns(rp);
            if read_alns.len() == 1 {
                for w in read_alns[0].windows().into_iter() {
                    window_depth[w as usize].0 += 1;
                }
                // const_term += aln_contrib * loc.ln_prob();
                continue;
            }

            for loc in read_alns.iter() {
                let var = model.add(var().obj(aln_contrib * loc.ln_prob()).bin());
                assgn_constr.push(Expr::var(&var));
                let ws = loc.windows();
                let inc: u32 = if ws[0] == ws[1] { 2 } else { 1 };
                for w in ws.into_iter() {
                    window_depth[w as usize].1 += inc;
                    window_depth_constrs[w as usize].push((f64::from(inc), Expr::var(&var)));
                    if inc == 2 {
                        break;
                    }
                }
                assgn_vars.push(var);
            }
            model.add(cons().expression(Expr::sum(assgn_constr.drain(..))).eq(1.0));
            assgn_constr.clear();
        }

        let mut depth_constr1 = Vec::new();
        for (w, mut depth_constr0) in window_depth_constrs.into_iter().enumerate() {
            let depth_distr = gt_alns.depth_distr(w);
            if depth_distr.is_trivial() {
                continue;
            }
            let (trivial_reads, non_trivial_reads) = window_depth[w];
            if non_trivial_reads == 0 {
                // const_term += depth_contrib * depth_distr.ln_pmf(trivial_reads);
                continue;
            }

            for depth_inc in 0..=non_trivial_reads {
                let var = model.add(
                    var().obj(depth_contrib * depth_distr.ln_prob(trivial_reads + depth_inc)).bin());
                depth_constr1.push(Expr::var(&var));
                if depth_inc > 0 {
                    depth_constr0.push((-f64::from(depth_inc), Expr::var(&var)));
                }
            }
            model.add(cons().expression(Expr::sum_weighted(depth_constr0.into_iter())).eq(0.0));
            model.add(cons().expression(Expr::sum(depth_constr1.drain(..))).eq(1.0));
            depth_constr1.clear();
        }
        (model, assgn_vars)
    }
}

impl super::Solver for ScipSolver {
    /// Distribute reads between several haplotypes in a best way.
    fn solve_nontrivial<'a>(
        &self,
        gt_alns: &'a GenotypeAlignments,
        _rng: &mut XoshiroRng
    ) -> crate::Result<ReadAssignment<'a>> {
        use russcip::WithSolutions;
        let (model, assgn_vars) = self.define_model(gt_alns);
        let solution = model.solve();
        let status = solution.status();
        if status != russcip::status::Status::Optimal {
            return Err(Error::solver("SCIP",
                format!("Model finished with non-optimal status {:?}", status)));
        }

        let sol = solution.best_sol().expect("Could not get an ILP solution");
        let mut i = 0;
        ReadAssignment::try_new(gt_alns, |locs| {
            let j = i + locs.len();
            let new_assgn = assgn_vars[i..j].iter().position(|var| {
                let v = sol.val(var);
                v >= 0.9999 && v <= 1.0001
            }).ok_or_else(|| Error::solver("Gurobi", "Cannot identify read assignment (output is not 0/1)"))?;
            i = j;
            Ok(new_assgn)
        })
    }

    fn describe_params(&self) -> String {
        String::new()
    }
}

impl super::SetParams for ScipSolver {
    fn set_param(&mut self, _key: &str, _val: &str) -> Result<(), super::ParamErr> {
        Err(super::ParamErr::Unknown)
    }
}

impl fmt::Display for ScipSolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SCIP")
    }
}
