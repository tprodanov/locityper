use std::fmt;
use highs::{RowProblem, Col, Sense, HighsModelStatus as Status};
use crate::{
    err::{Error, validate_param},
    model::assgn::{GenotypeAlignments, ReadAssignment},
    ext::{
        vec::F64Ext,
        rand::XoshiroRng,
    },
};

const HIGHS_NAME: &'static str = "HiGHS";

/// ILP Solver using HiGHS library.
#[derive(Clone)]
pub struct HighsSolver {
    /// Solver type: possible values are `choose`, `simplex` and `ipm`.
    mode: String,
}

impl Default for HighsSolver {
    fn default() -> Self {
        Self {
            mode: "simplex".to_owned(),
        }
    }
}

impl HighsSolver {
    pub fn set_mode(&mut self, mode: &str) -> crate::Result<()> {
        validate_param!(mode == "choose" || mode == "simplex" || mode == "ipm",
            "{} solver: unknown mode {:?}", HIGHS_NAME, mode);
        self.mode = mode.to_owned();
        Ok(())
    }

    fn define_model(&self, gt_alns: &GenotypeAlignments) -> RowProblem {
        let total_windows = gt_alns.total_windows();
        // Number of trivial and non-trivial reads mapped to a window.
        let mut window_depth = vec![(0_u32, 0_u32); total_windows];
        let mut window_depth_constrs: Vec<Vec<(Col, f64)>> = vec![Vec::new(); total_windows];

        let (depth_contrib, aln_contrib) = gt_alns.contributions();
        // let mut const_term = 0.0;
        let mut problem = RowProblem::default();
        let mut assgn_constr = Vec::new();
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
                let var = problem.add_integer_column(aln_contrib * loc.ln_prob(), 0.0..=1.0);
                assgn_constr.push((var, 1.0));
                let ws = loc.windows();
                let inc: u32 = if ws[0] == ws[1] { 2 } else { 1 };
                for w in ws.into_iter() {
                    window_depth[w as usize].1 += inc;
                    window_depth_constrs[w as usize].push((var, f64::from(inc)));
                    if inc == 2 {
                        break;
                    }
                }
            }
            problem.add_row(1.0..=1.0, &assgn_constr);
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
                let var = problem.add_integer_column(
                    depth_contrib * depth_distr.ln_prob(trivial_reads + depth_inc), 0.0..=1.0);
                depth_constr1.push((var, 1.0));
                if depth_inc > 0 {
                    depth_constr0.push((var, -f64::from(depth_inc)));
                }
            }
            problem.add_row(0.0..=0.0, &depth_constr0);
            problem.add_row(1.0..=1.0, &depth_constr1);
            depth_constr1.clear();
        }
        problem
    }
}

impl super::Solver for HighsSolver {
    fn name(&self) -> &'static str {
        "HiGHS ILP"
    }

    /// Distribute reads between several haplotypes in a best way.
    fn solve_nontrivial<'a>(
        &self,
        gt_alns: &'a GenotypeAlignments,
        _rng: &mut XoshiroRng
    ) -> crate::Result<ReadAssignment<'a>> {
        let problem = self.define_model(gt_alns);
        let mut model = problem.optimise(Sense::Maximise);
        model.set_option("parallel", "off");
        model.set_option("presolve", "on");
        model.set_option("solver", &self.mode as &str);
        model.make_quiet();

        let solved_model = model.solve();
        if solved_model.status() != Status::Optimal {
            return Err(Error::solver(HIGHS_NAME,
                format!("Model finished with non-optimal status {:?}", solved_model.status())));
        }

        let sol = solved_model.get_solution();
        let sol = sol.columns();
        let mut i = 0;
        Ok(ReadAssignment::new(gt_alns, |locs| {
            let j = i + locs.len();
            // Take argmax because HiGHS does not always output reasonable solutions,
            // this way we always have a read assignment.
            let new_assgn = F64Ext::argmax(&sol[i..j]).0;
            i = j;
            new_assgn
        }))
    }
}

impl super::SetParams for HighsSolver {
    fn set_param(&mut self, key: &str, val: &str) -> crate::Result<()> {
        if key == "mode" || key == "type" {
            self.set_mode(val)?;
        } else {
            log::error!("{} solver: unknown parameter {:?}", HIGHS_NAME, key);
        }
        Ok(())
    }
}

impl fmt::Display for HighsSolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}-{}", HIGHS_NAME, self.mode)
    }
}
