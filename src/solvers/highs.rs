use std::fmt;
use highs::{RowProblem, Col, Sense, HighsModelStatus as Status};
use crate::{
    Error,
    model::{
        windows::REG_WINDOW_SHIFT,
        assgn::ReadAssignment,
    },
    ext::vec::F64Ext,
    bg::ser::json_get,
    ext::rand::XoshiroRng,
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
    pub fn set_type(&mut self, mode: &str) -> &mut Self {
        assert!(mode == "choose" || mode == "simplex" || mode == "ipm",
            "Unknown {} solver type `{}`", HIGHS_NAME, mode);
        self.mode = mode.to_owned();
        self
    }

    fn define_model(&self, assignments: &mut ReadAssignment) -> RowProblem {
        let total_windows = assignments.contig_windows().total_windows() as usize;
        // Number of trivial and non-trivial reads mapped to a window.
        let mut window_depth = vec![(0_u32, 0_u32); total_windows];
        let mut window_depth_constrs: Vec<Vec<(Col, f64)>> = vec![Vec::new(); total_windows];

        // let mut const_term = 0.0;
        let mut problem = RowProblem::default();
        let mut assgn_constr = Vec::new();
        for rp in 0..assignments.total_reads() {
            let read_alns = assignments.possible_read_alns(rp);
            if read_alns.len() == 1 {
                let loc = &read_alns[0];
                let ((w1s, w1e), (w2s, w2e)) = loc.windows();
                for w in (w1s..w1e).chain(w2s..w2e) {
                    window_depth[w as usize].0 += 1;
                }
                // const_term += loc.ln_prob();
                continue;
            }

            for loc in read_alns.iter() {
                let var = problem.add_integer_column(loc.ln_prob(), 0.0..=1.0);
                assgn_constr.push((var, 1.0));
                let ((w1s, w1e), (w2s, w2e)) = loc.windows();
                for w in (w1s..w1e).chain(w2s..w2e) {
                    if w < REG_WINDOW_SHIFT {
                        window_depth[w as usize].0 += 1;
                    } else {
                        window_depth[w as usize].1 += 1;
                        window_depth_constrs[w as usize].push((var, 1.0));
                    }
                }
            }
            problem.add_row(1.0..=1.0, &assgn_constr);
            assgn_constr.clear();
        }

        let depth_contrib = assignments.depth_contrib();
        let mut depth_constr1 = Vec::new();
        for (w, mut depth_constr0) in window_depth_constrs.into_iter().enumerate() {
            let depth_distr = assignments.depth_distr(w);
            let (trivial_reads, non_trivial_reads) = window_depth[w];
            if non_trivial_reads == 0 {
                // const_term += depth_distr.ln_pmf(trivial_reads);
                continue;
            }

            for depth_inc in 0..=non_trivial_reads {
                let var = problem.add_integer_column(
                    depth_contrib * depth_distr.ln_pmf(trivial_reads + depth_inc), 0.0..=1.0);
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

    /// Query read assignments from the ILP solution, and set them to `assignments`.
    fn set_assignments(&self, assignments: &mut ReadAssignment, vals: &[f64]) -> f64 {
        let mut i = 0;
        assignments.init_assignments(|locs| {
            let j = i + locs.len();
            // Take argmax because HiGHS does not always output reasonable solutions,
            // this way we always have a read assignment.
            let new_assgn = F64Ext::argmax(&vals[i..j]).0;
            i = j;
            new_assgn
        })
    }
}

impl super::Solver for HighsSolver {
    /// Distribute reads between several haplotypes in a best way.
    fn solve(&self, assignments: &mut ReadAssignment, _rng: &mut XoshiroRng) -> Result<f64, Error> {
        let problem = self.define_model(assignments);
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
        let solution = solved_model.get_solution();

        let old_lik = assignments.likelihood();
        if old_lik.is_finite() {
            let old_assgns = assignments.read_assignments().to_vec();
            let new_lik = self.set_assignments(assignments, solution.columns());
            if new_lik < old_lik {
                // Previous solution was better.
                Ok(assignments.set_assignments(&old_assgns))
            } else {
                Ok(new_lik)
            }
        } else {
            Ok(self.set_assignments(assignments, solution.columns()))
        }
    }
}

impl super::SetParams for HighsSolver {
    fn set_params(&mut self, obj: &json::JsonValue) -> Result<(), Error> {
        json_get!(obj -> mode? (as_str));
        if let Some(mode) = mode {
            self.set_type(mode);
        }
        Ok(())
    }
}

impl fmt::Display for HighsSolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}-{}", HIGHS_NAME, self.mode)
    }
}
