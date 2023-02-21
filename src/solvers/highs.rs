use std::fmt;
use highs::{RowProblem, Model, Col, Sense, HighsModelStatus as Status};
use crate::model::assgn::{ReadAssignment, UNMAPPED_WINDOW};
use super::Solver;

pub struct HighsSolver {
    model: Option<Model>,
    #[allow(dead_code)]
    const_term: f64,
    assignments: ReadAssignment,
    solver_type: String,
}

impl HighsSolver {
    /// Creates a new HiGHS solver. Possible `solver_type` values: `choose`, `simplex`, or `ipm`.
    pub fn new(solver_type: String, mut assignments: ReadAssignment) -> Self {
        let contig_windows = assignments.contig_windows();
        let total_windows = contig_windows.total_windows() as usize;
        // Number of trivial and non-trivial reads mapped to a window.
        let mut window_depth = vec![(0_u32, 0_u32); total_windows];
        let mut window_depth_constrs: Vec<Vec<(Col, f64)>> = vec![Vec::new(); total_windows];

        let mut problem = RowProblem::default();
        let mut const_term = 0.0;
        let mut assgn_constr = Vec::new();
        for rp in 0..assignments.total_reads() {
            let read_alns = assignments.possible_read_alns(rp);
            if read_alns.len() == 1 {
                let loc = &read_alns[0];
                let (w1, w2) = contig_windows.get_pair_window_ixs(loc);
                window_depth[w1 as usize].0 += 1;
                window_depth[w2 as usize].0 += 1;
                const_term += loc.ln_prob();
                continue;
            }

            for loc in read_alns.iter() {
                let var = problem.add_integer_column(loc.ln_prob(), 0.0..=1.0);
                assgn_constr.push((var, 1.0));
                let (w1, w2) = contig_windows.get_pair_window_ixs(loc);
                let inc = if w1 == w2 { 2 } else { 1 };
                for &w in &[w1, w2] {
                    if w == UNMAPPED_WINDOW {
                        window_depth[w as usize].0 += inc;
                    } else {
                        window_depth[w as usize].1 += inc;
                        window_depth_constrs[w as usize].push((var, f64::from(inc)));
                    }
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
            let depth_distr = assignments.depth_distr(w);
            let (trivial_reads, non_trivial_reads) = window_depth[w];
            if non_trivial_reads == 0 {
                const_term += depth_distr.ln_pmf(trivial_reads);
                continue;
            }

            for depth_inc in 0..=non_trivial_reads {
                let var = problem.add_integer_column(depth_distr.ln_pmf(trivial_reads + depth_inc), 0.0..=1.0);
                depth_constr1.push((var, 1.0));
                if depth_inc > 0 {
                    depth_constr0.push((var, -f64::from(depth_inc)));
                }
            }
            problem.add_row(0.0..=0.0, &depth_constr0);
            problem.add_row(1.0..=1.0, &depth_constr1);
            depth_constr1.clear();
        }

        assignments.init_assignments(|_| 0);
        let mut model = problem.optimise(Sense::Maximise);
        model.set_option("parallel", "off");
        model.set_option("presolve", "on");
        model.set_option("solver", &solver_type as &str);
        model.make_quiet();

        Self {
            model: Some(model),
            const_term, assignments, solver_type,
        }
    }

    /// Query read assignments from the ILP solution, and set them in the `self.assignments`.
    fn set_assignments(&mut self, vals: &[f64]) {
        let mut i = 0;
        self.assignments.init_assignments(|locs| {
            let j = i + locs.len();
            // Take argmax because HiGHS does not always output reasonable solutions,
            // this way we always have a read assignment.
            let new_assgn = crate::algo::vec_ext::F64Ext::argmax(&vals[i..j]).0;
            i = j;
            new_assgn
        });
    }
}

impl Solver for HighsSolver {
    type Error = ();

    fn is_seedable() -> bool { false }

    fn set_seed(&mut self, _seed: u64) -> Result<(), Self::Error> {
        // Even so it is possible to set seed to Highs Solver, results are always the same.
        panic!("Cannot set seed to HiGHS");
    }

    /// Resets and initializes anew read assignments.
    fn reset(&mut self) -> Result<(), Self::Error> {
        Ok(())
    }

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> Result<(), Self::Error> {
        let solved_model = self.model.take().expect("Cannot run HighsSolver::step twice!").solve();
        if solved_model.status() != Status::Optimal {
            panic!("Highs model finished with non-optimal status {:?}", solved_model.status());
        }
        let solution = solved_model.get_solution();
        self.set_assignments(solution.columns());
        Ok(())
    }

    /// Returns true if the solver is finished.
    fn is_finished(&self) -> bool {
        self.model.is_none()
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

impl fmt::Display for HighsSolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "HiGHS({})", self.solver_type)
    }
}
