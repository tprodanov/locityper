use highs::{RowProblem, Col, Model, Sense, HighsModelStatus as Status};
use crate::model::assgn::{ReadAssignment, UNMAPPED_WINDOW};
use super::Solver;

pub struct HighsSolver {
    model: Option<Model>,
    #[allow(dead_code)]
    const_term: f64,
    assignments: ReadAssignment,
}

impl HighsSolver {
    pub fn new(assignments: ReadAssignment) -> Self {
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

        log::debug!("Try optimize");
        let mut model = problem.try_optimise(Sense::Maximise).unwrap();
        log::debug!("Constructed model");
        model.set_option("parallel", "off");
        Self {
            model: Some(model),
            const_term, assignments,
        }
    }

    /// Query read assignments from the ILP solution, and set them in the `self.assignments`.
    fn set_assignments(&mut self, vals: &[f64]) {
        let mut i = 0;
        self.assignments.init_assignments(|locs| {
            let j = i + locs.len();
            let new_assgn = vals[i..j].iter()
                .position(|&v| v >= 0.99999 && v <= 1.00001)
                .expect("Read has no assignment!");
            i = j;
            new_assgn
        });
    }
}

impl Solver for HighsSolver {
    type Error = ();

    fn is_seedable() -> bool { false }

    fn set_seed(&mut self, _seed: u64) -> Result<(), Self::Error> {
        panic!("Cannot set seed to the Highs solver!")
    }

    /// Resets and initializes anew read assignments.
    fn initialize(&mut self) -> Result<(), Self::Error> {
        assert!(self.model.is_some(), "Cannot run `HighsSolver::initialize` twice!");
        self.assignments.init_assignments(|_| 0);
        Ok(())
    }

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> Result<f64, Self::Error> {
        let solved_model = self.model.take().expect("Cannot run `HighsSolver::step` twice!").solve();
        if solved_model.status() != Status::Optimal {
            panic!("Highs model finished with non-optimal status {:?}", solved_model.status());
        }
        let old_lik = self.assignments.likelihood();
        let solution = solved_model.get_solution();
        self.set_assignments(solution.columns());
        let new_lik = self.assignments.likelihood();
        Ok(new_lik - old_lik)
    }

    /// Returns true if the solver is finished.
    fn is_finished(&self) -> bool {
        self.model.is_none()
        // unimplemented!()
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
