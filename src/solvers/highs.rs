use highs::{RowProblem, Col, Model, Sense, HighsModelStatus as Status};
use crate::model::assgn::{ReadAssignment, UNMAPPED_WINDOW, INIT_WSHIFT};
use super::Solver;

pub struct HighsSolver {
    model: Option<Model>,
    #[allow(dead_code)]
    const_term: f64,
    assignments: ReadAssignment,
}

impl HighsSolver {
    pub fn new(assignments: ReadAssignment) -> Self {
        let mut problem = RowProblem::new();

        let contig_windows = assignments.contig_windows();
        let total_windows = contig_windows.total_windows() as usize;
        let mut by_window_vars: Vec<Vec<Col>> = vec![Vec::new(); total_windows];
        let mut window_trivial_depth = vec![0; total_windows];

        let mut const_term = 0.0;
        let mut assgn_constr = Vec::new();
        for rp in 0..assignments.total_reads() {
            let read_locs = assignments.read_locs(rp);
            if read_locs.len() == 1 {
                let loc = &read_locs[0];
                let (w1, w2) = contig_windows.get_pair_window_ixs(loc);
                window_trivial_depth[w1 as usize] += 1;
                window_trivial_depth[w2 as usize] += 1;
                const_term += loc.ln_prob();
                continue;
            }

            assgn_constr.clear();
            for loc in read_locs.iter() {
                let var = problem.add_integer_column(loc.ln_prob(), 0.0_f64..=1.0);
                assgn_constr.push((var, 1.0));
                let (w1, w2) = contig_windows.get_pair_window_ixs(loc);
                if w1 != UNMAPPED_WINDOW {
                    by_window_vars[w1 as usize].push(var);
                }
                if w2 != UNMAPPED_WINDOW {
                    by_window_vars[w2 as usize].push(var);
                }
            }
            problem.add_row(1..=1, &assgn_constr);
        }

        let mut depth_constr1 = Vec::new();
        let mut depth_constr2 = Vec::new();
        for w in INIT_WSHIFT as usize..total_windows {
            let depth_distr = assignments.depth_distr(w);
            let potential_min_depth = window_trivial_depth[w] as u32;
            let potential_max_depth = potential_min_depth + by_window_vars[w].len() as u32;
            if potential_min_depth == potential_max_depth {
                const_term += depth_distr.ln_pmf(potential_min_depth);
                continue;
            }

            depth_constr1.clear();
            depth_constr2.clear();
            for &var in by_window_vars[w].iter() {
                depth_constr2.push((var, 1.0));
            }
            for depth in potential_min_depth..=potential_max_depth {
                let prob = depth_distr.ln_pmf(depth);
                let var = problem.add_integer_column(prob, 0..=1);
                depth_constr1.push((var, 1.0));
                if depth > potential_min_depth {
                    depth_constr2.push((var, -f64::from(depth - potential_min_depth)));
                }
            }
            problem.add_row(1..=1, &depth_constr1);
            problem.add_row(0..=0, &depth_constr2);
        }

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
            let new_assgn = vals[i..j].iter().position(|&v| v == 1.0).expect("Read has no assignment!");
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