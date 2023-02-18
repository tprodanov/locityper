use good_lp::{
    {ProblemVariables, variable, Variable, Expression, SolverModel, Solution},
    variable::UnsolvedProblem,
};
use crate::model::assgn::{ReadAssignment, UNMAPPED_WINDOW};
use super::Solver;

pub struct IlpSolver<S: SolverModel> {
    model: Option<S>,
    assignment_vars: Vec<Variable>,
    assignments: ReadAssignment,
}

impl<S: SolverModel> IlpSolver<S> {
    pub fn build<F>(assignments: ReadAssignment, builder: F) -> Self
    where F: Fn(UnsolvedProblem) -> S
    {
        let contig_windows = assignments.contig_windows();
        let total_windows = contig_windows.total_windows() as usize;
        let mut by_window_vars: Vec<Vec<Variable>> = vec![Vec::new(); total_windows];
        let mut window_trivial_depth = vec![0; total_windows];

        let mut variables = ProblemVariables::new();
        let mut objective = Expression::default();
        let mut assignment_vars = Vec::new();
        let mut all_constraints = Vec::new();
        for rp in 0..assignments.total_reads() {
            let read_alns = assignments.possible_read_alns(rp);
            if read_alns.len() == 1 {
                let loc = &read_alns[0];
                let (w1, w2) = contig_windows.get_pair_window_ixs(loc);
                window_trivial_depth[w1 as usize] += 1;
                window_trivial_depth[w2 as usize] += 1;
                objective += loc.ln_prob();
                continue;
            }

            let mut assgn_constr = Expression::default();
            for (j, loc) in read_alns.iter().enumerate() {
                let var = variables.add(variable().binary().name(format!("R{}_{}", rp, j)));
                objective.add_mul(loc.ln_prob(), var);
                assignment_vars.push(var);
                assgn_constr += var;
                let (w1, w2) = contig_windows.get_pair_window_ixs(loc);
                if w1 != UNMAPPED_WINDOW {
                    by_window_vars[w1 as usize].push(var);
                } else {
                    window_trivial_depth[w1 as usize] += 1;
                }
                if w2 != UNMAPPED_WINDOW {
                    by_window_vars[w2 as usize].push(var);
                } else {
                    window_trivial_depth[w2 as usize] += 1;
                }
            }
            all_constraints.push(assgn_constr.eq(1.0));
        }

        for w in 0..total_windows {
            let depth_distr = assignments.depth_distr(w);
            let potential_min_depth = window_trivial_depth[w] as u32;
            let potential_max_depth = potential_min_depth + by_window_vars[w].len() as u32;
            if potential_min_depth == potential_max_depth {
                objective += depth_distr.ln_pmf(potential_min_depth);
                continue;
            }

            let mut depth_constr1 = Expression::default();
            let mut depth_constr2 = Expression::default();
            for &var in by_window_vars[w].iter() {
                depth_constr2 += var;
            }
            for depth in potential_min_depth..=potential_max_depth {
                let var = variables.add(variable().binary().name(format!("D{}_{}", w, depth)));
                objective.add_mul(depth_distr.ln_pmf(depth), var);
                depth_constr1 += var;
                if depth > potential_min_depth {
                    depth_constr2.add_mul(-f64::from(depth - potential_min_depth), var);
                }
            }
            all_constraints.push(depth_constr1.eq(1.0));
            all_constraints.push(depth_constr2.eq(0.0));
        }

        let problem = variables.maximise(objective);
        let mut model = builder(problem);
        for constr in all_constraints.into_iter() {
            model.add_constraint(constr);
        }
        Self {
            model: Some(model),
            assignments,
            assignment_vars,
        }
    }

    /// Query read assignments from the ILP solution, and set them in the `self.assignments`.
    fn set_assignments(&mut self, sol: &S::Solution) {
        let vals: Vec<_> = self.assignment_vars.iter().map(|&var| sol.value(var)).collect();
        let mut i = 0;
        self.assignments.init_assignments(|locs| {
            let j = i + locs.len();
            let new_assgn = vals[i..j].iter()
                .position(|&v| v >= 0.99999 && v <= 1.00001)
                .expect("Read has no assignment!");
            i = j;
            new_assgn
        });
        assert_eq!(i, vals.len(), "Numbers of total read assignments do not match");
    }
}

impl<S: SolverModel> Solver for IlpSolver<S> {
    type Error = S::Error;

    fn is_seedable() -> bool { false }

    fn set_seed(&mut self, _seed: u64) -> Result<(), Self::Error> {
        panic!("Cannot set seed to the Highs solver!")
    }

    /// Resets and initializes anew read assignments.
    fn initialize(&mut self) -> Result<(), Self::Error> {
        assert!(self.model.is_some(), "Cannot run `IlpSolver::initialize` twice!");
        self.assignments.init_assignments(|_| 0);
        Ok(())
    }

    /// Perform one iteration, and return the likelihood improvement.
    fn step(&mut self) -> Result<f64, Self::Error> {
        let solution = self.model.take().expect("Cannot run `IlpSolver::step` twice!").solve()?;
        // if solved_model.status() != Status::Optimal {
        //     panic!("Highs model finished with non-optimal status {:?}", solved_model.status());
        // }
        let old_lik = self.assignments.likelihood();
        self.set_assignments(&solution);
        let new_lik = self.assignments.likelihood();
        Ok(new_lik - old_lik)
        // unimplemented!()
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

#[cfg(feature = "highs")]
pub type HighsSolver = IlpSolver<good_lp::solvers::highs::HighsProblem>;

#[cfg(feature = "highs")]
impl HighsSolver {
    pub fn new(assignments: ReadAssignment) -> Self {
        Self::build(assignments, good_lp::highs)
    }
}

#[cfg(feature = "cbc")]
pub type CbcSolver = IlpSolver<good_lp::solvers::coin_cbc::CoinCbcProblem>;

#[cfg(feature = "cbc")]
impl CbcSolver {
    pub fn new(assignments: ReadAssignment) -> Self {
        Self::build(assignments, good_lp::coin_cbc)
    }
}
