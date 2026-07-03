//! Greedy algorithm to identify a dominating set.

use crate::{
    ext::vec::BitArray,
};

/// Finds dominating set based on the adjacencies list.
/// If corresponding feature is enabled, tries SCIP ILP solver, otherwise uses a greedy solution.
pub fn find_dominating_set<'a>(
    n_vars: usize,
    constraints: impl Iterator<Item = &'a BitArray>,
) -> Vec<usize> {
    use russcip::prelude::*;
    let mut model = Model::default().hide_output().minimize();
    let mut output = Vec::new();
    let variables: Vec<_> = (0..n_vars).map(|_| model.add(var().bin().obj(1.0))).collect();

    for constraint in constraints {
        let mut constr_builder = cons();
        for i in constraint.true_indices() {
            constr_builder = constr_builder.coef(&variables[i], 1.0);
        }
        model.add(constr_builder.ge(1.0));
    }
    let solution = model.solve();
    let status = solution.status();
    if status != russcip::status::Status::Optimal {
        log::debug!("         Could not find an optimal solution (status: {:?})", status);
    }
    let best_sol = solution.best_sol().expect("Could not get an ILP solution");
    for (i, var) in variables.iter().enumerate() {
        if best_sol.val(var) > 0.5 {
            output.push(i);
        }
    }
    output
}
