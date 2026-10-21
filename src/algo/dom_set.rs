//! Greedy algorithm to identify a dominating set.

use crate::{
    ext::vec::BitArray,
};

/// Finds dominating set based on the adjacencies list.
pub fn find_dominating_set<'a>(
    n_vars: usize,
    constraints: impl Iterator<Item = &'a BitArray>,
) -> Vec<usize> {
    use russcip::prelude::*;
    let timer = std::time::Instant::now();
    let mut model = Model::default().hide_output().minimize();
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
    let mut output = Vec::new();
    for (i, var) in variables.iter().enumerate() {
        if best_sol.val(var) > 0.5 {
            output.push(i);
        }
    }
    log::debug!("    Found ILP solution in {}", crate::ext::fmt::Duration(timer.elapsed()));
    output
}

// /// Finds dominating set based on the adjacencies list.
// pub fn find_dominating_set_microlp<'a>(
//     n_vars: usize,
//     constraints: impl Iterator<Item = &'a BitArray>,
// ) -> crate::Result<Vec<usize>> {
//     let timer = std::time::Instant::now();
//     use microlp::{ComparisonOp, LinearExpr, OptimizationDirection, Problem};
//     let mut problem = Problem::new(OptimizationDirection::Minimize);
//     problem.set_time_limit(std::time::Duration::from_mins(5));
//     let variables: Vec<_> = (0..n_vars).map(|_| problem.add_binary_var(1.0)).collect();

//     for constraint in constraints {
//         let mut expression = LinearExpr::empty();
//         for i in constraint.true_indices() {
//             expression.add(variables[i], 1.0);
//         }
//         problem.add_constraint(expression, ComparisonOp::Ge, 1.0);
//     }
//     let solution = problem
//         .solve().map_err(|e| crate::Error::Solver("MicroLP", format!("{:?}", e)))?
//         .into_solution().map_err(|_|
//             crate::Error::Solver("MicroLP", "timie limit exceeded and no feasible solution found".to_string()))?;
//     if solution.status() != microlp::SolutionStatus::Optimal {
//         log::debug!("         Could not find an optimal solution (optimal solution feasible)");
//     }
//     let mut output = Vec::new();
//     for (i, &var) in variables.iter().enumerate() {
//         if solution.var_value(var) > 0.5 {
//             output.push(i);
//         }
//     }
//     log::debug!("    Found ILP solution in {}", crate::ext::fmt::Duration(timer.elapsed()));
//     Ok(output)
// }
