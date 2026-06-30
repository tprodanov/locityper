//! Greedy algorithm to identify a dominating set.

use std::{
    cmp::max,
};
use crate::{
    algo::IntSet,
};

pub type Adjacencies = Vec<Vec<usize>>;

/// Greedy algorithm for obtaining a dominant set, where on each iteration we examine a node with the largest number
/// of not dominated neighbors.
/// Runtime is O(E) since every node stops being white only once.
fn find_dominating_set_greedy(graph: &Adjacencies) -> Vec<usize> {
    #[repr(u8)]
    #[derive(Clone, Copy, PartialEq, Eq, Debug)]
    enum Color {
        White, // Not dominated
        Gray,  // Dominated
        Black, // In the dominating set
    }
    use Color::*;
    // For each node, store its color and degree = number of white neighbors.
    let mut colors_and_degrees: Vec<(Color, usize)> = Vec::with_capacity(graph.len());
    // Contains non-black nodes with index = degree.
    // Bucket 0 is always empty because these nodes are irrelevant.
    let mut buckets = Vec::new();

    let mut output = Vec::new();
    for (node, edges) in graph.iter().enumerate() {
        let d = edges.len();
        if d > 0 {
            colors_and_degrees.push((White, d));
            buckets.resize_with(max(buckets.len(), d + 1), IntSet::default);
            buckets[d].insert(node);
        } else {
            colors_and_degrees.push((Black, 0));
            output.push(node);
        }
    }

    let mut rem = graph.len() - output.len();
    while rem > 0 {
        // Removes first element in the set.
        let node_i = buckets.last_mut().expect("Buckets must not be empty")
            .extract_if(|_| true).next().expect("Last bucket must be non-empty");
        output.push(node_i);

        let i_cd = &mut colors_and_degrees[node_i];
        assert!(i_cd.0 != Black);
        // To disable bucket switching.
        i_cd.1 = 0;

        for &node_j in itertools::chain!(std::iter::once(&node_i), &graph[node_i]) {
            let j_cd = &mut colors_and_degrees[node_j];
            if j_cd.0 != White { continue };
            j_cd.0 = Gray;
            rem -= 1;
            for &node_k in &graph[node_j] {
                let k_cd = &mut colors_and_degrees[node_k];
                // This check includes black nodes.
                if k_cd.1 != 0 {
                    buckets[k_cd.1].remove(&node_k);
                    k_cd.1 -= 1;
                }
                if k_cd.1 != 0 {
                    // Keep 0 bucket empty.
                    buckets[k_cd.1].insert(node_k);
                }
            }
        }
        colors_and_degrees[node_i].0 = Black;

        while buckets.pop_if(|v| v.is_empty()).is_some() {
            // Repeat until we remove all empty buckets at the end.
        }
    }
    output
}

#[cfg(feature = "scip")]
fn find_dominating_set_ilp(graph: &Adjacencies) -> (Vec<usize>, russcip::status::Status) {
    use russcip::prelude::*;
    let mut model = Model::default().hide_output().minimize();
    let mut output = Vec::new();
    let variables: Vec<_> = graph.iter().enumerate().map(|(node, edges)| {
        if edges.is_empty() {
            output.push(node);
            None
        } else {
            Some(model.add(var().bin().obj(1.0)))
        }
    }).collect();
    for (node_var, edges) in itertools::izip!(&variables, graph) {
        let Some(node_var) = node_var else { continue };
        let mut constr_vars = Vec::with_capacity(edges.len() + 1);
        constr_vars.push(node_var);
        constr_vars.extend(edges.iter().map(|&j| variables[j].as_ref().expect("Neighbor node must have a variable")));
        model.add(cons().coefs(constr_vars, vec![1.0; edges.len() + 1]).ge(1.0));
    }
    let solution = model.solve();
    let status = solution.status();
    let best_sol = solution.best_sol().expect("Could not get an ILP solution");
    for (i, var) in variables.iter().enumerate() {
        if let Some(var) = var && best_sol.val(var) > 0.5 {
            output.push(i);
        }
    }
    (output, status)
}

fn validate_solution(graph: &Adjacencies, set: &[usize]) -> bool {
    let mut dominated = vec![false; graph.len()];
    for &i in set {
        dominated[i] = true;
        for &j in &graph[i] {
            dominated[j] = true;
        }
    }
    dominated.iter().all(|&d| d)
}

/// Finds dominating set based on the adjacencies list.
/// If corresponding feature is enabled, tries SCIP ILP solver, otherwise uses a greedy solution.
pub fn find_dominating_set(graph: &Adjacencies) -> Vec<usize> {
    #[cfg(feature = "scip")] {
        let (mut output, status) = find_dominating_set_ilp(graph);
        if status != russcip::status::Status::Optimal {
            let output2 = find_dominating_set_greedy(graph);
            if output2.len() <= output.len() || !validate_solution(graph, &output) {
                output = output2;
            }
        }
        output
    }
    #[cfg(not(feature = "scip"))] {
        find_dominating_set_greedy(graph)
    }
}
