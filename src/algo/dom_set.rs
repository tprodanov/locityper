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
pub fn find_dominating_set(graph: &Adjacencies) -> Vec<usize> {
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
