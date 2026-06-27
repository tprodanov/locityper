//! Greedy algorithm to identify a dominating set.

use std::{
    mem,
    cmp::{min, max, Reverse},
};
use rand::{Rng, RngExt};

pub type AdjacencyList = Vec<Vec<usize>>;

fn greedy_iteration(graph: &AdjacencyList, order: &[usize], covered: &mut [bool], output: &mut Vec<usize>) {
    output.clear();
    covered.fill(false);
    for &i in order {
        if mem::replace(&mut covered[i], true) { continue }
        output.push(i);
        for &j in &graph[i] {
            covered[j] = true;
        }
    }
}

pub fn banded_shuffle<T>(v: &mut [T], band: usize, rng: &mut impl Rng) {
    let n = v.len();
    for i in 0..n {
        let a = i.saturating_sub(band);
        let b = min(i + band + 1, n);
        let j = rng.random_range(a..b);
        if i != j {
            v.swap(i, j);
        }
    }
}

/// Find minimal dominating set.
/// First, nodes are sorted by their degree, and greedily included in the dominating set.
/// During subsequent iterations, ordered list of nodes is partially shuffled (up to shuffle_band * len positions),
/// after which greedy iteration is repeated.
pub fn find_dominating_set(
    graph: &AdjacencyList,
    n_iters: usize,
    shuffle_band: f64,
    rng: &mut impl Rng,
) -> Vec<usize> {
    let n = graph.len();
    let mut order: Vec<_> = (0..n).collect();
    order.sort_by_key(|&i| Reverse(graph[i].len()));

    let mut covered = vec![false; n];
    let mut output = Vec::with_capacity(n / 2);
    greedy_iteration(graph, &order, &mut covered, &mut output);
    if n_iters <= 1 { return output };

    let mut best = output.clone();
    let mut shuf_order = order.clone();
    let band = max(1, (shuffle_band * n as f64).round() as usize);
    for _ in 1..n_iters {
        shuf_order.copy_from_slice(&order);
        banded_shuffle(&mut shuf_order, band, rng);
        greedy_iteration(graph, &order, &mut covered, &mut output);
        if output.len() < best.len() {
            best.truncate(output.len());
            best.copy_from_slice(&output);
        }
    }
    best
}
