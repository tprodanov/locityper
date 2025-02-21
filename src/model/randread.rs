use rand::Rng;
use crate::model::locs::AllAlignments;

struct Bin {
    /// Start and end positions within the parent vector of indices.
    start: usize,
    end: usize,
    /// Cumulative adjusted count (number of reads in this bin, adjusted by the max weight).
    /// Later, values are divided by the total cumulative count, so that the last value = 1.
    cum_count: f64,
}

/// Structure, that allows to randomly obtain reads with different weights.
/// For speed, reads are distributed into bins with similar weights.
pub(crate) struct ReadSelector {
    /// Indices of non-trivial reads, reverse-ordered by weights.
    non_trivial_ixs: Vec<usize>,

    /// If all reads have the same weight, this array will be empty.
    bins: Vec<Bin>,
}

impl ReadSelector {
    /// all_alns.reads() are reverse-sorted by weight.
    /// non_trivial_ixs: indices of non-trivial reads, with more than one possible assignment to a given genotype.
    ///     must be strictly-increasing, meaning that resulting reads[ix] are also reverse-sorted by weights.
    pub(crate) fn new(non_trivial_ixs: Vec<usize>, all_alns: &AllAlignments) -> Self {
        // Max multiplicative factor between top and bottom weights in a single bin.
        // min_weight >= DROP * max_weight.
        const DROP: f64 = 0.7;
        let n = non_trivial_ixs.len();
        if n == 0 {
            return Self { non_trivial_ixs, bins: Vec::new() }
        }

        let reads = all_alns.reads();
        let w1 = reads[non_trivial_ixs[0]].weight();
        let wn = reads[non_trivial_ixs[n - 1]].weight();
        debug_assert!(wn <= w1, "Read weights must be ordered in decreasing order");
        debug_assert!(wn > 0.0, "All reads must have positive weights");

        let mut max_weight = w1;
        let mut thresh_weight = DROP * w1;
        if wn >= thresh_weight {
            return Self { non_trivial_ixs, bins: Vec::new() }
        }

        let mut start = 0;
        let mut cum_count = 0.0;
        let mut bins = Vec::new();
        for (i, &j) in non_trivial_ixs.iter().enumerate() {
            let w = reads[j].weight();
            if w < thresh_weight {
                let end = i;
                cum_count += (end - start) as f64 * max_weight;
                bins.push(Bin { start, end, cum_count });
                start = i;
                max_weight = w;
                thresh_weight = DROP * w;
            }
        }
        let end = n;
        cum_count += (end - start) as f64 * max_weight;
        bins.push(Bin { start, end, cum_count });
        for bin in &mut bins {
            bin.cum_count /= cum_count;
        }

        Self { non_trivial_ixs, bins }
    }

    /// Returns randomly selected read pair index.
    pub(crate) fn random(&self, rng: &mut impl Rng) -> usize {
        if self.bins.is_empty() {
            return self.non_trivial_ixs[rng.random_range(0..self.non_trivial_ixs.len())]
        }
        let r = rng.random::<f64>();
        for bin in &self.bins {
            if r <= bin.cum_count {
                return self.non_trivial_ixs[rng.random_range(bin.start..bin.end)];
            }
        }
        unreachable!("One of the bins must be selected at this point")
    }

    // /// Returns non-trivial indices.
    // #[inline(always)]
    // pub fn ixs(&self) -> &[usize] {
    //     &self.non_trivial_ixs
    // }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.non_trivial_ixs.len()
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.non_trivial_ixs.is_empty()
    }
}
