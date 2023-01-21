use statrs::distribution::Discrete;
use crate::{
    algo::{
        nbinom::{NBinom, UniformNBinom, CachedDistr},
    },
    bg::depth::GC_BINS,
};

/// Store read depth probabilities for values between 0 and 127 for each GC content.
const CACHE_SIZE: usize = 128;

/// Store cached depth distbruti
pub struct CachedDepthDistrs {
    /// Regular read depth distributions (first + second reads).
    regular: Vec<CachedDistr<NBinom>>,
    /// Read depth distributions at the boundary windows.
    boundary: Vec<CachedDistr<UniformNBinom>>,
}

/// Stores read depth values and distributions across a single contig.
pub struct ContigDepth<'a> {
    /// Current read depth at each window.
    depth: Vec<u32>,
    /// Read depth distribution at each window.
    distrs: Vec<Box<&'a dyn Discrete<u32, f64>>>,
}

impl<'a> ContigDepth<'a> {
    /// Creates a new `ContigDepth` with 0 read depth at all windows,
    /// and sets read depth distributions from a given vector of distributions.
    pub fn new(distrs: Vec<Box<&'a dyn Discrete<u32, f64>>>) -> Self {
        Self {
            depth: vec![0; distrs.len()],
            distrs,
        }
    }

    /// Calculates the total ln-probability of the read depth across all windows.
    pub fn sum_ln_prob(&self) -> f64 {
        self.distrs.iter()
            .zip(&self.depth)
            .map(|(distr, &depth)| distr.ln_pmf(depth))
            .sum()
    }

    /// Assuming that read depth in `window` will change by `depth_change`,
    /// calculates the difference between the new and the old ln-probabilities.
    /// Positive value means that the likelihood will improve.
    /// Does not actually update the read depth.
    pub fn lik_difference(&self, window: usize, depth_change: i32) -> f64 {
        let old_depth = self.depth[window];
        let new_depth = old_depth.checked_add_signed(depth_change).expect("Read depth became negative!");
        self.distrs[window].ln_pmf(new_depth) - self.distrs[window].ln_pmf(old_depth)
    }

    /// Updates read depth in `window` by `depth_change`.
    pub fn update(&mut self, window: usize, depth_change: i32) {
        self.depth[window] = self.depth[window].checked_add_signed(depth_change).expect("Read depth became negative!");
    }

    /// Adds one read to `window`.
    pub fn inc(&mut self, window: usize) {
        self.depth[window] += 1;
    }

    /// Removes one read from `window`.
    pub fn dec(&mut self, window: usize) {
        self.depth[window] -= 1;
    }
}

pub struct ReadAssignment {
    
}
