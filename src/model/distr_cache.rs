use std::{
    sync::Arc,
};
use crate::{
    math::{
        distr::{DiscretePmf, NBinom, LinearCache, BayesCalc},
    },
    bg::{self,
        depth::GC_BINS,
    },
};

/// Store read depth probabilities for values between 0 and 255 for each GC content.
const CACHE_SIZE: usize = 256;

/// Cached read depth distribution.
type CachedDistr = Arc<LinearCache<BayesCalc<NBinom, NBinom, 2>>>;

/// Window distribution.
pub struct WindowDistr {
    /// Window weight.
    weight: f64,
    /// Cached Bayesian calculator. If not provided - weight is too low, all probabilities will be 1.
    distr: Option<CachedDistr>,
}

impl WindowDistr {
    pub const TRIVIAL: Self = Self {
        weight: 0.0,
        distr: None,
    };

    /// Calculates weight ln-probability at read depth `k`.
    pub fn ln_prob(&self, k: u32) -> f64 {
        match &self.distr {
            Some(distr) => self.weight * distr.ln_pmf(k),
            None => 0.0,
        }
    }

    #[inline]
    pub fn weight(&self) -> f64 {
        self.weight
    }

    /// Cached Bayesian distribution, if weight is not too low.
    pub fn inner(&self) -> &Option<CachedDistr> {
        &self.distr
    }

    /// Returns true if this distribution always produces prob 1.
    pub fn is_trivial(&self) -> bool {
        self.distr.is_none()
    }
}

/// Collection of cached read depth distributions for each GC-content.
pub struct DistrCache(Vec<CachedDistr>);

impl DistrCache {
    pub fn new(bg_distr: &bg::BgDistr, alt_cn: (f64, f64)) -> Self {
        let depth = match bg_distr.opt_depth() {
            Some(depth) => depth,
            None => return Self(Vec::new()),
        };
        let mul_coef = if bg_distr.insert_distr().is_paired_end() { 2.0 } else { 1.0 };
        let mut cached_distrs = Vec::with_capacity(GC_BINS);
        for gc in 0..GC_BINS {
            let cn1_distr = depth.depth_distribution(gc.try_into().unwrap()).mul(mul_coef);
            let sub_distr = cn1_distr.mul(alt_cn.0);
            let super_distr = cn1_distr.mul(alt_cn.1);
            let bayes = LinearCache::new(BayesCalc::new(cn1_distr, [sub_distr, super_distr]), CACHE_SIZE);
            cached_distrs.push(Arc::new(bayes));
        }
        Self(cached_distrs)
    }

    /// Returns distribution corresponding to GC content `gc`, without accounting for window weight.
    pub fn get_inner_distribution(&self, gc: u8) -> &CachedDistr {
        &self.0[usize::from(gc)]
    }

    /// Returns a box to either `NBinom`, `WeightedDistr`, depending on GC-content and the weight of the window.
    pub fn get_distribution(&self, gc: u8, weight: f64) -> WindowDistr {
        if weight < 1e-7 {
            WindowDistr::TRIVIAL
        } else {
            WindowDistr {
                weight,
                distr: Some(Arc::clone(&self.0[usize::from(gc)])),
            }
        }
    }
}
