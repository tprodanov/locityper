use std::{
    sync::Arc,
};
use crate::{
    math::{
        distr::{DistrBox, DiscretePmf, NBinom, LinearCache, BayesCalc},
    },
    bg::{self,
        depth::GC_BINS,
    },
};

/// Proxy distribution, that produces 1.0 probability for all values.
#[derive(Clone, Copy, Debug)]
pub struct TrivialDistr;

impl DiscretePmf for TrivialDistr {
    fn ln_pmf(&self, _: u32) -> f64 { 0.0 }
}

/// Wrapper over another distribution, where all probabilities are raised to the power `weight`.
/// As a consequence, probabilities will not sum up to one.
#[derive(Clone, Debug)]
struct WeightedDistr<D> {
    weight: f64,
    inner: D,
}

impl<D: DiscretePmf> WeightedDistr<D> {
    /// Creates a new weighted distribution.
    fn new(inner: D, weight: f64) -> Self {
        assert!(weight >= 0.0 && weight <= 1.0, "Weight ({}) must be in [0, 1].", weight);
        Self { inner, weight }
    }
}

impl<D: DiscretePmf + Clone + Sync> DiscretePmf for WeightedDistr<D> {
    fn ln_pmf(&self, k: u32) -> f64 {
        self.inner.ln_pmf(k) * self.weight
    }
}

/// Store read depth probabilities for values between 0 and 255 for each GC content.
const CACHE_SIZE: usize = 256;

/// Read depth distribution for each window.
type WindowDistr = BayesCalc<NBinom, NBinom, 2>;

/// Cached read depth distribution.
type CachedDistr = Arc<LinearCache<WindowDistr>>;

/// Collection of cached read depth distributions for each GC-content.
pub struct DistrCache(Vec<CachedDistr>);

impl DistrCache {
    pub fn new(bg_distr: &bg::BgDistr, alt_cn: (f64, f64)) -> Self {
        let mul_coef = if bg_distr.insert_distr().is_paired_end() { 2.0 } else { 1.0 };
        let mut cache = Vec::with_capacity(GC_BINS);
        for gc in 0..GC_BINS {
            let cn1_distr = bg_distr.depth().depth_distribution(gc.try_into().unwrap()).mul(mul_coef);
            let sub_distr = cn1_distr.mul(alt_cn.0);
            let super_distr = cn1_distr.mul(alt_cn.1);
            let bayes = LinearCache::new(BayesCalc::new(cn1_distr, [sub_distr, super_distr]), CACHE_SIZE);
            cache.push(Arc::new(bayes));
        }
        Self(cache)
    }

    /// Returns a box to either `NBinom`, `WeightedDistr`, depending on GC-content and the weight of the window.
    pub fn get_distribution(&self, gc: u8, weight: f64) -> DistrBox {
        if weight < 1e-7 {
            Box::new(TrivialDistr)
        } else {
            let regular = Arc::clone(&self.0[usize::from(gc)]);
            if weight > 0.99999 {
                Box::new(regular)
            } else {
                Box::new(WeightedDistr::new(regular, weight))
            }
        }
    }
}
