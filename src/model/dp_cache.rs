use std::{
    sync::Arc,
};
use once_cell::sync::OnceCell;
use crate::{
    math::{
        distr::{DistrBox, DiscretePmf, NBinom, LinearCache, BayesCalc},
    },
    bg::{self,
        depth::GC_BINS,
    },
};

/// Fake proxy distribution, that has 1.0 probability for all values.
#[derive(Clone, Copy, Debug)]
pub struct AlwaysOneDistr;

impl DiscretePmf for AlwaysOneDistr {
    fn ln_pmf(&self, _: u32) -> f64 { 0.0 }
}

/// Wrapper over another distribution, where all probabilities are raised to the power `weight`.
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

type RegularDistr = BayesCalc<NBinom, NBinom, 2>;

/// Store cached depth distbrutions.
#[derive(Clone)]
pub struct CachedDepthDistrs<'a> {
    /// Background read depth distribution.
    bg_depth: &'a bg::depth::ReadDepth,
    /// Multiplication coefficient: assume that read depth is `mul_coef` * bg read depth.
    mul_coef: f64,

    /// Cached read depth distributions in windows with few common k-mers (one for each GC-content).
    cached: [OnceCell<Arc<LinearCache<RegularDistr>>>; GC_BINS],
    alt_cn: (f64, f64),
}

impl<'a> CachedDepthDistrs<'a> {
    /// Create a set of cached depth distributions.
    pub fn new(bg_distr: &'a bg::BgDistr, alt_cn: (f64, f64)) -> Self {
        const NBINOM_CELL: OnceCell<Arc<LinearCache<RegularDistr>>> = OnceCell::new();
        Self {
            bg_depth: bg_distr.depth(),
            mul_coef: if bg_distr.insert_distr().is_paired_end() { 2.0 } else { 1.0 },
            cached: [NBINOM_CELL; GC_BINS],
            alt_cn,
        }
    }

    /// Returns read depth distribution in regular windows at GC-content.
    fn regular_distr(&self, gc_content: u8) -> &Arc<LinearCache<RegularDistr>> {
        self.cached[usize::from(gc_content)].get_or_init(|| {
            // Probability at CN = 1.
            let cn1 = self.bg_depth.depth_distribution(gc_content).mul(self.mul_coef);
            // Probabilities at alternative CN values.
            let alt = [cn1.mul(self.alt_cn.0), cn1.mul(self.alt_cn.1)];
            let bayes = BayesCalc::new(cn1, alt);
            Arc::new(LinearCache::new(bayes, CACHE_SIZE))
        })
    }

    /// Returns a box to either `NBinom`, `WeightedDistr`, depending on GC-content and the weight of the window.
    pub fn get_distribution(&self, gc_content: u8, weight: f64) -> DistrBox {
        if weight < 0.00001 {
            Box::new(AlwaysOneDistr)
        } else {
            let regular = Arc::clone(self.regular_distr(gc_content));
            if weight > 0.99999 {
                Box::new(regular)
            } else {
                Box::new(WeightedDistr::new(regular, weight))
            }
        }
    }

    pub fn bg_depth(&self) -> &bg::depth::ReadDepth {
        &self.bg_depth
    }
}
