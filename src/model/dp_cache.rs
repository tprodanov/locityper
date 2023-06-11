use std::{
    sync::Arc,
};
use once_cell::sync::OnceCell;
use crate::{
    math::{
        distr::{DiscretePmf, WithQuantile, NBinom, Uniform, LinearCache, BayesCalc},
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

impl<D: DiscretePmf> DiscretePmf for WeightedDistr<D> {
    fn ln_pmf(&self, k: u32) -> f64 {
        self.inner.ln_pmf(k) * self.weight
    }
}

/// Windows with many common k-mers have a mixure of Neg.Binomial and Uniform distributions.
/// Size of the Uniform distribution is calculated as 2 x 0.99 quantile of the Neg.Binomial distribution.
const UNIFSIZE_QUANTILE: f64 = 0.99;
const UNIFSIZE_MULT: f64 = 2.0;

/// Store read depth probabilities for values between 0 and 255 for each GC content.
const CACHE_SIZE: usize = 256;

pub(super) type DistrBox = Box<dyn DiscretePmf>;

type RegularDistr = BayesCalc<NBinom>;

/// Store cached depth distbrutions.
#[derive(Clone)]
pub struct CachedDepthDistrs {
    /// Background read depth distribution.
    bg_depth: bg::depth::ReadDepth,
    /// Multiplication coefficient: assume that read depth is `mul_coef` * bg read depth.
    mul_coef: f64,

    /// Cached read depth distributions in windows with few common k-mers (one for each GC-content).
    cached: [OnceCell<Arc<LinearCache<RegularDistr>>>; GC_BINS],
    /// Uniform distribution size (one for each GC-content).
    unif_size: [OnceCell<u32>; GC_BINS],
    alt_cn: (f64, f64),
}

impl CachedDepthDistrs {
    /// Create a set of cached depth distributions.
    /// Assume that there are `mul_coef` as much reads, as in the background distribution.
    pub fn new(bg_distr: &bg::BgDistr, alt_cn: (f64, f64)) -> Self {
        const NBINOM_CELL: OnceCell<Arc<LinearCache<RegularDistr>>> = OnceCell::new();
        const U32_CELL: OnceCell<u32> = OnceCell::new();
        Self {
            bg_depth: bg_distr.depth().clone(),
            mul_coef: if bg_distr.insert_distr().is_paired_end() { 2.0 } else { 1.0 },
            cached: [NBINOM_CELL; GC_BINS],
            unif_size: [U32_CELL; GC_BINS],
            alt_cn,
        }
    }

    /// Returns a pointer to unmapped distribution (`AlwaysOneDistr`).
    pub fn unmapped_distr(&self) -> AlwaysOneDistr {
        AlwaysOneDistr
    }

    /// Returns a pointer to boundary distribution (`AlwaysOneDistr`).
    pub fn boundary_distr(&self) -> AlwaysOneDistr {
        AlwaysOneDistr
    }

    fn regular_nbinom(&self, gc_content: u8) -> NBinom {
        self.bg_depth.depth_distribution(gc_content).mul(self.mul_coef)
    }

    /// Returns read depth distribution in regular windows at GC-content.
    pub fn regular_distr(&self, gc_content: u8) -> &Arc<LinearCache<RegularDistr>> {
        self.cached[usize::from(gc_content)].get_or_init(|| {
            // Probability at CN = 1.
            let cn1 = self.regular_nbinom(gc_content);
            // Probabilities at alternative CN values.
            let alt = vec![cn1.mul(self.alt_cn.0), cn1.mul(self.alt_cn.1)];
            let bayes = BayesCalc::new(cn1, alt);
            Arc::new(LinearCache::new(bayes, CACHE_SIZE))
        })
    }

    /// Returns depth bound for the given GC-content (see `UNIFSIZE_QUANTILE`).
    fn uniform_size(&self, gc_content: u8) -> u32 {
        *self.unif_size[usize::from(gc_content)].get_or_init(||
            (UNIFSIZE_MULT * self.regular_nbinom(gc_content).quantile(UNIFSIZE_QUANTILE)) as u32)
    }

    /// Returns a box to either `NBinom`, `WeightedDistr`, depending on GC-content and the weight of the window.
    pub fn get_distribution(&self, gc_content: u8, weight: f64) -> DistrBox {
        if weight < 0.00001 {
            Box::new(Uniform::new(0, self.uniform_size(gc_content)))
        } else if weight > 0.99999 {
            Box::new(Arc::clone(&self.regular_distr(gc_content)))
        } else {
            let wdistr = WeightedDistr::new(Arc::clone(self.regular_distr(gc_content)), weight);
            Box::new(wdistr)
        }
    }
}
