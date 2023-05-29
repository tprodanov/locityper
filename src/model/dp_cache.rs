use std::{
    sync::Arc,
};
use once_cell::sync::OnceCell;
use crate::{
    math::{
        Ln,
        distr::{DiscretePmf, WithQuantile, NBinom, Uniform, LinearCache},
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

#[derive(Clone, Debug)]
struct WeightedDistr {
    ln_probs: Vec<f64>,
}

impl WeightedDistr {
    /// Creates a new distribution from all cached values in the `cached` distribution, raising them in the power of
    /// `weight` (ln-probabilities are multiplied by weight).
    ///
    /// Values out of the cache range are ignored.
    fn new(cached: &LinearCache<impl DiscretePmf>, weight: f64) -> Self {
        let mut ln_probs: Vec<_> = (0..cached.cache_size() as u32)
            .map(|k| cached.ln_pmf(k) * weight)
            .collect();
        let norm_fct = Ln::sum(&ln_probs);
        assert!(norm_fct.is_finite(), "Cannot normalize by an infinite factor.");
        ln_probs.iter_mut().for_each(|v| *v -= norm_fct);
        Self { ln_probs }
    }
}

impl DiscretePmf for WeightedDistr {
    fn ln_pmf(&self, k: u32) -> f64 {
        // This will produce an error if `k` is too large. This is intended behaviour.
        self.ln_probs[k as usize]
    }
}

/// Windows with many common k-mers have a mixure of Neg.Binomial and Uniform distributions.
/// Size of the Uniform distribution is calculated as 2 x 0.99 quantile of the Neg.Binomial distribution.
const UNIFSIZE_QUANTILE: f64 = 0.99;
const UNIFSIZE_MULT: f64 = 2.0;

/// Store read depth probabilities for values between 0 and 255 for each GC content.
const CACHE_SIZE: usize = 256;

pub(super) type DistrBox = Box<dyn DiscretePmf>;

type RegularDistr = NBinom;

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
}

impl CachedDepthDistrs {
    /// Create a set of cached depth distributions.
    /// Assume that there are `mul_coef` as much reads, as in the background distribution.
    pub fn new(bg_distr: &bg::BgDistr) -> Self {
        const NBINOM_CELL: OnceCell<Arc<LinearCache<RegularDistr>>> = OnceCell::new();
        const U32_CELL: OnceCell<u32> = OnceCell::new();
        Self {
            bg_depth: bg_distr.depth().clone(),
            mul_coef: if bg_distr.insert_distr().is_paired_end() { 2.0 } else { 1.0 },
            cached: [NBINOM_CELL; GC_BINS],
            unif_size: [U32_CELL; GC_BINS],
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
            let null = self.regular_nbinom(gc_content);
            // let alt = vec![null.mul(0.01), null.mul(2.0)];
            // let bayes = BayesCalc::new(null, alt);
            Arc::new(LinearCache::new(null, CACHE_SIZE))
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
            let wdistr = WeightedDistr::new(self.regular_distr(gc_content), weight);
            Box::new(wdistr)
        }
    }
}
