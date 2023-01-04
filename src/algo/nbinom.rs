use std::cell::Cell;
use statrs::{
    statistics::{Min, Max},
    distribution::{Discrete, DiscreteCDF},
    function::{
        beta::beta_reg,
        gamma::ln_gamma,
    },
};

/// Negative Binomial distribution,
/// similar to `statrs::distribution::NegativeBinomial`, but storing some precalculated values.
///
/// n: number of successes,
/// p: probability in one trial,
/// x: number of failures before n successes are achieved.
pub struct NBinom {
    n: f64,
    p: f64,
    // ln(1 - p)
    lnq: f64,
    // n * ln(p) - ln_gamma(n)
    lnpmf_const: f64,
}

impl NBinom {
    /// Creates a new Negative Binomial distribution.
    /// Conditions: n > 0, 0 <= p <= 1.
    pub fn new(n: f64, p: f64) -> NBinom {
        assert!(n > 0.0 && p >= 0.0 && p <= 1.0, "Incorrect Negative Binomial parameters n = {:?}, p = {:?}", n, p);
        NBinom {
            n, p,
            lnq: (-p).ln_1p(),
            lnpmf_const: n * p.ln() - ln_gamma(n),
        }
    }

    /// Create Negative Binomial distribution from mean `m` and variance `v`.
    /// Conditions: 0 < m < v.
    pub fn estimate(m: f64, v: f64) -> NBinom {
        NBinom::new(m * m / (v - m), m / v)
    }

    /// Create a new distribution where n is multiplied by `mul` and p stays the same.
    pub fn mul_n(&self, mul: f64) -> NBinom {
        NBinom::new(self.n * mul, self.p)
    }

    /// Creates a cached Negative Binomial distribution, with precomputed values from 0 to <0.999 quantile>.
    pub fn cached(self) -> CachedDistr<Self> {
        let q999 = self.inverse_cdf(0.999);
        CachedDistr::new(self, q999)
    }

    pub fn n(&self) -> f64 {
        self.n
    }

    pub fn p(&self) -> f64 {
        self.p
    }

    pub fn mean(&self) -> f64 {
        self.n * (1.0 - self.p) / self.p
    }

    pub fn variance(&self) -> f64 {
        self.n * (1.0 - self.p) / (self.p * self.p)
    }

    pub fn ln_pmf(&self, x: f64) -> f64 {
        self.lnpmf_const + ln_gamma(self.n + x) - ln_gamma(x + 1.0) + x * self.lnq
    }
}

macro_rules! nbinom_impl {
    (
        $int:ident
    ) => {
        impl Min<$int> for NBinom {
            fn min(&self) -> $int {
                0
            }
        }

        impl Max<$int> for NBinom {
            fn max(&self) -> $int {
                std::$int::MAX
            }
        }

        impl Discrete<$int, f64> for NBinom {
            fn pmf(&self, x: $int) -> f64 {
                self.ln_pmf(x as f64).exp()
            }

            fn ln_pmf(&self, x: $int) -> f64 {
                self.ln_pmf(x as f64)
            }
        }

        impl DiscreteCDF<$int, f64> for NBinom {
            fn cdf(&self, x: $int) -> f64 {
                beta_reg(self.n, (x + 1) as f64, self.p)
            }

            fn sf(&self, x: $int) -> f64 {
                beta_reg((x + 1) as f64, self.n, 1.0 - self.p)
            }
        }
    }
}

nbinom_impl!(u32);
nbinom_impl!(u64);

/// Distribution with cached ln_pmf values.
pub struct CachedDistr<D: Discrete<u64, f64>> {
    distr: D,
    cache: Vec<f64>,
}

impl<D: Discrete<u64, f64>> CachedDistr<D> {
    /// Creates the cached distribution.
    /// Precomputes ln_pmf for all integers in [0, n].
    pub fn new(distr: D, n: u64) -> Self {
        let cache = (0..=n).map(|i| distr.ln_pmf(i)).collect();
        CachedDistr { distr, cache }
    }

    pub fn ln_pmf(&self, k: u64) -> f64 {
        let i = k as usize;
        if i < self.cache.len() {
            self.cache[i]
        } else {
            self.distr.ln_pmf(k)
        }
    }
}
