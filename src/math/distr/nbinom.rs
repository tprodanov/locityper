use std::fmt::{self, Debug, Display, Formatter};
use statrs::function::{
    beta::beta_reg,
    gamma::ln_gamma,
};
use crate::bg::ser::{JsonSer, LoadError};
use super::{DiscretePmf, DiscreteCdf, WithMoments, LinearCache};

/// Negative Binomial distribution,
/// similar to `statrs::distribution::NegativeBinomial`, but storing some precalculated values.
///
/// n: number of successes,
/// p: probability in one trial,
/// x: number of failures before n successes are achieved.
#[derive(Clone)]
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
        assert!(v > m, "Negative Binomial: cannot estimate parameters for mean {:.3} and variance {:.3}", m, v);
        NBinom::new(m * m / (v - m), m / v)
    }

    /// Create a new distribution where n is multiplied by `coeff` and p stays the same.
    pub fn mul(&self, coeff: f64) -> NBinom {
        NBinom::new(self.n * coeff, self.p)
    }

    /// Createts a cached Negative Binomial distribution, storing ln_pmf values in 0..=n.
    pub fn cached(self, cache_size: usize) -> LinearCache<Self> {
        LinearCache::new(self, cache_size)
    }

    pub fn n(&self) -> f64 {
        self.n
    }

    pub fn p(&self) -> f64 {
        self.p
    }

    pub fn ln_pmf_f64(&self, x: f64) -> f64 {
        self.lnpmf_const + ln_gamma(self.n + x) - ln_gamma(x + 1.0) + x * self.lnq
    }
}

impl Debug for NBinom {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "N.Binom.(n = {:?}, p = {:?})", self.n, self.p)
    }
}

impl Display for NBinom {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "N.Binom.(n = {:.5}, p = {:.5})", self.n, self.p)
    }
}

impl JsonSer for NBinom {
    fn save(&self) -> json::JsonValue {
        json::object!{
            n: self.n,
            p: self.p,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        let n = obj["n"].as_f64().ok_or_else(|| LoadError(format!(
            "NBinom: Failed to parse '{}': missing or incorrect 'n' field!", obj)))?;
        let p = obj["p"].as_f64().ok_or_else(|| LoadError(format!(
            "NBinom: Failed to parse '{}': missing or incorrect 'p' field!", obj)))?;
        Ok(Self::new(n, p))
    }
}

impl DiscretePmf for NBinom {
    fn ln_pmf(&self, k: u32) -> f64 {
        self.ln_pmf_f64(f64::from(k))
    }
}

impl WithMoments for NBinom {
    fn mean(&self) -> f64 {
        self.n * (1.0 - self.p) / self.p
    }

    fn variance(&self) -> f64 {
        self.n * (1.0 - self.p) / (self.p * self.p)
    }
}

impl DiscreteCdf for NBinom {
    fn cdf(&self, k: u32) -> f64 {
        beta_reg(self.n, f64::from(k + 1), self.p)
    }

    fn sf(&self, k: u32) -> f64 {
        beta_reg(f64::from(k + 1), self.n, 1.0 - self.p)
    }
}
