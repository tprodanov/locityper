use std::fmt::{self, Debug, Display, Formatter};
use once_cell::unsync::OnceCell;
use statrs::{
    statistics::{Min, Max},
    distribution::Discrete,
    function::{
        beta::beta_reg,
        gamma::ln_gamma,
    },
};
use crate::bg::ser::{JsonSer, LoadError};

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
        NBinom::new(m * m / (v - m), m / v)
    }

    /// Create a new distribution where n is multiplied by `mul` and p stays the same.
    pub fn mul_n(&self, mul: f64) -> NBinom {
        NBinom::new(self.n * mul, self.p)
    }

    /// Creates a cached Negative Binomial distribution, storing ln_pmf values from 0 up to <0.999 quantile>.
    pub fn cached_q999(self) -> CachedDistr<Self> {
        let q999 = self.quantile(0.999).ceil() as usize;
        CachedDistr::new(self, q999 + 1)
    }

    /// Createts a cached Negative Binomial distribution, storing ln_pmf values in 0..=n.
    pub fn cached(self, cache_size: usize) -> CachedDistr<Self> {
        CachedDistr::new(self, cache_size)
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

    pub fn ln_pmf_f64(&self, x: f64) -> f64 {
        self.lnpmf_const + ln_gamma(self.n + x) - ln_gamma(x + 1.0) + x * self.lnq
    }

    pub fn cdf(&self, x: f64) -> f64 {
        beta_reg(self.n, x + 1.0, self.p)
    }

    pub fn sf(&self, x: f64) -> f64 {
        beta_reg(x + 1.0, self.n, 1.0 - self.p)
    }

    pub fn quantile(&self, q: f64) -> f64 {
        if q <= 0.0 {
            return 0.0;
        }
        if q >= 1.0 {
            return f64::INFINITY;
        }
        let mut low = 0;
        let mut high = (2.0 * self.mean()) as u64;
        while self.cdf(high as f64) < q {
            low = high;
            high *= 2;
        }
        while high >= low {
            let mid = (low + high) / 2;
            if self.cdf(mid as f64) >= q {
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }
        let k = high as f64;
        let cdf0 = self.cdf(k);
        assert!(cdf0 <= q);
        let cdf1 = self.cdf(k + 1.0);
        let diff = cdf1 - cdf0;
        match diff.total_cmp(&0.0) {
            std::cmp::Ordering::Equal => k,
            std::cmp::Ordering::Greater => {
                let r = (q - cdf0) / diff;
                k * (1.0 - r) + (k + 1.0) * r
            }
            _ => panic!("Assertion failed: CDF is non-increasing!"),
        }
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

impl Min<u64> for NBinom {
    fn min(&self) -> u64 {
        0
    }
}

impl Max<u64> for NBinom {
    fn max(&self) -> u64 {
        std::u64::MAX
    }
}

impl Discrete<u64, f64> for NBinom {
    fn pmf(&self, x: u64) -> f64 {
        self.ln_pmf_f64(x as f64).exp()
    }

    fn ln_pmf(&self, x: u64) -> f64 {
        self.ln_pmf_f64(x as f64)
    }
}

/// Distribution with cached ln_pmf values.
#[derive(Clone)]
pub struct CachedDistr<D: Discrete<u64, f64>> {
    distr: D,
    cache: Vec<OnceCell<f64>>,
}

impl<D: Discrete<u64, f64>> CachedDistr<D> {
    /// Creates the cached distribution.
    /// Caches ln_pmf values in `0..cache_size`.
    pub fn new(distr: D, cache_size: usize) -> Self {
        CachedDistr {
            distr,
            cache: vec![OnceCell::new(); cache_size],
        }
    }

    /// Get inner distrbution.
    pub fn distr(&self) -> &D {
        &self.distr
    }

    /// Returns ln(pmf(k)). Caches values in a certain range.
    pub fn ln_pmf(&self, k: u64) -> f64 {
        let i = k as usize;
        if i < self.cache.len() {
            *self.cache[i].get_or_init(|| self.distr.ln_pmf(k))
        } else {
            self.distr.ln_pmf(k)
        }
    }
}

impl<D: Discrete<u64, f64> + Debug> Debug for CachedDistr<D> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "Cached[{:?}, size={}]", self.distr, self.cache.len())
    }
}

impl<D: Discrete<u64, f64> + Display> Display for CachedDistr<D> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "Cached[{}, size={}]", self.distr, self.cache.len())
    }
}

impl<D: Discrete<u64, f64> + JsonSer> JsonSer for CachedDistr<D> {
    fn save(&self) -> json::JsonValue {
        json::object!{
            distr: self.distr.save(),
            size: self.cache.len(),
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        if !obj.has_key("distr") {
            return Err(LoadError(format!("CachedDistr: Failed to parse '{}': missing 'distr' field!", obj)));
        }
        let distr = D::load(&obj["distr"])?;
        let size = obj["size"].as_usize().ok_or_else(|| LoadError(format!(
            "CachedDistr: Failed to parse '{}': missing or incorrect 'size' field!", obj)))?;
        Ok(Self::new(distr, size))
    }
}
