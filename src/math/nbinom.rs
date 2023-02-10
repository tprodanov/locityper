use std::fmt::{self, Debug, Display, Formatter};
use once_cell::unsync::OnceCell;
use statrs::{
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
        assert!(v > m, "Negative Binomial: cannot estimate parameters for mean {:.3} and variance {:.3}", m, v);
        NBinom::new(m * m / (v - m), m / v)
    }

    /// Create a new distribution where n is multiplied by `coeff` and p stays the same.
    pub fn mul(&self, coeff: f64) -> NBinom {
        NBinom::new(self.n * coeff, self.p)
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

    /// Cumulative distribution function: P(X <= x).
    pub fn cdf(&self, x: f64) -> f64 {
        beta_reg(self.n, x + 1.0, self.p)
    }

    /// Survival function: P(X > x).
    pub fn sf(&self, x: f64) -> f64 {
        beta_reg(x + 1.0, self.n, 1.0 - self.p)
    }

    /// Inverse CDF function: returns approximate such `x` that `cdf(x) = q`.
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

impl Discrete<u32, f64> for NBinom {
    fn pmf(&self, x: u32) -> f64 {
        self.ln_pmf_f64(f64::from(x)).exp()
    }

    fn ln_pmf(&self, x: u32) -> f64 {
        self.ln_pmf_f64(f64::from(x))
    }
}

/// Chimeric distribution: discrete Uniform on the left, and Negative Binomial distribution on the right.
#[derive(Clone)]
pub struct UniformNBinom {
    partition: u32,
    uniform_ln_prob: f64,
    tail: NBinom,
    tail_ln_weight: f64,
}

impl UniformNBinom {
    /// Creates `UniformNBinom` from discrete uniform distribution in `[0, partition]` with `1 - tail_weight`,
    /// with Negative Binomial tail in `[partition+1, INF)` with weight `tail_weight`.
    ///
    /// If `tail_weight` is not given, select a weight that would give a smooth transition.
    pub fn new(partition: u32, tail: NBinom, tail_weight: Option<f64>) -> Self {
        let x = f64::from(partition);
        let tail_sf = tail.sf(x);
        // let tail_weight = tail_weight.unwrap_or(tail_sf);

        let tail_weight = if let Some(weight) = tail_weight {
            assert!(0.0 <= weight && weight <= 1.0, "UniformNBinom: weight ({}) must be in [0, 1]", weight);
            weight
        } else {
            let pmf = tail.pmf(partition);
            // Select weight, for which nbinom.pmf(x) ~= self.pmf(x)
            1.0 / ((x + 1.0) * pmf / tail_sf + 1.0)
        };

        let uniform_ln_prob = (-tail_weight).ln_1p() - x.ln_1p();
        let tail_ln_weight = tail_weight.ln() - tail_sf.ln();
        Self { partition, uniform_ln_prob, tail, tail_ln_weight }
    }
}

impl Discrete<u32, f64> for UniformNBinom {
    fn pmf(&self, x: u32) -> f64 {
        self.ln_pmf(x).exp()
    }

    fn ln_pmf(&self, x: u32) -> f64 {
        if x <= self.partition {
            self.uniform_ln_prob
        } else {
            self.tail_ln_weight + self.tail.ln_pmf(x)
        }
    }
}

/// Distribution with cached ln_pmf values.
#[derive(Clone)]
pub struct CachedDistr<D> {
    distr: D,
    cache: Vec<OnceCell<f64>>,
}

impl<D> CachedDistr<D> {
    /// Creates the cached distribution.
    /// Caches ln_pmf values in `0..cache_size`.
    pub fn new(distr: D, cache_size: usize) -> Self {
        CachedDistr {
            distr,
            cache: vec![OnceCell::new(); cache_size],
        }
    }

    /// Get inner distribution.
    pub fn distr(&self) -> &D {
        &self.distr
    }

    /// Returns the maximal cache size. Values 0..cache_size are stored once calculated the first time.
    pub fn cache_size(&self) -> usize {
        self.cache.len()
    }
}

impl<D: Discrete<u32, f64>> Discrete<u32, f64> for CachedDistr<D> {
    /// Returns ln(pmf(k)). Caches values in a certain range.
    fn ln_pmf(&self, k: u32) -> f64 {
        let i = k as usize;
        if i < self.cache.len() {
            *self.cache[i].get_or_init(|| self.distr.ln_pmf(k))
        } else {
            self.distr.ln_pmf(k)
        }
    }

    fn pmf(&self, k: u32) -> f64 {
        self.ln_pmf(k).exp()
    }
}

impl<D: Debug> Debug for CachedDistr<D> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "Cached[{:?}, size={}]", self.distr, self.cache.len())
    }
}

impl<D: Display> Display for CachedDistr<D> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "Cached[{}, size={}]", self.distr, self.cache.len())
    }
}

impl<D: JsonSer> JsonSer for CachedDistr<D> {
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
