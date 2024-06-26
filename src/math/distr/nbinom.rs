use std::fmt::{self, Debug, Display, Formatter};
use statrs::function::{
    beta::beta_reg,
    gamma::ln_gamma,
};
use argmin::{
    core::{State, Executor, CostFunction},
    solver::neldermead::NelderMead,
};
use crate::{
    err::error,
    bg::ser::JsonSer,
};
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
        assert!(0.0 < m && m < v, "Cannot estimate N.Binom. parameters from mean {:.3} and variance {:.3}", m, v);
        NBinom::new(m * m / (v - m), m / v)
    }

    /// Estimate Negative Binomial from mean `m` and variance `v`.
    /// Uses distribution, close to Poisson if variance is too low.
    pub fn estimate_corrected(m: f64, v: f64) -> NBinom {
        assert!(0.0 < m, "Cannot estimate N.Binom. parameters from mean {:.3} and variance {:.3}", m, v);
        const PMAX: f64 = 0.99999;
        let p = m / v;
        if p > PMAX {
            if 0.9 * v < m {
                log::warn!("    Cannot not accurately fit Negative Binomial to mean {:.3} and variance {:.3}", m, v);
            }
            NBinom::new(PMAX * m / (1.0 - PMAX), PMAX)
        } else {
            NBinom::new(m * m / (v - m), p)
        }
    }

    /// Create a new distribution where n is multiplied by `coeff` and p stays the same.
    pub fn mul(&self, coeff: f64) -> NBinom {
        NBinom::new(self.n * coeff, self.p)
    }

    /// Createts a cached Negative Binomial distribution, storing ln_pmf values in 0..=n.
    pub fn cached(self, cache_size: usize) -> LinearCache<Self> {
        LinearCache::new(self, cache_size)
    }

    /// Returns the mode of the distribution: value, which produces the highest probability.
    pub fn mode(&self) -> u32 {
        ((self.n - 1.0) * (1.0 - self.p) / self.p).floor().max(0.0) as u32
    }

    pub fn n(&self) -> f64 {
        self.n
    }

    pub fn p(&self) -> f64 {
        self.p
    }

    /// Assuming a situation, where each random value is then used in a separate Binomial trial with `p = rate`,
    /// new random variable also follows Negative Binomial.
    /// This distribution is returned by this function.
    pub fn binomial_subsample(&self, rate: f64) -> Self {
        Self::new(self.n, self.p / (self.p + rate - self.p * rate))
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

    fn load(obj: &json::JsonValue) -> crate::Result<Self> {
        let n = obj["n"].as_f64().ok_or_else(|| error!(JsonLoad,
            "NBinom: Failed to parse '{}': missing or incorrect 'n' field!", obj))?;
        let p = obj["p"].as_f64().ok_or_else(|| error!(JsonLoad,
            "NBinom: Failed to parse '{}': missing or incorrect 'p' field!", obj))?;
        Ok(Self::new(n, p))
    }
}

impl DiscretePmf for NBinom {
    fn ln_pmf(&self, k: u32) -> f64 {
        let x = f64::from(k);
        self.lnpmf_const + ln_gamma(self.n + x) - ln_gamma(x + 1.0) + x * self.lnq
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

struct NBinomProblem {
    /// Observed mean and variance.
    sample_mean: f64,
    sample_var: f64,
    /// Subsampling rate.
    rate: f64,
    /// Regularization factor.
    lambda: f64,
}

impl NBinomProblem {
    fn new(sample_mean: f64, sample_var: f64, rate: f64, lambda: f64) -> Self {
        Self { sample_mean, sample_var, rate, lambda }
    }

    /// See here https://math.stackexchange.com/questions/4700260,
    /// estimating Negative Binomial parameters from Binomial subsampling.
    fn mean(&self, n: f64, p: f64) -> f64 {
        self.rate * n * (1.0 - p) / p
    }

    fn var(&self, n: f64, p: f64) -> f64 {
        self.rate * n * (1.0 - p) * (p + self.rate - p * self.rate) / (p * p)
    }
}

impl CostFunction for NBinomProblem {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        const LARGE_NUM: f64 = 1e30;
        let n = param[0];
        let p = param[1];
        if n <= 0.0 || p <= 0.0 || p >= 1.0 {
            Ok(LARGE_NUM)
        } else {
            let mean_err = self.mean(n, p) - self.sample_mean;
            let var_err = self.var(n, p) - self.sample_var;
            Ok(mean_err * mean_err + var_err * var_err + self.lambda * n)
        }
    }
}

/// L1-regularized parameter estimator for Negative Binomial.
pub struct RegularizedEstimator {
    /// Subsampling rate (1 by default = no subsampling).
    rate: f64,
    /// Regularization factor (1e-6 by default).
    lambda: f64,
}

impl Default for RegularizedEstimator {
    fn default() -> Self {
        Self {
            rate: 1.0,
            lambda: 1e-5,
        }
    }
}

impl RegularizedEstimator {
    pub fn set_subsampling_rate(&mut self, rate: f64) -> &mut Self {
        assert!(rate > 0.0 && rate <= 1.0, "Subsampling rate ({}) must be within (0, 1].", rate);
        self.rate = rate;
        self
    }

    pub fn set_lambda(&mut self, lambda: f64) -> &mut Self {
        assert!(lambda >= 0.0, "Regularization factor ({}) must be non-negative.", lambda);
        self.lambda = lambda;
        self
    }

    /// Estimate Negative Binomial distribution using numerical method Nelder-Mead.
    pub fn estimate(&self, sample_mean: f64, sample_var: f64) -> NBinom {
        let start_points = vec![
            vec![10.0, 0.3],
            vec![20.0, 0.7],
            vec![30.0, 0.3],
        ];
        let solver = NelderMead::new(start_points)
            .with_sd_tolerance(1e-6).unwrap();
        let problem = NBinomProblem::new(sample_mean, sample_var, self.rate, self.lambda);
        let solution = Executor::new(problem, solver)
            .run()
            .expect("Nelder-Mead finished with an error");
        let param = solution.state.get_param().expect("Nelder-Mead failed to find appropriate N.Binom. parameters");
        NBinom::new(param[0], param[1])
    }
}