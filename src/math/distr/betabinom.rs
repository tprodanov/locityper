use std::{
    fmt,
};
use statrs::function::beta::ln_beta;
use argmin::{
    core::{State, Executor, CostFunction},
    solver::neldermead::NelderMead,
};
use crate::math::Ln;

/// Beta-binomial distribution, see https://en.wikipedia.org/wiki/Beta-binomial_distribution
#[derive(Clone)]
pub struct BetaBinomial {
    alpha: f64,
    beta: f64,
    /// Constant term: ln B(a, b).
    ln_beta_ab: f64,
}

impl BetaBinomial {
    pub fn new(alpha: f64, beta: f64) -> Self {
        assert!(alpha > 0.0, "Alpha ({}) must be positive", alpha);
        assert!(beta > 0.0, "Beta ({}) must be positive", beta);
        Self {
            alpha, beta,
            ln_beta_ab: ln_beta(alpha, beta),
        }
    }

    /// Returns alpha parameter.
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// Returns beta parameter.
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// Calculates variable part of PMF(k, n; alpha, beta), without the constant terms.
    #[inline]
    fn ln_pmf_inner(&self, k: f64, n: f64) -> f64 {
        -ln_beta(n - k + 1.0, k + 1.0) + ln_beta(k + self.alpha, n - k + self.beta)
    }

    /// Calculates ln-PMF(k, n).
    pub fn ln_pmf(&self, k: u32, n: u32) -> f64 {
        let n = f64::from(n);
        self.ln_pmf_inner(f64::from(k), n) - (n + 1.0).ln() - self.ln_beta_ab
    }

    /// Returns largest possible `k`, such that `CDF(k) <= cdf`.
    pub fn inv_cdf(&self, n: u32, cdf: f64) -> u32 {
        let m = f64::from(n);
        let const_term = -(m + 1.0).ln() - self.ln_beta_ab;

        // PDF(0). k = 0 will be always applicable, even if it is over p-val threshold.
        let mut ln_cdf = -ln_beta(m + 1.0, 1.0) + ln_beta(self.alpha, m + self.beta) + const_term;
        for i in 0..n {
            let k = f64::from(i + 1);
            ln_cdf = Ln::add(ln_cdf, self.ln_pmf_inner(k, m) + const_term);
            // Next CDF is out of bounds, return i.
            if ln_cdf.exp() > cdf {
                return i;
            }
        }
        n
    }

    /// Returns two inverse-CDF values, `argmin CDF(k) <= cdf1` and `argmin CDF(k) <= cdf2`.
    /// `cdf1` must not be greater than `cdf2`.
    pub fn inv_cdf2(&self, n: u32, cdf1: f64, cdf2: f64) -> (u32, u32) {
        let m = f64::from(n);
        let const_term = -(m + 1.0).ln() - self.ln_beta_ab;

        // PDF(0). k = 0 will be always applicable, even if it is over p-val threshold.
        let mut ln_cdf = -ln_beta(m + 1.0, 1.0) + ln_beta(self.alpha, m + self.beta) + const_term;
        let mut k1 = n;
        for i in 0..n {
            let k = f64::from(i + 1);
            ln_cdf = Ln::add(ln_cdf, self.ln_pmf_inner(k, m) + const_term);
            // Next CDF is out of bounds, return i.
            if ln_cdf.exp() > cdf1 {
                k1 = i;
                break;
            }
        }
        if ln_cdf.exp() > cdf2 {
            return (k1, k1);
        }
        for i in (k1 + 1)..n {
            let k = f64::from(i + 1);
            ln_cdf = Ln::add(ln_cdf, self.ln_pmf_inner(k, m) + const_term);
            // Next CDF is out of bounds, return i.
            if ln_cdf.exp() > cdf2 {
                return (k1, i);
            }
        }
        (k1, n)
    }

    /// Use Nelder-Mead to find max-likelihood `alpha` & `beta` for an observed vector of tuples `(k, n, weight)`.
    pub fn max_lik_estimate(observations: &[(u32, u32, f64)]) -> Self {
        let start_points = vec![
            vec![0.7,  50.0],
            vec![0.3, 100.0],
            vec![0.5,  10.0],
        ];
        let solver = NelderMead::new(start_points)
            .with_sd_tolerance(1e-6).unwrap();
        let solution = Executor::new(MLEProblem { observations }, solver)
            .run()
            .expect("Nelder-Mead finished with an error");
        let param = solution.state.get_param().expect("Nelder-Mead failed to find appropriate N.Binom. parameters");
        Self::new(param[0], param[1])
    }
}

struct MLEProblem<'a> {
    observations: &'a [(u32, u32, f64)],
}

impl<'a> CostFunction for MLEProblem<'a> {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        const LARGE_NUM: f64 = 1e30;
        const PARAM_LIMIT: f64 = 100000.0;
        let alpha = param[0];
        let beta = param[1];
        if alpha <= 0.0 || beta <= 0.0 || alpha >= PARAM_LIMIT || beta >= PARAM_LIMIT {
            Ok(LARGE_NUM)
        } else {
            let beta_binom = BetaBinomial::new(alpha, beta);
            Ok(-self.observations.iter()
                .map(|&(k, n, weight)| weight * beta_binom.ln_pmf(k, n))
                .sum::<f64>())
        }
    }
}

impl fmt::Debug for BetaBinomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "BetaBinom({}, {})", self.alpha, self.beta)
    }
}
