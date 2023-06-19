use crate::math::Ln;
use super::DiscretePmf;

/// Bayesian probability calculator for the main hypothesis and several alternative hypotheses.
/// Here, all priors are equal to 1.
pub struct BayesCalc<T, U, const N: usize> {
    null_hypoth: T,
    alternative: [U; N],
}

impl<T, U, const N: usize> BayesCalc<T, U, N> {
    pub fn new(null_hypoth: T, alternative: [U; N]) -> Self {
        Self { null_hypoth, alternative }
    }
}

impl<T: DiscretePmf, U: DiscretePmf, const N: usize> DiscretePmf for BayesCalc<T, U, N> {
    fn ln_pmf(&self, k: u32) -> f64 {
        let null_prob = self.null_hypoth.ln_pmf(k);
        let mut sum_prob = null_prob;
        for distr in self.alternative.iter() {
            sum_prob = Ln::add(sum_prob, distr.ln_pmf(k));
        }
        null_prob - sum_prob
    }
}
