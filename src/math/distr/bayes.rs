use crate::math::Ln;
use super::DiscretePmf;

/// Bayesian probability calculator for the main hypothesis and several alternative hypotheses.
/// Here, all priors are equal to 1.
#[derive(Clone)]
pub struct BayesCalc<T, U, const N: usize> {
    null_hypoth: T,
    alternative: [U; N],
}

impl<T, U, const N: usize> BayesCalc<T, U, N> {
    pub fn new(null_hypoth: T, alternative: [U; N]) -> Self {
        Self { null_hypoth, alternative }
    }
}

impl<T: DiscretePmf + Clone, U: DiscretePmf + Clone, const N: usize> DiscretePmf for BayesCalc<T, U, N> {
    fn ln_pmf(&self, k: u32) -> f64 {
        let null_prob = self.null_hypoth.ln_pmf(k);
        let mut probs = [0.0; N];
        for (prob, distr) in probs.iter_mut().zip(self.alternative.iter()) {
            *prob = distr.ln_pmf(k);
        }
        let sum_prob = Ln::sum_init(&probs, null_prob);
        null_prob - sum_prob
    }
}
