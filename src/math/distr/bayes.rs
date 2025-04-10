use crate::math::Ln;
use super::DiscretePmf;

pub const N_ALTS: usize = 16;

/// Bayesian probability calculator for the main hypothesis and several alternative hypotheses.
/// Here, all priors are equal to 1.
#[derive(Clone)]
pub struct BayesCalc<T, U> {
    null_distr: T,
    alternatives: Vec<U>,
}

impl<T, U> BayesCalc<T, U> {
    pub fn new(null_distr: T, alternatives: Vec<U>) -> Self {
        assert!(alternatives.len() < N_ALTS);
        Self { null_distr, alternatives }
    }

    /// Returns distribution according to the null hypothesis.
    pub fn null_distr(&self) -> &T {
        &self.null_distr
    }
}

impl<T: DiscretePmf + Clone, U: DiscretePmf + Clone> DiscretePmf for BayesCalc<T, U> {
    fn ln_pmf(&self, k: u32) -> f64 {
        let null_prob = self.null_distr.ln_pmf(k);
        let mut probs = [0.0; N_ALTS];
        for (prob, distr) in probs.iter_mut().zip(self.alternatives.iter()) {
            *prob = distr.ln_pmf(k);
        }
        let sum_prob = Ln::sum_init(&probs[..self.alternatives.len()], null_prob);
        null_prob - sum_prob
    }
}
