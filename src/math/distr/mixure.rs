use crate::math::Ln;
use super::DiscretePmf;

/// Mixure of two distributions with different weights.
#[derive(Clone)]
pub struct Mixure<T, U> {
    /// First distribution.
    distr1: T,
    /// ln-weight of the first distribution.
    lnw1: f64,
    /// Second distribution.
    distr2: U,
    /// ln-weight of the second distribution.
    lnw2: f64,
}

impl<T, U> Mixure<T, U> {
    /// Creates a new Mixure distribution from two distributions and the weight of the first one.
    pub fn new(distr1: T, weight1: f64, distr2: U) -> Self {
        assert!(weight1 > 0.0 && weight1 < 1.0,
            "Cannot create Mixure distribution: weight ({:.5}) must be within (0, 1)", weight1);
        Self {
            distr1, distr2,
            lnw1: weight1.ln(),
            lnw2: (-weight1).ln_1p(),
        }
    }
}

impl<T: DiscretePmf, U: DiscretePmf> DiscretePmf for Mixure<T, U> {
    fn ln_pmf(&self, k: u32) -> f64 {
        Ln::add(self.lnw1 + self.distr1.ln_pmf(k), self.lnw2 + self.distr2.ln_pmf(k))
    }
}