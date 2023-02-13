use std::fmt::Debug;
use super::{DiscretePmf, DiscreteCdf};

/// Chimeric distribution: one distribution before partition, and another after.
#[derive(Clone)]
pub struct Chimeric<T, U> {
    partition: u32,
    distr1: T,
    lnw1: f64,
    distr2: U,
    lnw2: f64,
}

impl<T, U> Chimeric<T, U>
where T: DiscretePmf + DiscreteCdf + Debug,
      U: DiscretePmf + DiscreteCdf + Debug,
{
    /// Creates a new distribution such that for `k <= partition`, distribution is `distr1`,
    /// and for `k > partition`, distribution is `distr2`.
    ///
    /// Weights are set such that `distr1.pmf(partition) / distr2.pmf(partition) = ratio`.
    pub fn new(partition: u32, distr1: T, distr2: U, ratio: f64) -> Self {
        let val1 = distr1.ln_pmf(partition).exp();
        let val2 = ratio * distr2.ln_pmf(partition).exp();
        assert!(val1 > 0.0 && val2 > 0.0,
            "Cannot create Chimeric distribution ({:?}, {:?}, partition={}), coeffient or PMF is 0",
            distr1, distr2, partition);
        let cdf1 = distr1.cdf(partition);
        let sf2 = distr2.sf(partition);

        let weight1 = val2 / (cdf1 * val2 + sf2 * val1);
        let weight2 = val1 * weight1 / val2;
        Self {
            partition, distr1, distr2,
            lnw1: weight1.ln(),
            lnw2: weight2.ln(),
        }
    }

    /// Creates a new Chimeric distribution such that PMFs at the partition point are equal.
    pub fn new_smooth(partition: u32, distr1: T, distr2: U) -> Self {
        Self::new(partition, distr1, distr2, 1.0)
    }
}

impl<T: DiscretePmf, U: DiscretePmf> DiscretePmf for Chimeric<T, U> {
    fn ln_pmf(&self, k: u32) -> f64 {
        if k <= self.partition {
            self.lnw1 + self.distr1.ln_pmf(k)
        } else {
            self.lnw2 + self.distr2.ln_pmf(k)
        }
    }
}
