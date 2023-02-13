use std::ops::Deref;

pub mod uniform;
pub mod nbinom;
pub mod mnom;
pub mod mixure;
pub mod chimeric;
pub mod lincache;

pub use uniform::Uniform;
pub use nbinom::NBinom;
pub use mnom::Multinomial;
pub use mixure::Mixure;
pub use chimeric::Chimeric;
pub use lincache::LinearCache;

/// Discrete distribution with a single argument.
pub trait DiscretePmf {
    /// Returns ln-probability of `k`.
    fn ln_pmf(&self, k: u32) -> f64;
}

/// Distribution has mean and variance values.
pub trait WithMoments {
    fn mean(&self) -> f64;

    fn variance(&self) -> f64;
}

/// Discrete distribution with CDF, SF and quantile functions.
pub trait DiscreteCdf: WithMoments {
    /// Cumulative distribution function: P(X <= k).
    fn cdf(&self, k: u32) -> f64;

    /// Survival function: P(X > k).
    fn sf(&self, k: u32) -> f64 {
        1.0 - self.cdf(k)
    }

    /// Inverse CDF function: returns approximate such `x` that `cdf(x) = q`.
    fn quantile(&self, q: f64) -> f64 {
        if q <= 0.0 {
            return 0.0;
        }
        if q >= 1.0 {
            return f64::INFINITY;
        }
        let mut low = 0;
        let mut high = (2.0 * self.mean()) as u32;
        while self.cdf(high) < q {
            low = high;
            high *= 2;
        }
        while high >= low {
            let mid = (low + high) / 2;
            if self.cdf(mid) >= q {
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }
        let cdf0 = self.cdf(high);
        assert!(cdf0 <= q);
        let cdf1 = self.cdf(high + 1);
        let diff = cdf1 - cdf0;
        let x = f64::from(high);
        match diff.total_cmp(&0.0) {
            std::cmp::Ordering::Equal => x,
            std::cmp::Ordering::Greater => {
                let r = (q - cdf0) / diff;
                x * (1.0 - r) + (x + 1.0) * r
            }
            _ => panic!("Assertion failed: CDF is non-increasing!"),
        }
    }
}

impl<T: DiscretePmf + ?Sized> DiscretePmf for std::rc::Rc<T> {
    fn ln_pmf(&self, k: u32) -> f64 {
        self.deref().ln_pmf(k)
    }
}