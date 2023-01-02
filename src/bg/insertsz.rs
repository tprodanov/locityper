//! Traits and structures related to insert size (distance between read mates).

use std::cmp::min;
use std::io::{Write, Result};
use statrs::distribution::{NegativeBinomial, Discrete};

use crate::algo::vec_ext::*;

/// Trait for insert size distribution.
pub trait Distr {
    /// Maximal insert size. Do not allow read mates with distance over `max_size`.
    fn max_size(&self) -> u32;

    /// Ln-probability of the insert size.
    fn ln_prob(&self, sz: u32) -> f64;

    /// Save distribution into a file/stream.
    fn save<W: Write>(&self, f: W) -> Result<()>;
}

/// Negative Binomial insert size.
#[derive(Debug, Clone)]
pub struct NegBinom {
    max_size: u32,
    distr: NegativeBinomial,
}

impl NegBinom {
    /// Creates the Neg. Binom. insert size distribution from an iterator of insert sizes.
    pub fn create<T, I>(insert_sizes: I) -> Self
    where T: Into<f64>,
          I: Iterator<Item = T>,
    {
        const QUANTILE: f64 = 0.99;
        const QUANT_MULT: f64 = 2.0;
        // Calculate max_size from input values as 2.0 * <99-th quantile>.
        // This is needed to remove read mates that were mapped to the same chromosome but very far apart.

        let mut insert_sizes: Vec<f64> = insert_sizes.map(T::into).collect();
        insert_sizes.sort();
        let n = insert_sizes.len();
        let i = min(n - 1, (n as f64 * QUANTILE) as usize);
        let max_size = QUANT_MULT * insert_sizes[i];
        // Find index after the limiting value.
        let m = insert_sizes.binary_search_right_at(&max_size, i, n);
        let lim_insert_sizes = &insert_sizes[..m];

        let mean = lim_insert_sizes.mean();
        // Increase variance, if less-equal than mean.
        let var = lim_insert_sizes.variance(Some(mean)).max(1.000001 * mean);
        let r = mean * mean / (var - mean);
        let p = 1.0 - mean / var;

        NegBinom {
            max_size: max_size.ceil() as u32,
            distr: NegativeBinomial::new(r, p).unwrap(),
        }
    }
}

impl Distr for NegBinom {
    fn max_size(&self) -> u32 {
        self.max_size
    }

    fn ln_prob(&self, sz: u32) -> f64 {
        self.distr.ln_pmf(sz as u64)
    }

    fn save<W: Write>(&self, mut f: W) -> Result<()> {
        f.write(b"Fragment size distribution: Negative Binomial\n")?;
        write!(f, "max_size: {}\nr: {:.?}\np: {:?}\n", self.max_size, self.distr.r(), self.distr.p())
    }
}