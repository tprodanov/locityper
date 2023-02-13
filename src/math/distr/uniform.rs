use std::{
    cmp::{min, max},
    fmt,
};
use super::{DiscretePmf, WithMoments, DiscreteCdf};

/// Discrete uniform distribution with values in [a, b] (both inclusive).
#[derive(Clone)]
pub struct Uniform {
    a: u32,
    b: u32,
    lnp: f64,
}

impl Uniform {
    /// Creates a new discrete uniform distribution with values in [a, b] (both inclusive).
    pub fn new(a: u32, b: u32) -> Self {
        Self {
            a, b,
            lnp: -f64::from(b - a + 1).ln(),
        }
    }
}

impl DiscretePmf for Uniform {
    fn ln_pmf(&self, k: u32) -> f64 {
        if k < self.a || k > self.b {
            f64::NEG_INFINITY
        } else {
            self.lnp
        }
    }
}

impl WithMoments for Uniform {
    fn mean(&self) -> f64 {
        0.5 * (f64::from(self.a) + f64::from(self.b))
    }

    fn variance(&self) -> f64 {
        let n = f64::from(self.b - self.a + 1);
        (n * n - 1.0) / 12.0
    }
}

impl DiscreteCdf for Uniform {
    fn cdf(&self, k: u32) -> f64 {
        f64::from((min(self.b, k) + 1).saturating_sub(self.a)) / f64::from(self.b - self.a + 1)
    }

    fn sf(&self, k: u32) -> f64 {
        f64::from(self.b.saturating_sub(max(k, self.a.saturating_sub(1)))) / f64::from(self.b - self.a + 1)
    }
}

impl fmt::Debug for Uniform {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Uniform{{{}, {}}}", self.a, self.b)
    }
}

impl fmt::Display for Uniform {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}