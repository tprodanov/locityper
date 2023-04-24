use std::fmt;
use statrs::function::gamma::ln_gamma;
use once_cell::sync::OnceCell;
use crate::{
    Error,
    math::Ln,
    bg::ser::JsonSer,
};

/// Calculates ln-factorial. Stores up to 1024 values in a cache.
fn ln_factorial(n: u32) -> f64 {
    const CACHE_SIZE: usize = 1024;
    const EMPTY_CELL: OnceCell<f64> = OnceCell::new();
    static FACTORIAL_CACHE: [OnceCell<f64>; CACHE_SIZE] = [EMPTY_CELL; CACHE_SIZE];

    let i = n as usize;
    if i < CACHE_SIZE {
        *FACTORIAL_CACHE[i].get_or_init(|| ln_gamma(f64::from(n + 1)))
    } else {
        ln_gamma(f64::from(n + 1))
    }
}

/// Multinomial distribution. `n` is not fixed and is selected based on the input counts.
#[derive(Clone, Debug, Default)]
pub struct Multinomial {
    lnp: Vec<f64>,
}

impl Multinomial {
    /// Creates Multinomial distribution based on probabilities (may not sum up to one).
    pub fn new(p: &[f64]) -> Self {
        let mut lnp: Vec<f64> = p.iter().copied().map(f64::ln).collect();
        let s = Ln::sum(&lnp);
        assert!(s.is_finite(), "Cannot create Multinomial distribution from probabilities ({:?})", p);
        lnp.iter_mut().for_each(|v| *v -= s);
        Multinomial { lnp }
    }

    /// Returns ln-probabilities of all events.
    pub fn ln_probabilities(&self) -> &[f64] {
        &self.lnp
    }

    /// Calculates ln-probability of the counts `ks` assuming that `sum(ks)` is fixed.
    pub fn ln_pmf(&self, ks: &[u32]) -> f64 {
        assert_eq!(self.lnp.len(), ks.len(),
            "Multinomial distribution: number of probabilities does not match the number of classes");
        let mut res = 0.0;
        let mut n = 0;
        for (&lnp, &k) in self.lnp.iter().zip(ks) {
            n += k;
            res += f64::from(k) * lnp - ln_factorial(k);
        }
        res + ln_factorial(n)
    }
}

impl fmt::Display for Multinomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Multinomial(p = {:?})", self.lnp.iter().copied().map(f64::exp).collect::<Vec<_>>())
    }
}

impl JsonSer for Multinomial {
    fn save(&self) -> json::JsonValue {
        json::object!{
            lnp: &self.lnp as &[f64],
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        let mut lnp = Vec::new();
        if let json::JsonValue::Array(arr) = &obj["lnp"] {
            for val in arr.into_iter() {
                let valf = val.as_f64().ok_or_else(||
                    Error::JsonLoad(format!("Failed to parse '{}': array '{}' has non-float element", obj, "lnp")))?;
                lnp.push(valf);
            }
            Ok(Self { lnp })
        } else {
            Err(Error::JsonLoad(format!("Failed to parse '{}': missing or incorrect array '{}'", obj, "lnp")))
        }
    }
}
