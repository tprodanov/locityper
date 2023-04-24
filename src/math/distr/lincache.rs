use std::fmt::{self, Debug, Display, Formatter};
use once_cell::unsync::OnceCell;
use crate::{
    Error,
    bg::ser::JsonSer,
};
use super::{DiscretePmf, DiscreteCdf, WithMoments};

/// Distribution with a fixed number of cached `ln_pmf` values.
#[derive(Clone)]
pub struct LinearCache<D> {
    distr: D,
    cache: Vec<OnceCell<f64>>,
}

impl<D> LinearCache<D> {
    /// Creates the cached distribution.
    /// Caches ln_pmf values in `0..cache_size`.
    pub fn new(distr: D, cache_size: usize) -> Self {
        Self {
            distr,
            cache: vec![OnceCell::new(); cache_size],
        }
    }

    /// Get inner distribution.
    pub fn distr(&self) -> &D {
        &self.distr
    }

    /// Returns the maximal cache size. Values 0..cache_size are stored once calculated the first time.
    pub fn cache_size(&self) -> usize {
        self.cache.len()
    }
}

impl<D: DiscretePmf> DiscretePmf for LinearCache<D> {
    /// Returns ln(pmf(k)). Caches values in a certain range.
    fn ln_pmf(&self, k: u32) -> f64 {
        let i = k as usize;
        if i < self.cache.len() {
            *self.cache[i].get_or_init(|| self.distr.ln_pmf(k))
        } else {
            self.distr.ln_pmf(k)
        }
    }
}

impl<D: DiscreteCdf> DiscreteCdf for LinearCache<D> {
    fn cdf(&self, k: u32) -> f64 {
        self.distr.cdf(k)
    }

    fn sf(&self, k: u32) -> f64 {
        self.distr.sf(k)
    }
}

impl<D: WithMoments> WithMoments for LinearCache<D> {
    fn mean(&self) -> f64 {
        self.distr.mean()
    }

    fn variance(&self) -> f64 {
        self.distr.variance()
    }
}

impl<D: Debug> Debug for LinearCache<D> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "Cached[{:?}, size={}]", self.distr, self.cache.len())
    }
}

impl<D: Display> Display for LinearCache<D> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "Cached[{}, size={}]", self.distr, self.cache.len())
    }
}

impl<D: JsonSer> JsonSer for LinearCache<D> {
    fn save(&self) -> json::JsonValue {
        json::object!{
            distr: self.distr.save(),
            size: self.cache.len(),
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        if !obj.has_key("distr") {
            return Err(Error::JsonLoad(format!("LinearCache: Failed to parse '{}': missing 'distr' field!", obj)));
        }
        let distr = D::load(&obj["distr"])?;
        let size = obj["size"].as_usize().ok_or_else(|| Error::JsonLoad(format!(
            "LinearCache: Failed to parse '{}': missing or incorrect 'size' field!", obj)))?;
        Ok(Self::new(distr, size))
    }
}
