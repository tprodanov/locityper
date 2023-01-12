//! Traits and structures related to insert size (distance between read mates).

use htslib::bam::record::{Record, Aux};

use crate::{
    algo::{
        vec_ext::*,
        nbinom::{NBinom, CachedDistr},
    },
    bg::ser::{JsonSer, LoadError},
};

/// Trait for insert size distribution.
pub trait InsertDistr {
    /// Maximal insert size. Do not allow read mates with distance over `max_size`.
    fn max_size(&self) -> u32;

    /// Ln-probability of the insert size.
    fn ln_prob(&self, sz: u32) -> f64;
}

/// Get an absolute value of the insert size.
/// Returns None if the read is unmapped, unpaired, two mates are on diff chromosomes, or
/// any of the mates have low mapping quality (< 30).
pub fn get_insert_size(record: &Record) -> Option<u32> {
    fn mate_mapq(record: &Record) -> u8 {
        match record.aux(b"MQ") {
            Ok(Aux::U8(val)) => val,
            Ok(Aux::I8(val)) => val as u8,
            Ok(Aux::I16(val)) => val as u8,
            Ok(Aux::U16(val)) => val as u8,
            Ok(Aux::I32(val)) => val as u8,
            Ok(Aux::U32(val)) => val as u8,
            Ok(_) => panic!("BAM record tag MQ has non-integer value!"),
            Err(_) => panic!("BAM record does not have a MQ tag!"),
        }
    }
    const MIN_MAPQ: u8 = 30;
    // Both mates mapped, primary alignments.
    if (record.flags() & 3980) == 0 && record.tid() == record.mtid() &&
            record.mapq() >= MIN_MAPQ && mate_mapq(record) >= MIN_MAPQ {
        Some(record.insert_size().abs() as u32)
    } else {
        None
    }
}

/// Negative Binomial insert size.
#[derive(Debug, Clone)]
pub struct InsertNegBinom {
    max_size: u32,
    distr: CachedDistr<NBinom>,
}

impl InsertNegBinom {
    /// Creates the Neg. Binom. insert size distribution from an iterator of insert sizes.
    pub fn estimate<T, I>(insert_sizes: I) -> Self
    where T: Into<f64>,
          I: Iterator<Item = T>,
    {
        log::info!("    Estimating insert size distribution");
        const QUANTILE: f64 = 0.99;
        const QUANT_MULT: f64 = 2.0;
        // Calculate max_size from input values as 2.0 * <99-th quantile>.
        // This is needed to remove read mates that were mapped to the same chromosome but very far apart.

        let mut insert_sizes: Vec<f64> = insert_sizes.map(T::into).collect();
        insert_sizes.sort();
        let max_size = QUANT_MULT * insert_sizes.quantile_sorted(QUANTILE);
        // Find index after the limiting value.
        let m = insert_sizes.binary_search_right(&max_size);
        let lim_insert_sizes = &insert_sizes[..m];
        let max_size = max_size.ceil() as u32;

        let mean = lim_insert_sizes.mean();
        // Increase variance, if less-equal than mean.
        let var = lim_insert_sizes.variance(Some(mean)).max(1.000001 * mean);
        log::info!("        Insert size mean = {:.1},  st.dev. = {:.1}", mean, var.sqrt());
        log::info!("        Treat reads with insert size > {} as unpaired", max_size);
        Self {
            max_size,
            distr: NBinom::estimate(mean, var).cached(max_size as usize),
        }
    }
}

impl InsertDistr for InsertNegBinom {
    fn max_size(&self) -> u32 {
        self.max_size
    }

    fn ln_prob(&self, sz: u32) -> f64 {
        self.distr.ln_pmf(sz as u64)
    }
}

impl JsonSer for InsertNegBinom {
    fn save(&self) -> json::JsonValue {
        json::object!{
            max_size: self.max_size,
            n: self.distr.distr().n(),
            p: self.distr.distr().p(),
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        let max_size = obj["max_size"].as_usize().ok_or_else(|| LoadError(format!(
            "NBinom: Failed to parse '{}': missing or incorrect 'max_size' field!", obj)))?;
        let n = obj["n"].as_f64().ok_or_else(|| LoadError(format!(
            "NBinom: Failed to parse '{}': missing or incorrect 'n' field!", obj)))?;
        let p = obj["p"].as_f64().ok_or_else(|| LoadError(format!(
            "NBinom: Failed to parse '{}': missing or incorrect 'p' field!", obj)))?;
        Ok(Self {
            max_size: max_size as u32,
            distr: NBinom::new(n, p).cached(max_size),
        })
    }
}