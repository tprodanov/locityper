//! Traits and structures related to insert size (distance between read mates).

use htslib::bam::record::{Record, Aux};
use crate::{
    algo::{
        vec_ext::{VecExt, F64Ext},

        bisect,
    },
    bg::ser::{JsonSer, LoadError},
    math::distr::{DiscretePmf, NBinom, LinearCache},
};

/// Trait for insert size distribution.
pub trait InsertDistr {
    /// Ln-probability of the insert size. `same_orient` is true if FF/RR, false if FR/RF.
    fn ln_prob(&self, sz: u32, same_orient: bool) -> f64;

    /// Maximum insert size. Over this size, all pairs are deemed unpaired.
    fn max_size(&self) -> u32;
}

/// Get MAPQ of the mate record.
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

/// Returns true if FF/RR orientation, false if FR/RF.
#[inline]
fn pair_orientation(record: &Record) -> bool {
    let flag = record.flags();
    (flag & 0x10) == (flag & 0x20)
}

/// Negative Binomial insert size.
#[derive(Debug, Clone)]
pub struct InsertNegBinom {
    max_size: u32,
    /// Log-probabilities of (FR/RF) orientation, and of (RR/FF) orientation.
    orient_probs: [f64; 2],
    distr: LinearCache<NBinom>,
}

impl InsertNegBinom {
    /// Creates the Neg. Binom. insert size distribution from an iterator of insert sizes.
    pub fn estimate<'a>(records: impl Iterator<Item = &'a Record>) -> Self {
        const MIN_MAPQ: u8 = 30;
        const QUANTILE: f64 = 0.99;
        const QUANT_MULT: f64 = 3.0;
        // Calculate max_size from input values as 3.0 * <99-th quantile>.
        // This is needed to remove read mates that were mapped to the same chromosome but very far apart.

        log::info!("    Estimating insert size distribution");
        let mut insert_sizes = Vec::<f64>::new();
        let mut orient_counts = [0_u64; 2];
        for record in records {
            // 3980 - ignore reads with any mate unmapped, ignore sec./supp. alignments, ignore second mates.
            if (record.flags() & 3980) == 0 && record.tid() == record.mtid() &&
                    record.mapq() >= MIN_MAPQ && mate_mapq(record) >= MIN_MAPQ {
                insert_sizes.push(record.insert_size().abs() as f64);
                orient_counts[pair_orientation(record) as usize] += 1;
            }
        }
        let total = (orient_counts[0] + orient_counts[1]) as f64;
        log::info!("        Analyzed {} read pairs", total);
        log::info!("        FR/RF: {} ({:.3}%)   FF/RR: {} ({:.3}%)",
            orient_counts[0], 100.0 * orient_counts[0] as f64 / total,
            orient_counts[1], 100.0 * orient_counts[1] as f64 / total);
        // Use ln1p in order to fix 0-probabilities in case of low number of reads.
        let orient_probs = [
            (orient_counts[0] as f64).ln_1p() - total.ln_1p(),
            (orient_counts[1] as f64).ln_1p() - total.ln_1p()];

        VecExt::sort(&mut insert_sizes);
        let max_size = QUANT_MULT * F64Ext::quantile_sorted(&insert_sizes, QUANTILE);
        // Find index after the limiting value.
        let m = bisect::right(&insert_sizes, &max_size);
        let lim_insert_sizes = &insert_sizes[..m];
        let max_size = max_size.ceil() as u32;

        let mean = F64Ext::mean(lim_insert_sizes);
        // Increase variance, if less-equal than mean.
        let var = F64Ext::variance(lim_insert_sizes, Some(mean)).max(1.000001 * mean);
        log::info!("        Insert size mean = {:.1},  st.dev. = {:.1}", mean, var.sqrt());
        log::info!("        Treat reads with insert size > {} as unpaired", max_size);
        Self {
            max_size, orient_probs,
            distr: NBinom::estimate(mean, var).cached(max_size as usize + 1),
        }
    }
}

impl InsertDistr for InsertNegBinom {
    fn ln_prob(&self, sz: u32, same_orient: bool) -> f64 {
        if sz > self.max_size {
            f64::NEG_INFINITY
        } else {
            self.distr.ln_pmf(sz) + self.orient_probs[same_orient as usize]
        }
    }

    fn max_size(&self) -> u32 {
        self.max_size
    }
}

impl JsonSer for InsertNegBinom {
    fn save(&self) -> json::JsonValue {
        json::object!{
            max_size: self.max_size,
            fr_prob: self.orient_probs[0],
            ff_prob: self.orient_probs[1],
            n: self.distr.distr().n(),
            p: self.distr.distr().p(),
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        let max_size = obj["max_size"].as_usize().ok_or_else(|| LoadError(format!(
            "InsertNegBinom: Failed to parse '{}': missing or incorrect 'max_size' field!", obj)))?;
        let n = obj["n"].as_f64().ok_or_else(|| LoadError(format!(
            "InsertNegBinom: Failed to parse '{}': missing or incorrect 'n' field!", obj)))?;
        let p = obj["p"].as_f64().ok_or_else(|| LoadError(format!(
            "InsertNegBinom: Failed to parse '{}': missing or incorrect 'p' field!", obj)))?;
        let orient_probs = [
            obj["fr_prob"].as_f64().ok_or_else(|| LoadError(format!(
                "InsertNegBinom: Failed to parse '{}': missing or incorrect 'fr_prob' field!", obj)))?,
            obj["ff_prob"].as_f64().ok_or_else(|| LoadError(format!(
                "InsertNegBinom: Failed to parse '{}': missing or incorrect 'ff_prob' field!", obj)))?,
        ];
        Ok(Self {
            max_size: u32::try_from(max_size).unwrap(),
            orient_probs,
            distr: NBinom::new(n, p).cached(max_size),
        })
    }
}