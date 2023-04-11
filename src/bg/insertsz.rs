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

/// Group read in pairs.
/// Only full pairs (both mates are present in some specific region) are stored.
pub struct ReadMateGrouping<'a> {
    pairs: Vec<(&'a Record, &'a Record)>,
}

/// Read is paired and is first mate.
const FIRST_AND_PAIRED: u16 = 65;
/// Read is unmapped, mate is unmapped, or alignment is secondary or supplementary, or does not pass some checks.
const BAD_READ: u16 = 3852;

/// If possible, puts a next possible first mate to `rec1`, and returns true.
/// Input `rec1` can also be used, if it is appropriate.
fn next_first_mate<'a>(rec1: &mut &'a Record, records: &mut impl Iterator<Item = &'a Record>) -> bool {
    // No flags from `BAD_READ`, all flags from `FIRST_AND_PAIRED`.
    if (rec1.flags() & (BAD_READ | FIRST_AND_PAIRED)) == FIRST_AND_PAIRED {
        return true;
    }
    while let Some(rec) = records.next() {
        if (rec.flags() & (BAD_READ | FIRST_AND_PAIRED)) == FIRST_AND_PAIRED {
            *rec1 = rec;
            return true;
        }
    }
    false
}

impl<'a> ReadMateGrouping<'a> {
    /// Group read mates from an unsorted BAM file, possibly filtered.
    /// This means that consecutive reads are either from the same mate,
    /// or one of the mates is missing (because it was discarded previously).
    ///
    /// Panics, if some of the reads are supplementary/secondary, or any of the mates are unmapped.
    pub fn from_unsorted_bam(mut records: impl Iterator<Item = &'a Record>, max_reads: Option<usize>) -> Self {
        let mut rec1 = match records.next() {
            Some(rec) => rec,
            None => return Self {
                pairs: Vec::new(),
            },
        };

        let mut pairs = if let Some(n) = max_reads { Vec::with_capacity(n / 2) } else { Vec::new() };
        loop {
            if !next_first_mate(&mut rec1, &mut records) {
                break;
            }
            let rec2 = match records.next() {
                Some(rec) => rec,
                None => break,
            };
            assert_eq!(rec2.flags() & BAD_READ, 0,
                "ReadMateGrouping failed: read {} has flag {}", String::from_utf8_lossy(rec2.qname()), rec2.flags());
            if rec2.is_last_in_template() && rec1.qname() == rec2.qname() {
                pairs.push((rec1, rec2));
                match records.next() {
                    Some(rec) => rec1 = rec,
                    None => break,
                }
            } else {
                rec1 = rec2;
            }
        }
        Self { pairs }
    }

    pub fn len(&self) -> usize {
        self.pairs.len()
    }
}

/// Insert size calculation parameters.
#[derive(Clone, Debug)]
pub struct InsertSizeParams {
    // Calculate max insert size from input reads as `<quantile_mult> * <quantile>-th insert size quantile`.
    // This is needed to remove read mates that were mapped to the same chromosome but very far apart.
    pub quantile: f64,
    pub quantile_mult: f64,
}

impl Default for InsertSizeParams {
    fn default() -> Self {
        Self {
            quantile: 0.99,
            quantile_mult: 3.0,
        }
    }
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
pub struct InsertDistr {
    max_size: u32,
    /// Log-probabilities of (FR/RF) orientation, and of (RR/FF) orientation.
    orient_probs: [f64; 2],
    distr: LinearCache<NBinom>,
}

/// Counts reads with insert size over 1Mb as certainly unpaired.
const MAX_REASONABLE_INSERT: f64 = 1e6;

impl InsertDistr {
    /// Creates the Neg. Binom. insert size distribution from an iterator of insert sizes.
    pub fn estimate<'a>(read_pairs: &ReadMateGrouping<'a>, params: &InsertSizeParams) -> Self {
        log::info!("    Estimating insert size distribution");
        let mut insert_sizes = Vec::<f64>::new();
        let mut orient_counts = [0_u64; 2];
        for (first, second) in read_pairs.pairs.iter() {
            let insert_size = first.insert_size().abs() as f64;
            if first.tid() == second.tid() && insert_size < MAX_REASONABLE_INSERT {
                insert_sizes.push(insert_size);
                orient_counts[pair_orientation(first) as usize] += 1;
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
        let max_size = params.quantile_mult * F64Ext::quantile_sorted(&insert_sizes, params.quantile);
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

    /// Ln-probability of the insert size. `same_orient` is true if FF/RR, false if FR/RF.
    pub fn ln_prob(&self, sz: u32, same_orient: bool) -> f64 {
        if sz > self.max_size {
            f64::NEG_INFINITY
        } else {
            self.distr.ln_pmf(sz) + self.orient_probs[same_orient as usize]
        }
    }

    /// Maximum insert size. Over this size, all pairs are deemed unpaired.
    pub fn max_size(&self) -> u32 {
        self.max_size
    }
}

impl JsonSer for InsertDistr {
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
            "InsertDistr: Failed to parse '{}': missing or incorrect 'max_size' field!", obj)))?;
        let n = obj["n"].as_f64().ok_or_else(|| LoadError(format!(
            "InsertDistr: Failed to parse '{}': missing or incorrect 'n' field!", obj)))?;
        let p = obj["p"].as_f64().ok_or_else(|| LoadError(format!(
            "InsertDistr: Failed to parse '{}': missing or incorrect 'p' field!", obj)))?;
        let orient_probs = [
            obj["fr_prob"].as_f64().ok_or_else(|| LoadError(format!(
                "InsertDistr: Failed to parse '{}': missing or incorrect 'fr_prob' field!", obj)))?,
            obj["ff_prob"].as_f64().ok_or_else(|| LoadError(format!(
                "InsertDistr: Failed to parse '{}': missing or incorrect 'ff_prob' field!", obj)))?,
        ];
        Ok(Self {
            max_size: u32::try_from(max_size).unwrap_or(u32::MAX / 2),
            orient_probs,
            distr: NBinom::new(n, p).cached(max_size),
        })
    }
}
