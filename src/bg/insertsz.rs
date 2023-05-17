//! Traits and structures related to insert size (distance between read mates).

use std::ops::Deref;
use fnv::FnvHashMap;
use crate::{
    Error,
    algo::bisect,
    ext::vec::{VecExt, F64Ext},
    bg::ser::{JsonSer, json_get},
    math::distr::{DiscretePmf, NBinom, LinearCache},
    seq::aln::Alignment,
};

/// Group read alignments in pairs, and returns indices of the pair in the original vector.
/// Only full pairs (both mates are present in some specific region) are stored.
///
/// Input: vector of primary alignments.
/// There can be at most one first and one second mate for each name.
pub fn group_mates<'a>(alns: &'a [Alignment]) -> Result<Vec<(usize, usize)>, Error> {
    let mut pairs: FnvHashMap<&'a [u8], [Option<usize>; 2]> =
        FnvHashMap::with_capacity_and_hasher(alns.len() / 2, Default::default());
    for (i, aln) in alns.iter().enumerate() {
        let ix = aln.read_end().ix();
        // Insert alignment and return Error, if there already was an alignment there.
        if pairs.entry(aln.name()).or_default()[ix].replace(i).is_some() {
            return Err(Error::InvalidData(format!("Read {} has several {} mates",
                aln.name_utf8(), if ix == 0 { "first" } else { "second" })));
        }
    }
    Ok(pairs.into_values().filter_map(|[opt_i, opt_j]| opt_i.zip(opt_j)).collect())
}

/// Allow read-pair orientation, if at least 5% pairs support it.
const ORIENT_THRESH: f64 = 0.05;

/// Insert size distribution.
#[derive(Debug, Clone)]
pub struct InsertDistr {
    /// Maximum allowed insert size. Higher values wil automatically produce -inf.
    max_size: u32,
    /// Are different read-pair orientation allowed? (FR/RF v. RR/FF)
    orient_allowed: [bool; 2],
    /// Cached insert size distribution.
    distr: Option<LinearCache<NBinom>>,
    /// Highest insert size ln-probability, achievable at this distribution.
    mode_prob: f64,
}

/// Counts reads with insert size over 500 kb as certainly unpaired.
const MAX_REASONABLE_INSERT: u32 = 500_000;

impl InsertDistr {
    pub fn undefined() -> Self {
        Self {
            max_size: 0,
            orient_allowed: [false; 2],
            distr: None,
            mode_prob: f64::NAN,
        }
    }

    /// Creates the Neg. Binom. insert size distribution from an iterator of insert sizes.
    ///
    // Calculate max insert size from input reads as `<quantile_mult> * <quantile>-th insert size quantile`.
    // This is needed to remove read mates that were mapped to the same chromosome but very far apart.
    pub fn estimate(
        alignments: &[Alignment],
        pair_ixs: &[(usize, usize)],
        params: &super::Params,
    ) -> Result<Self, Error>
    {
        if pair_ixs.is_empty() {
            return Err(Error::InvalidData("BAM records are supposed to be paired!".to_owned()));
        } else if pair_ixs.len() < 1000 {
            return Err(Error::InvalidData("Not enough paired reads to calculate insert size distribution".to_owned()));
        }
        log::info!("Estimating insert size distribution");
        let mut insert_sizes = Vec::<f64>::new();
        let mut orient_counts = [0_u64; 2];
        for (i, j) in pair_ixs.iter() {
            let first = &alignments[*i].deref();
            let second = &alignments[*j].deref();
            let insert_size = first.insert_size(second);
            if insert_size < MAX_REASONABLE_INSERT {
                insert_sizes.push(f64::from(insert_size));
                orient_counts[usize::from(first.pair_orientation(second))] += 1;
            }
        }

        // First: FR/RF, Second: FF/RR.
        let total = (orient_counts[0] + orient_counts[1]) as f64;
        let orient_frac = [orient_counts[0] as f64 / total, orient_counts[1] as f64 / total];
        let orient_allowed = [orient_frac[0] >= ORIENT_THRESH, orient_frac[1] >= ORIENT_THRESH];
        log::info!("    Analyzed {} read pairs", total);
        for i in 0..2 {
            log::info!("    {}: {:8} ({:.3}%),  {}",
                if i == 0 { "FR/RF" } else { "FF/RR" },
                orient_counts[i], orient_frac[i],
                if orient_allowed[i] { "allowed" } else { "forbidden" });
        }

        VecExt::sort(&mut insert_sizes);
        let max_size = params.ins_quantile_mult * F64Ext::quantile_sorted(&insert_sizes, params.ins_quantile);
        // Find index after the limiting value.
        let m = bisect::right(&insert_sizes, &max_size);
        let lim_insert_sizes = &insert_sizes[..m];
        let max_size = max_size.ceil() as u32;

        let mean = F64Ext::mean(lim_insert_sizes);
        // Increase variance, if less-equal than mean.
        let var = F64Ext::variance(lim_insert_sizes, Some(mean)).max(1.000001 * mean);
        log::info!("    Insert size mean = {:.1},  st.dev. = {:.1}", mean, var.sqrt());
        log::info!("    Treat reads with insert size > {} as unpaired", max_size);
        let distr = NBinom::estimate(mean, var).cached(max_size as usize + 1);
        let mode_prob = distr.ln_pmf(distr.inner().mode());

        Ok(Self {
            max_size, orient_allowed, mode_prob,
            distr: Some(distr),
        })
    }

    /// Ln-probability of the insert size. `same_orient` is true if FF/RR, false if FR/RF.
    pub fn ln_prob(&self, sz: u32, same_orient: bool) -> f64 {
        if sz <= self.max_size && self.orient_allowed[usize::from(same_orient)] {
            self.distr.as_ref().unwrap().ln_pmf(sz)
        } else {
            f64::NEG_INFINITY
        }
    }

    /// Maximum insert size. Over this size, all pairs are considered unpaired.
    pub fn max_size(&self) -> u32 {
        self.max_size
    }

    /// Returns true if the reads are paired-end, false if single-end.
    pub fn is_paired_end(&self) -> bool {
        self.distr.is_some()
    }

    /// Maximum achievable insert size ln-probability.
    pub fn mode_prob(&self) -> f64 {
        self.mode_prob
    }
}

impl JsonSer for InsertDistr {
    fn save(&self) -> json::JsonValue {
        if let Some(distr) = &self.distr {
            json::object!{
                max_size: self.max_size,
                fr_allowed: self.orient_allowed[0],
                ff_allowed: self.orient_allowed[1],
                n: distr.inner().n(),
                p: distr.inner().p(),
            }
        } else {
            json::object!{
                max_size: self.max_size,
            }
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> max_size (as_usize));
        if max_size == 0 {
            return Ok(Self::undefined());
        }
        json_get!(obj -> n (as_f64), p (as_f64), fr_allowed (as_bool), ff_allowed (as_bool));
        let distr = NBinom::new(n, p).cached(max_size + 1);
        let mode_prob = distr.ln_pmf(distr.inner().mode());
        Ok(Self {
            max_size: u32::try_from(max_size).unwrap_or(u32::MAX / 2),
            orient_allowed: [fr_allowed, ff_allowed],
            distr: Some(distr),
            mode_prob,
        })
    }
}
