//! Traits and structures related to insert size (distance between read mates).

use std::{
    path::Path,
    ops::Deref,
    io::Write,
};
use fnv::FnvHashMap;
use nohash::IntMap;
use crate::{
    Error,
    algo::bisect,
    ext::{self, vec::{VecExt, F64Ext}},
    bg::ser::{JsonSer, json_get},
    math::distr::{
        DiscretePmf, LinearCache, WithMoments, WithQuantile,
        nbinom::{NBinom, RegularizedEstimator},
    },
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
    /// Are different read-pair orientation allowed? (FR/RF v. RR/FF)
    orient_allowed: [bool; 2],
    /// Cached insert size distribution.
    distr: Option<LinearCache<NBinom>>,
    /// Minimum and maximum allowed insert sizes.
    conf_interval: (u32, u32),
}

impl InsertDistr {
    pub fn undefined() -> Self {
        Self {
            orient_allowed: [false; 2],
            distr: None,
            conf_interval: (0, 0),
        }
    }

    /// Creates the Neg. Binom. insert size distribution from an iterator of insert sizes.
    ///
    /// Calculate max insert size from input reads as `<quantile_mult> * <quantile>-th insert size quantile`.
    /// This is needed to remove read mates that were mapped to the same chromosome but very far apart.
    ///
    /// Write debug information if `out_dir` is Some.
    pub fn estimate(
        alignments: &[Alignment],
        pair_ixs: &[(usize, usize)],
        params: &super::Params,
        out_dir: Option<&Path>,
    ) -> Result<Self, Error>
    {
        // Counts reads with insert size over 500 kb as certainly unpaired.
        const MAX_REASONABLE_INSERT: u32 = 500_000;
        // Filter out insert sizes that are larger than `<INS_QUANTILE_MULT> * <INS_QUANTILE>-th quantile`.
        // This is needed to remove read mates that were mapped to the same chromosome but very far apart.
        const INS_QUANTILE: f64 = 0.99;
        const INS_QUANTILE_MULT: f64 = 3.0;

        if pair_ixs.is_empty() {
            return Err(Error::InvalidData("BAM records are supposed to be paired!".to_owned()));
        } else if pair_ixs.len() < 1000 {
            return Err(Error::InvalidData("Not enough paired reads to calculate insert size distribution".to_owned()));
        }
        log::info!("Estimating insert size distribution from {} read pairs", pair_ixs.len());
        let mut insert_sizes = Vec::<f64>::new();
        let mut orient_counts = [0_u64; 2];
        let mut histogram = IntMap::<u32, u32>::default();
        for (i, j) in pair_ixs.iter() {
            let first = &alignments[*i].deref();
            let second = &alignments[*j].deref();
            let insert_size = first.insert_size(second);
            if insert_size < MAX_REASONABLE_INSERT {
                insert_sizes.push(f64::from(insert_size));
                histogram.entry(insert_size)
                    .and_modify(|counter| *counter = counter.saturating_add(1))
                    .or_insert(1);
                orient_counts[usize::from(first.pair_orientation(second))] += 1;
            }
        }

        if let Some(dir) = out_dir {
            let mut dbg_writer = ext::sys::create_gzip(&dir.join("insert_sizes.csv.gz"))?;
            writeln!(dbg_writer, "size\tcount")?;
            let mut histogram: Vec<(u32, u32)> = histogram.into_iter().collect();
            histogram.sort_unstable();
            for (size, count) in histogram.into_iter() {
                writeln!(dbg_writer, "{}\t{}", size, count)?;
            }
        }

        // First: FR/RF, Second: FF/RR.
        let total = (orient_counts[0] + orient_counts[1]) as f64;
        let orient_frac = [orient_counts[0] as f64 / total, orient_counts[1] as f64 / total];
        let orient_allowed = [orient_frac[0] >= ORIENT_THRESH, orient_frac[1] >= ORIENT_THRESH];
        for i in 0..2 {
            log::info!("    {}: {:8} ({:.3}%),  {}",
                if i == 0 { "FR/RF" } else { "FF/RR" },
                orient_counts[i], orient_frac[i],
                if orient_allowed[i] { "allowed" } else { "forbidden" });
        }

        VecExt::sort(&mut insert_sizes);
        let max_size = INS_QUANTILE_MULT * F64Ext::quantile_sorted(&insert_sizes, INS_QUANTILE);
        // Find index after the limiting value.
        let m = bisect::right(&insert_sizes, &max_size);
        let lim_insert_sizes = &insert_sizes[..m];

        let mean = F64Ext::mean(lim_insert_sizes);
        // Increase variance, if less-equal than mean.
        let var = F64Ext::variance(lim_insert_sizes, Some(mean));
        let distr = RegularizedEstimator::default().estimate(mean, var);
        log::info!("    Insert size mean = {:.1},  st.dev. = {:.1}", distr.mean(), distr.variance().sqrt());

        let quantile = 0.5 * (1.0 - params.ins_conf_level);
        let min_size = (distr.quantile(quantile) - 1e-8).floor().max(0.0) as u32;
        let max_size = (distr.quantile(1.0 - quantile) + 1e-8).ceil() as u32;
        log::info!("    Allowed insert size: [{}, {}]  ({}%-confidence interval)",
            min_size, max_size, 100.0 * params.ins_conf_level);

        Ok(Self {
            orient_allowed,
            distr: Some(distr.cached(max_size as usize + 1)),
            conf_interval: (min_size, max_size),
        })
    }

    /// Ln-probability of the insert size. `same_orient` is true if FF/RR, false if FR/RF.
    pub fn ln_prob(&self, sz: u32, same_orient: bool) -> f64 {
        if self.conf_interval.0 <= sz && sz <= self.conf_interval.1 && self.orient_allowed[usize::from(same_orient)] {
            self.distr.as_ref().unwrap().ln_pmf(sz)
        } else {
            f64::NEG_INFINITY
        }
    }

    /// Maximum insert size. Over this size, all pairs are considered unpaired.
    pub fn confidence_interval(&self) -> (u32, u32) {
        self.conf_interval
    }

    /// Returns true if the reads are paired-end, false if single-end.
    pub fn is_paired_end(&self) -> bool {
        self.distr.is_some()
    }
}

impl JsonSer for InsertDistr {
    fn save(&self) -> json::JsonValue {
        if let Some(distr) = &self.distr {
            json::object!{
                fr_allowed: self.orient_allowed[0],
                ff_allowed: self.orient_allowed[1],
                n: distr.inner().n(),
                p: distr.inner().p(),
                min_size: self.conf_interval.0,
                max_size: self.conf_interval.1,
            }
        } else {
            json::object!{
                max_size: 0,
            }
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> max_size (as_u32));
        if max_size == 0 {
            return Ok(Self::undefined());
        }
        json_get!(obj -> n (as_f64), p (as_f64), fr_allowed (as_bool), ff_allowed (as_bool), min_size (as_u32));

        let distr = NBinom::new(n, p).cached(max_size as usize + 1);
        Ok(Self {
            orient_allowed: [fr_allowed, ff_allowed],
            distr: Some(distr),
            conf_interval: (min_size, max_size),
        })
    }
}
