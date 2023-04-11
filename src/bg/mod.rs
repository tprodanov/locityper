pub mod err_prof;
pub mod insertsz;
pub mod depth;
pub mod ser;

use std::io;
use htslib::bam::Record;
use crate::{
    seq::{
        Interval,
        kmers::KmerCounts,
        cigar::CigarItem,
    },
    err::{Error, validate_param},
};
use {
    depth::{ReadDepth, ReadDepthParams},
    insertsz::InsertDistr,
    err_prof::ErrorProfile,
    ser::{LoadError},
};
pub use ser::JsonSer;

/// Parameters for background distributions estimation.
#[derive(Debug, Clone)]
pub struct Params {
    /// Background read depth parameters.
    pub depth: ReadDepthParams,

    /// When calculating insert size distributions and read error profiles,
    /// ignore reads with `clipping > max_clipping * read_len`.
    pub max_clipping: f64,

    // Calculate max insert size from input reads as `<ins_quantile_mult> * <ins_quantile>-th insert size quantile`.
    // This is needed to remove read mates that were mapped to the same chromosome but very far apart.
    pub ins_quantile: f64,
    pub ins_quantile_mult: f64,

    /// Similar as for insert size quantiles, calculate minimal possible alignment probability
    /// as `<err_quantile_mult> & <err_quantile>` (in log-space).
    pub err_quantile: f64,
    pub err_quantile_mult: f64,
    /// Error probability multiplier: multiply read error probabilities (mismatches, insertions, deletions, clipping),
    /// by this value. This will soften overly harsh read alignment penalties.
    pub err_prob_mult: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            depth: Default::default(),

            max_clipping: 0.02,

            ins_quantile: 0.99,
            ins_quantile_mult: 3.0,

            err_quantile: 0.01,
            err_quantile_mult: 3.0,
            err_prob_mult: 1.2,
        }
    }
}

impl Params {
    /// Validate all parameter values.
    pub fn validate(&self) -> Result<(), Error> {
        // TODO: self.depth.validate();

        validate_param!(0.0 <= self.max_clipping && self.max_clipping <= 1.0,
            "Max clipping ({:.5}) must be within [0, 1]", self.max_clipping);

        validate_param!(0.5 <= self.ins_quantile && self.ins_quantile <= 1.0,
            "Insert size quantile ({:.5}) must be within [0.5, 1]", self.ins_quantile);
        validate_param!(self.ins_quantile_mult >= 1.0,
            "Insert size quantile multiplier ({:.5}) must be at least 1", self.ins_quantile_mult);

        validate_param!(0.0 <= self.err_quantile && self.err_quantile <= 0.5,
            "Alignment likelihood quantile ({:.5}) must be within [0, 0.5]", self.err_quantile);
        validate_param!(self.err_quantile_mult >= 1.0,
            "Alignment likelihood quantile multiplier ({:.5}) must be at least 1", self.err_quantile_mult);
        validate_param!(0.2 <= self.err_prob_mult,
            "Error prob. multiplier ({:.5}) should not be too low", self.err_prob_mult);
        Ok(())
    }
}

/// Returns true if the alignment has no soft or hard clipping.
fn is_full_aln(record: &Record) -> bool {
    let raw_cigar = record.raw_cigar();
    let n = raw_cigar.len();
    CigarItem::from_u32(raw_cigar[0]).operation().consumes_ref() &&
        CigarItem::from_u32(raw_cigar[n - 1]).operation().consumes_ref()
}

/// Various background distributions, including
/// - read depth distribution,
/// - insert size distribution,
/// - error profile.
pub struct BgDistr {
    depth: ReadDepth,
    insert_sz: InsertDistr,
    err_prof: ErrorProfile,
}

impl BgDistr {
    // /// Estimates read depth, insert size and error profile given a slice of BAM records.
    // pub fn estimate(
    //     records: &[Record],
    //     interval: &Interval,
    //     ref_seq: &[u8],
    //     kmer_counts: &KmerCounts,
    //     params: &Params,
    // ) -> io::Result<Self> {
    //     unimplemented!()
    //     // log::info!("Estimating background parameters");
    //     // log::debug!("    Use {} reads on {} bp interval", records.len(), interval.len());

    //     // let full_alns: Vec<&Record> = records.iter()
    //     //     .filter(|record| is_full_aln(record)).collect();
    //     // let rec_grouping = insertsz::ReadMateGrouping::from_unsorted_bam(
    //     //     full_alns.iter().map(|record| *record), Some(full_alns.len()));
    //     // let insert_sz = InsertDistr::estimate(&rec_grouping, &params.insert_size);
    //     // log::debug!("insert size: {:?}", insert_sz);
    //     // let err_prof = ErrorProfile::estimate(records.iter(), params.max_clipping, params.err_prob_mult);
    //     // let depth = ReadDepth::estimate(interval, &ref_seq, records.iter(), &params.depth, insert_sz.max_size());
    //     // Ok(Self { depth, insert_sz, err_prof })
    // }

    pub fn depth(&self) -> &ReadDepth {
        &self.depth
    }

    pub fn insert_size(&self) -> &InsertDistr {
        &self.insert_sz
    }

    pub fn error_profile(&self) -> &ErrorProfile {
        &self.err_prof
    }
}

impl JsonSer for BgDistr {
    fn save(&self) -> json::JsonValue {
        json::object!{
            bg_depth: self.depth.save(),
            insert_size: self.insert_sz.save(),
            error_profile: self.err_prof.save(),
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        if obj.has_key("bg_depth") && obj.has_key("insert_size") && obj.has_key("error_profile") {
            Ok(Self {
                depth: ReadDepth::load(&obj["bg_depth"])?,
                insert_sz: InsertDistr::load(&obj["insert_size"])?,
                err_prof: ErrorProfile::load(&obj["error_profile"])?,
            })
        } else {
            Err(LoadError(format!(
                "BgDistr: Failed to parse '{}': missing 'bg_depth', 'insert_size' or 'error_profile' keys!", obj)))
        }
    }
}
