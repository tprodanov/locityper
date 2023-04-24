pub mod err_prof;
pub mod insertsz;
pub mod depth;
pub mod ser;

use htslib::bam::Record;
use crate::err::{Error, validate_param};
use {
    depth::{ReadDepth, ReadDepthParams},
    insertsz::InsertDistr,
    err_prof::ErrorProfile,
};
pub use ser::JsonSer;

/// Parameters for background distributions estimation.
#[derive(Debug, Clone)]
pub struct Params {
    /// Background read depth parameters.
    pub depth: ReadDepthParams,

    /// Minimal mapping quality.
    pub min_mapq: u8,
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
    pub err_rate_mult: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            depth: Default::default(),

            min_mapq: 20,
            max_clipping: 0.02,

            ins_quantile: 0.99,
            ins_quantile_mult: 3.0,

            err_quantile: 0.01,
            err_quantile_mult: 2.0,
            err_rate_mult: 1.2,
        }
    }
}

impl Params {
    /// Validate all parameter values.
    pub fn validate(&self) -> Result<(), Error> {
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
        validate_param!(0.2 <= self.err_rate_mult,
            "Error rate multiplier ({:.5}) should not be too low", self.err_rate_mult);

        self.depth.validate()
    }
}

/// Returns true for reads alignments that are
/// - Primary and pass all checks,
/// - Are unpaired, OR
/// - Are paired and within an appropriate insert size.
fn read_unpaired_or_proper_pair(record: &Record, max_insert_size: i64) -> bool {
    (record.flags() & 3844) == 0
        && (!record.is_paired()
            || (!record.is_mate_unmapped()
                && record.tid() == record.mtid()
                && record.insert_size().abs() <= max_insert_size))
}

/// Various background distributions, including
/// - read depth distribution,
/// - insert size distribution,
/// - error profile.
pub struct BgDistr {
    insert_distr: InsertDistr,
    err_prof: ErrorProfile,
    depth: ReadDepth,
}

impl BgDistr {
    pub fn new(insert_distr: InsertDistr, err_prof: ErrorProfile, depth: ReadDepth) -> Self {
        Self { insert_distr, err_prof, depth }
    }

    pub fn depth(&self) -> &ReadDepth {
        &self.depth
    }

    pub fn insert_distr(&self) -> &InsertDistr {
        &self.insert_distr
    }

    pub fn error_profile(&self) -> &ErrorProfile {
        &self.err_prof
    }
}

impl JsonSer for BgDistr {
    fn save(&self) -> json::JsonValue {
        json::object!{
            insert_distr: self.insert_distr.save(),
            error_profile: self.err_prof.save(),
            bg_depth: self.depth.save(),
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        if obj.has_key("bg_depth") && obj.has_key("insert_distr") && obj.has_key("error_profile") {
            Ok(Self {
                insert_distr: InsertDistr::load(&obj["insert_distr"])?,
                err_prof: ErrorProfile::load(&obj["error_profile"])?,
                depth: ReadDepth::load(&obj["bg_depth"])?,
            })
        } else {
            Err(Error::JsonLoad(format!(
                "BgDistr: Failed to parse '{}': missing 'bg_depth', 'insert_distr' or 'error_profile' keys!", obj)))
        }
    }
}
