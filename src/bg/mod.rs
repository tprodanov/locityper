pub mod err_prof;
pub mod insertsz;
pub mod depth;
pub mod ser;

use crate::err::{Error, validate_param};
pub use {
    depth::{ReadDepth, ReadDepthParams},
    insertsz::InsertDistr,
    err_prof::ErrorProfile,
    ser::JsonSer,
};
use ser::json_get;

/// Parameters for background distributions estimation.
#[derive(Debug, Clone)]
pub struct Params {
    /// Background read depth parameters.
    pub depth: ReadDepthParams,

    // Calculate max insert size from input reads as `<ins_quantile_mult> * <ins_quantile>-th insert size quantile`.
    // This is needed to remove read mates that were mapped to the same chromosome but very far apart.
    pub ins_quantile: f64,
    pub ins_quantile_mult: f64,

    /// Error probability multiplier: multiply read error probabilities (mismatches, insertions, deletions, clipping),
    /// by this value. This will soften overly harsh read alignment penalties.
    pub err_rate_mult: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            depth: Default::default(),
            ins_quantile: 0.99,
            ins_quantile_mult: 3.0,
            err_rate_mult: 1.0,
        }
    }
}

impl Params {
    /// Validate all parameter values.
    pub fn validate(&self) -> Result<(), Error> {
        validate_param!(0.5 <= self.ins_quantile && self.ins_quantile <= 1.0,
            "Insert size quantile ({:.5}) must be within [0.5, 1]", self.ins_quantile);
        validate_param!(self.ins_quantile_mult >= 1.0,
            "Insert size quantile multiplier ({:.5}) must be at least 1", self.ins_quantile_mult);
        validate_param!(0.2 <= self.err_rate_mult,
            "Error rate multiplier ({:.5}) should not be too low", self.err_rate_mult);
        self.depth.validate()
    }
}

/// Various background distributions.
pub struct BgDistr {
    /// Mean read length.
    read_len: f64,
    /// Insert size distribution,
    insert_distr: InsertDistr,
    /// Error profile.
    err_prof: ErrorProfile,
    /// Read depth distribution,
    depth: ReadDepth,
}

impl BgDistr {
    pub fn new(read_len: f64, insert_distr: InsertDistr, err_prof: ErrorProfile, depth: ReadDepth) -> Self {
        Self { read_len, insert_distr, err_prof, depth }
    }

    pub fn mean_read_len(&self) -> f64 {
        self.read_len
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
            read_len: self.read_len,
            insert_distr: self.insert_distr.save(),
            error_profile: self.err_prof.save(),
            bg_depth: self.depth.save(),
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> read_len (as_f64));
        if obj.has_key("bg_depth") && obj.has_key("insert_distr") && obj.has_key("error_profile") {
            Ok(Self {
                read_len,
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
