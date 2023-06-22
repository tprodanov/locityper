pub mod err_prof;
pub mod insertsz;
pub mod depth;
pub mod ser;

use std::{
    fmt,
    str::FromStr,
};
use crate::err::{Error, validate_param};
use self::ser::json_get;
pub use {
    depth::{ReadDepth, ReadDepthParams},
    insertsz::InsertDistr,
    err_prof::ErrorProfile,
    ser::JsonSer,
};

/// Parameters for background distributions estimation.
#[derive(Debug, Clone)]
pub struct Params {
    /// Background read depth parameters.
    pub depth: ReadDepthParams,

    /// Confidence level for the insert size.
    pub ins_conf_level: f64,
    /// Confidence level for the edit distance.
    pub err_conf_level: f64,
    /// Error probability multiplier: multiply read error probabilities (mismatches, insertions, deletions, clipping),
    /// by this value. This will soften overly harsh read alignment penalties.
    pub err_rate_mult: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            depth: Default::default(),
            ins_conf_level: 0.999,
            err_conf_level: 0.99,
            err_rate_mult: 1.0,
        }
    }
}

impl Params {
    /// Validate all parameter values.
    pub fn validate(&self) -> Result<(), Error> {
        validate_param!(0.2 <= self.err_rate_mult,
            "Error rate multiplier ({:.5}) should not be too low", self.err_rate_mult);
        for &conf_lvl in &[self.ins_conf_level, self.err_conf_level] {
            validate_param!(0.5 < conf_lvl && conf_lvl < 1.0,
                "Confidence level ({}) must be in (0.5, 1).", conf_lvl);
        }
        self.depth.validate()
    }
}

/// Various background distributions.
pub struct BgDistr {
    seq_info: SequencingInfo,
    /// Insert size distribution,
    insert_distr: InsertDistr,
    /// Error profile.
    err_prof: ErrorProfile,
    /// Read depth distribution,
    depth: ReadDepth,
}

impl BgDistr {
    pub fn new(
        seq_info: SequencingInfo,
        insert_distr: InsertDistr,
        err_prof: ErrorProfile,
        depth: ReadDepth
    ) -> Self
    {
        Self { seq_info, insert_distr, err_prof, depth }
    }

    /// Access sequencing information (read length and technology).
    pub fn seq_info(&self) -> &SequencingInfo {
        &self.seq_info
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
            seq_info: self.seq_info.save(),
            insert_distr: self.insert_distr.save(),
            error_profile: self.err_prof.save(),
            bg_depth: self.depth.save(),
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        if obj.has_key("seq_info") && obj.has_key("bg_depth")
                && obj.has_key("insert_distr") && obj.has_key("error_profile") {
            Ok(Self {
                seq_info: SequencingInfo::load(&obj["seq_info"])?,
                insert_distr: InsertDistr::load(&obj["insert_distr"])?,
                err_prof: ErrorProfile::load(&obj["error_profile"])?,
                depth: ReadDepth::load(&obj["bg_depth"])?,
            })
        } else {
            Err(Error::JsonLoad(format!(
                "BgDistr: Failed to parse '{}': missing 'seq_info', 'bg_depth', \
                    'insert_distr' or 'error_profile' keys!", obj)))
        }
    }
}

/// Sequencing technology.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Technology {
    Illumina,
    HiFi,
    PacBio,
    Nanopore,
}

impl Technology {
    pub fn to_str(&self) -> &'static str {
        match self {
            Self::Illumina => "illumina",
            Self::HiFi => "hifi",
            Self::PacBio => "pacbio",
            Self::Nanopore => "nanopore",
        }
    }

    pub fn minimap_preset(&self) -> &'static str {
        match self {
            Self::Illumina => {
                log::warn!("Using Minimap2 for short reads!");
                "sr"
            },
            Self::HiFi => "map-hifi",
            Self::PacBio => "map-pb",
            Self::Nanopore => "map-ont",
        }
    }
}

impl FromStr for Technology {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match &s.to_lowercase() as &str {
            "illumina" | "sr" => Ok(Self::Illumina),
            "hifi" => Ok(Self::HiFi),
            "pacbio" | "pb" => Ok(Self::PacBio),
            "nanopore" | "ont" => Ok(Self::Nanopore),
            _ => Err(format!("Unknown technology {:?}", s)),
        }
    }
}

impl fmt::Display for Technology {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(self.to_str())
    }
}

impl fmt::Debug for Technology {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(match self {
            Self::Illumina => "'Illumina'",
            Self::HiFi => "'PacBio HiFi'",
            Self::PacBio => "'Pacbio CLR'",
            Self::Nanopore => "'Oxford Nanopore'",
        })
    }
}

#[derive(Debug, Clone)]
pub struct SequencingInfo {
    /// Mean read length.
    read_len: f64,
    /// Sequencing technology.
    technology: Technology,
}

impl SequencingInfo {
    pub fn new(read_len: f64, technology: Technology) -> Self {
        if read_len > 400.0 && technology == Technology::Illumina {
            log::error!("Unusual mean read length ({:.0}) for the {:?} sequencing technology",
                read_len, technology);
        }
        Self { read_len, technology }
    }

    pub fn mean_read_len(&self) -> f64 {
        self.read_len
    }

    pub fn technology(&self) -> Technology {
        self.technology
    }
}

impl JsonSer for SequencingInfo {
    fn save(&self) -> json::JsonValue {
        json::object!{
            read_len: self.read_len,
            technology: self.technology.to_str(),
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> read_len (as_f64), technology (as_str));
        let technology = Technology::from_str(technology).map_err(|e| Error::JsonLoad(e))?;
        Ok(Self { read_len, technology })
    }
}
