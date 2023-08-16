pub mod err_prof;
pub mod insertsz;
pub mod depth;
pub mod windows;
pub mod ser;

use std::{
    fmt,
    str::FromStr,
};
use crate::{
    math,
    err::{Error, validate_param},
};
use self::ser::json_get;
pub use {
    depth::{ReadDepth, ReadDepthParams},
    insertsz::InsertDistr,
    err_prof::ErrorProfile,
    ser::JsonSer,
    windows::Windows,
};

/// Parameters for background distributions estimation.
#[derive(Debug, Clone)]
pub struct Params {
    /// Background read depth parameters.
    pub depth: ReadDepthParams,

    /// p-value threshold for the insert size.
    pub insert_pval: f64,
    /// p-value threhsold for the edit distance.
    pub edit_pval: f64,
    /// Error probability multiplier: multiply read error probabilities (mismatches, insertions, deletions, clipping),
    /// by this value. This will soften overly harsh read alignment penalties.
    pub err_rate_mult: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            depth: Default::default(),
            insert_pval: 0.001,
            edit_pval: err_prof::DEF_EDIT_PVAL.0,
            err_rate_mult: 1.0,
        }
    }
}

impl Params {
    /// Validate all parameter values.
    pub fn validate(&self) -> Result<(), Error> {
        validate_param!(0.2 <= self.err_rate_mult,
            "Error rate multiplier ({:.5}) should not be too low", self.err_rate_mult);
        for &pval in &[self.insert_pval, self.edit_pval] {
            validate_param!(0.0 < pval && pval < 0.5,
                "p-value threshold ({}) must be in (0, 0.5)", pval);
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

    pub fn set_edit_pvals(&mut self, edit_pvals: (f64, f64)) {
        self.err_prof.set_edit_pvals(edit_pvals);
        let read_len = math::round_signif(self.seq_info.mean_read_len(), 2).round() as u32;
        let (good_dist, passable_dist) = self.err_prof.allowed_edit_dist(read_len);
        log::info!("Edit distances for read length {}: {} (good) and {} (passable)",
            read_len, good_dist, passable_dist);
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
    pub fn to_str(self) -> &'static str {
        match self {
            Self::Illumina => "illumina",
            Self::HiFi => "hifi",
            Self::PacBio => "pacbio",
            Self::Nanopore => "nanopore",
        }
    }

    pub fn long_name(self) -> &'static str {
        match self {
            Self::Illumina => "Illumina",
            Self::HiFi => "PacBio HiFi",
            Self::PacBio => "Pacbio CLR",
            Self::Nanopore => "Oxford Nanopore",
        }
    }

    pub fn minimap_preset(self) -> &'static str {
        match self {
            Self::Illumina => unreachable!("Minimap2 is not used for short reads"),
            Self::HiFi => "map-hifi",
            Self::PacBio => "map-pb",
            Self::Nanopore => "map-ont",
        }
    }

    /// Returns true if the technology exhibits different depth values at different GC-contents.
    pub fn has_gc_bias(self) -> bool {
        self == Self::Illumina
    }

    /// Returns ranges of expected sequencing lengths for various sequencing technologies.
    pub fn expect_mean_length(self) -> (f64, f64) {
        match self {
            Self::Illumina => (100.0, 400.0),
            Self::HiFi => (5_000.0, 30_000.0),
            Self::PacBio => (5_000.0, 150_000.0),
            Self::Nanopore => (5_000.0, 500_000.0),
        }
    }

    pub fn paired_end_allowed(self) -> bool {
        self == Self::Illumina
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

#[derive(Clone)]
pub struct SequencingInfo {
    /// Mean read length.
    read_len: f64,
    /// Sequencing technology.
    technology: Technology,
}

impl SequencingInfo {
    /// Creates sequencing info. Does not return Error if `warn`, even if the read length is suspicious.
    pub fn new(read_len: f64, technology: Technology, warn: bool) -> Result<Self, Error> {
        let (min_mean_len, max_mean_len) = technology.expect_mean_length();
        if read_len < min_mean_len || read_len > max_mean_len {
            if warn {
                log::warn!("Unusual mean read length ({:.0}) for the {} sequencing technology",
                    read_len, technology.long_name());
            } else {
                return Err(Error::InvalidInput(format!(
                    "Unusual mean read length ({:.0}) for the {} sequencing technology",
                    read_len, technology.long_name())));
            }
        }
        Ok(Self { read_len, technology })
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
