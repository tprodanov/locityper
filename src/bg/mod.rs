pub mod err_prof;
pub mod insertsz;
pub mod depth;
pub mod windows;
pub mod ser;

pub use {
    depth::{ReadDepth, ReadDepthParams},
    insertsz::InsertDistr,
    err_prof::ErrorProfile,
    ser::JsonSer,
    windows::Windows,
};

use std::{
    fmt::{self, Write},
    str::FromStr,
    path::Path,
};
use crate::{
    math, ext,
    err::{Error, validate_param, add_path},
};
use self::ser::json_get;

/// Parameters for background distributions estimation.
#[derive(Debug, Clone)]
pub struct Params {
    /// Background read depth parameters.
    pub depth: ReadDepthParams,

    /// p-value threshold for the insert size.
    pub insert_pval: f64,
    /// p-value threhsold for the edit distance.
    pub edit_pval: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            depth: Default::default(),
            insert_pval: 0.001,
            edit_pval: err_prof::DEF_EDIT_PVAL.0,
        }
    }
}

impl Params {
    /// Validate all parameter values.
    pub fn validate(&self) -> Result<(), Error> {
        for &pval in &[self.insert_pval, self.edit_pval] {
            validate_param!(0.0 < pval && pval < 0.5,
                "p-value threshold ({}) must be in (0, 0.5)", pval);
        }
        self.depth.validate()
    }
}

/// Various background distributions.
#[derive(Clone)]
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

    pub fn load_from(path: &Path) -> Result<Self, Error> {
        let mut stream = ext::sys::open(&path)?;
        let mut s = String::new();
        stream.read_to_string(&mut s).map_err(add_path!(path))?;
        BgDistr::load(&json::parse(&s)?)
    }

    /// Access sequencing information (read length and technology).
    pub fn seq_info(&self) -> &SequencingInfo {
        &self.seq_info
    }

    pub fn set_seq_info(&mut self, seq_info: SequencingInfo) {
        self.seq_info = seq_info;
    }

    pub fn depth(&self) -> &ReadDepth {
        &self.depth
    }

    pub fn depth_mut(&mut self) -> &mut ReadDepth {
        &mut self.depth
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

    /// Describe values for each technology in format
    /// "illumina: X, hifi: Y, ..."
    pub fn describe_values<D, F>(f: F) -> String
    where D: fmt::Display,
          F: Fn(Self) -> D,
    {
        // NOTE: In nightly version, can use `std::mem::variant_count`.
        const TECHS: [Technology; 4] = [Technology::Illumina, Technology::HiFi,
            Technology::PacBio, Technology::Nanopore];
        let mut s = String::new();
        for (i, &tech) in TECHS.iter().enumerate() {
            if i > 0 {
                s.push_str(", ");
            }
            write!(s, "{}: {}", tech, f(tech)).unwrap();
        }
        s
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

    /// Returns minimizer matching fraction for different technologies.
    pub fn default_matches_frac(self) -> f64 {
        match self {
            Self::Illumina => 0.3,
            Self::HiFi => 0.2,
            Self::PacBio | Self::Nanopore => 0.1,
        }
    }

    /// Returns true if mean read lengths in two datasets are similar enough.
    pub fn is_read_len_similar(self, len1: f64, len2: f64) -> bool {
        if self == Self::Illumina {
            (len1 - len2).abs() < 3.0
        } else {
            // Read length does not differ by more than 20%.
            (len1 - len2).abs() / len1.min(len2) < 0.2
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

#[derive(Clone)]
pub struct SequencingInfo {
    /// Mean read length.
    read_len: f64,
    /// Sequencing technology.
    technology: Technology,
    /// Total number of reads/read pairs in the input file.
    total_reads: Option<u64>,
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
        Ok(Self {
            read_len, technology,
            total_reads: None,
        })
    }

    pub fn mean_read_len(&self) -> f64 {
        self.read_len
    }

    pub fn technology(&self) -> Technology {
        self.technology
    }

    pub fn set_total_reads(&mut self, count: u64) {
        self.total_reads = Some(count);
    }

    pub fn total_reads(&self) -> Option<u64> {
        self.total_reads
    }
}

impl JsonSer for SequencingInfo {
    fn save(&self) -> json::JsonValue {
        json::object!{
            read_len: self.read_len,
            technology: self.technology.to_str(),
            total_reads: self.total_reads,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> read_len (as_f64), technology (as_str), total_reads? (as_u64));
        let technology = Technology::from_str(technology).map_err(|e| Error::JsonLoad(e))?;
        Ok(Self { read_len, technology, total_reads })
    }
}
