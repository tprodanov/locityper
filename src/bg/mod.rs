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
    fmt,
    str::FromStr,
    path::Path,
};
use crate::{
    ext,
    err::{Error, validate_param, add_path, error},
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
            edit_pval: 0.01,
        }
    }
}

impl Params {
    /// Validate all parameter values.
    pub fn validate(&self) -> crate::Result<()> {
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
    depth: Option<ReadDepth>,
}

impl BgDistr {
    pub fn new(
        seq_info: SequencingInfo,
        insert_distr: InsertDistr,
        err_prof: ErrorProfile,
        depth: Option<ReadDepth>,
    ) -> Self
    {
        Self { seq_info, insert_distr, err_prof, depth }
    }

    pub fn load_from(path: &Path, success_file: Option<&Path>) -> crate::Result<Self> {
        if let Some(filename) = success_file {
            if !filename.exists() {
                log::warn!("File {} does not exist, possibly preprocessing was not completed",
                   ext::fmt::path(filename));
            }
        }
        let mut stream = ext::sys::open(&path)?;
        let mut s = String::new();
        stream.read_to_string(&mut s).map_err(add_path!(path))?;
        BgDistr::load(&json::parse(&s)?)
    }

    /// Describe distributions to log.
    pub fn describe(&self) {
        self.err_prof.describe();
        self.insert_distr.describe();
        match &self.depth {
            Some(depth) => depth.describe(self.insert_distr.is_paired_end()),
            None => log::warn!("Background read depth was not estimated"),
        };
    }

    /// Access sequencing information (read length and technology).
    #[inline(always)]
    pub fn seq_info(&self) -> &SequencingInfo {
        &self.seq_info
    }

    #[inline(always)]
    pub fn set_seq_info(&mut self, seq_info: SequencingInfo) {
        self.seq_info = seq_info;
    }

    #[inline(always)]
    pub fn has_read_depth(&self) -> bool {
        self.depth.is_some()
    }

    #[inline(always)]
    pub fn depth(&self) -> &ReadDepth {
        self.depth.as_ref().expect("Read depth must be defined")
    }

    #[inline(always)]
    pub fn opt_depth(&self) -> Option<&ReadDepth> {
        self.depth.as_ref()
    }

    #[inline(always)]
    pub fn depth_mut(&mut self) -> &mut ReadDepth {
        self.depth.as_mut().expect("Read depth must be defined")
    }

    #[inline(always)]
    pub fn insert_distr(&self) -> &InsertDistr {
        &self.insert_distr
    }

    #[inline(always)]
    pub fn error_profile(&self) -> &ErrorProfile {
        &self.err_prof
    }
}

impl JsonSer for BgDistr {
    fn save(&self) -> json::JsonValue {
        let mut obj = json::object!{
            seq_info: self.seq_info.save(),
            insert_distr: self.insert_distr.save(),
            error_profile: self.err_prof.save(),
        };
        if let Some(depth) = &self.depth {
            obj["bg_depth"] = depth.save();
        }
        obj
    }

    fn load(obj: &json::JsonValue) -> crate::Result<Self> {
        let depth = if obj.has_key("bg_depth") {
            Some(ReadDepth::load(&obj["bg_depth"])?)
        } else { None };

        if obj.has_key("seq_info") && obj.has_key("insert_distr") && obj.has_key("error_profile") {
            Ok(Self {
                seq_info: SequencingInfo::load(&obj["seq_info"])?,
                insert_distr: InsertDistr::load(&obj["insert_distr"])?,
                err_prof: ErrorProfile::load(&obj["error_profile"])?,
                depth,
            })
        } else {
            Err(error!(JsonLoad, "BgDistr: Failed to parse '{}': missing 'seq_info', \
                'insert_distr' or 'error_profile' keys!", obj))
        }
    }
}

/// Sequencing technology.
#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Technology {
    Illumina,
    HiFi,
    PacBio,
    Nanopore,
}

// NOTE: In nightly version, can use `std::mem::variant_count`.
pub const TECHNOLOGIES: [Technology; 4] = [
    Technology::Illumina,
    Technology::HiFi,
    Technology::PacBio,
    Technology::Nanopore
];

impl Technology {
    pub fn to_str(self) -> &'static str {
        match self {
            Self::Illumina => "illumina",
            Self::HiFi => "hifi",
            Self::PacBio => "pacbio",
            Self::Nanopore => "ont",
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

    /// Returns minimizer matching fraction for different technologies.
    pub fn default_match_frac(self, is_paired_end: bool) -> f64 {
        match (self, is_paired_end) {
            (Self::Illumina, false) => 0.7,
            (Self::Illumina, true)  => 0.5,
            (_, false) => 0.5,
            (_, true) => unreachable!("Paired-end long reads are not supported"),
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
    /// Full input file size.
    file_size: Option<u64>,
}

impl SequencingInfo {
    /// Creates sequencing info. Does not return Error if `warn`, even if the read length is suspicious.
    pub fn new(read_len: f64, technology: Technology, warn: bool) -> crate::Result<Self> {
        let (min_mean_len, max_mean_len) = technology.expect_mean_length();
        if read_len < min_mean_len || read_len > max_mean_len {
            if warn {
                log::warn!("Unusual mean read length ({:.0}) for the {} sequencing technology",
                    read_len, technology.long_name());
            } else {
                return Err(error!(InvalidInput,
                    "Unusual mean read length ({:.0}) for the {} sequencing technology.\n\
                    Please specify technology using `--tech` argument (`--tech illumina` overrides this error)",
                    read_len, technology.long_name()));
            }
        }
        Ok(Self {
            read_len, technology,
            total_reads: None,
            file_size: None,
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

    pub fn set_file_size(&mut self, size: u64) {
        self.file_size = Some(size);
    }

    pub fn file_size(&self) -> Option<u64> {
        self.file_size
    }
}

impl JsonSer for SequencingInfo {
    fn save(&self) -> json::JsonValue {
        json::object!{
            read_len: self.read_len,
            technology: self.technology.to_str(),
            total_reads: self.total_reads,
            file_size: self.file_size,
        }
    }

    fn load(obj: &json::JsonValue) -> crate::Result<Self> {
        json_get!(obj => read_len (as_f64), technology (as_str), total_reads? (as_u64), file_size? (as_u64));
        let technology = Technology::from_str(technology).map_err(|e| Error::JsonLoad(e))?;
        Ok(Self { read_len, technology, total_reads, file_size })
    }
}
