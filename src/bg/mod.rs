pub mod err_prof;
pub mod insertsz;
pub mod depth;
pub mod ser;

use std::io;
use htslib::bam::Record;
use crate::seq::{
    Interval,
    kmers::KmerCounts,
    cigar::CigarItem,
};
use {
    depth::{ReadDepth, ReadDepthParams},
    insertsz::{InsertNegBinom, InsertDistr},
    err_prof::ErrorProfile,
    ser::{LoadError},
};
pub use ser::JsonSer;

/// Ignore reads with MAPQ < 20.
pub const MIN_MAPQ: u8 = 20;

/// Parameters for background distributions estimation.
#[derive(Debug, Clone)]
pub struct Params {
    /// Background read depth parameters.
    pub depth: ReadDepthParams,

    /// Insert size calculation parameters.
    pub insert_size: insertsz::InsertSizeParams,

    /// Error probability multiplier: multiply read error probabilities (mismatches, insertions, deletions, clipping),
    /// by this value. This will soften overly harsh read alignment penalties.
    pub err_prob_mult: f64,
    /// When calculating read error profiles, ignore reads with `clipping > max_clipping * read_len`.
    pub max_clipping: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            depth: Default::default(),
            insert_size: Default::default(),
            err_prob_mult: 1.2,
            max_clipping: 0.02,
        }
    }
}

impl Params {
    /// Checks all parameter values.
    pub fn check(&self) {
        self.depth.check();
        assert!(0.2 <= self.err_prob_mult, "Error prob. multiplier ({:.5}) should not be too low", self.err_prob_mult);
        assert!(0.0 <= self.max_clipping && self.max_clipping <= 1.0, "Max clipping ({:.5}) must be within [0, 1]",
            self.max_clipping);
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
    insert_sz: InsertNegBinom,
    err_prof: ErrorProfile,
}

impl BgDistr {
    /// Estimates read depth, insert size and error profile given a slice of BAM records.
    pub fn estimate(
        records: &[Record],
        interval: &Interval,
        ref_seq: &[u8],
        kmer_counts: &KmerCounts,
        params: &Params,
    ) -> io::Result<Self> {
        log::info!("Estimating background parameters");
        log::debug!("    Use {} reads on {} bp interval", records.len(), interval.len());

        let full_alns: Vec<&Record> = records.iter()
            .filter(|record| is_full_aln(record)).collect();
        let rec_grouping = insertsz::ReadMateGrouping::from_unsorted_bam(
            full_alns.iter().map(|record| *record), Some(full_alns.len()));
        let insert_sz = InsertNegBinom::estimate(&rec_grouping, &params.insert_size);
        log::debug!("insert size: {:?}", insert_sz);
        let err_prof = ErrorProfile::estimate(records.iter(), params.max_clipping, params.err_prob_mult);
        let depth = ReadDepth::estimate(interval, &ref_seq, records.iter(), &params.depth, insert_sz.max_size());
        Ok(Self { depth, insert_sz, err_prof })
    }

    pub fn depth(&self) -> &ReadDepth {
        &self.depth
    }

    pub fn insert_size(&self) -> &InsertNegBinom {
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
                insert_sz: InsertNegBinom::load(&obj["insert_size"])?,
                err_prof: ErrorProfile::load(&obj["error_profile"])?,
            })
        } else {
            Err(LoadError(format!(
                "BgDistr: Failed to parse '{}': missing 'bg_depth', 'insert_size' or 'error_profile' keys!", obj)))
        }
    }
}
