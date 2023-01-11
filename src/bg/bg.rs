use std::io::{self, Read, Seek};
use htslib::bam::Record;
use bio::io::fasta::IndexedReader;
use crate::{
    bg::{
        depth::{ReadDepth, ReadDepthParams},
        insertsz::{InsertNegBinom, InsertDistr, get_insert_size},
        err_prof::ErrorProfile,
        ser::{JsonSer, LoadError},
    },
    seq::interv::Interval,
};

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
    pub fn estimate<R: Read + Seek>(
            records: &[Record],
            interval: &Interval,
            fasta: &mut IndexedReader<R>,
            params: &ReadDepthParams)
            -> io::Result<Self>
    {
        // TODO: Consider using lighter dependency (such as faimm).
        fasta.fetch(interval.contig_name(), interval.start() as u64, interval.end() as u64)?;
        let mut ref_seq = Vec::new();
        fasta.read(&mut ref_seq)?;

        let insert_sz = InsertNegBinom::estimate(records.iter().filter_map(get_insert_size));
        let err_prof = ErrorProfile::estimate(records.iter());
        let depth = ReadDepth::estimate(interval, &ref_seq, records.iter(), params, insert_sz.max_size());
        Ok(Self { depth, insert_sz, err_prof })
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
