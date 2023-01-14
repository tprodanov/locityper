use std::rc::Rc;
use std::fmt;
use htslib::bam::Record;
use intmap::{IntMap, Entry};
use crate::{
    seq::{
        contigs::ContigNames,
        // interv::Interval,
        // cigar::Cigar,
        aln::{READ_ENDS, ReadEnd, Alignment},
    },
    bg::err_prof::ErrorProfile,
    algo::hash::fnv1a,
};

/// Extension over alignment: store alignment, alignment probability and the read end.
#[derive(Clone, Debug)]
struct ExtAln {
    aln: Alignment,
    aln_prob: f64,
    read_end: ReadEnd,
}

impl ExtAln {
    /// Creates a new alignment extension from a htslib `Record`.
    fn from_record<E: ErrorProfile>(contigs: Rc<ContigNames>, record: &Record, err_prof: &E) -> Self {
        let aln = Alignment::from_record(contigs, record);
        let read_end = ReadEnd::from_record(record);
        let aln_prob = err_prof.ln_prob(aln.cigar(), read_end);
        Self { aln, aln_prob, read_end }
    }

    /// Get the index: contig_id * 2 + read_end.
    fn ix(&self) -> usize {
        self.aln.ref_interval().contig_id().ix() * READ_ENDS + self.read_end.ix()
    }
}

impl fmt::Display for ExtAln {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ExtAln[{}, {}, prob={:.2}]", self.read_end, self.aln, self.aln_prob)
    }
}

/// Preliminary read locations for a single read pair.
/// Store all potential alignments in a single vector.
#[derive(Clone, Debug)]
struct PrelimReadLocations(Vec<ExtAln>);

impl PrelimReadLocations {
    /// Creates an empty vector of preliminary locations.
    fn new() -> Self {
        Self(Vec::new())
    }

    /// Push a new `ExtAln`.
    #[inline]
    fn push(&mut self, ext_aln: ExtAln) {
        self.0.push(ext_aln);
    }

    /// Sort alignments first by contig, and then by the read end.
    #[inline]
    fn sort(&mut self) {
        self.0.sort_by_key(|ext_aln| ext_aln.ix())
        // self.0.sort_by(|a, b| b.aln_prob.total_cmp(&a.aln_prob))
    }

    fn debug(&mut self) {
        self.sort();
        for aln in self.0.iter() {
            log::debug!("        {}", aln);
        }
    }
}

/// Preliminary read locations for multiple reads.
pub struct AllLocations<'a, E: ErrorProfile> {
    contigs: Rc<ContigNames>,
    err_prof: &'a E,
    reads: IntMap<PrelimReadLocations>,
}

impl<'a, E: ErrorProfile> AllLocations<'a, E> {
    /// Creates a new, empty, `AllLocations`.
    pub fn new(contigs: Rc<ContigNames>, err_prof: &'a E) -> Self {
        Self {
            contigs, err_prof,
            reads: IntMap::new(),
        }
    }

    /// Create a new `ExtAln` from the record, and push it into IntMap.
    pub fn push(&mut self, record: &Record) {
        let hash = fnv1a(record.qname());
        let ext_aln = ExtAln::from_record(Rc::clone(&self.contigs), record, self.err_prof);
        log::debug!("    {}     {}", ext_aln, hash);
        match self.reads.entry(hash) {
            Entry::Occupied(mut entry) => entry.get_mut().push(ext_aln),
            Entry::Vacant(entry) => entry.insert(PrelimReadLocations::new()).push(ext_aln),
        }
    }

    /// Push multiple records in the `AllLocations`.
    pub fn extend<'b>(&mut self, records: impl Iterator<Item = &'b Record>) {
        for record in records {
            self.push(record);
        }
    }

    pub fn debug(&mut self) {
        log::debug!("Debug read locations!");
        for (key, val) in self.reads.iter_mut() {
            log::debug!("    Read {}", key);
            val.debug();
        }
    }
}