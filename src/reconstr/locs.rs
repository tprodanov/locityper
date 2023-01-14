use std::rc::Rc;
use std::fmt;
use htslib::bam::Record;
use intmap::{IntMap, Entry};
use crate::{
    seq::{
        contigs::ContigNames,
        interv::Interval,
        cigar::Cigar,
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
}

impl fmt::Display for ExtAln {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ExtAln[{}, {}, prob={:.2}]", self.read_end, self.aln, self.aln_prob)
    }
}

/// Preliminary read locations for a single read pair.
/// `self.0[contig_id * 2 + read_end]` - Vector of potential alignment at a certain contig and read end.
struct PrelimReadLocations(Vec<Vec<ExtAln>>);

impl PrelimReadLocations {
    /// Creates an empty vector of preliminary locations.
    fn new(n_contigs: usize) -> Self {
        Self(vec![Vec::new(); READ_ENDS * n_contigs])
    }

    /// Push a new `ExtAln`.
    fn push(&mut self, ext_aln: ExtAln) {
        let ix = ext_aln.aln.ref_interval().contig_id().ix() * READ_ENDS + ext_aln.read_end.ix();
        self.0[ix].push(ext_aln);
    }
}

/// Preliminary read locations for multiple reads.
struct AllLocations<E: ErrorProfile> {
    contigs: Rc<ContigNames>,
    err_prof: E,
    reads: IntMap<PrelimReadLocations>,
}

impl<E: ErrorProfile> AllLocations<E> {
    /// Creates a new, empty, `AllLocations`.
    fn new(contigs: Rc<ContigNames>, err_prof: E) -> Self {
        Self {
            contigs, err_prof,
            reads: IntMap::new(),
        }
    }

    /// Create a new `ExtAln` from the record, and push it into IntMap.
    fn push(&mut self, record: &Record) {
        let hash = fnv1a(record.qname());
        let ext_aln = ExtAln::from_record(Rc::clone(&self.contigs), record, &self.err_prof);
        match self.reads.entry(hash) {
            Entry::Occupied(mut entry) => entry.get_mut().push(ext_aln),
            Entry::Vacant(entry) => entry.insert(PrelimReadLocations::new(self.contigs.len())).push(ext_aln),
        }
    }

    /// Push multiple records in the `AllLocations`.
    fn extend<'a>(&mut self, records: impl Iterator<Item = &'a Record>) {
        for record in records {
            self.push(record);
        }
    }
}