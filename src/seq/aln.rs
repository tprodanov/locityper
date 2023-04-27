use std::{
    fmt,
    sync::Arc,
};
use htslib::bam::{
    Record,
    ext::BamRecordExtensions,
};
use crate::{
    seq::{
        Interval, ContigId, ContigNames,
        cigar::Cigar,
    },
};

/// Newtype over strand: false = negative, true = positive.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Strand(bool);

impl Strand {
    #[inline]
    pub fn new(forward: bool) -> Self {
        Self(forward)
    }

    #[inline]
    pub fn from_record(record: &Record) -> Self {
        Self(!record.is_reverse())
    }

    #[inline]
    pub fn is_forward(self) -> bool {
        self.0
    }
}

impl fmt::Debug for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.0 {
            write!(f, "+")
        } else {
            write!(f, "-")
        }
    }
}

/// Read-end: first or second.
#[derive(Clone, Copy, PartialOrd, Ord, PartialEq, Eq)]
pub enum ReadEnd {
    First,
    Second,
}

impl ReadEnd {
    /// Gets record read-end information.
    #[inline]
    pub fn from_record(record: &Record) -> ReadEnd {
        if record.is_last_in_template() {
            ReadEnd::Second
        } else {
            ReadEnd::First
        }
    }

    /// Converts read-end into integer (First => 0, Second => 1).
    pub const fn as_u16(self) -> u16 {
        match self {
            ReadEnd::First => 0,
            ReadEnd::Second => 1,
        }
    }

    /// Converts read-end into integer (First => 0, Second => 1).
    pub const fn ix(self) -> usize {
        match self {
            ReadEnd::First => 0,
            ReadEnd::Second => 1,
        }
    }
}

impl fmt::Debug for ReadEnd {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl fmt::Display for ReadEnd {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ReadEnd::First => write!(f, "Mate1"),
            ReadEnd::Second => write!(f, "Mate2"),
        }
    }
}

/// In total, there can be only two read ends.
pub const READ_ENDS: usize = 2;

/// Read alignment.
/// Stores reference interval, strand, extended CIGAR and read-end information.
/// Does not store read name.
#[derive(Clone)]
pub struct Alignment {
    ref_interval: Interval,
    strand: Strand,
    /// Extended CIGAR with X/= instead of M.
    cigar: Cigar,
}

impl Alignment {
    /// Creates a new Alignment from the record.
    pub fn from_record(record: &Record, contigs: Arc<ContigNames>) -> Self {
        let cigar = Cigar::infer_ext_cigar_md(record, ());
        let contig_id = ContigId::new(record.tid());
        let start = u32::try_from(record.pos()).unwrap();
        let ref_interval = if cigar.is_empty() {
            assert!(record.is_unmapped(), "Read is mapped, but has empty CIGAR!");
            Interval::new_empty(contigs, contig_id, start)
        } else {
            Interval::new(contigs, contig_id, start, start + cigar.ref_len())
        };
        debug_assert_eq!(cigar.ref_len() as i64, record.reference_end() - record.pos());
        Self {
            ref_interval, cigar,
            strand: Strand::from_record(record),
        }
    }

    /// Returns reference interval.
    pub fn ref_interval(&self) -> &Interval {
        &self.ref_interval
    }

    /// Consumes alignments and returns reference interval.
    pub fn take_interval(self) -> Interval {
        self.ref_interval
    }

    /// Returns contig id of the alignment.
    pub fn contig_id(&self) -> ContigId {
        self.ref_interval.contig_id()
    }

    /// Returns extended alignment CIGAR.
    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    /// Returns alignment strand.
    pub fn strand(&self) -> Strand {
        self.strand
    }
}

impl fmt::Debug for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Aln({:?}:{:?}, {:?})", self.ref_interval, self.strand, self.cigar)
    }
}

impl fmt::Display for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Aln({}:{}, {})", self.ref_interval, self.strand, self.cigar)
    }
}
