use std::{
    fmt,
    sync::Arc,
    ops::{Deref, DerefMut},
    borrow::Cow,
};
use htslib::bam::Record;
use crate::{
    bg::InsertDistr,
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
        match self {
            ReadEnd::First => write!(f, "Mate1"),
            ReadEnd::Second => write!(f, "Mate2"),
        }
    }
}

impl fmt::Display for ReadEnd {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ReadEnd::First => write!(f, "1"),
            ReadEnd::Second => write!(f, "2"),
        }
    }
}

/// Returns reference interval for the alignment. Panics if the record is unmapped.
pub fn ref_interval(record: &Record, cigar: &Cigar, contigs: Arc<ContigNames>) -> Interval {
    assert!(!record.is_unmapped(),
        "Cannot get alignment interval for an unmapped record {}", String::from_utf8_lossy(record.qname()));
    let contig_id = ContigId::new(record.tid());
    let start = u32::try_from(record.pos()).unwrap();
    Interval::new(contigs, contig_id, start, start + cigar.ref_len())
}

/// Light alignment information. Use `Alignemnt` to store alignment name and CIGAR as well.
#[derive(Clone)]
pub struct LightAlignment {
    interval: Interval,
    strand: Strand,
    read_end: ReadEnd,
    ln_prob: f64,
}

impl LightAlignment {
    /// Creates a new alignment information from a bam record.
    pub fn new(
        record: &Record,
        cigar: &Cigar,
        read_end: ReadEnd,
        contigs: Arc<ContigNames>,
        ln_prob: f64,
    ) -> Self
    {
        Self {
            interval: ref_interval(record, cigar, contigs),
            strand: Strand::from_record(record),
            ln_prob, read_end,
        }
    }

    /// Get contig id of the alignment.
    pub fn contig_id(&self) -> ContigId {
        self.interval.contig_id()
    }

    /// Reference interval, to which the record is aligned.
    pub fn interval(&self) -> &Interval {
        &self.interval
    }

    /// Read end of the alignment.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Read end of the alignment.
    pub fn read_end(&self) -> ReadEnd {
        self.read_end
    }

    /// Ln-probability of the alignment.
    pub fn ln_prob(&self) -> f64 {
        self.ln_prob
    }

    /// Sets ln-probability of the alignment.
    pub fn set_ln_prob(&mut self, ln_prob: f64) {
        self.ln_prob = ln_prob;
    }

    /// Returns sorted key: first sort by contig id, then by read end.
    pub fn sort_key(&self) -> u16 {
        // Total number of contigs is checked in `ContigNames::new`.
        (self.interval.contig_id().get() << 1) | (self.read_end.as_u16())
    }

    /// Returns insert size (distance between smallest start and largest end).
    /// Panics if two mates are not on the same strand.
    #[inline]
    pub fn insert_size(&self, mate: &LightAlignment) -> u32 {
        self.interval.furthest_distance(&mate.interval)
            .expect("Alignments must be on the same contig")
    }

    /// Returns false for `FR/RF` and true for `FF/RR`.
    #[inline]
    pub fn pair_orientation(&self, mate: &LightAlignment) -> bool {
        self.strand == mate.strand
    }

    /// Returns insert size probability.
    pub fn insert_size_prob(&self, mate: &LightAlignment, insert_distr: &InsertDistr) -> f64 {
        insert_distr.ln_prob(self.insert_size(mate), self.pair_orientation(mate))
    }

    /// Returns probability of two alignments: probability of first * probability of second * insert size probability.
    pub fn paired_prob(&self, mate: &LightAlignment, insert_distr: &InsertDistr) -> f64 {
        self.ln_prob + mate.ln_prob + self.insert_size_prob(mate, insert_distr)
    }
}

/// Full alignment information: light alignment + name + CIGAR.
#[derive(Clone)]
pub struct Alignment {
    name: Vec<u8>,
    cigar: Cigar,
    info: LightAlignment,
}

impl Alignment {
    pub fn new(
        record: &Record,
        cigar: Cigar,
        read_end: ReadEnd,
        contigs: Arc<ContigNames>,
        ln_prob: f64,
    ) -> Self
    {
        Self {
            name: record.qname().to_vec(),
            info: LightAlignment::new(record, &cigar, read_end, contigs, ln_prob),
            cigar,
        }
    }

    /// Returns record name as non-decoded bytes.
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    /// Converts name into UTF-8.
    pub fn name_utf8(&self) -> Cow<'_, str> {
        String::from_utf8_lossy(&self.name)
    }

    /// Returns alignment CIGAR.
    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }
}

impl Deref for Alignment {
    type Target = LightAlignment;

    fn deref(&self) -> &LightAlignment {
        &self.info
    }
}

impl DerefMut for Alignment {
    fn deref_mut(&mut self) -> &mut LightAlignment {
        &mut self.info
    }
}
