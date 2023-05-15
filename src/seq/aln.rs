use std::{
    fmt,
    sync::Arc,
};
use htslib::bam::Record;
use crate::{
    seq::{
        Interval, ContigId, ContigNames,
        cigar::Cigar,
    },
    bg::{
        err_prof::ErrorProfile,
        insertsz::InsertDistr,
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

/// Alignment location, strand, read-end, and alignment ln-probability.
/// CIGAR is not stored for now, can be added later, if needed.
pub struct Alignment {
    interval: Interval,
    strand: Strand,
    read_end: ReadEnd,
    ln_prob: f64,
}

impl Alignment {
    /// Creates a new alignment extension from a htslib `Record`.
    pub fn from_record(
        record: &Record,
        cigar: &Cigar,
        read_end: ReadEnd,
        contigs: Arc<ContigNames>,
        err_prof: &ErrorProfile,
    ) -> Self
    {
        Self {
            ln_prob: err_prof.ln_prob(cigar),
            strand: Strand::from_record(record),
            interval: ref_interval(record, cigar, contigs),
            read_end,
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

    /// Returns sorted key: first sort by contig id, then by read end.
    pub fn sort_key(&self) -> u16 {
        // Total number of contigs is checked in `ContigNames::new`.
        (self.interval.contig_id().get() << 1) | (self.read_end.as_u16())
    }

    /// Returns probability of two alignments: probability of first * probability of second * insert size probability.
    pub fn paired_prob(&self, mate: &Alignment, insert_distr: &InsertDistr) -> f64 {
        let insert_size = self.interval.furthest_distance(&mate.interval)
            .expect("Alignments must be on the same contig");
        self.ln_prob + mate.ln_prob + insert_distr.ln_prob(insert_size, self.strand == mate.strand)
    }
}
