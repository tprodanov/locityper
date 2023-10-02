use std::{
    fmt,
    sync::Arc,
    ops::{Deref, DerefMut},
    borrow::Cow,
    cmp::{min, max},
};
use htslib::bam::Record;
use crate::{
    bg::{InsertDistr, err_prof::OperCounts},
    seq::{
        Interval, ContigId, ContigNames,
        cigar::{Cigar, Operation},
    },
};

/// Newtype over strand: false = negative, true = positive.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    #[inline]
    pub fn from_record(record: &Record) -> Self {
        if record.is_reverse() { Self::Reverse } else { Self::Forward }
    }
}

impl fmt::Debug for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::Forward => f.write_str("+"),
            Self::Reverse => f.write_str("-"),
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

/// Light alignment information. Use `Alignemnt` to store alignment name and CIGAR as well.
#[derive(Clone)]
pub struct Alignment {
    interval: Interval,
    cigar: Cigar,
    strand: Strand,
    read_end: ReadEnd,
    ln_prob: f64,
}

impl Alignment {
    /// Creates a new alignment information from a bam record.
    pub fn new(
        record: &Record,
        cigar: Cigar,
        read_end: ReadEnd,
        contigs: Arc<ContigNames>,
        ln_prob: f64,
    ) -> Self
    {
        assert!(!record.is_unmapped(),
            "Cannot get alignment interval for an unmapped record {}", String::from_utf8_lossy(record.qname()));
        let contig_id = ContigId::new(record.tid());
        let start = u32::try_from(record.pos()).unwrap();
        let interval = Interval::new(contigs, contig_id, start, start + cigar.ref_len());
        Self {
            interval,
            strand: Strand::from_record(record),
            cigar, ln_prob, read_end,
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

    /// Returns alignment CIGAR.
    pub fn cigar(&self) -> &Cigar {
        &self.cigar
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
    pub fn insert_size(&self, mate: &Alignment) -> u32 {
        self.interval.furthest_distance(&mate.interval)
            .expect("Alignments must be on the same contig")
    }

    /// Returns false for `FR/RF` and true for `FF/RR`.
    #[inline]
    pub fn pair_orientation(&self, mate: &Alignment) -> bool {
        self.strand == mate.strand
    }

    /// Returns probability of two alignments: probability of first * probability of second * insert size probability.
    pub fn paired_prob(&self, mate: &Alignment, insert_distr: &InsertDistr) -> f64 {
        self.ln_prob + mate.ln_prob + insert_distr.ln_prob(self.insert_size(mate), self.pair_orientation(mate))
    }

    /// Counts operations in the alignment, excluding operations outside the boundary of the `region`.
    pub fn count_region_operations(&self, region: &Interval) -> OperCounts<u32> {
        debug_assert_eq!(self.interval.contig_id(), region.contig_id());
        let region_start = region.start();
        let region_end = region.end();

        let mut counts = OperCounts::default();
        let mut rpos = self.interval.start();
        let mut first = true;
        for item in self.cigar.iter() {
            let oplen = item.len();
            match item.operation() {
                Operation::Equal => {
                    // Equal, Diff, Del: same clause, different operations.
                    counts.matches += min(rpos + oplen, region_end).saturating_sub(max(rpos, region_start));
                    rpos += oplen;
                }
                Operation::Diff => {
                    counts.mismatches += min(rpos + oplen, region_end).saturating_sub(max(rpos, region_start));
                    rpos += oplen;
                }
                Operation::Del => {
                    counts.deletions += min(rpos + oplen, region_end).saturating_sub(max(rpos, region_start));
                    rpos += oplen;
                }
                Operation::Ins => {
                    // Only add insertion if it is within boundaries.
                    counts.insertions += if region_start <= rpos && rpos < region_end { oplen } else { 0 };
                }
                Operation::Soft => {
                    counts.clipping += if first {
                        // If first operation of the CIGAR,
                        // limit clipping to the distance between the start of the region and current position.
                        min(oplen, rpos.saturating_sub(region_start))
                    } else {
                        // Otherwise, we should be at the end of the CIGAR,
                        // limit clipping to the distance between curr position and the end.
                        min(oplen, region_end.saturating_sub(rpos))
                    };
                }
                op => panic!("Unsupported CIGAR operation {}", op),
            }
            first = false;
        }
        counts
    }

    /// Counts operations in the alignment, excluding operations outside the boundary `0..contig_end`.
    /// In contrast to `count_region_operations`, this function assumes that there is no alignment
    /// out of the contig boundaries.
    pub fn count_region_operations_fast(&self, contig_end: u32) -> OperCounts<u32> {
        let mut counts = OperCounts::default();
        let mut first = true;
        for item in self.cigar.iter() {
            let oplen = item.len();
            match item.operation() {
                Operation::Equal => counts.matches += oplen,
                Operation::Diff => counts.mismatches += oplen,
                Operation::Del => counts.deletions += oplen,
                Operation::Ins => counts.insertions += oplen,
                Operation::Soft => {
                    counts.clipping += if first {
                        // If first operation of the CIGAR,
                        // limit clipping to the distance between contig start (0) and alignment start.
                        min(oplen, self.interval.start())
                    } else {
                        // Otherwise, we should be at the end of the CIGAR,
                        // limit clipping to the distance distance between curr position and contig end.
                        min(oplen, contig_end.saturating_sub(self.interval.end()))
                    };
                }
                op => panic!("Unsupported CIGAR operation {}", op),
            }
            first = false;
        }
        counts
    }
}

/// Full alignment information: light alignment + name + CIGAR.
#[derive(Clone)]
pub struct NamedAlignment {
    name: Vec<u8>,
    aln: Alignment,
}

impl NamedAlignment {
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
            aln: Alignment::new(record, cigar, read_end, contigs, ln_prob),
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

    pub fn take_light_aln(self) -> Alignment {
        self.aln
    }
}

impl Deref for NamedAlignment {
    type Target = Alignment;

    fn deref(&self) -> &Alignment {
        &self.aln
    }
}

impl DerefMut for NamedAlignment {
    fn deref_mut(&mut self) -> &mut Alignment {
        &mut self.aln
    }
}
