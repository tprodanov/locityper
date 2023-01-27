use std::rc::Rc;
use std::fmt;
use std::cmp::{min, max, Ordering};
// use num_format::{Locale, WriteFormatted};
use crate::seq::contigs::{ContigId, ContigNames};

/// Genomic interval.
#[derive(Clone)]
pub struct Interval {
    contig_names: Rc<ContigNames>,
    contig_id: ContigId,
    start: u32,
    end: u32,
}

impl Interval {
    /// Create a new interval from contig names, contig id, start and end.
    pub fn new(contig_names: Rc<ContigNames>, contig_id: ContigId, start: u32, end: u32) -> Self {
        assert!(contig_id.ix() < contig_names.len(),
            "Cannot create interval with id {}, when there are {} contigs in total",
            contig_id.ix(), contig_names.len());
        assert!(start < end, "Cannot create interval {}:{}-{}", contig_id, start, end);
        Self { contig_names, contig_id, start, end }
    }

    /// Creates an empty interval from the contig id and interval start.
    pub fn new_empty(contig_names: Rc<ContigNames>, contig_id: ContigId, start: u32) -> Self {
        assert!(contig_id.ix() < contig_names.len(),
            "Cannot create interval with id {}, when there are {} contigs in total",
            contig_id.ix(), contig_names.len());
        Self {
            contig_names, contig_id, start,
            end: start,
        }
    }

    /// Parse interval name from string "name:start-end", where start is 1-based, inclusive.
    pub fn parse(contig_names: Rc<ContigNames>, s: &str) -> Result<Self, String> {
        if let Some((name, coord)) = s.split_once(':') {
            if let Some((start, end)) = coord.split_once('-') {
                if let Some(contig_id) = contig_names.get_id(name) {
                    let start: u32 = start.parse().map_err(|e: std::num::ParseIntError| e.to_string())?;
                    let end: u32 = end.parse().map_err(|e: std::num::ParseIntError| e.to_string())?;
                    return Ok(Self::new(contig_names, contig_id, start - 1, end));
                }
                return Err(format!("Unknown contig '{}'!", name));
            }
        }
        Err(format!("Cannot parse Interval '{}'", s))
    }

    /// Contig id.
    #[inline]
    pub fn contig_id(&self) -> ContigId {
        self.contig_id
    }

    /// Contig name.
    #[inline]
    pub fn contig_name(&self) -> &str {
        self.contig_names.name(self.contig_id)
    }

    pub fn contigs(&self) -> &Rc<ContigNames> {
        &self.contig_names
    }

    /// Interval start.
    #[inline]
    pub fn start(&self) -> u32 {
        self.start
    }

    /// Interval end.
    #[inline]
    pub fn end(&self) -> u32 {
        self.end
    }

    /// Returns the middle of the interval.
    pub fn middle(&self) -> u32 {
        (self.start + self.end) / 2
    }

    /// Get interval length.
    #[inline]
    pub fn len(&self) -> u32 {
        self.end - self.start
    }

    /// Distance between closest points on two intervals (0 if overlap).
    /// Returns None if on different contigs.
    pub fn distance(&self, other: &Interval) -> Option<u32> {
        if self.contig_id == other.contig_id {
            Some(max(
                (self.start + 1).saturating_sub(other.end),
                (other.start + 1).saturating_sub(self.end),
            ))
        } else {
            None
        }
    }

    /// Distance between furthest points on two intervals.
    /// Returns None if on different contigs.
    pub fn furthest_distance(&self, other: &Interval) -> Option<u32> {
        if self.contig_id == other.contig_id {
            Some(max(self.end, other.end) - min(self.start, other.start))
        } else {
            None
        }
    }

    /// Returns true if two intervals overlap by at least one base-pair.
    pub fn overlaps(&self, other: &Interval) -> bool {
        self.contig_id == other.contig_id && self.start < other.end && other.start < self.end
    }

    // /// Convert interval into string with comma separator.
    // pub fn comma_string(&self) -> String {
    //     let mut s = String::new();
    //     write!(s, "{}:", self.contig_name()).unwrap();
    //     s.write_formatted(&(self.start + 1), &Locale::en).unwrap();
    //     s.write_str("-").unwrap();
    //     s.write_formatted(&(self.end), &Locale::en).unwrap();
    //     s
    // }

    /// Convert interval into bed string (tab separator, 0-indexing).
    pub fn bed_string(&self) -> String {
        format!("{}\t{}\t{}", self.contig_name(), self.start, self.end)
    }
}

impl fmt::Debug for Interval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}={}:{}-{}", self.contig_id, self.contig_name(), self.start + 1, self.end)
    }
}

impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.contig_name(), self.start + 1, self.end)
    }
}

impl PartialEq for Interval {
    fn eq(&self, other: &Self) -> bool {
        Rc::ptr_eq(&self.contig_names, &other.contig_names) && self.contig_id == other.contig_id
            && self.start == other.start && self.end == other.end
    }
}

impl Eq for Interval {}

impl Ord for Interval {
    fn cmp(&self, other: &Self) -> Ordering {
        debug_assert!(Rc::ptr_eq(&self.contig_names, &other.contig_names),
            "Cannot compare intervals from different contig sets!");
        (self.contig_id, self.start, self.end).cmp(&(other.contig_id, other.start, other.end))
    }
}

impl PartialOrd for Interval {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
