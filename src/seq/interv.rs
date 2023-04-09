use std::{
    fmt,
    io::{self, Read, Seek},
    rc::Rc,
    cmp::{min, max, Ordering},
};
use regex::Regex;
use lazy_static::lazy_static;
use const_format::formatcp;
use bio::io::fasta;
use crate::Error;
use super::{ContigId, ContigNames};

/// SAM file specification (https://samtools.github.io/hts-specs/SAMv1.pdf) allows
/// many possible contig names: `[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*`.
///
/// However, here we provide a bit more limiting pattern for contig names.
const CONTIG_PATTERN: &'static str = r"[0-9A-Za-z][0-9A-Za-z+._|~=@^-]*";

/// Interval: `contig:start-end`.
const INTERVAL_PATTERN: &'static str = formatcp!("^({}):([0-9]+)-([0-9]+)$", CONTIG_PATTERN);

/// Name of the interval (almost the same as `CONTIG_PATTERN`, but includes `:`).
const NAME_PATTERN: &'static str = r"[0-9A-Za-z][0-9A-Za-z:+._|~=@^-]*";

/// Optionally named interval: `contig:start-end[@name]`.
const NAMED_INTERVAL_PATTERN: &'static str = formatcp!("^({}):([0-9]+)-([0-9]+)(={})?$", CONTIG_PATTERN, NAME_PATTERN);

/// Genomic interval.
#[derive(Clone)]
pub struct Interval {
    contigs: Rc<ContigNames>,
    contig_id: ContigId,
    start: u32,
    end: u32,
}

impl Interval {
    /// Create a new interval from contig names, contig id, start and end.
    pub fn new(contigs: Rc<ContigNames>, contig_id: ContigId, start: u32, end: u32) -> Self {
        assert!(contig_id.ix() < contigs.len(),
            "Cannot create interval with contig id {} (total {} contigs)", contig_id.ix(), contigs.len());
        assert!(start < end, "Cannot create an empty interval {}:{}-{}",
            contigs.get_name(contig_id), start + 1, end);
        Self { contigs, contig_id, start, end }
    }

    /// Creates an empty interval from the contig id and interval start.
    pub fn new_empty(contigs: Rc<ContigNames>, contig_id: ContigId, start: u32) -> Self {
        assert!(contig_id.ix() < contigs.len(),
            "Cannot create interval with id {}, when there are {} contigs in total",
            contig_id.ix(), contigs.len());
        Self {
            contigs, contig_id, start,
            end: start,
        }
    }

    /// Creates a new interval, covering the full contig.
    pub fn full_contig(contigs: Rc<ContigNames>, contig_id: ContigId) -> Self {
        let end = contigs.get_len(contig_id);
        Self::new(contigs, contig_id, 0, end)
    }

    /// Parses interval from string "name:start-end", where start is 1-based, inclusive.
    pub fn parse(s: &str, contigs: &Rc<ContigNames>) -> Result<Self, Error> {
        lazy_static! {
            static ref RE: Regex = Regex::new(INTERVAL_PATTERN).unwrap();
        }
        let captures = RE.captures(s).ok_or_else(|| Error::ParsingError(format!("Cannot parse interval '{}'", s)))?;
        let contig_id = contigs.try_get_id(&captures[1])
            .ok_or_else(|| Error::ParsingError(format!("Cannot parse interval '{}': unknown contig", s)))?;
        let start: u32 = captures[2].parse()
            .map_err(|_| Error::ParsingError(format!("Cannot parse interval '{}'", s)))?;
        let end: u32 = captures[3].parse()
            .map_err(|_| Error::ParsingError(format!("Cannot parse interval '{}'", s)))?;
        Ok(Self::new(Rc::clone(contigs), contig_id, start - 1, end))
    }

    /// Parses interval from iterator over strings. Moves iterator by three positions (chrom, start, end).
    /// Start is 0-based, inclusive, end is 0-based, non-inclusive.
    pub fn parse_bed<'a, I>(split: &mut I, contigs: &Rc<ContigNames>) -> Result<Self, Error>
    where I: Iterator<Item = &'a str>,
    {
        let contig_name = split.next()
            .ok_or_else(|| Error::ParsingError("Cannot parse BED line, not enough columns".to_string()))?;
        let contig_id = contigs.try_get_id(contig_name)
            .ok_or_else(|| Error::ParsingError(format!("Cannot parse interval '{}': unknown contig", contig_name)))?;
        let start = split.next()
            .ok_or_else(|| Error::ParsingError("Cannot parse BED line, not enough columns".to_string()))?
            .parse::<u32>().map_err(|e| Error::ParsingError(format!("Cannot parse BED line: {}", e)))?;
        let end = split.next()
            .ok_or_else(|| Error::ParsingError("Cannot parse BED line, not enough columns".to_string()))?
            .parse::<u32>().map_err(|e| Error::ParsingError(format!("Cannot parse BED line: {}", e)))?;
        Ok(Self::new(Rc::clone(contigs), contig_id, start, end))
    }

    /// Contig id.
    #[inline]
    pub fn contig_id(&self) -> ContigId {
        self.contig_id
    }

    /// Contig name.
    #[inline]
    pub fn contig_name(&self) -> &str {
        self.contigs.get_name(self.contig_id)
    }

    pub fn contigs(&self) -> &Rc<ContigNames> {
        &self.contigs
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

    /// Returns start and end of the interval.
    pub fn range(&self) -> (u32, u32) {
        (self.start, self.end)
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
    pub fn bed_fmt<'a>(&'a self) -> BedFormat<'a> {
        BedFormat(self)
    }

    /// Expand the interval by `left` and `right` bp to the left and to the right.
    /// Limit new start to 0 and new end to the contig length.
    pub fn expand(&self, left: u32, right: u32) -> Self {
        Self {
            contigs: Rc::clone(&self.contigs),
            contig_id: self.contig_id,
            start: self.start.saturating_sub(left),
            end: min(self.end + right, self.contigs.get_len(self.contig_id)),
        }
    }

    /// Fetches sequence of the interval from an indexed fasta reader.
    pub fn fetch_seq<R: Read + Seek>(&self, fasta: &mut fasta::IndexedReader<R>) -> io::Result<Vec<u8>> {
        fasta.fetch(self.contig_name(), u64::from(self.start), u64::from(self.end))?;
        let mut seq = Vec::with_capacity(self.len() as usize);
        fasta.read(&mut seq)?;
        crate::seq::standardize(&mut seq);
        Ok(seq)
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
        Rc::ptr_eq(&self.contigs, &other.contigs) && self.contig_id == other.contig_id
            && self.start == other.start && self.end == other.end
    }
}

impl Eq for Interval {}

impl Ord for Interval {
    fn cmp(&self, other: &Self) -> Ordering {
        debug_assert!(Rc::ptr_eq(&self.contigs, &other.contigs),
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

#[doc(hidden)]
pub struct BedFormat<'a>(&'a Interval);

impl<'a> fmt::Display for BedFormat<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\t{}\t{}", self.0.contig_name(), self.0.start, self.0.end)
    }
}

/// Stores an interval with its name.
#[derive(Clone)]
pub struct NamedInterval {
    interval: Interval,
    /// Name of the interval, must satisfy the `NAME_PATTERN` regular expression.
    name: String,
    /// True if `name` is was set explicitely (not equal to the interval itself).
    explicit_name: bool,
}

impl NamedInterval {
    /// Create a named interval.
    /// If name is not provided, set it to `contig:start-end`.
    pub fn new(interval: Interval, name: Option<&str>) -> Result<Self, Error> {
        let explicit_name = name.is_some();
        let name = name.map(Into::into).unwrap_or_else(|| interval.to_string());

        lazy_static! {
            static ref NAME_RE: Regex = Regex::new(formatcp!("^{}$", NAME_PATTERN)).unwrap();
        }
        if !NAME_RE.is_match(&name) {
            Err(Error::ParsingError(format!("Interval name '{}' contains forbidden symbols", name)))
        } else {
            Ok(Self { name, interval, explicit_name })
        }
    }

    /// Parse interval name from string "contig:start-end[@name]".
    /// If name is not set, set it to `contig:start-end`.
    pub fn parse(s: &str, contigs: &Rc<ContigNames>) -> Result<Self, Error> {
        lazy_static! {
            static ref RE: Regex = Regex::new(NAMED_INTERVAL_PATTERN).unwrap();
        }
        let captures = RE.captures(s).ok_or_else(|| Error::ParsingError(format!("Cannot parse interval '{}'", s)))?;
        let contig_id = contigs.try_get_id(&captures[1])
            .ok_or_else(|| Error::ParsingError(format!("Cannot parse interval '{}': unknown contig", s)))?;
        let start: u32 = captures[2].parse()
            .map_err(|_| Error::ParsingError(format!("Cannot parse interval '{}'", s)))?;
        let end: u32 = captures[3].parse()
            .map_err(|_| Error::ParsingError(format!("Cannot parse interval '{}'", s)))?;
        let interval = Interval::new(Rc::clone(contigs), contig_id, start - 1, end);
        NamedInterval::new(interval, captures.get(4).map(|m| &m.as_str()[1..]))
    }

    /// Parses interval from iterator over strings. Moves iterator by three positions (chrom, start, end).
    /// Start is 0-based, inclusive, end is 0-based, non-inclusive.
    pub fn parse_bed<'a, I>(split: &mut I, contigs: &Rc<ContigNames>) -> Result<Self, Error>
    where I: Iterator<Item = &'a str>,
    {
        let interval = Interval::parse_bed(split, contigs)?;
        Self::new(interval, split.next())
    }

    /// Return interval.
    pub fn interval(&self) -> &Interval {
        &self.interval
    }

    /// Returns interval name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Creates a new interval with updated start and end.
    /// Name remains the same if it was set explicitely, otherwise it is created from the new interval positions.
    pub fn with_new_range(&self, new_start: u32, new_end: u32) -> Self {
        let new_interval = Interval::new(
            Rc::clone(self.interval.contigs()), self.interval.contig_id(), new_start, new_end);
        Self::new(new_interval, if self.explicit_name { Some(&self.name) } else { None }).unwrap()
    }
}

impl fmt::Debug for NamedInterval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.explicit_name {
            write!(f, "{:?}@{}", self.interval, self.name)
        } else {
            fmt::Debug::fmt(&self.interval, f)
        }
    }
}

impl fmt::Display for NamedInterval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.explicit_name {
            write!(f, "{} ({})", self.name, self.interval)
        } else {
            fmt::Display::fmt(&self.interval, f)
        }
    }
}

/// Split intervals (start, end, depth) into disjoint intervals with fixed depth.
/// Depth is summed for overlapping intervals.
/// Only output intervals for which `keep(depth)` is true.
pub(crate) fn split_disjoint<P>(intervals: &[(u32, u32, i32)], mut keep: P) -> Vec<(u32, u32, i32)>
where P: FnMut(i32) -> bool,
{
    let mut res = Vec::new();
    if intervals.len() == 0 {
        return res;
    }
    let mut endpoints = Vec::with_capacity(intervals.len() * 2);
    // Starts of interval: positive depth.
    endpoints.extend(intervals.iter().map(|(start, _end, depth)| (*start, *depth)));
    // Ends of interval: negative depth.
    endpoints.extend(intervals.iter().map(|(_start, end, depth)| (*end, -depth)));
    endpoints.sort_unstable();

    let mut curr_start = endpoints[0].0;
    let mut pushed_end = curr_start;
    let mut curr_depth = 0;
    let mut pushed_depth = 0;

    for (curr_end, depth) in endpoints {
        // Add interval (curr_start, curr_end, curr_depth), if needed.
        // and then update curr_depth += depth.
        if curr_start < curr_end && keep(curr_depth) {
            if pushed_end == curr_start && pushed_depth == curr_depth {
                res.last_mut().unwrap().1 = curr_end;
            } else {
                res.push((curr_start, curr_end, curr_depth));
                pushed_depth = curr_depth;
            }
            pushed_end = curr_end;
        }
        curr_depth += depth;
        curr_start = curr_end;
    }
    // No need to do anything after the last endpoint.
    debug_assert_eq!(curr_depth, 0);
    res
}
