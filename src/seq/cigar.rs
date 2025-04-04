//! Various functions and structures related to sequence alignment.

use std::{
    fmt::{self, Write},
    cmp::min,
    ops::Index,
};
use htslib::bam::record;
use crate::ext::vec::VecOrNone;

/// Subset of CIGAR operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Operation {
    Match,
    Equal,
    Diff,
    Soft,
    Hard,
    Ins,
    Del,
}

impl Operation {
    /// Does the Cigar operation consume reference sequence?
    #[inline]
    pub const fn consumes_ref(self) -> bool {
        match self {
            Operation::Match | Operation::Equal | Operation::Diff | Operation::Del => true,
            Operation::Soft | Operation::Hard | Operation::Ins => false,
        }
    }

    /// Does the Cigar operation consume query sequence?
    #[inline]
    pub const fn consumes_query(self) -> bool {
        match self {
            Operation::Match | Operation::Equal | Operation::Diff | Operation::Ins | Operation::Soft => true,
            Operation::Hard | Operation::Del => false,
        }
    }

    /// Does the Cigar operation consume both reference and query sequences?
    #[inline]
    pub const fn consumes_both(self) -> bool {
        match self {
            Operation::Match | Operation::Equal | Operation::Diff => true,
            Operation::Soft | Operation::Hard | Operation::Ins | Operation::Del => false,
        }
    }

    /// Convert operation into char.
    pub const fn to_char(self) -> char {
        match self {
            Operation::Match => 'M',
            Operation::Equal => '=',
            Operation::Diff => 'X',
            Operation::Ins => 'I',
            Operation::Del => 'D',
            Operation::Soft => 'S',
            Operation::Hard => 'H',
        }
    }

    /// Get operation from a raw u32 value.
    pub const fn from_u32(val: u32) -> Operation {
        match val {
            0 => Operation::Match,
            1 => Operation::Ins,
            2 => Operation::Del,
            // 3 => RefSkip,
            4 => Operation::Soft,
            5 => Operation::Hard,
            // 6 => Padding,
            7 => Operation::Equal,
            8 => Operation::Diff,
            _ => panic!("Unsupported cigar operation"),
        }
    }

    /// Convert operation into an index (values in 0..RAW_OPERATIONS).
    pub const fn ix(self) -> usize {
        match self {
            Operation::Match => 0,
            Operation::Ins => 1,
            Operation::Del => 2,
            // RefSkip,
            Operation::Soft => 4,
            Operation::Hard => 5,
            // Padding,
            Operation::Equal => 7,
            Operation::Diff => 8,
        }
    }

    /// inverts reference and query.
    pub const fn invert(self) -> Self {
        match self {
            Operation::Ins | Operation::Soft => Operation::Del,
            Operation::Del => Operation::Ins,

            Operation::Match => Operation::Match,
            Operation::Equal => Operation::Equal,
            Operation::Diff => Operation::Diff,

            Operation::Hard => panic!("Cigar::invert not supported with Hard clipping"),
        }
    }
}

impl fmt::Display for Operation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_char(self.to_char())
    }
}

/// Tuple operation + length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarItem {
    op: Operation,
    len: u32,
}

impl CigarItem {
    /// Creates a new tuple.
    #[inline(always)]
    pub fn new(op: Operation, len: u32) -> Self {
        CigarItem { op, len }
    }

    /// Get operation.
    #[inline(always)]
    pub fn operation(&self) -> Operation {
        self.op
    }

    /// Get length.
    #[inline(always)]
    pub fn len(&self) -> u32 {
        self.len
    }

    /// Creates a new tuple from a raw u32 value.
    #[inline(always)]
    pub fn from_u32(val: u32) -> Self {
        Self {
            op: Operation::from_u32(val & 0b1111),
            len: val >> 4,
        }
    }

    pub const fn to_htslib(self) -> record::Cigar {
        match self.op {
            Operation::Match => record::Cigar::Match(self.len),
            Operation::Ins => record::Cigar::Ins(self.len),
            Operation::Del => record::Cigar::Del(self.len),
            Operation::Soft => record::Cigar::SoftClip(self.len),
            Operation::Hard => record::Cigar::HardClip(self.len),
            Operation::Equal => record::Cigar::Equal(self.len),
            Operation::Diff => record::Cigar::Diff(self.len),
        }
    }

    pub fn invert(&self) -> Self {
        Self {
            op: self.op.invert(),
            len: self.len,
        }
    }
}

impl fmt::Display for CigarItem {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}", self.len, self.op)
    }
}

/// Wrapper over vector of `CigarItem`.
#[derive(Default, Clone)]
pub struct Cigar {
    tuples: Vec<CigarItem>,
    /// Reference length.
    rlen: u32,
    /// Query length.
    qlen: u32,
}

impl Cigar {
    /// Create an empty CIGAR.
    pub fn new() -> Cigar {
        Cigar::default()
    }

    /// Creates cigar of `len` matches (=).
    pub fn new_full_match(len: u32) -> Self {
        assert!(len > 0);
        Self {
            tuples: vec![CigarItem::new(Operation::Equal, len)],
            rlen: len,
            qlen: len,
        }
    }

    /// Creates cigar of `len` mismatches (X).
    pub fn new_full_mismatch(len: u32) -> Self {
        assert!(len > 0);
        Self {
            tuples: vec![CigarItem::new(Operation::Diff, len)],
            rlen: len,
            qlen: len,
        }
    }

    /// Converts raw CIGAR to this struct.
    pub fn from_raw(raw_cigar: &[u32]) -> Self {
        let mut res = Cigar::new();
        for &val in raw_cigar.iter() {
            res.push_item_unchecked(CigarItem::from_u32(val));
        }
        res
    }

    /// Returns true if the record has Hard clipping. Must not be empty.
    pub fn has_hard_clipping(&self) -> bool {
        self.tuples[0].op == Operation::Hard || self.tuples.last().unwrap().op == Operation::Hard
    }

    /// Replace hard clipping with soft.
    pub fn hard_to_soft(&mut self) {
        let first = &mut self.tuples[0];
        if first.op == Operation::Hard {
            first.op = Operation::Soft;
            self.qlen += first.len;
        }
        let last = self.tuples.last_mut().unwrap();
        if last.op == Operation::Hard {
            last.op = Operation::Soft;
            self.qlen += last.len;
        }
    }

    /// inverts reference and query (for example, `10M3I5M -> 10M3D15M`).
    pub fn invert(&self) -> Self {
        Self {
            tuples: self.tuples.iter().map(|item| CigarItem::new(item.op.invert(), item.len)).collect(),
            rlen: self.qlen,
            qlen: self.rlen,
        }
    }

    /// Length of the reference sequence.
    pub fn ref_len(&self) -> u32 {
        self.rlen
    }

    /// Length of the query sequence.
    pub fn query_len(&self) -> u32 {
        self.qlen
    }

    pub fn push_item_unchecked(&mut self, item: CigarItem) {
        if item.op.consumes_ref() {
            self.rlen += item.len;
        }
        if item.op.consumes_query() {
            self.qlen += item.len;
        }
        self.tuples.push(item);
    }

    /// Push a new `CigarItem`, does not merge with the latest entry.
    #[inline]
    pub fn push_unchecked(&mut self, op: Operation, len: u32) {
        self.push_item_unchecked(CigarItem::new(op, len));
    }

    /// Push a new entry, merges with the latest entry, if relevant.
    pub fn push_checked(&mut self, op: Operation, len: u32) {
        if op.consumes_ref() {
            self.rlen += len;
        }
        if op.consumes_query() {
            self.qlen += len;
        }
        match self.tuples.last_mut() {
            Some(item) if item.op == op => item.len += len,
            _ => self.tuples.push(CigarItem::new(op, len)),
        }
    }

    pub fn extend(&mut self, other: &Cigar) {
        if other.tuples.is_empty() {
            return;
        }
        let n = self.tuples.len();
        let oth_first = other.tuples[0];
        if n > 0 && self.tuples[n - 1].op == oth_first.op {
            self.tuples[n - 1].len += oth_first.len;
            self.tuples.extend_from_slice(&other.tuples[1..]);
        } else {
            self.tuples.extend_from_slice(&other.tuples);
        };
        self.rlen += other.rlen;
        self.qlen += other.qlen;
    }

    /// Returns iterator over CIGAR items.
    pub fn iter<'a>(&'a self) -> std::slice::Iter<'a, CigarItem> {
        self.tuples.iter()
    }

    /// Returns CIGAR length.
    #[inline]
    pub fn len(&self) -> usize {
        self.tuples.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.tuples.is_empty()
    }

    #[inline]
    pub fn shrink_to_fit(&mut self) {
        self.tuples.shrink_to_fit()
    }

    /// Returns true if the CIGAR does not contain M operation.
    pub fn is_extended(&self) -> bool {
        self.tuples.iter().all(|item| item.op != Operation::Match)
    }

    /// Infer extended CIGAR from a region with already known reference sequence
    /// (`ref_seq`, shifted by `ref_seq_shift`).
    ///
    /// If the record goes outside of the reference sequence, return None.
    pub fn infer_ext_cigar(rec: &record::Record, ref_seq: &[u8], ref_seq_shift: u32) -> Option<Cigar> {
        let ref_seq_end = ref_seq_shift + ref_seq.len() as u32;
        let ref_start = u32::try_from(rec.pos()).unwrap();
        let query_seq = rec.seq();
        let mut cigar = Cigar::new();

        for &val in rec.raw_cigar() {
            let item = CigarItem::from_u32(val);
            if item.op != Operation::Match {
                cigar.push_item_unchecked(item);
                continue;
            }

            if ref_start + cigar.rlen < ref_seq_shift || ref_start + cigar.rlen + item.len >= ref_seq_end {
                // Read goes beyond the boundary of the reference sequence, try to infer extended CIGAR using MD tag.
                return None;
            }

            let ref_shift = (ref_start + cigar.rlen - ref_seq_shift) as usize;
            let query_shift = cigar.qlen as usize;
            let mut curr_len = 0;
            let mut curr_equal = true;
            for i in 0..item.len as usize {
                let now_equal = ref_seq[ref_shift + i] == query_seq[query_shift + i];
                if now_equal != curr_equal && curr_len > 0 {
                    cigar.push_unchecked(if curr_equal { Operation::Equal } else { Operation::Diff }, curr_len);
                    curr_len = 0;
                }
                curr_equal = now_equal;
                curr_len += 1;
            }
            cigar.push_unchecked(if curr_equal { Operation::Equal } else { Operation::Diff }, curr_len);
        }
        let in_query_len = query_seq.len() as u32;
        // If the query sequence is missing, in_query_len will be 0.
        if in_query_len > 0 {
            assert_eq!(cigar.qlen, in_query_len,
                "Failed to convert CIGAR for read {:?}", String::from_utf8_lossy(rec.qname()));
        }
        Some(cigar)
    }

    /// Create an extended CIGAR from short CIGAR and MD string. Returns Cigar.
    /// Second argument: either `()`, or `&mut Vec<u8>`.
    pub fn infer_ext_cigar_md<V>(rec: &record::Record, mut ref_seq: V) -> Cigar
    where V: VecOrNone<u8>
    {
        let md_str = match rec.aux(b"MD") {
            Ok(record::Aux::String(s)) => s,
            _ => panic!("Cannot create extended CIGAR: record {} has no MD tag",
                String::from_utf8_lossy(rec.qname())),
        };
        let raw_md = md_str.as_bytes();
        let mut data = ExtCigarData {
            new_cigar: Cigar::new(),
            query_seq: rec.seq(),
            ref_seq: &mut ref_seq,
            md_str, raw_md,
            md_entries: parse_md(raw_md),
            md_ix: 0,
            md_shift: 0,
        };

        for &val in rec.raw_cigar() {
            data.process_op(CigarItem::from_u32(val));
        }
        debug_assert_eq!(data.md_ix, data.md_entries.len(), "Failed to parse MD tag {:?}", data.md_str);
        if let Some(ref_len) = data.ref_seq.try_len() {
            assert_eq!(ref_len as u32, data.new_cigar.ref_len(), "Failed to parse MD tag {:?}", md_str);
        }
        let in_query_len = data.query_seq.len() as u32;
        // If the query sequence is missing, in_query_len will be 0.
        if in_query_len > 0 {
            assert_eq!(data.new_cigar.query_len(), in_query_len, "Failed to parse MD tag {:?}", md_str);
        }
        data.new_cigar.shrink_to_fit();
        data.new_cigar
    }

    /// Returns soft clipping on the left and right.
    pub fn soft_clipping(&self) -> (u32, u32) {
        assert!(!self.is_empty(), "Cannot calculate soft clipping on an empty CIGAR!");
        let first = &self.tuples[0];
        let last = &self.tuples[self.len() - 1];
        (
            if first.op == Operation::Soft { first.len() } else { 0 },
            if last.op == Operation::Soft { last.len() } else { 0 },
        )
    }

    /// Calculates "true" clipping: total query length before the first match (=), and after the last match.
    pub fn true_clipping(&self) -> (u32, u32) {
        let mut left = 0;
        for item in self.tuples.iter() {
            if item.op == Operation::Equal {
                break;
            } else if item.op.consumes_query() {
                left += item.len;
            }
        }
        assert!(left < self.qlen, "Error in true_clipping({}): there are no matches in the CIGAR.", self);

        let mut right = 0;
        for item in self.tuples.iter().rev() {
            if item.op == Operation::Equal {
                break;
            } else if item.op.consumes_query() {
                right += item.len;
            }
        }
        (left, right)
    }

    /// Returns the number of matches and the total size of the alignment.
    pub fn frac_matches(&self) -> (u32, u32) {
        let mut nmatches = 0;
        let mut total_size = 0;
        for item in self.iter() {
            total_size += item.len;
            nmatches += if item.op == Operation::Equal { item.len } else { 0 };
        }
        (nmatches, total_size)
    }
}

impl PartialEq for Cigar {
    fn eq(&self, oth: &Self) -> bool {
        self.tuples.eq(&oth.tuples)
    }
}

impl Eq for Cigar {}

impl fmt::Debug for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_empty() {
            write!(f, "*")
        } else {
            let mut first = true;
            for tup in self.iter() {
                if first {
                    first = false;
                } else {
                    f.write_char(' ')?;
                }
                write!(f, "{}", tup)?;
            }
            Ok(())
        }
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_empty() {
            write!(f, "*")
        } else {
            for tup in self.iter() {
                write!(f, "{}", tup)?;
            }
            Ok(())
        }
    }
}

impl Index<usize> for Cigar {
    type Output = CigarItem;

    fn index(&self, i: usize) -> &CigarItem {
        &self.tuples.index(i)
    }
}

/// Entry in a MD tag.
#[derive(Clone, Debug)]
enum MdEntry {
    /// Match(length): reference matches with the query for `length` nucleotides,
    Match(u32),
    /// Mismatch(start, end): there is a mismatch, see reference sequence in md_tag[start..end],
    Mismatch(u32, u32),
    /// Deletion(start, end): there is a deletion, see reference sequence in md_tag[start..end].
    Deletion(u32, u32),
}

/// Parse MD tag, returns Vector of `MdEntry`.
fn parse_md(md: &[u8]) -> Vec<MdEntry> {
    const MATCH: u8 = 0;
    const MISM: u8 = 1;
    const DEL: u8 = 2;
    let mut status = MATCH;
    // Value is used as length for match, start index for mismatch and deletion.
    let mut value: u32 = 0;
    let mut md_entries = Vec::new();
    for (i, &byte) in md.iter().enumerate() {
        let i = i as u32;
        if b'0' <= byte && byte <= b'9' {
            if status >= MISM {
                debug_assert!(value < i);
                md_entries.push(
                    if status == MISM { MdEntry::Mismatch(value, i) } else { MdEntry::Deletion(value, i) });
                status = MATCH;
                value = u32::from(byte - b'0');
            } else {
                value = 10 * value + u32::from(byte - b'0');
            }
        }

        else if byte == b'^' {
            debug_assert!(status == MATCH);
            if value > 0 {
                md_entries.push(MdEntry::Match(value));
            }
            status = DEL;
            value = i + 1;
        }

        else if status == MATCH {
            if value > 0 {
                md_entries.push(MdEntry::Match(value));
            }
            status = MISM;
            value = i;
        }
    }
    debug_assert!(status == MATCH);
    if value > 0 {
        md_entries.push(MdEntry::Match(value));
    }
    md_entries
}

/// Temporary data structure, used to reduce the number of function parameters.
struct ExtCigarData<'a, V: VecOrNone<u8>> {
    query_seq: record::Seq<'a>,
    new_cigar: Cigar,
    ref_seq: &'a mut V,
    md_str: &'a str,
    raw_md: &'a [u8],
    md_entries: Vec<MdEntry>,
    md_ix: usize,
    md_shift: u32,
}

impl<'a, V: VecOrNone<u8>> ExtCigarData<'a, V> {
    fn process_op(&mut self, tup: CigarItem) {
        if tup.op.consumes_both() {
            let mut op_len = tup.len;
            while op_len != 0 {
                match &self.md_entries[self.md_ix] {
                    MdEntry::Match(length) => self.process_match_match(&mut op_len, *length),
                    MdEntry::Mismatch(start, end) => self.process_match_mismatch(&mut op_len, *start, *end),
                    MdEntry::Deletion(_, _) => unreachable!(
                        "Failed to parse MD tag {:?}: operation M coincides with deletion", self.md_str),
                }
            }

        } else if tup.op.consumes_query() {
            self.new_cigar.push_item_unchecked(tup);

        } else {
            if let MdEntry::Deletion(start, end) = self.md_entries[self.md_ix] {
                debug_assert_eq!(self.md_shift, 0);
                debug_assert_eq!(end - start, tup.len);
                if !V::IS_SINK {
                    for i in start..end {
                        self.ref_seq.push(self.raw_md[i as usize]);
                    }
                }
                self.new_cigar.push_item_unchecked(tup);
                self.md_ix += 1;
            }
        }
    }

    fn process_match_match(&mut self, op_len: &mut u32, md_len: u32) {
        let rem_md_len = md_len - self.md_shift;
        let pos_inc = min(*op_len, rem_md_len);
        if !V::IS_SINK {
            for qpos in self.new_cigar.qlen..self.new_cigar.qlen + pos_inc {
                self.ref_seq.push(self.query_seq[qpos as usize]);
            }
        }

        if rem_md_len > pos_inc {
            self.md_shift += pos_inc;
        } else {
            self.md_ix += 1;
            self.md_shift = 0;
        }
        self.new_cigar.push_checked(Operation::Equal, pos_inc);
        *op_len -= pos_inc;
    }

    fn process_match_mismatch(&mut self, op_len: &mut u32, md_start: u32, md_end: u32) {
        let rem_md_len = md_end - md_start - self.md_shift;
        let pos_inc = min(*op_len, rem_md_len);
        if !V::IS_SINK {
            for pos in md_start + self.md_shift..md_start + self.md_shift + pos_inc {
                self.ref_seq.push(self.raw_md[pos as usize]);
            }
        }

        if rem_md_len > pos_inc {
            self.md_shift += pos_inc;
        } else {
            self.md_ix += 1;
            self.md_shift = 0;
        }
        self.new_cigar.push_checked(Operation::Diff, pos_inc);
        *op_len -= pos_inc;
    }
}

/// Extract Soft/Hard clipping size from a raw CIGAR.
/// This function only considers left-most and right-most operations, so 10H10S10M would produce 10, not 20.
fn raw_clipping(raw_cigar: &[u32]) -> u32 {
    let n = raw_cigar.len();
    let first = CigarItem::from_u32(raw_cigar[0]);
    let mut clipping = u32::from(!first.op.consumes_ref()) * first.len;
    if n > 0 {
        let last = CigarItem::from_u32(raw_cigar[raw_cigar.len() - 1]);
        clipping += u32::from(!last.op.consumes_ref()) * last.len;
    }
    clipping
}

/// Returns Soft/Hard clipping size divided by the read length.
/// This function only considers left-most and right-most CIGAR operations,
/// and does not account for Hard clipping operations when calculating read length.
/// Therefore `10H10S10M10S` would produce `20 / 30`, and not `30 / 40`.
pub fn clipping_rate(record: &record::Record) -> f64 {
    let clipping = raw_clipping(record.raw_cigar());
    if clipping == 0 {
        0.0
    } else {
        f64::from(clipping) / record.seq_len() as f64
    }
}
