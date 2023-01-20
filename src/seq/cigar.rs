//! Various functions and structures related to sequence alignment.

use std::fmt::{self, Write};
use std::cmp::min;
use std::ops::Index;

use htslib::bam::record;

/// Subset of CIGAR operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operation {
    Match,
    Equal,
    Diff,
    Soft,
    Ins,
    Del,
}

impl Operation {
    /// Does the Cigar operation consume reference sequence?
    pub const fn consumes_ref(self) -> bool {
        match self {
            Operation::Match | Operation::Equal | Operation::Diff | Operation::Del => true,
            _ => false,
        }
    }

    /// Does the Cigar operation consume query sequence?
    pub const fn consumes_query(self) -> bool {
        match self {
            Operation::Match | Operation::Equal | Operation::Diff | Operation::Ins | Operation::Soft => true,
            _ => false,
        }
    }

    /// Does the Cigar operation consume both reference and query sequences?
    pub const fn consumes_both(self) -> bool {
        match self {
            Operation::Match | Operation::Equal | Operation::Diff => true,
            _ => false,
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
            // 5 => Hard,
            // 6 => Padding,
            7 => Operation::Equal,
            8 => Operation::Diff,
            _ => panic!("Unexpected cigar operation"),
        }
    }
}

impl fmt::Display for Operation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_char(self.to_char())
    }
}

/// Tuple operation + length.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CigarItem {
    op: Operation,
    len: u32,
}

impl CigarItem {
    /// Creates a new tuple.
    pub fn new(op: Operation, len: u32) -> Self {
        CigarItem { op, len }
    }

    /// Get operation.
    pub fn operation(&self) -> Operation {
        self.op
    }

    /// Get length.
    pub fn len(&self) -> u32 {
        self.len
    }

    /// Creates a new tuple from a raw u32 value.
    pub fn from_u32(val: u32) -> Self {
        Self {
            op: Operation::from_u32(val & 0b1111),
            len: val >> 4,
        }
    }

    /// Creates a new tuple from HTSLIB format.
    pub fn from_htslib(op: record::Cigar) -> Self {
        match op {
            record::Cigar::Match(len) => CigarItem::new(Operation::Match, len),
            record::Cigar::Ins(len) => CigarItem::new(Operation::Ins, len),
            record::Cigar::Del(len) => CigarItem::new(Operation::Del, len),
            record::Cigar::SoftClip(len) => CigarItem::new(Operation::Soft, len),
            record::Cigar::Equal(len) => CigarItem::new(Operation::Equal, len),
            record::Cigar::Diff(len) => CigarItem::new(Operation::Diff, len),
            _ => panic!("Unexpected CIGAR operation {}", op.char()),
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
    #[inline]
    pub fn new() -> Cigar {
        Cigar::default()
    }

    pub fn from_raw(raw_cigar: &[u32]) -> Cigar {
        let mut res = Cigar::new();
        for &val in raw_cigar.iter() {
            res.push(CigarItem::from_u32(val));
        }
        res
    }

    /// Length of the reference sequence.
    pub fn ref_len(&self) -> u32 {
        self.rlen
    }

    /// Length of the query sequence.
    pub fn query_len(&self) -> u32 {
        self.qlen
    }

    /// Push a new `CigarItem`, does not merge with the latest entry.
    pub fn push(&mut self, item: CigarItem) {
        if item.op.consumes_ref() {
            self.rlen += item.len;
        }
        if item.op.consumes_query() {
            self.qlen += item.len;
        }
        self.tuples.push(item);
    }

    /// Push a new entry, merges with the latest entry, if relevant.
    pub fn checked_push(&mut self, op: Operation, len: u32) {
        let n = self.tuples.len();
        if n == 0 || self.tuples[n - 1].op != op {
            self.push(CigarItem::new(op, len));
        } else {
            if op.consumes_ref() {
                self.rlen += len;
            }
            if op.consumes_query() {
                self.qlen += len;
            }
            self.tuples[n - 1].len += len;
        }
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

    /// Create an extended CIGAR from short CIGAR and MD string. Returns Cigar.
    pub fn infer_ext_cigar(rec: &record::Record) -> Cigar {
        let md_str = if let Ok(record::Aux::String(s)) = rec.aux(b"MD") {
            s
        } else {
            panic!("Cannot create extended CIGAR: record {} either has no MD tag, or MD tag has incorrect type",
                String::from_utf8_lossy(rec.qname()))
        };
        ExtCigarData::<()>::process(rec.seq(), rec.raw_cigar(), md_str).0
    }

    /// Create an extended CIGAR from short CIGAR and MD string, as well as the reference sequence.
    /// Returns pair (Cigar, Vec<u8>).
    pub fn infer_ext_cigar_seq(rec: &record::Record) -> (Cigar, Vec<u8>) {
        let md_str = if let Ok(record::Aux::String(s)) = rec.aux(b"MD") {
            s
        } else {
            panic!("Cannot create extended CIGAR: record {} either has no MD tag, or MD tag has incorrect type",
                String::from_utf8_lossy(rec.qname()))
        };
        ExtCigarData::process(rec.seq(), rec.raw_cigar(), md_str)
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

/// Trait that stores either Vec<u8>, or ().
pub trait SeqOrNone : std::fmt::Debug {
    /// Creates a new object.
    fn new() -> Self;

    /// Returns true if it is possible to push nucletides in the object (false for ()).
    fn push_possible() -> bool;

    /// Pushes a new nucleotide.
    fn push(&mut self, nt: u8);

    /// Complete the object creation (for example, shrinks the vector).
    fn complete(&mut self);

    /// Returns length, if available.
    fn len(&self) -> Option<usize>;
}

impl SeqOrNone for () {
    fn new() -> Self {}
    fn push_possible() -> bool { false }
    fn push(&mut self, _nt: u8) {}
    fn complete(&mut self) {}
    fn len(&self) -> Option<usize> { None }
}

impl SeqOrNone for Vec<u8> {
    fn new() -> Self {
        Vec::new()
    }

    fn push_possible() -> bool { true }

    fn push(&mut self, nt: u8) {
        self.push(nt);
    }

    fn complete(&mut self) {
        self.shrink_to_fit();
    }

    fn len(&self) -> Option<usize> {
        Some(self.len())
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
struct ExtCigarData<'a, S: SeqOrNone> {
    query_seq: record::Seq<'a>,
    new_cigar: Cigar,
    ref_seq: S,

    md_str: &'a str,
    raw_md: &'a [u8],
    md_entries: Vec<MdEntry>,
    md_ix: usize,
    md_shift: u32,
}

impl<'a, S: SeqOrNone> ExtCigarData<'a, S> {
    fn process(query_seq: record::Seq<'a>, raw_cigar: &[u32], md_str: &'a str) -> (Cigar, S) {
        let raw_md = md_str.as_bytes();
        let mut data = Self {
            query_seq,
            new_cigar: Cigar::new(),
            ref_seq: S::new(),
            md_str,
            raw_md,
            md_entries: parse_md(raw_md),
            md_ix: 0,
            md_shift: 0,
        };

        for &val in raw_cigar {
            data.process_op(CigarItem::from_u32(val));
        }
        debug_assert_eq!(data.md_ix, data.md_entries.len(), "Failed to parse MD tag \"{}\"", data.md_str);

        if let Some(ref_len) = data.ref_seq.len() {
            assert_eq!(ref_len as u32, data.new_cigar.ref_len(), "Failed to parse MD tag \"{}\"", md_str);
        }
        let in_query_len = query_seq.len() as u32;
        // If the query sequence is missing, in_query_len will be 0.
        if in_query_len > 0 {
            assert_eq!(data.new_cigar.query_len(), in_query_len, "Failed to parse MD tag \"{}\"", md_str);
        }

        data.new_cigar.shrink_to_fit();
        data.ref_seq.complete();
        (data.new_cigar, data.ref_seq)
    }

    fn process_op(&mut self, tup: CigarItem) {
        if tup.op.consumes_both() {
            let mut op_len = tup.len;
            while op_len != 0 {
                match &self.md_entries[self.md_ix] {
                    MdEntry::Match(length) => self.process_match_match(&mut op_len, *length),
                    MdEntry::Mismatch(start, end) => self.process_match_mismatch(&mut op_len, *start, *end),
                    MdEntry::Deletion(_, _) => unreachable!(
                        "Failed to parse MD tag \"{}\": operation M coincides with deletion", self.md_str),
                }
            }

        } else if tup.op.consumes_query() {
            self.new_cigar.push(CigarItem::new(tup.op, tup.len));

        } else {
            if let MdEntry::Deletion(start, end) = self.md_entries[self.md_ix] {
                debug_assert_eq!(self.md_shift, 0);
                debug_assert_eq!(end - start, tup.len);
                if S::push_possible() {
                    for i in start..end {
                        self.ref_seq.push(self.raw_md[i as usize]);
                    }
                }
                self.new_cigar.push(CigarItem::new(tup.op, tup.len));
                self.md_ix += 1;
            }
        }
    }

    fn process_match_match(&mut self, op_len: &mut u32, md_len: u32) {
        let rem_md_len = md_len - self.md_shift;
        let pos_inc = min(*op_len, rem_md_len);
        if S::push_possible() {
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
        self.new_cigar.checked_push(Operation::Equal, pos_inc);
        *op_len -= pos_inc;
    }

    fn process_match_mismatch(&mut self, op_len: &mut u32, md_start: u32, md_end: u32) {
        let rem_md_len = md_end - md_start - self.md_shift;
        let pos_inc = min(*op_len, rem_md_len);
        if S::push_possible() {
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
        self.new_cigar.checked_push(Operation::Diff, pos_inc);
        *op_len -= pos_inc;
    }
}
