//! Various functions and structures related to sequence alignment.

use std::{
    fmt::{self, Write},
    cmp::{min, max},
    ops::Index,
};
use htslib::bam::record;
use crate::{
    err::error,
    ext::vec::VecOrNone,
    seq::wfa::{self, Aligner},
    algo::bisect,
    math::{RoundDiv, bool_mask, ifelse0},
};

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Consumes {
    Query,
    Ref,
    Both,
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

    /// Returns tuple `(consumes_query(), consumes_ref())` in one go.
    #[inline]
    pub const fn consumes_query_ref(self) -> (bool, bool) {
        match self {
            Operation::Match | Operation::Equal | Operation::Diff => (true, true),
            Operation::Soft | Operation::Ins => (true, false),
            Operation::Del => (false, true),
            Operation::Hard => (false, false),
        }
    }

    /// Returns whether this operation consumes query, reference or both. Panics on hard clipping.
    #[inline]
    pub const fn consumes(self) -> Consumes {
        match self {
            Operation::Match | Operation::Equal | Operation::Diff => Consumes::Both,
            Operation::Soft | Operation::Ins => Consumes::Query,
            Operation::Del => Consumes::Ref,
            Operation::Hard => panic!("Unexpected hard clipping"),
        }
    }

    pub fn from_char(c: u8) -> crate::Result<Self> {
        match c {
            b'M' => Ok(Self::Match),
            b'=' => Ok(Self::Equal),
            b'X' => Ok(Self::Diff),
            b'I' => Ok(Self::Ins),
            b'D' => Ok(Self::Del),
            b'S' => Ok(Self::Soft),
            b'H' => Ok(Self::Hard),
            b'N' | b'P' => Err(error!(RuntimeError, "CIGAR operations N and P are not supported")),
            _ => Err(error!(InvalidData, "Unexpected CIGAR operation {}", c as char)),
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

    /// Returns true if insertion and deletion go right after another.
    pub const fn gap_conflict(self, next: Operation) -> bool {
        match (self, next) {
            (Operation::Ins, Operation::Del) | (Operation::Del, Operation::Ins) => true,
            _ => false
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

    pub fn clear(&mut self) {
        self.tuples.clear();
        self.rlen = 0;
        self.qlen = 0;
    }

    /// Parses string representation of CIGAR.
    pub fn from_str(s: &[u8]) -> crate::Result<Self> {
        let mut cigar = Self::new();
        let mut len = 0;
        for &c in s {
            match c {
                b'0' ..= b'9' => len = 10 * len + u32::from(c - b'0'),
                _ => {
                    cigar.push_unchecked(Operation::from_char(c)?, len);
                    len = 0;
                }
            }
        }
        Ok(cigar)
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
        let Some(first) = self.tuples.first_mut() else { return };
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
        let (cons_query, cons_ref) = item.op.consumes_query_ref();
        self.qlen += ifelse0(cons_query, item.len);
        self.rlen += ifelse0(cons_ref, item.len);
        self.tuples.push(item);
    }

    /// Push a new `CigarItem`, does not merge with the latest entry.
    #[inline]
    pub fn push_unchecked(&mut self, op: Operation, len: u32) {
        self.push_item_unchecked(CigarItem::new(op, len));
    }

    /// Push a new entry, merges with the latest entry, if relevant.
    pub fn push_checked(&mut self, op: Operation, len: u32) {
        let (cons_query, cons_ref) = op.consumes_query_ref();
        self.qlen += ifelse0(cons_query, len);
        self.rlen += ifelse0(cons_ref, len);
        match self.tuples.last_mut() {
            Some(item) if item.op == op => item.len += len,
            _ => self.tuples.push(CigarItem::new(op, len)),
        }
    }

    /// Removes last entry in this CIGAR.
    pub fn pop(&mut self) -> Option<CigarItem> {
        let item = self.tuples.pop()?;
        let (cons_query, cons_ref) = item.op.consumes_query_ref();
        self.qlen -= ifelse0(cons_query, item.len);
        self.rlen -= ifelse0(cons_ref, item.len);
        Some(item)
    }

    /// Removes last entry in the CIGAR if predicate returns true.
    pub fn pop_if(&mut self, predicate: impl FnOnce(&CigarItem) -> bool) -> Option<CigarItem> {
        let item = self.tuples.pop_if(|item| predicate(item))?;
        let (cons_query, cons_ref) = item.op.consumes_query_ref();
        self.qlen -= ifelse0(cons_query, item.len);
        self.rlen -= ifelse0(cons_ref, item.len);
        Some(item)
    }

    #[inline(always)]
    pub fn last(&self) -> Option<&CigarItem> {
        self.tuples.last()
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
    pub fn infer_ext_cigar_from_md<V>(rec: &record::Record, mut ref_seq: V) -> Cigar
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

    #[inline(always)]
    pub fn full_sequence_match(&self) -> bool {
        self.tuples.len() == 1 && self.tuples[0].op == Operation::Equal
    }

    /// Returns soft clipping on the left and right.
    pub fn soft_clipping(&self) -> (u32, u32) {
        assert!(!self.is_empty(), "Cannot calculate soft clipping on an empty CIGAR!");
        let first = self.tuples.first().unwrap();
        let last = self.tuples.last().unwrap();
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

    /// Replace Insertion at the boundary into soft clipping.
    #[inline]
    pub fn boundary_ins_to_soft(&mut self) {
        let n = self.len();
        assert!(n > 0);
        let first = unsafe { self.tuples.get_unchecked_mut(0) };
        first.op = if first.op == Operation::Ins { Operation::Soft } else { first.op };
        let last = unsafe { self.tuples.get_unchecked_mut(n - 1) };
        last.op = if last.op == Operation::Ins { Operation::Soft } else { last.op };
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

    /// Compares `=` and `X` entries between query and reference sequences
    /// and returns true if there are no discrepancies.
    #[allow(unused)]
    pub fn validate(&self, query_seq: &[u8], ref_seq: &[u8]) -> bool {
        if query_seq.len() as u32 != self.qlen || ref_seq.len() as u32 != self.rlen {
            return false;
        }
        let mut qpos = 0;
        let mut rpos = 0;
        for item in &self.tuples {
            let (cons_query, cons_ref) = item.op.consumes_query_ref();
            if cons_query && cons_ref && item.op != Operation::Match {
                let subseq1 = &query_seq[qpos as usize..(qpos + item.len) as usize];
                let subseq2 = &ref_seq[rpos as usize..(rpos + item.len) as usize];
                let must_match = item.op == Operation::Equal;
                if itertools::izip!(subseq1, subseq2).any(|(&c1, &c2)| (c1 == c2) != must_match) {
                    return false;
                }
            }
            qpos += ifelse0(cons_query, item.len);
            rpos += ifelse0(cons_ref, item.len);
        }
        true
    }

    /// Returns maximum edit distance, observed across any window of the given size.
    /// Returns full alignment edit distance if window is larger than the alignment size.
    pub fn max_local_edit(&self, window: u32) -> u32 {
        // l_ = window start, r_ = window end.
        // For both positions, store `rem` - remainder of the current CIGAR item and
        // `mask`, which is 0 for "=" operation and 11..111b for everything else.
        let mut r_iter = self.tuples.iter();
        let mut r_rem;
        let mut r_mask;
        let mut window_rem = window;
        let mut edit = 0;
        loop {
            match r_iter.next() {
                Some(item) if item.len > window_rem => {
                    r_rem = item.len - window_rem;
                    r_mask = bool_mask(item.op != Operation::Equal);
                    edit += r_mask & window_rem;
                    break;
                }
                Some(item) => {
                    edit += bool_mask(item.op != Operation::Equal) & item.len;
                    window_rem -= item.len;
                }
                // Reached the end of alignment before we reached the end of window.
                None => return edit,
            }
        }
        let mut max_edit = edit;

        let mut l_iter = self.tuples.iter();
        let &CigarItem { op: tmp_op, len: mut l_rem } = l_iter.next().expect("CIGAR must be non-empty");
        let mut l_mask = bool_mask(tmp_op != Operation::Equal);

        // Simultaneously go through the left and right iterator, subtract left edit, add right edit.
        loop {
            let shift = min(l_rem, r_rem);
            edit = edit + (r_mask & shift) - (l_mask & shift);
            max_edit = max(max_edit, edit);

            if shift == r_rem {
                let Some(item) = r_iter.next() else { break };
                r_rem = item.len;
                r_mask = bool_mask(item.op != Operation::Equal);
            } else {
                r_rem -= shift;
            }
            if shift == l_rem {
                let item = l_iter.next().expect("Left iterator could not overtake the right one");
                l_rem = item.len;
                l_mask = bool_mask(item.op != Operation::Equal);
            } else {
                l_rem -= shift;
            }
        }
        max_edit
    }

    /// Finds locally similar moving windows (size = `window`, step = `step`)
    /// whose edit distances do not exceed `max_edit`.
    ///
    /// Appends `indices` with the corresponding window indices (vector is not cleared).
    pub fn locally_similar<const IN_QUERY: bool>(
        &self,
        window: u32,
        step: u32,
        max_edit: u32,
        indices: &mut Vec<usize>,
    ) {
        indices.clear();
        // ...1 = variables at the window start, ...2 = variables at the window end.
        let mut pos2 = 0;
        let mut iter2 = self.tuples.iter();
        let mut edit = 0;
        // Store `rem` - remainder of the current CIGAR item and last operation `op`.
        let (mut rem2, mut op2) = loop {
            if let Some(item) = iter2.next() {
                let moves_pos = if IN_QUERY { item.op.consumes_query() } else { item.op.consumes_ref() };
                if moves_pos {
                    let window_rem = window - pos2;
                    let shift = min(item.len, window_rem);
                    edit += ifelse0(item.op != Operation::Equal, shift);
                    pos2 += shift;
                    if item.len > window_rem {
                        break (item.len - window_rem, item.op);
                    }
                } else {
                    edit += item.len;
                }
            } else {
                // Reached the end of alignment.
                if edit <= max_edit {
                    indices.push(0);
                }
                return;
            }
        };

        let mut iter1 = self.tuples.iter();
        let &CigarItem { op: mut op1, len: mut rem1 } = iter1.next().expect("CIGAR must be non-empty");
        let mut pos1 = 0;
        // Simultaneously go through the left and right iterator, subtract left edit, add right edit.
        loop {
            let shift = min(rem1, rem2);
            let moves_pos1 = if IN_QUERY { op1.consumes_query() } else { op1.consumes_ref() };
            let moves_pos2 = if IN_QUERY { op2.consumes_query() } else { op2.consumes_ref() };
            // Move either both positions, or only update CIGAR entry that does not move position.
            let upd1 = (moves_pos1 == moves_pos2) || !moves_pos1;
            let upd2 = (moves_pos1 == moves_pos2) || !moves_pos2;

            if moves_pos1 && moves_pos2 {
                let mask1 = bool_mask(op1 != Operation::Equal);
                let mask2 = bool_mask(op2 != Operation::Equal);
                // Smallest i * step >= pos1.
                let next_saved_pos = pos1 + i32::rem_euclid(-(pos1 as i32), step as i32) as u32;
                for pos in (next_saved_pos..pos1 + shift).step_by(step as usize) {
                    let curr_shift = pos - pos1;
                    if edit + (mask2 & curr_shift) - (mask1 & curr_shift) <= max_edit {
                        indices.push((pos / step) as usize);
                    }
                }
                pos1 += shift;
                pos2 += shift;
                edit = edit + (mask2 & shift) - (mask1 & shift);
            } else {
                // Do not move left/right positions, only update these CIGAR entries that do not move positions.
                edit = edit + (bool_mask(upd2) & shift) - (bool_mask(upd1) & shift);
            }

            if upd2 {
                if shift == rem2 {
                    let Some(item) = iter2.next() else { break };
                    rem2 = item.len;
                    op2 = item.op;
                } else {
                    rem2 -= shift;
                }
            }
            if upd1 {
                if shift == rem1 {
                    let item = iter1.next().expect("Left iterator could not overtake the right one");
                    rem1 = item.len;
                    op1 = item.op;
                } else {
                    rem1 -= shift;
                }
            }
        }
        if edit <= max_edit {
            indices.push(pos1.fast_ceil_div(step) as usize);
        }
        assert!(pos1 + window == pos2);
        assert!(pos2 == if IN_QUERY { self.qlen } else { self.rlen });
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
        self.tuples.index(i)
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
    let mut clipping = ifelse0(!first.op.consumes_ref(), first.len);
    if n > 0 {
        let last = CigarItem::from_u32(raw_cigar[raw_cigar.len() - 1]);
        clipping += ifelse0(!last.op.consumes_ref(), last.len);
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

/// To speed up binary search, store cigar index for each position = i << STEP_PWR. With pwr = 8, step is 256 bp.
pub const SPARSE_STEP_PWR: u32 = 8;
const SPARSE_MASK: u32 = (1 << SPARSE_STEP_PWR) - 1;

#[inline(always)]
fn update_sparse_index(
    sparse_index: &mut Vec<(u32, u32)>,
    cigar_ix: u32,
    len: u32,
    pos1: u32,
    pos2: u32,
    consumes_other: bool,
) {
    for i in sparse_index.len() as u32..=((pos1 + len - 1) >> SPARSE_STEP_PWR) {
        let sparse_pos1 = i << SPARSE_STEP_PWR;
        let sparse_pos2 = pos2 + u32::from(consumes_other) * (sparse_pos1 - pos1);
        sparse_index.push((cigar_ix, sparse_pos2));
    }
}

/// Search direction, either from query to reference, or reference to query.
pub trait CigarDirection : Copy + fmt::Debug {
    /// In an array of two values [qpos, rpos], returns index of the qpos.
    /// In q_to_r direction, this is 0, otherwise 1.
    fn query_ix(self) -> usize;

    /// In an array of two values [qpos, rpos], returns index of the rpos.
    /// In q_to_r direction, this is 1, otherwise 0.
    fn ref_ix(self) -> usize;

    /// Returns query position out of two positions.
    /// In q_to_r direction, returns the first value, otherwise the second.
    fn query_pos(self, pos: [u32; 2]) -> u32;

    /// Returns reference position out of two positions.
    /// In q_to_r direction, returns the second value, otherwise the first.
    fn ref_pos(self, pos: [u32; 2]) -> u32;

    /// In q_to_r direction, returns operation itself, otherwise its inversion.
    fn operation(self, operation: Operation) -> Operation;

    /// False for q_to_r, true otherwise.
    fn invert(self) -> bool;

    #[inline(always)]
    fn query_len(self, cigar: &Cigar) -> u32 {
        if self.invert() { cigar.rlen } else { cigar.qlen }
    }

    #[inline(always)]
    fn ref_len(self, cigar: &Cigar) -> u32 {
        if self.invert() { cigar.qlen } else { cigar.rlen }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct QueryToRef;

impl CigarDirection for QueryToRef {
    #[inline(always)]
    fn query_ix(self) -> usize { 0 }

    #[inline(always)]
    fn ref_ix(self) -> usize { 1 }

    #[inline(always)]
    fn query_pos(self, [qpos, _]: [u32; 2]) -> u32 { qpos }

    #[inline(always)]
    fn ref_pos(self, [_, rpos]: [u32; 2]) -> u32 { rpos }

    #[inline(always)]
    fn operation(self, operation: Operation) -> Operation { operation }

    #[inline(always)]
    fn invert(self) -> bool { false }
}

#[derive(Clone, Copy, Debug)]
pub struct RefToQuery;

impl CigarDirection for RefToQuery {
    #[inline(always)]
    fn query_ix(self) -> usize { 1 }

    #[inline(always)]
    fn ref_ix(self) -> usize { 0 }

    #[inline(always)]
    fn query_pos(self, [_, rpos]: [u32; 2]) -> u32 { rpos }

    #[inline(always)]
    fn ref_pos(self, [qpos, _]: [u32; 2]) -> u32 { qpos }

    #[inline(always)]
    fn operation(self, operation: Operation) -> Operation { operation.invert() }

    #[inline(always)]
    fn invert(self) -> bool { true}
}

/// Structure, returned by `CigarIndex::find_approx_position`.
#[derive(Clone, Copy, Debug)]
pub struct ApproxPosition {
    pub min_cigar_ix: u32,
    pub max_cigar_ix: u32,
    pub approx_pos: u32,
}

/// Structure, returned by `CigarIndex::find_cigar_offset`.
#[derive(Clone, Copy, Debug, Default)]
pub struct CigarOffset {
    pub cigar_ix: usize,
    pub qpos_at_ix: u32,
    pub rpos_at_ix: u32,
}

/// Extension over CIGAR that answers queries "convert range start..end to another sequence" in log(n) time,
/// where n is the number of operations in this CIGAR.
#[derive(Clone)]
pub struct CigarIndex {
    /// Length = cigar.len()
    /// For each iteration, stores query and reference position at the start of the element.
    positions: Vec<[u32; 2]>,
    /// Outer index: 0 = query -> ref, 1 = ref -> query.
    /// Inner vector: for each position = i * STEP store cigar index and corresponding position in the other sequence.
    /// As the final value, each vector contains cigar length - 1 and length of the other sequence.
    sparse_index: [Vec<(u32, u32)>; 2],
}

impl CigarIndex {
    pub fn new(cigar: &Cigar) -> Self {
        let mut positions = Vec::with_capacity(cigar.len());
        let mut qpos = 0;
        let mut rpos = 0;
        let mut sparse_index_qr = Vec::with_capacity(2 + (cigar.qlen.saturating_sub(1) >> SPARSE_STEP_PWR) as usize);
        let mut sparse_index_rq = Vec::with_capacity(2 + (cigar.rlen.saturating_sub(1) >> SPARSE_STEP_PWR) as usize);
        for (cigar_ix, &CigarItem { len, op }) in cigar.iter().enumerate() {
            let cigar_ix = cigar_ix as u32;
            positions.push([qpos, rpos]);
            let (cons_query, cons_ref) = op.consumes_query_ref();
            let old_qpos = qpos;
            if cons_query {
                update_sparse_index(&mut sparse_index_qr, cigar_ix, len, qpos, rpos, cons_ref);
                qpos += len;
            }
            if cons_ref {
                update_sparse_index(&mut sparse_index_rq, cigar_ix, len, rpos, old_qpos, cons_query);
                rpos += len
            }
        }
        sparse_index_qr.push((cigar.len().strict_sub(1) as u32, cigar.rlen));
        sparse_index_rq.push((cigar.len().strict_sub(1) as u32, cigar.qlen));
        Self {
            positions,
            sparse_index: [sparse_index_qr, sparse_index_rq],
        }
    }

    /// Based on query position, returns approximate reference position.
    #[inline]
    pub fn find_approx_position(&self, qpos: u32, direction: impl CigarDirection) -> ApproxPosition {
        let sparse_index = &self.sparse_index[direction.query_ix()];
        let i = qpos >> SPARSE_STEP_PWR;
        assert!(i as usize + 1 <= sparse_index.len());
        let (min_cigar_ix, rpos1) = unsafe { *sparse_index.get_unchecked(i as usize) };
        let (max_cigar_ix, rpos2) = unsafe { *sparse_index.get_unchecked(i as usize + 1) };
        // Same as rpos1 + (qpos - qpos1) * (rpos2 - rpos1) / (qpos2 - qpos1), because:
        //     qpos & SPARSE_MASK = qpos - qpos1
        //     X >> SPARSE_STEP_PWR = X / (qpos2 - qpos1)
        let approx_pos = rpos1 + (((qpos & SPARSE_MASK) * (rpos2 - rpos1)) >> SPARSE_STEP_PWR);
        ApproxPosition { min_cigar_ix, max_cigar_ix, approx_pos }
    }

    /// Refines approximate position, finding specific CIGAR index, query and ref positions at that point.
    #[inline]
    pub fn find_cigar_offset(
        &self,
        qpos: u32,
        approx_pos: ApproxPosition,
        direction: impl CigarDirection,
    ) -> CigarOffset {
        let cigar_ix = if approx_pos.min_cigar_ix == approx_pos.max_cigar_ix {
            approx_pos.min_cigar_ix as usize
        } else {
            bisect::right_by_at(&self.positions, |&pos| direction.query_pos(pos).cmp(&qpos),
                approx_pos.min_cigar_ix as usize, approx_pos.max_cigar_ix as usize + 1).strict_sub(1)
        };
        let pos = self.positions[cigar_ix];
        CigarOffset {
            cigar_ix,
            qpos_at_ix: direction.query_pos(pos),
            rpos_at_ix: direction.ref_pos(pos),
        }
    }
}

impl Cigar {
    /// Realign potentially incorrect parts of the alignment.
    fn optimize(
        &mut self,
        ref_seq: &[u8],
        query_seq: &[u8],
        aligner: &Aligner,
        max_gap: impl wfa::Threshold,
        anchor_size: u32,
    ) {
        // Index after the last match (= or X), and the corresponding query and reference positions.
        let mut i = 0;
        let mut qpos1 = 0;
        let mut rpos1 = 0;
        // Will contain 0b01 if encountered deletion and 0b10 if encountered insertion.
        let mut flag = 0_u8;
        // Actual query and reference position at the current operation.
        let mut qpos2 = 0;
        let mut rpos2 = 0;
        let mut lazy_new_cigar = None::<Cigar>;

        for (j, &CigarItem { op, len }) in self.iter().enumerate() {
            let (cons_query, cons_ref) = op.consumes_query_ref();
            if cons_query && cons_ref && len >= anchor_size {
                let qshift = qpos2 - qpos1;
                let rshift = rpos2 - rpos1;
                // There was both a deletion and an insertion after the last anchor.
                if flag == 0b11 && !max_gap.under(qshift) && !max_gap.under(rshift) {
                    let new_cigar = lazy_new_cigar.get_or_insert_with(|| Cigar {
                        tuples: self.tuples[..i].to_vec(),
                        qlen: qpos1,
                        rlen: rpos1,
                    });
                    aligner.smart_align(ref_seq, rpos1, rpos2, query_seq, qpos1, qpos2, (), new_cigar);
                    i = j;
                }

                qpos2 += len;
                rpos2 += len;
                qpos1 = qpos2;
                rpos1 = rpos2;
                flag = 0;
                if let Some(new_cigar) = &mut lazy_new_cigar {
                    new_cigar.tuples.extend_from_slice(&self.tuples[i..j]);
                    new_cigar.push_checked(op, len);
                    new_cigar.qlen = qpos2;
                    new_cigar.rlen = rpos2;
                }
                i = j + 1;
            } else {
                qpos2 += ifelse0(cons_query, len);
                rpos2 += ifelse0(cons_ref, len);
                flag |= u8::from(!cons_query) | (u8::from(!cons_ref) << 1);
            }
        }

        let qshift = qpos2 - qpos1;
        let rshift = rpos2 - rpos1;
        if flag == 0b11 && !max_gap.under(qshift) && !max_gap.under(rshift) {
            // There were >=2 gap operations in a row at the end of the cigar.
            let new_cigar = lazy_new_cigar.get_or_insert_with(|| Cigar {
                tuples: self.tuples[..i].to_vec(),
                qlen: qpos1,
                rlen: rpos1,
            });
            aligner.smart_align(ref_seq, rpos1, rpos2, query_seq, qpos1, qpos2, (), new_cigar);
            i = self.len();
        }
        if let Some(mut new_cigar) = lazy_new_cigar {
            new_cigar.tuples.extend_from_slice(&self.tuples[i..]);
            self.tuples = new_cigar.tuples;
        }
    }

    /// Given two alignments (cigar_ij and cigar_jk), where sequence i starts at start_i compared to sequence j,
    /// computes new alignment from i to k.
    ///
    /// Returns new CIGAR and start in position k.
    ///
    /// Based on `direction_ij`, i can be query or reference compared to k;
    /// `direction_jk`, j can be query or reference compared to k.
    ///
    /// Anchor size represents the minimal size of the `=` operation that will act as an anchor.
    fn transfer_alignment<const GLOBAL_ALN: bool>(
        start_j: u32,
        offset: CigarOffset,
        cigar_ij: &Cigar,
        direction_ij: impl CigarDirection,
        cigar_jk: &Cigar,
        direction_jk: impl CigarDirection,
        seq_i: &[u8],
        seq_k: &[u8],
        aligner: &Aligner,
        max_gap: impl wfa::Threshold,
        anchor_size: u32,
    ) -> (u32, Cigar) {
        use Operation::Equal;
        let mut cigar_jk_iter = cigar_jk.tuples[offset.cigar_ix..].iter();
        let tmp = cigar_jk_iter.next().unwrap();
        let mut op2 = direction_jk.operation(tmp.op);
        let init_shift = start_j - offset.qpos_at_ix;
        let mut len2 = tmp.len;
        let mut rem2 = len2 - init_shift;

        // Alignment starting position in the new haplotype. Could change when soft clipping is examined.
        let mut start_k = offset.rpos_at_ix + ifelse0(op2.consumes_ref(), init_shift);
        let len_k = seq_k.len() as u32;

        if GLOBAL_ALN {
            if cigar_ij.full_sequence_match() {
                return (0, if direction_jk.invert() { cigar_jk.invert() } else { cigar_jk.clone() });
            } else if cigar_jk.full_sequence_match() {
                return (0, if direction_ij.invert() { cigar_ij.invert() } else { cigar_ij.clone() });
            }
        } else {
            // If read is completely enclosed in a "=" entry + padding, simply copy the alignment.
            const FULL_MATCH_PADDING: u32 = 3;
            debug_assert!(!direction_ij.invert());
            if op2 == Equal && init_shift >= FULL_MATCH_PADDING && rem2 >= cigar_ij.rlen + FULL_MATCH_PADDING {
                return (start_k, cigar_ij.clone());
            }
        }

        let mut cigar_ij_iter = cigar_ij.iter();
        let tmp = cigar_ij_iter.next().expect("CIGAR is empty");
        let mut len1 = tmp.len;
        let mut rem1 = len1;
        let mut op1 = direction_ij.operation(tmp.op);
        // For sequences i (1) and k (2), `pos` is the current position in both cigars.
        // `last` is the end of the sequence added to the new cigar.
        // So, gap `last..pos` is not yet aligned.
        let mut last1 = 0;
        let mut pos1 = 0;
        let mut last2 = start_k;
        let mut pos2 = start_k;

        // When aligning padding, compare read tail v. haplotype sequence of length(tail) + PADDING.
        const CLIP_PADDING: u32 = 3;
        // Anchor margin represents the size of the `=` operation to the left of the anchor
        // when the second operation is not `=`.
        const ANCHOR_MARGIN: u32 = 5;
        let mut new_cigar = Cigar::new();
        loop {
            // Copy operation from one of the CIGARs if one is "=" of sufficient length.
            let add_operation = match (op1 == Equal, op2 == Equal) {
                (true, true) if min(rem1, rem2) >= anchor_size => Some(Equal),
                (true, false) if rem1 >= anchor_size && len1 - rem1 >= ANCHOR_MARGIN => Some(op2),
                (false, true) if rem2 >= anchor_size && len2 - rem2 >= ANCHOR_MARGIN => Some(op1),
                _ => None,
            };
            if add_operation.is_some() {
                if !GLOBAL_ALN && last1 == 0 && pos1 > 0 {
                    aligner.align_ends::<true>(seq_k, last2.saturating_sub(pos1 + CLIP_PADDING), pos2,
                        seq_i, last1, pos1, &mut new_cigar);
                    start_k = start_k + pos2 - last2 - new_cigar.rlen;
                } else {
                    aligner.smart_align(seq_k, last2, pos2, seq_i, last1, pos1, max_gap, &mut new_cigar);
                }
            }
            // Based on the two operations, calculate how much positions and CIGAR shifts should be updated.
            let shift = double_cigar_move_and_shift(op1, op2, &mut pos1, &mut rem1, &mut pos2, &mut rem2);
            if let Some(op) = add_operation {
                new_cigar.push_checked(op, shift);
                last1 = pos1;
                last2 = pos2;
            }

            if rem1 == 0 {
                let Some(tmp) = cigar_ij_iter.next() else { break };
                len1 = tmp.len;
                rem1 = len1;
                op1 = direction_ij.operation(tmp.op);
            }
            if rem2 == 0 {
                let Some(tmp) = cigar_jk_iter.next() else { break };
                len2 = tmp.len;
                rem2 = len2;
                op2 = direction_jk.operation(tmp.op);
            }
        }
        let len_i = seq_i.len() as u32;
        if last1 != len_i {
            if GLOBAL_ALN {
                aligner.smart_align(seq_k, last2, len_k, seq_i, last1, len_i, max_gap, &mut new_cigar);
            } else {
                aligner.align_ends::<false>(seq_k, last2, min(len_k, last2 + len_i - last1 + CLIP_PADDING),
                    seq_i, last1, len_i, &mut new_cigar);
            }
        }
        assert_eq!(len_i, new_cigar.qlen);
        if GLOBAL_ALN {
            const MAX_OPTIMIZATION_GAP: u32 = 1000;
            const OPTIMIZATION_ANCHOR: u32 = 51;
            new_cigar.optimize(seq_k, seq_i, aligner, MAX_OPTIMIZATION_GAP, OPTIMIZATION_ANCHOR);
        } else {
            new_cigar.boundary_ins_to_soft();
        }
        (start_k, new_cigar)
    }

    #[inline]
    pub fn transfer_read_alignment(
        &self,
        qpos: u32,
        offset: CigarOffset,
        read_cigar: &Cigar,
        read_seq: &[u8],
        ref_seq: &[u8],
        aligner: &Aligner,
        direction: impl CigarDirection,
    ) -> (u32, Cigar) {
        const ANCHOR_SIZE: u32 = 5;
        Self::transfer_alignment::<false>(
            qpos, offset, read_cigar, QueryToRef, self, direction, read_seq, ref_seq, aligner, (), ANCHOR_SIZE)
    }

    /// Using alignments between haplotypes i, j and j, k; find new alignment between i (query) and k (reference).
    /// In alignment i-j j is reference if `j_ref_ij`;
    /// similarly, in alignment j-k k is reference ij `k_ref_jk`.
    pub fn find_transitive_alignment(
        cigar_ij: &Cigar,
        j_ref_ij: bool,
        cigar_jk: &Cigar,
        k_ref_jk: bool,
        seq_i: &[u8],
        seq_k: &[u8],
        aligner: &Aligner,
        max_gap: u32,
        anchor_size: u32,
    ) -> Cigar {
        let offset0 = CigarOffset::default();
        let (pos, cigar) = match (j_ref_ij, k_ref_jk) {
            (true, true) => Self::transfer_alignment::<true>(
                0, offset0, cigar_ij, QueryToRef, cigar_jk, QueryToRef, seq_i, seq_k, aligner, max_gap, anchor_size),
            (true, false) => Self::transfer_alignment::<true>(
                0, offset0, cigar_ij, QueryToRef, cigar_jk, RefToQuery, seq_i, seq_k, aligner, max_gap, anchor_size),
            (false, true) => Self::transfer_alignment::<true>(
                0, offset0, cigar_ij, RefToQuery, cigar_jk, QueryToRef, seq_i, seq_k, aligner, max_gap, anchor_size),
            (false, false) => Self::transfer_alignment::<true>(
                0, offset0, cigar_ij, RefToQuery, cigar_jk, RefToQuery, seq_i, seq_k, aligner, max_gap, anchor_size),
        };
        assert_eq!(pos, 0);
        assert_eq!(cigar.rlen, seq_k.len() as u32);
        assert!(cigar.validate(seq_i, seq_k)); // [TODO] REMOVE
        cigar
    }
}

/// Depending on the two parallel CIGAR operations (read v. hapQ) and (hapQ v. hapT)
/// depermine whether read position and hapT position should change,
/// and whether read and hapT CIGAR iterator should move.
/// Returns shift size.
#[inline(always)]
fn double_cigar_move_and_shift(
    op1: Operation,
    op2: Operation,
    pos1: &mut u32,
    rem1: &mut u32,
    pos2: &mut u32,
    rem2: &mut u32,
) -> u32 {
    // 1. Read CIGAR shifts, unless *both*
    //    a. deletion in hapQ relative to hapQ,
    //    b. no insertion in the read rel to hapQ.
    // 2. Read moves, whenever
    //    a. read CIGAR shifts,
    //    b. there is no deletion in the read rel to hapQ.
    // Analogous ideas for the hapT shift and move.
    let (read_moves, read_cigar_shifts, hap_moves, hap_cigar_shifts) = match (op1.consumes(), op2.consumes()) {
        (Consumes::Both,  Consumes::Both) => (true,  true,  true,  true),
        // Insertion in read.
        (Consumes::Query, Consumes::Both) => (true,  true,  false, false),
        // Deletion in read.
        (Consumes::Ref,   Consumes::Both) => (false, true,  true,  true),

        // Insertion in hapQ rel to hapT.
        (Consumes::Both,  Consumes::Query) => (true,  true,  false, true),
        // Insertion in read & insertion in hapQ, which is now ignored.
        (Consumes::Query, Consumes::Query) => (true,  true,  false, false),
        // Deletion in the read and insertion in hapQ. Shift both CIGARs.
        (Consumes::Ref,   Consumes::Query) => (false, true,  false, true),

        // Deletion in hapQ, shift and move hapT.
        (Consumes::Both,  Consumes::Ref) => (false, false, true,  true),
        // Twice deletion in hapQ, move everything.
        (Consumes::Query, Consumes::Ref) => (true,  true,  true,  true),
        // Deletion in read, deletion in hapQ. Shift and move hapT.
        (Consumes::Ref,   Consumes::Ref) => (false, false, true,  true),
    };

    // If both CIGAR shifts, take minimal of rems, otherwise, take one of the rems directly.
    let shift = if read_cigar_shifts && (!hap_cigar_shifts || *rem1 <= *rem2) { *rem1 } else { *rem2 };
    *pos1 += ifelse0(read_moves, shift);
    *rem1 -= ifelse0(read_cigar_shifts, shift);
    *pos2 += ifelse0(hap_moves, shift);
    *rem2 -= ifelse0(hap_cigar_shifts, shift);
    shift
}
