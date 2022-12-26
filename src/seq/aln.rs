//! Various functions and structures related to sequence alignment.

use std::fmt::Write;

use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::bam::record::Cigar;

/// Does the Cigar operation consume reference sequence?
const fn consumes_ref(op: Cigar) -> bool {
    match op {
        Cigar::Match(_) | Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::Equal(_) | Cigar::Diff(_) => true,
        _ => false,
    }
}

/// Does the Cigar operation consume query sequence?
const fn consumes_query(op: Cigar) -> bool {
    match op {
        Cigar::Match(_) | Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::Equal(_) | Cigar::Diff(_) => true,
        _ => false,
    }
}

/// Does the Cigar operation consume both reference and query sequences?
const fn consumes_both(op: Cigar) -> bool {
    match op {
        Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => true,
        _ => false,
    }
}

pub trait StrOrNone {
    fn new() -> Self;

    fn append(&mut self, ch: char);

    fn extend(&mut self, s: &str);
}

impl StrOrNone for () {
    fn new() -> Self { }

    fn append(&mut self, _: char) { }

    fn extend(&mut self, _: &str) { }
}

impl StrOrNone for String {
    fn new() -> Self {
        String::new()
    }

    fn append(&mut self, ch: char) {
        self.write_char(ch).unwrap()
    }

    fn extend(&mut self, s: &str) {
        self.write_str(s).unwrap()
    }
}

enum MdEntry<'a> {
    Match(u16),
    Mismatch(&'a str),
    Deletion(&'a str),
}

fn parse_md<'a>(md: &'a str) -> Vec<MdEntry<'a>> {
    lazy_static! {
        static ref RE_NUMBER: Regex = Regex::new("[0-9]+").unwrap();
    }
    let md_bytes = md.as_bytes();
    let mut res = Vec::new();
    let mut last_end = 0;
    for m in RE_NUMBER.find_iter(md) {
        if last_end > 0 {
            let dist = last_end - m.start();
            assert!(dist > 0, "Failed to parse MD tag {}", md);
            if md_bytes[last_end] == b'^' {
                assert!(dist > 1, "Failed to parse MD tag {}", md);
                res.push(MdEntry::Deletion(&md[last_end + 1..m.start()]));
            } else {
                res.push(MdEntry::Mismatch(&md[last_end..m.start()]));
            }
        } else {
            assert_eq!(m.start(), 0, "Failed to parse MD tag {}", md);
        }
        let length = m.as_str().parse().unwrap();
        if length > 0 {
            res.push(MdEntry::Match(length));
        }
        last_end = m.end();
    }
    assert_eq!(last_end, md.len(), "Failed to parse MD tag {}", md);
    res
}

pub fn to_extended_cigar<Out: StrOrNone>(cigar: impl Iterator<Item = Cigar>, md: &str)
        -> (Vec<Operation>, Out) {
    let mut qpos = 0;
    let mut rpos = 0;
    let mut new_cigar = Vec::new();
    let mut ref_seq = Out::new();

    let md_entries = parse_md(md);
    let mut md_ix = 0;
    let mut md_shift = 0;

    for op in cigar {
        if consumes_both(op) {
            
        }
    }

    (new_cigar, ref_seq)
}
