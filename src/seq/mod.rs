pub mod contigs;
pub mod interv;
pub mod cigar;
pub mod aln;
pub mod kmers;
pub mod panvcf;
pub mod fastx;
pub mod recruit;
pub mod wfa;
pub mod compl;
pub mod counts;
#[cfg(feature = "align")]
pub mod dist;
pub mod div;

pub use interv::{Interval, NamedInterval};
pub use contigs::{ContigId, ContigNames, ContigSet};
pub use fastx::{write_fasta, write_multiline_fasta};

use std::io::{Read, Seek};
use bio::io::fasta;
use crate::err::{error, add_path};

/// Make nucleotide sequence standard: only letters A,C,G,T,N.
pub fn standardize(seq: &mut [u8]) -> Result<(), u8> {
    for nt in seq.iter_mut() {
        *nt = match &nt {
            b'A' | b'a' => b'A',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            b'T' | b't' => b'T',
            b'N' | b'R' | b'Y' | b'K' | b'M' | b'S' | b'W' | b'B' | b'D' | b'H' | b'V'
                | b'n' | b'r' | b'y' | b'k' | b'm' | b's' | b'w' | b'b' | b'd' | b'h' | b'v' => b'N',
            _ => return Err(*nt),
        };
    }
    Ok(())
}

/// Does the sequence include unknown (N) nucleotides?
pub fn has_n(seq: &[u8]) -> bool {
    seq.iter().any(|&nt| nt == b'N')
}

/// Count the number of C,G nucleotides in the sequence.
pub fn gc_count(seq: &[u8]) -> u32 {
    seq.iter().fold(0_u32, |acc, &nt| acc + u32::from(nt == b'C' || nt == b'G'))
}

/// Calculate GC-content (between 0 and 100).
pub fn gc_content(seq: &[u8]) -> f64 {
    100.0 * f64::from(gc_count(seq)) / seq.len() as f64
}

/// Finds runs of Ns in the sequence and returns pairs (start, end).
pub fn n_runs(seq: &[u8]) -> Vec<(u32, u32)> {
    let mut start = 0;
    let mut n_run = false;
    let mut runs = Vec::new();
    for (i, &nt) in seq.iter().enumerate() {
        if nt == b'N' && !n_run {
            start = i as u32;
            n_run = true;
        } else if nt != b'N' && n_run {
            runs.push((start, i as u32));
            n_run = false;
        }
    }
    if n_run {
        runs.push((start, seq.len() as u32));
    }
    runs
}

/// Sequence with its name.
#[derive(Clone)]
pub struct NamedSeq {
    name: String,
    seq: Vec<u8>,
}

impl NamedSeq {
    /// Constructs new named sequence.
    pub fn new(name: String, seq: Vec<u8>) -> Self {
        Self { name, seq }
    }

    /// Returns reference to the name.
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn name_mut(&mut self) -> &mut String {
        &mut self.name
    }

    /// Returns reference to the sequence.
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    /// Returns mutable sequence reference.
    pub fn seq_mut(&mut self) -> &mut Vec<u8> {
        &mut self.seq
    }

    /// Consumes this object and returns name.
    pub fn take_name(self) -> String {
        self.name
    }

    /// Consumes this object and returns sequence.
    pub fn take_seq(self) -> Vec<u8> {
        self.seq
    }

    /// Sequence length.
    pub fn len(&self) -> u32 {
        self.seq.len() as u32
    }
}

#[inline]
pub fn complement_nt(nt: u8) -> u8 {
    match nt {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        b'N' | b'R' | b'Y' | b'K' | b'M' | b'S' | b'W' | b'B' | b'D' | b'H' | b'V'
            | b'n' | b'r' | b'y' | b'k' | b'm' | b's' | b'w' | b'b' | b'd' | b'h' | b'v' => b'N',
        _ => panic!("Unknown nucleotide {} ({})", char::from(nt), nt),
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().copied().map(complement_nt).collect()
}

/// Fetches sequence of the interval from an indexed fasta reader.
pub fn fetch_seq(
    fasta: &mut fasta::IndexedReader<impl Read + Seek>,
    contig: &str,
    start: u64,
    end: u64,
) -> crate::Result<Vec<u8>> {
    assert!(start < end, "Cannot fetch sequence for an empty interval {}:{}-{}", contig, start + 1, end);
    fasta.fetch(contig, start.into(), end.into()).map_err(add_path!(!))?;
    let len = (end - start) as usize;
    let mut seq = Vec::with_capacity(len);
    fasta.read(&mut seq).map_err(add_path!(!))?;
    if seq.len() != len {
        return Err(error!(RuntimeError, "Fetched sequence of incorrect size ({}) for interval {}:{}-{} ({} bp)",
        seq.len(), contig, start + 1, end, len));
    }

    crate::seq::standardize(&mut seq).map_err(|nt| error!(InvalidData,
        "Unknown nucleotide `{}` ({}) in {}:{}-{}", char::from(nt), nt, contig, start + 1, end))?;
    Ok(seq)
}
