pub mod contigs;
pub mod interv;
pub mod cigar;
pub mod aln;
pub mod compl;
pub mod kmers;
#[cfg(feature = "devel")]
pub mod dist;

pub use interv::Interval;
pub use contigs::{ContigId, ContigNames};

use std::{
    cmp::min,
    io::{self, Write},
};

/// Make nucleotide sequence standard: only letters A,C,G,T,N.
pub fn standardize(seq: &mut [u8]) {
    for nt in seq.iter_mut() {
        *nt = match &nt {
            b'A' | b'a' => b'A',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            b'T' | b't' => b'T',
            b'n' | b'R' | b'Y' | b'K' | b'M' | b'S' | b'W' | b'B' | b'D' | b'H' | b'V'
                 | b'r' | b'y' | b'k' | b'm' | b's' | b'w' | b'b' | b'd' | b'h' | b'v'
                 => b'N',
            _ => panic!("Unknown nucleotide {}", *nt as char),
        };
    }
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

/// Write a single sequence to the FASTA file.
/// Use this function instead of `bio::fasta::Writer` as the latter
/// writes the sequence into a single line, without splitting.
pub fn write_fasta<W: Write>(f: &mut W, name: &[u8], seq: &[u8]) -> io::Result<()> {
    const WIDTH: usize = 120;
    f.write_all(b">")?;
    f.write_all(name)?;
    f.write_all(b"\n")?;

    let n = seq.len();
    for i in (0..n).step_by(WIDTH) {
        f.write_all(&seq[i..min(i + WIDTH, n)])?;
        f.write_all(b"\n")?;
    }
    Ok(())
}