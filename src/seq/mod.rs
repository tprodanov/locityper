pub mod contigs;
pub mod interv;
pub mod cigar;
pub mod aln;
#[cfg(feature = "devel")]
pub mod compl;
pub mod kmers;
pub mod panvcf;
pub mod dist;
pub mod fastx;
pub mod recruit;
pub mod wfa;

pub use interv::{Interval, NamedInterval};
pub use contigs::{ContigId, ContigNames, ContigSet};
pub use fastx::write_fasta;

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
