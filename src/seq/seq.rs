
/// Count the number of C,G nucleotides in the sequence.
pub fn gc_count(seq: &[u8]) -> u32 {
    seq.iter().fold(0_u32, |acc, &nt| acc + (nt == b'C' || nt == b'G') as u32)
}

/// Calculate GC-content (between 0 and 100).
pub fn gc_content(seq: &[u8]) -> f64 {
    100.0 * gc_count(seq) as f64 / seq.len() as f64
}
