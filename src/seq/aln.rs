//! Various functions and structures related to sequence alignment.

/// Alignment metrics.
pub struct Metrics {
    pub matches: u32,
    pub mismatches: u32,
    pub insertions: u32,
    pub deletions: u32,
}
