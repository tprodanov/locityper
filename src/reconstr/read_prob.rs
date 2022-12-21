use rust_htslib::bam::record::Record;

use crate::algo::math::Phred;

/// Interface for read probabilities.
pub trait ReadProb {
    /// Get natural log probability of a read alignment.
    fn read_log_prob(record: &Record) -> f64;
}

/// Simple `ReadProb` implementation that simply takes probability from mapping quality.
pub struct MapqReadProb;

impl MapqReadProb {
    pub fn new() -> Self {
        Self {}
    }
}

impl ReadProb for MapqReadProb {
    fn read_log_prob(record: &Record) -> f64 {
        Phred::to_lprob(record.mapq() as f64)
    }
}