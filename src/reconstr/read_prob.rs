use rust_htslib::bam::record::Record;
// use statrs::distribution::{Binomial, Discrete};

use crate::algo::math::Phred;

/// Interface for read probabilities.
pub trait ReadProb {
    /// Get natural log probability of a read alignment.
    fn read_log_prob(record: &Record) -> f64;
}

/// Simple `ReadProb` implementation that simply takes probability from mapping quality.
#[derive(Debug, Clone, Copy)]
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

// /// Simple `ReadProb` implementation, which equally penalizes all errors.
// #[derive(Debug, Clone)]
// pub struct SingleErrorReadProb {q
//     error_prob: f64,
// }

// impl SingleErrorReadProb {
//     pub fn new(error_prob: f64) -> Self {
//         Self { error_prob }
//     }
// }

// impl Default for SingleErrorReadProb {
//     fn default() -> Self {
//         Self {
//             error_prob: 0.01,
//         }
//     }
// }

// impl ReadProb for MapqReadProb {
//     fn read_log_prob(record: &Record) -> f64 {
//         let ext_cigar = cigar::get_ext_cigar
//     }
// }
