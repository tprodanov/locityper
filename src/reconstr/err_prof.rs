use crate::seq::cigar::{Cigar, Operation};

use rust_htslib::bam::record::Record;

const MATCH: usize = 0;
const MISM: usize = 1;
const INS: usize = 2;
const DEL: usize = 3;
const SOFT: usize = 4;
const N_OPS: usize = 5;
const N_OPS_SQ: usize = N_OPS * N_OPS;

/// Converts CIGAR operation into an integer (ignore unnecessary operations).
const fn cigar_op_to_int(op: Operation) -> usize {
    match op {
        Operation::Equal => MATCH,
        Operation::Diff => MISM,
        Operation::Ins => INS,
        Operation::Del => DEL,
        Operation::Soft => SOFT,
        _ => panic!("Unexpected operation in an extended CIGAR!"),
    }
}

/// Private structure for calculating error profiles.
struct ErrorCounts {
    start_counts: [u64; N_OPS],
    trans_counts: [u64; N_OPS_SQ],
    end_counts: [u64; N_OPS],
}

impl ErrorCounts {
    /// Create new error counts, start from ones.
    fn new() -> Self {
        Self {
            start_counts: [1; N_OPS],
            trans_counts: [1; N_OPS_SQ],
            end_counts: [1; N_OPS],
        }
    }

    /// Create mate-error profile from read counts.
    fn to_profile(&self) -> MateErrorProfile {
        let n_reads_ln = (self.start_counts.iter().cloned().sum::<u64>() as f64).ln();
        let n_trans_ln = (self.trans_counts.iter().cloned().sum::<u64>() as f64).ln();
        let mut prof = MateErrorProfile {
            start_probs: [0.0; N_OPS],
            trans_probs: [0.0; N_OPS_SQ],
            end_probs: [0.0; N_OPS],
        };
        for i in 0..N_OPS {
            prof.start_probs[i] = (self.start_counts[i] as f64).ln() - n_reads_ln;
            prof.end_probs[i] = (self.end_counts[i] as f64).ln() - n_reads_ln;
        }
        for i in 0..N_OPS_SQ {
            prof.trans_probs[i] = (self.trans_counts[i] as f64).ln() - n_trans_ln;
        }
        prof
    }

    /// Updates mate error profile from an extended CIGAR.
    fn update(&mut self, ext_cigar: &Cigar) {
        let mut iter = ext_cigar.iter();
        let first = iter.next().expect("Cannot calculate error counts from an empty CIGAR!");
        let mut prev = cigar_op_to_int(first.operation());
        self.start_counts[prev] += 1;
        self.trans_counts[prev * N_OPS + prev] += (first.len() - 1) as u64;

        let mut curr;
        for tup in iter {
            curr = cigar_op_to_int(tup.operation());
            self.trans_counts[prev * N_OPS + curr] += 1;
            self.trans_counts[curr * N_OPS + curr] += (tup.len() - 1) as u64;
            prev = curr;
        }
        self.end_counts[prev] += 1;
    }
}

/// Single mate error profile. All probabilities are in log-space.
struct MateErrorProfile {
    start_probs: [f64; N_OPS],
    trans_probs: [f64; N_OPS_SQ],
    end_probs: [f64; N_OPS],
}

impl MateErrorProfile {
    pub fn ln_prob(&self, ext_cigar: &Cigar) -> f64 {
        let mut iter = ext_cigar.iter();
        let first = iter.next().expect("Cannot calculate error counts from an empty CIGAR!");
        let mut prev = cigar_op_to_int(first.operation());
        let mut prob = self.start_probs[prev]
            + (first.len() - 1) as f64 * self.trans_probs[prev * N_OPS + prev];

        let mut curr;
        for tup in iter {
            curr = cigar_op_to_int(tup.operation());
            prob += self.trans_probs[prev * N_OPS + curr]
                + (tup.len() - 1) as f64 * self.trans_probs[curr * N_OPS + curr];
            prev = curr;
        }
        prob + self.end_probs[prev]
    }
}

/// Error profile.
/// Takes into account transition probabilities between CIGAR operations (match->match), etc.
pub struct ErrorProfile {
    prof1: MateErrorProfile,
    prof2: MateErrorProfile,
}

impl ErrorProfile {
    /// Create error profile from the iterator over records.
    pub fn create(records: impl Iterator<Item = Record>) -> ErrorProfile {
        let mut counts1 = ErrorCounts::new();
        let mut counts2 = ErrorCounts::new();
        for record in records {
            let ext_cigar = Cigar::infer_ext_cigar(&record);
            if record.is_last_in_template() {
                counts2.update(&ext_cigar)
            } else {
                counts1.update(&ext_cigar)
            }
        }
        ErrorProfile {
            prof1: counts1.to_profile(),
            prof2: counts2.to_profile(),
        }
    }

    /// Get log-probability of the record alignment.
    pub fn ln_prob(&self, record: &Record) -> f64 {
        let ext_cigar = Cigar::infer_ext_cigar(record);
        if record.is_last_in_template() {
            self.prof2.ln_prob(&ext_cigar)
        } else {
            self.prof1.ln_prob(&ext_cigar)
        }
    }
}
