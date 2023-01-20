use std::fmt;
use htslib::bam::record::Record;
use crate::{
    seq::cigar::{Cigar, Operation},
    seq::aln::{READ_ENDS, ReadEnd},
    bg::ser::{JsonSer, LoadError, parse_f64_arr},
};

/// Trait to calculate probabilities of an alignment.
pub trait ErrorProfile : Clone + fmt::Debug + JsonSer {
    /// Get log-probability of the record alignment.
    fn ln_prob(&self, ext_cigar: &Cigar, read_end: ReadEnd) -> f64;
}

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
        self.trans_counts[prev * N_OPS + prev] += u64::from(first.len() - 1);

        let mut curr;
        for item in iter {
            curr = cigar_op_to_int(item.operation());
            self.trans_counts[prev * N_OPS + curr] += 1;
            self.trans_counts[curr * N_OPS + curr] += u64::from(item.len() - 1);
            prev = curr;
        }
        self.end_counts[prev] += 1;
    }
}

/// Single mate error profile. All probabilities are in log-space.
#[derive(Clone, Debug, Default)]
struct MateErrorProfile {
    start_probs: [f64; N_OPS],
    trans_probs: [f64; N_OPS_SQ],
    end_probs: [f64; N_OPS],
}

impl MateErrorProfile {
    fn ln_prob(&self, ext_cigar: &Cigar) -> f64 {
        let mut iter = ext_cigar.iter();
        let first = iter.next().expect("Cannot calculate error counts from an empty CIGAR!");
        let mut prev = cigar_op_to_int(first.operation());
        let mut prob = self.start_probs[prev]
            + f64::from(first.len() - 1) * self.trans_probs[prev * N_OPS + prev];

        let mut curr;
        for item in iter {
            curr = cigar_op_to_int(item.operation());
            prob += self.trans_probs[prev * N_OPS + curr]
                + f64::from(item.len() - 1) * self.trans_probs[curr * N_OPS + curr];
            prev = curr;
        }
        prob + self.end_probs[prev]
    }
}

impl JsonSer for MateErrorProfile {
    fn save(&self) -> json::JsonValue {
        json::object!{
            start: &self.start_probs as &[f64],
            trans: &self.trans_probs as &[f64],
            end: &self.end_probs as &[f64],
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        let mut res = Self::default();
        parse_f64_arr(obj, "start", &mut res.start_probs)?;
        parse_f64_arr(obj, "trans", &mut res.trans_probs)?;
        parse_f64_arr(obj, "end", &mut res.end_probs)?;
        Ok(res)
    }
}

/// Error profile.
/// Takes into account transition probabilities between CIGAR operations (match->match), etc.
#[derive(Debug, Clone)]
pub struct TransErrorProfile {
    profiles: [MateErrorProfile; READ_ENDS],
}

impl TransErrorProfile {
    /// Create error profile from the iterator over records.
    pub fn estimate<'a>(records: impl Iterator<Item = &'a Record>) -> TransErrorProfile {
        log::info!("    Estimating read error profiles");
        let mut counts1 = ErrorCounts::new();
        let mut counts2 = ErrorCounts::new();
        for record in records {
            // Unmapped or secondary.
            if (record.flags() & 3844) != 0 {
                continue;
            }
            let ext_cigar = Cigar::infer_ext_cigar(&record);

            if record.is_last_in_template() {
                counts2.update(&ext_cigar)
            } else {
                counts1.update(&ext_cigar)
            }
        }
        TransErrorProfile {
            profiles: [counts1.to_profile(), counts2.to_profile()],
        }
    }
}

impl JsonSer for TransErrorProfile {
    fn save(&self) -> json::JsonValue {
        json::object!{
            prof1: self.profiles[0].save(),
            prof2: self.profiles[1].save(),
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        if obj.has_key("prof1") && obj.has_key("prof2") {
            Ok(Self {
                profiles: [MateErrorProfile::load(&obj["prof1"])?, MateErrorProfile::load(&obj["prof2"])?],
            })
        } else {
            Err(LoadError(format!("TransErrorProfile: Failed to parse '{}': missing 'prof1' or 'prof2' keys!", obj)))
        }
    }
}

impl ErrorProfile for TransErrorProfile {
    fn ln_prob(&self, ext_cigar: &Cigar, read_end: ReadEnd) -> f64 {
        self.profiles[read_end.ix()].ln_prob(ext_cigar)
    }
}
