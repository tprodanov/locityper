use std::{
    ops::{Add, AddAssign},
    fmt,
};
use crate::{
    Error,
    math::Ln,
    seq::cigar::{Operation, Cigar, RAW_OPERATIONS},
    math::distr::{DiscretePmf, NBinom, Multinomial},
    bg::ser::{JsonSer, json_get},
};

/// Counts of five operation types (=, X, I, D, S).
#[derive(Clone, Default, Debug)]
struct OpCounts<T> {
    matches: T,
    mismatches: T,
    insertions: T,
    deletions: T,
    clipping: T,
}

impl<T> OpCounts<T>
{
    /// Counts operations in a CIGAR. CIGAR must not contain "M" operation.
    fn calculate(cigar: &Cigar) -> Self
    where u32: TryInto<T>,
          <u32 as TryInto<T>>::Error: fmt::Debug,
          T: Eq + Add<Output = T> + Copy + fmt::Debug,
    {
        let mut counts = [0_u32; RAW_OPERATIONS];
        let mut sum_len = 0;
        for item in cigar.iter() {
            counts[item.operation().ix()] += item.len();
            sum_len += item.len();
        }
        let res = Self {
            matches: counts[Operation::Equal.ix()].try_into().unwrap(),
            mismatches: counts[Operation::Diff.ix()].try_into().unwrap(),
            insertions: counts[Operation::Ins.ix()].try_into().unwrap(),
            deletions: counts[Operation::Del.ix()].try_into().unwrap(),
            clipping: counts[Operation::Soft.ix()].try_into().unwrap(),
        };
        assert_eq!(res.matches + res.mismatches + res.insertions + res.deletions + res.clipping,
            sum_len.try_into().unwrap(), "Cigar {} contains unexpected operations!", cigar);
        res
    }
}

impl OpCounts<u64> {
    /// Converts operation counts into error profile.
    /// Error probabilities (mismatches, insertions & deletions) are increased by `err_rate_mult`
    /// to account for possible mutations in the data.
    ///
    /// Clipping is ignored.
    fn get_profile(&self, err_rate_mult: f64) -> ErrorProfile {
        // Clipping is ignored.
        let sum_len = (self.matches + self.mismatches + self.insertions + self.deletions) as f64;
        let mism_prob = self.mismatches as f64 / sum_len;
        let corr_mism_prob = mism_prob * err_rate_mult;
        let ins_prob = self.insertions as f64 / sum_len;
        let corr_ins_prob = ins_prob * err_rate_mult;
        let del_prob = self.deletions as f64 / sum_len;
        let corr_del_prob = del_prob * err_rate_mult;
        let corr_match_prob = 1.0 - corr_mism_prob - corr_ins_prob - corr_del_prob;

        log::info!("    {:12} matches    ({:.6} -> {:.6})",
            self.matches, self.matches as f64 / sum_len, corr_match_prob);
        log::info!("    {:12} mismatches ({:.6} -> {:.6})",
            self.mismatches, mism_prob, corr_mism_prob);
        log::info!("    {:12} insertions ({:.6} -> {:.6})",
            self.insertions, ins_prob, corr_ins_prob);
        log::info!("    {:12} deletions  ({:.6} -> {:.6})",
            self.deletions, del_prob, corr_del_prob);
        ErrorProfile::new(corr_match_prob, corr_mism_prob, corr_ins_prob, corr_del_prob)
    }
}

impl<T> fmt::Display for OpCounts<T>
where T: Add<Output = T> + Copy + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Matches: {}, Mism: {}, Ins: {}, Del: {}, Clip: {}",
            self.matches, self.mismatches, self.insertions, self.deletions, self.clipping)
    }
}

impl<T, U> AddAssign<&OpCounts<U>> for OpCounts<T>
where U: TryInto<T> + Copy,
      U::Error: fmt::Debug,
      T: AddAssign,
{
    /// Sums operation counts with other operation counts.
    fn add_assign(&mut self, oth: &OpCounts<U>) {
        self.matches += oth.matches.try_into().unwrap();
        self.mismatches += oth.mismatches.try_into().unwrap();
        self.insertions += oth.insertions.try_into().unwrap();
        self.clipping += oth.clipping.try_into().unwrap();
        self.deletions += oth.deletions.try_into().unwrap();
    }
}

/// Read error profile.
/// Basically, each operation is penalized according to its probability in the background data.
/// This ensures that an alignment with fewer errors receives higher probability than the more errorneous alignment,
/// even at high error probabilities (< 50% is still required).
///
/// Note that two alignments of the same read can be directly compared between each other,
/// however, absolute probability values are not relevant, and cannot be directly compared.
#[derive(Clone)]
pub struct ErrorProfile {
    /// ln-probabilities of different operations, summing up to one.
    match_lik: f64,
    mism_lik: f64,
    ins_lik: f64,
    del_lik: f64,
    /// Clipping probability is selected as max(mismatch_lik, insert_lik).
    clip_lik: f64,
}

impl ErrorProfile {
    pub fn new(match_prob: f64, mism_prob: f64, ins_prob: f64, del_prob: f64) -> Self {
        assert!(match_prob > 0.5, "Match probability ({:.5}) must be over 50%", match_prob);
        assert!((match_prob + mism_prob + ins_prob + del_prob - 1.0).abs() < 1e-8,
            "Error probabilities do not sum to one");
        Self {
            match_lik: match_prob.ln(),
            mism_lik: mism_prob.ln(),
            ins_lik: ins_prob.ln(),
            del_lik: del_prob.ln(),
            clip_lik: mism_prob.max(ins_prob).ln(),
        }
    }

    /// Create error profile from the iterator over CIGARs (with extended operations X/=).
    pub fn estimate<'a>(
        count: usize,
        cigars: impl Iterator<Item = &'a Cigar>,
        err_rate_mult: f64,
    ) -> ErrorProfile
    {
        log::info!("Estimating read error profiles from {} reads", count);
        let mut prof_builder = OpCounts::<u64>::default();
        for cigar in cigars {
            prof_builder += &OpCounts::<u32>::calculate(cigar);
        }
        prof_builder.get_profile(err_rate_mult)
    }

    /// Returns alignment ln-probability.
    pub fn ln_prob(&self, cigar: &Cigar) -> f64 {
        let read_prof = OpCounts::<u32>::calculate(cigar);
        self.match_lik * f64::from(read_prof.matches)
            + self.mism_lik * f64::from(read_prof.mismatches)
            + self.ins_lik * f64::from(read_prof.insertions)
            + self.del_lik * f64::from(read_prof.deletions)
            + self.clip_lik * f64::from(read_prof.clipping)
    }

}

impl fmt::Debug for ErrorProfile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "ErrorProfile {{ M: {}, X: {}, I: {}, D: {} }}",
            self.match_lik, self.mism_lik, self.ins_lik, self.del_lik)
    }
}

impl JsonSer for ErrorProfile {
    fn save(&self) -> json::JsonValue {
        json::object!{
            m: self.match_lik,
            x: self.mism_lik,
            i: self.ins_lik,
            d: self.del_lik,
            s: self.clip_lik,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> m (as_f64), x (as_f64), i (as_f64), d (as_f64), s (as_f64));
        Ok(Self {
            match_lik: m,
            mism_lik: x,
            ins_lik: i,
            del_lik: d,
            clip_lik: s,
        })
    }
}
