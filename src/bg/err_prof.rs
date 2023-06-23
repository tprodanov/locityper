use std::{
    ops::{Add, AddAssign, Sub},
    fmt,
    cell::RefCell,
    borrow::Borrow,
    collections::HashMap,
};
use nohash::IntMap;
use crate::{
    Error,
    seq::{
        aln::Alignment,
        cigar::{Operation, Cigar, RAW_OPERATIONS},
    },
    math::distr::BetaBinomial,
    bg::ser::{JsonSer, json_get},
};

/// Counts of five operation types (=, X, I, D, S).
#[derive(Clone, Default, Debug)]
pub struct OpCounts<T> {
    matches: T,
    mismatches: T,
    insertions: T,
    deletions: T,
    clipping: T,
}

impl OpCounts<u32> {
    /// Counts operations in a CIGAR. CIGAR must not contain "M" operation.
    pub fn from_cigar(cigar: &Cigar) -> Self {
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
            sum_len, "Cigar {} contains unexpected operations!", cigar);
        res
    }
}

impl<T: Add<Output = T> + Sub<Output = T> + Copy> OpCounts<T> {
    /// Calculates edit distance based on read length in three operations, not four.
    pub fn edit_distance(&self, read_len: T) -> T {
        read_len - self.matches + self.deletions
    }
}

impl OpCounts<u64> {
    /// Converts operation counts into error profile.
    /// Error probabilities (mismatches, insertions & deletions) are increased by `err_rate_mult`
    /// to account for possible mutations in the data.
    ///
    /// Returns (match_prob, mism_prob, ins_prob, del_prob).
    pub fn get_op_probs(&self, err_rate_mult: f64) -> (f64, f64, f64, f64) {
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
        (corr_match_prob, corr_mism_prob, corr_ins_prob, corr_del_prob)
    }
}

impl<T: fmt::Display> fmt::Display for OpCounts<T> {
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
    ln_match: f64,
    ln_mism: f64,
    ln_ins: f64,
    ln_del: f64,
    /// Clipping probability is selected as max(misln_match, insert_lik).
    ln_clip: f64,

    /// Distribution of edit distance (n = read length).
    edit_dist_distr: BetaBinomial,
    /// Maximum edit distance that is in the one-sided Beta-Binomial confidence interval.
    max_edit_dist: RefCell<IntMap<u32, u32>>,
    /// Confidence level.
    conf_lvl: f64,
}

impl ErrorProfile {
    /// Create error profile from the iterator over CIGARs (with extended operations X/=).
    pub fn estimate<'a>(
        alns: &[impl Borrow<Alignment>],
        mean_read_len: f64,
        params: &super::Params,
    ) -> ErrorProfile
    {
        log::info!("Estimating read error profiles from {} reads", alns.len());
        let mut prof_builder = OpCounts::<u64>::default();
        // (edit_dist, read_len) -> count.
        let mut edit_distances = HashMap::<(u32, u32), u64>::new();

        for aln in alns.iter() {
            let cigar = aln.borrow().cigar();
            let counts = OpCounts::from_cigar(cigar);
            prof_builder += &counts;

            let read_len = cigar.query_len();
            let edit_dist = counts.edit_distance(read_len);
            edit_distances.entry((edit_dist, read_len))
                .and_modify(|counter| *counter += 1)
                .or_insert(1);
        }

        let (match_prob, mism_prob, ins_prob, del_prob) = prof_builder.get_op_probs(params.err_rate_mult);
        assert!(match_prob > 0.5, "Match probability ({:.5}) must be over 50%", match_prob);
        assert!((match_prob + mism_prob + ins_prob + del_prob - 1.0).abs() < 1e-8,
            "Error probabilities do not sum to one");

        let edit_distances: Vec<_> = edit_distances.into_iter()
            .map(|((k, n), count)| (k, n, count as f64))
            .collect();
        let err_prof = Self {
            ln_match: match_prob.ln(),
            ln_mism: mism_prob.ln(),
            ln_ins: ins_prob.ln(),
            ln_del: del_prob.ln(),
            ln_clip: mism_prob.max(ins_prob).ln(),
            edit_dist_distr: BetaBinomial::max_lik_estimate(&edit_distances),
            max_edit_dist: RefCell::default(),
            conf_lvl: params.err_conf_level,
        };

        let read_len = mean_read_len.round() as u32;
        log::info!("    Maximum allowed edit distance: {} (mean read length {}, {}%-confidence interval)",
            err_prof.allowed_edit_dist(read_len), read_len, 100.0 * err_prof.conf_lvl);
        err_prof
    }

    /// Returns ln-probability for operation counts.
    pub fn ln_prob<T>(&self, counts: &OpCounts<T>) -> f64
    where T: Copy + TryInto<f64>,
          <T as TryInto<f64>>::Error: std::fmt::Debug,
    {
        self.ln_match * counts.matches.try_into().unwrap()
            + self.ln_mism * counts.mismatches.try_into().unwrap()
            + self.ln_ins * counts.insertions.try_into().unwrap()
            + self.ln_del * counts.deletions.try_into().unwrap()
            + self.ln_clip * counts.clipping.try_into().unwrap()
    }

    /// Returns the maximum allowed edit distance for the given read length.
    /// Values are cached for future use.
    pub fn allowed_edit_dist(&self, read_len: u32) -> u32 {
        let mut cache = self.max_edit_dist.borrow_mut();
        *cache.entry(read_len).or_insert_with(||
            self.edit_dist_distr.confidence_right_border(read_len, self.conf_lvl))
    }

    /// Returns true if alignment edit distance is not too high.
    pub fn aln_passes(&self, cigar: &Cigar) -> bool {
        let read_prof = OpCounts::from_cigar(cigar);
        let read_len = cigar.query_len();
        let edit_dist = read_prof.edit_distance(read_len);
        edit_dist <= self.allowed_edit_dist(read_len)
    }
}

impl fmt::Debug for ErrorProfile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "ErrorProfile {{ M: {}, X: {}, I: {}, D: {} }}",
            self.ln_match, self.ln_mism, self.ln_ins, self.ln_del)
    }
}

impl JsonSer for ErrorProfile {
    fn save(&self) -> json::JsonValue {
        json::object!{
            ln_match: self.ln_match,
            ln_mism: self.ln_mism,
            ln_ins: self.ln_ins,
            ln_del: self.ln_del,
            ln_clip: self.ln_clip,

            alpha: self.edit_dist_distr.alpha(),
            beta: self.edit_dist_distr.beta(),
            conf_lvl: self.conf_lvl,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> ln_match (as_f64), ln_mism (as_f64), ln_ins (as_f64), ln_del (as_f64), ln_clip (as_f64),
            alpha (as_f64), beta (as_f64), conf_lvl (as_f64));
        Ok(Self {
            ln_match, ln_mism, ln_ins, ln_del, ln_clip, conf_lvl,
            edit_dist_distr: BetaBinomial::new(alpha, beta),
            max_edit_dist: RefCell::default(),
        })
    }
}
