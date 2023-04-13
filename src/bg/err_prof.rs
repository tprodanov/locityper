use std::{
    ops::{Add, AddAssign},
    fmt,
};
use htslib::bam::Record;
use crate::{
    algo::vec_ext::{VecExt, F64Ext},
    seq::cigar::{Operation, Cigar, RAW_OPERATIONS},
    math::{
        Ln,
        distr::{NBinom, Multinomial},
    },
    bg::ser::{JsonSer, LoadError},
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
        let read_len = self.matches + self.mismatches + self.insertions;
        let read_lenf = read_len as f64;
        let mism_prob = self.mismatches as f64 / read_lenf;
        let corr_mism_prob = mism_prob * err_rate_mult;
        let ins_prob = self.insertions as f64 / read_lenf;
        let corr_ins_prob = ins_prob * err_rate_mult;
        let del_prob = self.deletions as f64 / (read_lenf + self.deletions as f64);
        let corr_del_prob = del_prob * err_rate_mult;
        let corr_match_prob = 1.0 - corr_mism_prob - corr_ins_prob;
        assert!(corr_match_prob >= 0.5, "Match probability ({:.5}) canont be less than 50%", corr_match_prob);
        assert!(corr_del_prob < 0.5, "Deletion probability ({:.5}) cannot be over 50%", corr_del_prob);

        log::info!("    {:12} matches    ({:.6} -> {:.6})",
            self.matches, self.matches as f64 / read_lenf, corr_match_prob);
        log::info!("    {:12} mismatches ({:.6} -> {:.6})",
            self.mismatches, mism_prob, corr_mism_prob);
        log::info!("    {:12} insertions ({:.6} -> {:.6})",
            self.insertions, ins_prob, corr_ins_prob);
        log::info!("    {:12} deletions  ({:.6} -> {:.6})",
            self.deletions, del_prob, corr_del_prob);

        ErrorProfile {
            // Probabilities are normalized in Multinomial::new.
            op_distr: Multinomial::new(&[corr_match_prob, corr_ins_prob, corr_del_prob]),
            no_del_prob: 1.0 - corr_del_prob,
            min_ln_prob: f64::NEG_INFINITY,
        }
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

/// Read error profile, characterised by two distributions:
/// * Multinomial(P(matches), P(mismatches), P(insertions)),
/// * Neg.Binomial(p = no_deletion, n = read_length):
///       where success = extend read length by one; failure = skip one base-pair in the reference.
#[derive(Clone)]
pub struct ErrorProfile {
    /// Distribution of matches, mismatches and insertions.
    op_distr: Multinomial,
    /// Probability of success (no deletion) per base-pair in the read sequence.
    /// Used in Neg.Binomial distribution.
    no_del_prob: f64,
    /// Minimum allowed alignment likelihood (in ln-space).
    min_ln_prob: f64,
}

impl ErrorProfile {
    /// Create error profile from the iterator over records.
    pub fn estimate<'a, I, F>(
        records: I,
        cigar_getter: F,
        max_insert_size: i64,
        params: &super::Params,
    ) -> ErrorProfile
    where I: Iterator<Item = &'a Record>,
          F: Fn(&Record) -> Cigar,
    {
        log::info!("Estimating read error profiles");
        let mut n_reads = 0;
        let mut prof_builder = OpCounts::<u64>::default();
        let mut cigars = Vec::new();
        for record in records {
            if super::read_unpaired_or_proper_pair(record, max_insert_size) {
                n_reads += 1;
                let cigar = cigar_getter(record);
                prof_builder += &OpCounts::<u32>::calculate(&cigar);
                cigars.push(cigar);
            }
        }
        log::info!("    Used {} reads", n_reads);
        let mut prof = prof_builder.get_profile(params.err_rate_mult);

        let mut probs: Vec<_> = cigars.iter().map(|cigar| prof.ln_prob(cigar)).collect();
        VecExt::sort(&mut probs);
        let min_ln_prob = params.err_quantile_mult * F64Ext::quantile_sorted(&probs, params.err_quantile);
        prof.min_ln_prob = min_ln_prob;
        log::info!("    Minimum alignment likelihood: 10^({:.2})", Ln::to_log10(min_ln_prob));
        prof
    }

    /// Returns minimum allowed alignment ln-probability.
    pub fn min_aln_ln_prob(&self) -> f64 {
        self.min_ln_prob
    }

    /// Returns alignment ln-probability.
    pub fn ln_prob(&self, cigar: &Cigar) -> f64 {
        let read_prof = OpCounts::<u32>::calculate(cigar);
        NBinom::new(f64::from(cigar.query_len()), self.no_del_prob).ln_pmf_f64(f64::from(read_prof.deletions))
            + self.op_distr.ln_pmf(
                &[read_prof.matches, read_prof.mismatches + read_prof.clipping, read_prof.insertions])
    }
}

impl fmt::Debug for ErrorProfile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "ErrorProfile {{ (M+S, X, I): {:?}; D: NBinom(p = {:.6}) }}", self.op_distr, self.no_del_prob)
    }
}

impl JsonSer for ErrorProfile {
    fn save(&self) -> json::JsonValue {
        json::object!{
            op_distr: self.op_distr.save(),
            no_del_prob: self.no_del_prob,
            min_ln_prob: self.min_ln_prob,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        if !obj.has_key("op_distr") {
            return Err(LoadError(format!(
                "ErrorProfile: Failed to parse '{}': missing 'op_distr' field!", obj)))?;
        }
        let op_distr = Multinomial::load(&obj["op_distr"])?;
        let no_del_prob = obj["no_del_prob"].as_f64().ok_or_else(|| LoadError(format!(
            "ErrorProfile: Failed to parse '{}': missing or incorrect 'no_del_prob' field!", obj)))?;
        let min_ln_prob = obj["min_ln_prob"].as_f64().ok_or_else(|| LoadError(format!(
            "ErrorProfile: Failed to parse '{}': missing or incorrect 'min_ln_prob' field!", obj)))?;
        Ok(Self { op_distr, no_del_prob, min_ln_prob })
    }
}
