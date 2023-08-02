use std::{
    fmt,
    io::Write,
    ops::{Add, AddAssign},
    cell::RefCell,
    borrow::Borrow,
    collections::HashMap,
    path::Path,
};
use nohash::IntMap;
use crate::{
    err::{Error, add_path},
    seq::{
        Interval,
        aln::Alignment,
    },
    math::{self, distr::BetaBinomial},
    bg::ser::{JsonSer, json_get},
    ext,
};

/// Counts of five operation types (=, X, I, D, S).
#[derive(Clone, Default, Debug)]
pub struct OperCounts<T> {
    pub matches: T,
    pub mismatches: T,
    pub insertions: T,
    pub deletions: T,
    pub clipping: T,
}

impl<T: Add<Output = T> + Copy> OperCounts<T> {
    /// Returns pair (edit distance, read length).
    pub fn edit_and_read_len(&self) -> (T, T) {
        let common = self.mismatches + self.insertions + self.clipping;
        (common + self.deletions, common + self.matches)
    }
}

impl OperCounts<u64> {
    /// Calculates ln-probabilities of various error types.
    /// Error probabilities (mismatches, insertions & deletions) are increased by `err_rate_mult`
    /// to account for possible mutations in the data.
    ///
    /// Returned clipping probability is set to the maximum of mismatch and insert probabilities.
    pub fn to_ln_probs(&self, err_rate_mult: f64) -> OperCounts<f64> {
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
        assert!(corr_match_prob > 0.5, "Match probability ({:.5}) must be over 50%", corr_match_prob);
        assert!((corr_match_prob + corr_mism_prob + corr_ins_prob + corr_del_prob - 1.0).abs() < 1e-8,
            "Error probabilities do not sum to one");
        OperCounts {
            matches: corr_match_prob.ln(),
            mismatches: corr_mism_prob.ln(),
            insertions: corr_ins_prob.ln(),
            deletions: corr_del_prob.ln(),
            clipping: corr_ins_prob.max(corr_mism_prob).ln(),
        }
    }
}

impl<T: fmt::Display> fmt::Display for OperCounts<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Matches: {}, Mism: {}, Ins: {}, Del: {}, Clip: {}",
            self.matches, self.mismatches, self.insertions, self.deletions, self.clipping)
    }
}

impl<T, U> AddAssign<&OperCounts<U>> for OperCounts<T>
where U: TryInto<T> + Copy,
      U::Error: fmt::Debug,
      T: AddAssign,
{
    /// Sums operation counts with other operation counts.
    fn add_assign(&mut self, oth: &OperCounts<U>) {
        self.matches += oth.matches.try_into().unwrap();
        self.mismatches += oth.mismatches.try_into().unwrap();
        self.insertions += oth.insertions.try_into().unwrap();
        self.clipping += oth.clipping.try_into().unwrap();
        self.deletions += oth.deletions.try_into().unwrap();
    }
}

pub const DEF_EDIT_PVAL: (f64, f64) = (0.01, 0.001);

/// Read error profile.
/// Basically, each operation is penalized according to its probability in the background data.
/// This ensures that an alignment with fewer errors receives higher probability than the more errorneous alignment,
/// even at high error probabilities (< 50% is still required).
///
/// Note that two alignments of the same read can be directly compared between each other,
/// however, absolute probability values are not relevant, and cannot be directly compared.
#[derive(Clone)]
pub struct ErrorProfile {
    /// ln-probabilities of different operations.
    oper_probs: OperCounts<f64>,
    /// Distribution of edit distance (n = read length).
    edit_distr: BetaBinomial,

    /// Two CDF levels (`1 - pval`), specifying edit distance thresholds:
    /// one for good alignments, second for passable.
    edit_cdf: (f64, f64),
    /// For each read length, store two maximum allowed edit distances (good and passable).
    max_edit_dist: RefCell<IntMap<u32, (u32, u32)>>,
}

impl ErrorProfile {
    /// Create error profile from the iterator over CIGARs (with extended operations X/=).
    /// If `out_dir` is Some, write debug information.
    pub fn estimate<'a>(
        alns: &[impl Borrow<Alignment>],
        region: &Interval,
        windows: &super::Windows,
        mean_read_len: f64,
        params: &super::Params,
        out_dir: Option<&Path>,
    ) -> Result<Self, Error>
    {
        log::info!("Estimating read error profiles from {} reads", alns.len());
        let mut total_counts = OperCounts::<u64>::default();
        // HashMap, values: (edit_dist, read_len), keys: number of appearances.
        let mut edit_distances = HashMap::<(u32, u32), u64>::new();
        for aln in alns.iter() {
            let aln = aln.borrow();
            if !windows.keep_window(aln.interval().middle()) {
                continue;
            }
            let counts = aln.count_region_operations(region);
            total_counts += &counts;
            let (edit_dist, read_len) = counts.edit_and_read_len();
            edit_distances.entry((edit_dist, read_len))
                .and_modify(|counter| *counter += 1)
                .or_insert(1);
        }

        let oper_probs = total_counts.to_ln_probs(params.err_rate_mult);
        let edit_distances: Vec<_> = edit_distances.into_iter()
            .map(|((k, n), count)| (k, n, count as f64))
            .collect();
        let edit_distr = BetaBinomial::max_lik_estimate(&edit_distances);

        if let Some(dir) = out_dir {
            let dbg_filename = dir.join("edit_dist.csv.gz");
            let mut dbg_writer = ext::sys::create_gzip(&dbg_filename)?;
            writeln!(dbg_writer, "edit\tsize\tcount").map_err(add_path!(dbg_filename))?;
            for (k, n, count) in edit_distances.iter() {
                writeln!(dbg_writer, "{}\t{}\t{}", k, n, count).map_err(add_path!(dbg_filename))?;
            }
            writeln!(dbg_writer, "# alpha = {}, beta = {}", edit_distr.alpha(), edit_distr.beta())
                .map_err(add_path!(dbg_filename))?;
        }

        let err_prof = Self {
            oper_probs,
            edit_distr,
            max_edit_dist: RefCell::default(),
            // Set cdf2 = cdf1 so that we do not have to go far.
            edit_cdf: (1.0 - params.edit_pval, 1.0 - params.edit_pval),
        };

        let read_len = math::round_signif(mean_read_len, 2).round() as u32;
        log::info!("    Maximum allowed edit distance: {} (for read length {}, {}%-confidence interval)",
            err_prof.allowed_edit_dist(read_len).0, read_len, 100.0 * err_prof.edit_cdf.0);
        Ok(err_prof)
    }

    pub fn set_edit_pvals(&mut self, (pval1, pval2): (f64, f64)) {
        self.edit_cdf = (1.0 - pval1, 1.0 - pval2);
    }

    /// Returns ln-probability for operation counts.
    pub fn ln_prob<T>(&self, counts: &OperCounts<T>) -> f64
    where T: Copy + TryInto<f64>,
          <T as TryInto<f64>>::Error: std::fmt::Debug,
    {
        self.oper_probs.matches * counts.matches.try_into().unwrap()
            + self.oper_probs.mismatches * counts.mismatches.try_into().unwrap()
            + self.oper_probs.insertions * counts.insertions.try_into().unwrap()
            + self.oper_probs.deletions * counts.deletions.try_into().unwrap()
            + self.oper_probs.clipping * counts.clipping.try_into().unwrap()
    }

    /// Returns the maximum allowed edit distance for the given read length.
    /// Values are cached for future use.
    pub fn allowed_edit_dist(&self, read_len: u32) -> (u32, u32) {
        let mut cache = self.max_edit_dist.borrow_mut();
        *cache.entry(read_len).or_insert_with(||
            self.edit_distr.inv_cdf2(read_len, self.edit_cdf.0, self.edit_cdf.1))
    }
}

impl JsonSer for ErrorProfile {
    fn save(&self) -> json::JsonValue {
        json::object!{
            matches: self.oper_probs.matches,
            mismatches: self.oper_probs.mismatches,
            insertions: self.oper_probs.insertions,
            deletions: self.oper_probs.deletions,
            clipping: self.oper_probs.clipping,

            alpha: self.edit_distr.alpha(),
            beta: self.edit_distr.beta(),
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> matches (as_f64), mismatches (as_f64), insertions (as_f64),
            deletions (as_f64), clipping (as_f64),
            alpha (as_f64), beta (as_f64));
        Ok(Self {
            oper_probs: OperCounts { matches, mismatches, insertions, deletions, clipping },
            edit_distr: BetaBinomial::new(alpha, beta),
            max_edit_dist: RefCell::default(),
            edit_cdf: (1.0 - DEF_EDIT_PVAL.0, 1.0 - DEF_EDIT_PVAL.1),
        })
    }
}
