use std::{
    fmt,
    cmp::{min, max},
    io::Write,
    ops::AddAssign,
    cell::RefCell,
    borrow::Borrow,
    path::Path,
};
use crate::{
    err::{add_path},
    seq::{
        aln::{NamedAlignment, OpCounter},
    },
    math::{self, distr::BetaBinomial},
    bg::ser::{JsonSer, json_get},
    algo::IntMap,
    ext,
};

/// Edit distance and read length.
/// Both values can be limited by the length of region of interest.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct EditDist {
    edit: u32,
    read_len: u32,
}

impl EditDist {
    #[inline]
    pub fn edit(&self) -> u32 {
        self.edit
    }

    #[inline]
    pub fn read_len(&self) -> u32 {
        self.read_len
    }
}

impl fmt::Display for EditDist {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}/{}", self.edit, self.read_len)
    }
}

impl std::hash::Hash for EditDist {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        hasher.write_u64((u64::from(self.edit) << 32) | u64::from(self.read_len))
    }
}

impl nohash::IsEnabled for EditDist {}

/// Counts of five operation types (=, X, I, D, S).
#[derive(Clone, Default, Debug)]
pub struct OperCounts<T> {
    pub matches: T,
    pub mismatches: T,
    pub insertions: T,
    pub deletions: T,
    pub clipping: T,
}

impl OperCounts<u32> {
    /// Returns edit distance and read length.
    pub fn edit_distance(&self) -> EditDist {
        let common = self.mismatches + self.insertions + self.clipping;
        EditDist {
            edit: common + self.deletions,
            read_len: common + self.matches,
        }
    }
}

impl OperCounts<u64> {
    /// Calculates ln-probabilities of various error types.
    /// Error probabilities (mismatches, insertions & deletions)
    /// to account for possible mutations in the data.
    ///
    /// Returned clipping probability is set to the maximum of mismatch and insert probabilities.
    pub fn to_ln_probs(&self) -> OperCounts<f64> {
        // Clipping is ignored.
        let sum_len = (self.matches + self.mismatches + self.insertions + self.deletions) as f64;
        let mism_prob = self.mismatches as f64 / sum_len;
        let ins_prob = self.insertions as f64 / sum_len;
        let del_prob = self.deletions as f64 / sum_len;
        let match_prob = 1.0 - mism_prob - ins_prob - del_prob;

        log::info!("    {:12} matches    ({:.6})", self.matches, match_prob);
        log::info!("    {:12} mismatches ({:.6})", self.mismatches, mism_prob);
        log::info!("    {:12} insertions ({:.6})", self.insertions, ins_prob);
        log::info!("    {:12} deletions  ({:.6})", self.deletions, del_prob);
        log::info!("    {:12} clippings  ({:.6})", self.clipping, self.clipping as f64 / sum_len);
        assert!(match_prob > 0.5, "Match probability ({:.5}) must be over 50%", match_prob);
        assert!((match_prob + mism_prob + ins_prob + del_prob - 1.0).abs() < 1e-8,
            "Error probabilities do not sum to one");
        OperCounts {
            matches: match_prob.ln(),
            mismatches: mism_prob.ln(),
            insertions: ins_prob.ln(),
            deletions: del_prob.ln(),
            clipping: ins_prob.max(mism_prob).ln(),
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
}

impl ErrorProfile {
    /// Create error profile from the iterator over CIGARs (with extended operations X/=).
    /// If `out_dir` is Some, write debug information.
    pub fn estimate<'a>(
        alns: &[impl Borrow<NamedAlignment>],
        op_counter: &OpCounter,
        windows: Option<&super::Windows>,
        out_dir: Option<&Path>,
    ) -> crate::Result<Self>
    {
        log::info!("Estimating read error profiles from {} reads", alns.len());
        let mut total_counts = OperCounts::<u64>::default();
        // HashMap, keys: (edit_dist, read_len), values: number of appearances.
        let mut edit_distances = IntMap::<EditDist, u64>::default();
        for aln in alns {
            let aln = aln.borrow();
            if !windows.map(|ws| ws.keep_window(aln.interval().middle())).unwrap_or(true) {
                continue;
            }
            let counts = op_counter.count(aln);
            total_counts += &counts;
            edit_distances.entry(counts.edit_distance())
                .and_modify(|counter| *counter += 1)
                .or_insert(1);
        }

        let oper_probs = total_counts.to_ln_probs();
        // Use `min(edit, read_len)` as due to deletions it is possible to get bigger edit distance than read length.
        let triples: Vec<_> = edit_distances.iter()
            .map(|(&EditDist { edit, read_len }, &count)| (min(edit, read_len), read_len, count as f64))
            .collect();

        // When finding Beta-Binomial distribution, add Uniform distribution with multiplier
        // UNIF_NOMINATOR / n_reads.
        const UNIF_NOMINATOR: f64 = 3.0;
        let unif_coef = f64::min(UNIF_NOMINATOR / alns.len() as f64, 0.1);
        let edit_distr = BetaBinomial::max_lik_estimate(&triples, unif_coef);

        if let Some(dir) = out_dir {
            let dbg_filename = dir.join("edit_dist.csv.gz");
            let mut dbg_writer = ext::sys::create_gzip(&dbg_filename)?;
            writeln!(dbg_writer, "edit\tsize\tcount").map_err(add_path!(dbg_filename))?;
            for (edit_dist, count) in edit_distances.iter() {
                writeln!(dbg_writer, "{}\t{}\t{}", edit_dist.edit, edit_dist.read_len, count)
                    .map_err(add_path!(dbg_filename))?;
            }
            writeln!(dbg_writer, "# alpha = {}, beta = {}", edit_distr.alpha(), edit_distr.beta())
                .map_err(add_path!(dbg_filename))?;
        }

        Ok(Self { oper_probs, edit_distr })
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

    fn load(obj: &json::JsonValue) -> crate::Result<Self> {
        json_get!(obj => matches (as_f64), mismatches (as_f64), insertions (as_f64),
            deletions (as_f64), clipping (as_f64),
            alpha (as_f64), beta (as_f64));
        Ok(Self {
            oper_probs: OperCounts { matches, mismatches, insertions, deletions, clipping },
            edit_distr: BetaBinomial::new(alpha, beta),
        })
    }
}

/// For each read length, stores ONE edit distance (see EditDistCache).
pub struct SingleEditDistCache {
    /// Distribution of edit distance (n = read length).
    edit_distr: BetaBinomial,
    /// 1 - pvalue.
    edit_cdf: f64,
    /// For each read length, store maximum allowed edit distance.
    cache: RefCell<IntMap<u32, u32>>,
}

impl SingleEditDistCache {
    pub fn new(err_prof: &ErrorProfile, pval: f64) -> Self {
        Self {
            edit_distr: err_prof.edit_distr.clone(),
            cache: RefCell::default(),
            edit_cdf: 1.0 - pval,
        }
    }

    pub fn get(&self, read_len: u32) -> u32 {
        *self.cache.borrow_mut().entry(read_len).or_insert_with(||
            self.edit_distr.inv_cdf(read_len, self.edit_cdf))
    }

    pub fn print_log(&self, mean_read_len: f64) {
        let read_len = math::round_signif(mean_read_len, 2).round() as u32;
        log::info!("    Maximum allowed edit distance: {} (for read length {}, {}%-confidence interval)",
            self.get(read_len), read_len, 100.0 * self.edit_cdf);
    }
}

/// For each read length, stores two edit distances: good and passable.
pub struct EditDistCache {
    /// Distribution of edit distance (n = read length).
    edit_distr: BetaBinomial,
    /// 1 - pvalue.
    edit_cdf: (f64, f64),
    /// Adds this number to the good distance (to allow for differences between alleles).
    good_dist_add: u32,
    /// For each read length, store two maximum allowed edit distances (good and passable).
    cache: RefCell<IntMap<u32, (u32, u32)>>,
}

impl EditDistCache {
    pub fn new(err_prof: &ErrorProfile, good_dist_add: u32, (pval1, pval2): (f64, f64)) -> Self {
        assert!(pval1 >= pval2);
        Self {
            edit_distr: err_prof.edit_distr.clone(),
            cache: RefCell::default(),
            good_dist_add,
            edit_cdf: (1.0 - pval1, 1.0 - pval2),
        }
    }

    /// Returns the maximum allowed edit distance for the given read length.
    /// Values are cached for future use.
    pub fn get(&self, read_len: u32) -> (u32, u32) {
        *self.cache.borrow_mut().entry(read_len).or_insert_with(|| {
            let (mut dist1, dist2) = self.edit_distr.inv_cdf2(read_len, self.edit_cdf.0, self.edit_cdf.1);
            dist1 += self.good_dist_add;
            (dist1, max(dist1, dist2))
        })
    }

    pub fn print_log(&self, mean_read_len: f64) {
        let read_len = math::round_signif(mean_read_len, 2).round() as u32;
        let (dist1, dist2) = self.get(read_len);
        log::info!("    Edit distances for read length {}: {} (good) and {} (passable)",
            read_len, dist1, dist2);
    }
}
