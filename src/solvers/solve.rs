//! Series of solvers.
//! Each next solver is expected to take more time, therefore is called on a smaller subset of haplotypes.

use std::{
    thread, fmt,
    cmp::{min, max},
    io::Write,
    time::{Instant, Duration},
    path::Path,
    sync::{
        Arc,
        mpsc::{self, Sender, Receiver},
    },
    ops::Deref,
};
use rand::{
    Rng,
    prelude::SliceRandom,
};
use crate::{
    ext::{
        self,
        TriangleMatrix,
        sys::BrFile,
        vec::{F64Ext, IterExt},
        rand::XoshiroRng,
        fmt::PrettyUsize,
        thread_files::{ThreadFile, ThreadFiles},
    },
    bg::{BgDistr, ReadDepth},
    math::{
        self, Ln, Phred, RoundDiv,
        // distr::WithMoments,
    },
    err::{error, validate_param, add_path},
    seq::{
        contigs::{ContigNames, Genotype},
        recruit::UPDATE_SECS,
    },
    model::{
        self,
        Params as AssgnParams,
        locs::AllAlignments,
        assgn::{GenotypeAlignments, ReadAssignment},
        windows::ContigInfos,
        distr_cache::DistrCache,
    },
    command::DebugLvl,
};
use super::Solver;

trait GenotypingScore {
    fn as_f64(&self) -> f64;
}

impl GenotypingScore for f64 {
    #[inline]
    fn as_f64(&self) -> f64 { *self }
}

impl GenotypingScore for Option<Prediction> {
    #[inline]
    fn as_f64(&self) -> f64 {
        self.as_ref().expect("Unavailable prediction").lik_mean
    }
}

struct CutRange {
    min_size: usize,
    max_size: usize,
    best_score: f64,
    worst_score: f64,
}

/// Indices of the currently contended genotypes.
pub struct GenotypeIxs {
    /// Retained indices (can be freely rearranged).
    retained: Vec<usize>,
    /// Index of a genotype that cannot be discarded.
    /// The same index must be inside the `retained` vector.
    fixed: Option<usize>,
}

impl GenotypeIxs {
    /// Removes fixed genotype from `ixs`, if needed.
    pub fn new(retained: Vec<usize>, fixed: Option<usize>) -> Self {
        Self { retained, fixed }
    }

    /// Iterate over all retained ixs + the fixed index.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = usize> {
        self.retained.iter().copied()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.retained.len()
    }

    /// Sorts non-fixed indices by score (highest first).
    #[inline]
    fn sort_indices(&mut self, scores: &[impl GenotypingScore]) {
        self.retained.sort_unstable_by(|&i, &j| scores[j].as_f64().total_cmp(&scores[i].as_f64()));
    }

    /// Sorts indices, writes debug messages, and returns `CutRange` containing min_size and max_size,
    /// as well as total best and total worst scores.
    fn identify_cut_range(
        &mut self,
        scores: &[impl GenotypingScore],
        genotypes: &[GenotypeInfo],
        mut threads: usize,
        mut min_size: usize,
        mut max_size: usize,
        verbose: bool,
    ) -> CutRange {
        self.sort_indices(scores);
        let n = self.retained.len();
        // Safety: we just checked the vector size.
        let first_ix = unsafe { *self.retained.get_unchecked(0) };
        let last_ix = unsafe { *self.retained.get_unchecked(n - 1) };
        let best_score = scores[first_ix].as_f64();
        let worst_score = scores[last_ix].as_f64();
        if verbose {
            log::debug!("        Worst: {} ({:.0}), best: {} ({:.0})",
                genotypes[last_ix].genotype, Ln::to_log10(worst_score),
                genotypes[first_ix].genotype, Ln::to_log10(best_score));
        }

        max_size = min(n, max_size);
        threads = threads.saturating_sub(usize::from(self.fixed.is_some()));
        min_size = max(min_size, threads).min(max_size).min(n);
        CutRange { min_size, max_size, best_score, worst_score }
    }

    /// Truncates indices at size `m` (not including fixed genotype), and writes a debug message.
    fn truncate(
        &mut self,
        mut m: usize,
        cut_range: &CutRange,
        scores: &[impl GenotypingScore],
        verbose: bool,
    ) {
        let n = self.retained.len();
        let last_score = scores[self.retained[m - 1]].as_f64();
        let opt_j = self.fixed
            .and_then(|fixed_ix| self.retained[m..].iter().position(|&ix| ix == fixed_ix))
            .map(|j| m + j);

        if verbose {
            log::debug!("        Keep {}/{} genotypes ({:.1}%), last retained: {:.0} ({:.1}% of full range){}",
                m, n, 100.0 * m as f32 / n as f32, Ln::to_log10(last_score),
                100.0 * (last_score - cut_range.worst_score) / (cut_range.best_score - cut_range.worst_score).max(1.0),
                if opt_j.is_some() { " +1 reference genotype" } else { "" });
        }

        if let Some(j) = opt_j {
            self.retained.swap(j, m);
            m += 1;
        }
        self.retained.truncate(m);
    }
}

/// Filter genotypes based on read alignments alone (without accounting for read depth).
fn run_filter(
    contigs: &ContigNames,
    genotypes: &[GenotypeInfo],
    ixs: &mut GenotypeIxs,
    all_alns: &AllAlignments,
    mut writer: Option<&mut impl Write>,
    params: &AssgnParams,
    out_size: usize,
    threads: usize,
) -> crate::Result<()>
{
    describe_stage(None);
    let contig_ids: Vec<_> = contigs.ids().collect();
    let best_aln_matrix = all_alns.best_aln_matrix(&contig_ids);
    // Vector (genotype index, score).
    let mut gt_best_probs = vec![0.0_f64; all_alns.reads().len()];
    let mut scores = vec![f64::NEG_INFINITY; genotypes.len()];
    for i in ixs.iter() {
        let info = &genotypes[i];
        let gt_ids = info.genotype.ids();
        gt_best_probs.copy_from_slice(&best_aln_matrix[gt_ids[0].ix()]);
        for id in gt_ids[1..].iter() {
            gt_best_probs.iter_mut().zip(&best_aln_matrix[id.ix()])
                .for_each(|(best_prob, &curr_prob)| *best_prob = best_prob.max(curr_prob));
        }
        let score = info.prior + gt_best_probs.iter().sum::<f64>();
        if let Some(ref mut w) = writer {
            writeln!(w, "0\t{}\t{:.3}", info.genotype, Ln::to_log10(score)).map_err(add_path!(!))?;
        }
        scores[i] = score;
    }
    let cut_range = ixs.identify_cut_range(&scores, genotypes, threads, out_size, 2 * out_size, true);
    let thresh = cut_range.best_score - params.filt_diff;
    let m = if cut_range.min_size >= cut_range.max_size && thresh <= cut_range.worst_score {
        cut_range.min_size
    } else {
        cut_range.min_size
            + ixs.retained[cut_range.min_size..cut_range.max_size].partition_point(|&i| scores[i] >= thresh)
    };
    ixs.truncate(m, &cut_range, &scores, true);
    Ok(())
}

fn describe_stage(stage: Option<&Stage>) {
    let Some(stage) = stage else {
        log::info!("*** Preliminary filtering");
        return;
    };

    const N_NUMS: usize = 10;
    const LATIN_NUMS: [&'static str; N_NUMS] = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"];
    let params = stage.solver.describe_params();
    log::info!("*** Stage {:>3}. {:<20}   [a={}{}{}]",
        if stage.ix < N_NUMS { LATIN_NUMS[stage.ix] } else { "..." },
        stage.solver.to_string(), stage.attempts, if params.is_empty() { "" } else { "," }, params);
}

struct Stage {
    ix: usize,
    solver: Box<dyn Solver>,
    /// Number of attempts per genotype.
    attempts: u16,
    /// Take approximately this number of top genotypes.
    in_size: usize,
}

impl Stage {
    /// Parse stage from string, with format
    /// SOLVER[:PARAM=VALUE,PARAM=VALUE,...]
    pub fn parse(ix: usize, s: &str) -> crate::Result<Self> {
        let (solver_name, params) = s.split_once(':').unwrap_or_else(|| (s, ""));
        let mut solver: Box<dyn Solver> = match solver_name {
            "greedy" => Box::new(super::Greedy::default()),
            "anneal" | "simanneal" | "annealing" | "simannealing"
                => Box::new(super::SimAnneal::default()),
            "highs" => {
                #[cfg(feature = "highs")]
                { Box::new(super::HighsSolver::default()) }
                #[cfg(not(feature = "highs"))]
                { return Err(error!(RuntimeError,
                    "HiGHS solver is disabled. Please recompile with `highs` feature.")) }
            }
            "gurobi" => {
                #[cfg(feature = "gurobi")]
                { Box::new(super::GurobiSolver::default()) }
                #[cfg(not(feature = "gurobi"))]
                { return Err(error!(RuntimeError,
                    "Gurobi solver is disabled. Please recompile with `gurobi` feature.")) }
            }
            _ => return Err(error!(InvalidInput, "Unknown solver {:?}", solver_name)),
        };

        let mut in_size = 1000;
        let mut attempts = 20;
        if !params.is_empty() {
            for kv in params.split(',') {
                let (key, val) = kv.split_once("=")
                    .ok_or_else(|| error!(InvalidInput, "Could not parse solver definition `{}`", s))?;
                match key {
                    "i" | "input" | "in-size" => in_size = val.parse::<PrettyUsize>()
                        .map_err(|_| error!(InvalidInput, "Could not parse `{}` in solver definition `{}`", kv, s))?
                        .0,
                    "a" | "attempts" => attempts = val.parse()
                        .map_err(|_| error!(InvalidInput, "Could not parse `{}` in solver definition `{}`", kv, s))?,
                    _ => {
                        match solver.set_param(key, val) {
                            Err(super::ParamErr::Unknown) =>
                                log::error!("Unknown parameter `{}` for solver `{}` in `{}`", key, solver_name, s),
                            Err(super::ParamErr::Parse) => return Err(error!(InvalidInput,
                                "Could not parse `{}` in `{}` (possibly, incorrect type)", kv, s)),
                            Err(super::ParamErr::Invalid(e)) => return Err(error!(InvalidInput,
                                "Invalid value of `{}` in `{}`: {}", kv, s, e)),
                            Ok(()) => {}
                        }
                    }
                }
            }
        }
        validate_param!(attempts > 0, "At least one attempt is required for each stage (`{}`)", s);
        validate_param!(in_size > 0, "At least one input genotype is required for each stage (`{}`)", s);
        Ok(Self { ix, solver, attempts, in_size })
    }
}

/// Solver scheme.
/// Consists of multiple solvers, each executed on (possibly) smaller subset of genotypes.
pub struct Scheme {
    stages: Vec<Stage>,
}

pub const DEFAULT_STAGES: &'static str = "-S greedy:i=5k,a=1 -S anneal:i=20,a=20";

impl Default for Scheme {
    fn default() -> Self {
        Self { stages: vec![
            Stage {
                ix: 0,
                solver: Box::new(super::Greedy::default()),
                attempts: 1,
                in_size: 5000,
            },
            Stage {
                ix: 1,
                solver: Box::new(super::SimAnneal::default()),
                attempts: 20,
                in_size: 20,
            },
        ]}
    }
}

impl Scheme {
    pub fn parse(solvers: &[String]) -> crate::Result<Self> {
        if solvers.is_empty() {
            return Ok(Self::default());
        }

        let stages = solvers.iter().enumerate()
            .map(|(i, s)| Stage::parse(i, s))
            .collect::<Result<Vec<_>, _>>()?;
        debug_assert!(!stages.is_empty());
        for i in 0..stages.len() - 1 {
            if stages[i].in_size <= stages[i + 1].in_size {
                log::warn!(
                    "Stage {} is redundant (input size {} not greater than next input size {}). May be skipped",
                    i + 1, stages[i].in_size, stages[i + 1].in_size);
            }
        }
        Ok(Self { stages })
    }
}

/// Genotype and associated information.
pub struct GenotypeInfo {
    pub genotype: Genotype,
    pub prior: f64,
}

impl GenotypeInfo {
    #[inline]
    pub fn new(genotype: Genotype, prior: f64) -> Self {
        Self { genotype, prior }
    }
}

/// Various structures, needed for sample genotyping.
pub struct Data {
    /// Solving scheme.
    pub scheme: Arc<Scheme>,
    pub contigs: Arc<ContigNames>,
    pub contig_distances: Option<TriangleMatrix<Option<u32>>>,
    pub true_edit_distances: bool,
    /// All read alignments, groupped by contigs.
    pub all_alns: AllAlignments,
    pub contig_infos: Arc<ContigInfos>,
    pub distr_cache: Arc<DistrCache>,

    /// Genotypes and their priors.
    pub genotypes: Vec<GenotypeInfo>,
    /// Index of a fixed genotype.
    pub fixed_gt_ix: Option<usize>,

    pub assgn_params: AssgnParams,
    pub debug: DebugLvl,
    pub threads: usize,
    pub is_paired_end: bool,
}

// const ALNS_CSV_HEADER: &'static str = "read_hash\tcontig\tpos1\tpos2\tlik\tselected";

// /// Write reads and their assignments to a CSV file in the following format (tab-separated):
// /// `prefix  read_hash  aln1  aln2  w1  w2  log10-prob  selected`
// fn write_alns(
//     f: &mut impl io::Write,
//     prefix: &str,
//     gt_alns: &GenotypeAlignments,
//     all_alns: &AllAlignments,
//     assgn_counts: &[u16],
//     attempts: u16,
// ) -> io::Result<()> {
//     let attempts = f32::from(attempts);
//     let mut counts_iter = assgn_counts.iter();
//     for (rp, groupped_alns) in all_alns.reads().iter().enumerate() {
//         let hash = groupped_alns.read_data().name_hash();
//         for curr_windows in gt_alns.possible_read_alns(rp) {
//             write!(f, "{}{}\t", prefix, hash)?;
//             match curr_windows.parent() {
//                 None => write!(f, "*\t*\t*\t")?,
//                 Some((contig_ix, aln_pair)) => {
//                     write!(f, "{}\t", contig_ix + 1)?;
//                     match aln_pair.ix1() {
//                         Some(i) => write!(f, "{}\t", groupped_alns.ith_aln(i).interval().start() + 1)?,
//                         None => write!(f, "*\t")?,
//                     };
//                     match aln_pair.ix2() {
//                         Some(i) => write!(f, "{}\t", groupped_alns.ith_aln(i).interval().start() + 1)?,
//                         None => write!(f, "*\t")?,
//                     };
//                 }
//             }
//             let prob = Ln::to_log10(curr_windows.ln_prob());
//             let count = *counts_iter.next().expect("Not enough assignment counts");
//             writeln!(f, "{:.3}\t{:.2}", prob, f32::from(count) / attempts)?;
//         }
//     }
//     assert!(counts_iter.next().is_none(), "Too many assignment counts");
//     Ok(())
// }

/// Returns probability that the first mean is larger than the second.
/// If variances are well defined, returns ln(p-value) of the t-test.
/// Otherwise, returns normalized first log-probability.
fn compare_two_likelihoods(pred1: &Prediction, pred2: &Prediction) -> f64 {
    const DIFF_VAR: bool = false;
    // Simple normalization.
    let simple_norm = pred1.lik_mean - Ln::add(pred1.lik_mean, pred2.lik_mean);
    if pred1.lik_var.is_normal() && pred2.lik_var.is_normal() {
        let t_pval = if pred1.attempts == pred2.attempts {
            math::unpaired_onesided_t_test::<DIFF_VAR>(
                pred1.lik_mean, pred1.lik_var, pred2.lik_mean, pred2.lik_var, f64::from(pred1.attempts))
        } else {
            math::unpaired_onesided_t_test_diffsizes::<DIFF_VAR>(
                pred1.lik_mean, pred1.lik_var, pred2.lik_mean, pred2.lik_var,
                f64::from(pred1.attempts), f64::from(pred2.attempts))
        };
        f64::max(simple_norm, t_pval.ln())
    } else {
        simple_norm
    }
}

/// Calculates distance between two genotypes as minimum sum distance between permutations of contigs within genotypes.
fn genotype_distance(gt1: &Genotype, gt2: &Genotype, distances: &TriangleMatrix<Option<u32>>) -> Option<u32> {
    let mut min_dist = u32::MAX;
    ext::vec::gen_permutations(gt1.ids(), |perm_ids1| {
        let mut dist = 0;
        for (i, j) in itertools::izip!(perm_ids1, gt2.ids()) {
            if i != j {
                match *distances.get_symmetric(i.ix(), j.ix()) {
                    Some(d) => dist += d,
                    None => {
                        dist = u32::MAX;
                        break;
                    }
                }
            }
        }
        min_dist = min(min_dist, dist);
    });
    (min_dist != u32::MAX).then_some(min_dist)
}

/// Prediction for a single genotype.
#[derive(Clone)]
pub struct Prediction {
    lik_mean: f64,
    lik_var: f64,
    attempts: u16,
    assgn_counts: Option<Vec<u16>>,
}

impl Prediction {
    fn new(lik_mean: f64, lik_var: f64, attempts: u16, assgn_counts: Vec<u16>) -> Self {
        Self {
            lik_mean, lik_var, attempts,
            assgn_counts: Some(assgn_counts),
        }
    }

    /// Clone without assignment counts.
    fn cheap_clone(&self) -> Self {
        Self {
            lik_mean: self.lik_mean,
            lik_var: self.lik_var,
            attempts: self.attempts,
            assgn_counts: None,
        }
    }

    /// Number of assignment attempts.
    #[allow(unused)]
    pub fn attempts(&self) -> u16 {
        self.attempts
    }

    /// How many times each read assignment was observed.
    #[allow(unused)]
    pub fn assgn_counts(&self) -> &[u16] {
        self.assgn_counts.as_ref().expect("Assignment counts undefined")
    }
}

/// Compare up to 100 genotypes with the best one using
const SOPHISTICATED_COUNT: usize = 100;

struct Predictions {
    predictions: Vec<Option<Prediction>>,
    /// Remaining indices after the latest iteration.
    ixs: GenotypeIxs,
}

impl Predictions {
    fn new(n_gts: usize, fixed: Option<usize>) -> Self {
        Self {
            predictions: vec![None; n_gts],
            ixs: GenotypeIxs::new((0..n_gts).collect(), fixed),
        }
    }

    /// Discard improbable genotypes. At least `min_size` genotypes will be preserved,
    /// genotypes between min_size and max_size will be compared with the best genotype and kept only if necessary.
    ///
    /// If out_ln_probs is provided, writes comparison log-likelihoods there.
    fn discard_improbable_genotypes(
        &mut self,
        genotypes: &[GenotypeInfo],
        prob_thresh: f64,
        min_size: usize,
        mut max_size: usize,
        threads: usize,
        verbose: bool,
        mut out_ln_probs: Option<&mut Vec<f64>>,
    ) {
        let n = self.ixs.len();
        if prob_thresh == f64::NEG_INFINITY || min_size >= n {
            if let Some(ref mut v) = out_ln_probs {
                let best = self.predictions[self.ixs.retained[0]].as_ref().expect("Unavailable prediction");
                v.push(0.0);
                v.extend(self.ixs.retained[1..].iter().map(|&ix| compare_two_likelihoods(
                    self.predictions[ix].as_ref().expect("Unavailable prediction"), best)));
            }
            return;
        }

        max_size = max(min_size, max_size);
        let cut_range = self.ixs.identify_cut_range(&self.predictions, genotypes, threads, min_size, max_size, verbose);
        let mut cut_pos = cut_range.min_size;

        let best = self.predictions[self.ixs.retained[0]].as_ref().expect("Unavailable prediction").cheap_clone();
        if let Some(ref mut v) = out_ln_probs {
            v.push(0.0);
            v.extend(self.ixs.retained[1..cut_range.min_size].iter().map(|&ix| compare_two_likelihoods(
                self.predictions[ix].as_ref().expect("Unavailable prediction"), &best)));
        }

        // Stop when encountered > 5 unlikely genotypes.
        const MAX_DROPPED: u8 = 5;
        let mut dropped_counter = 0..MAX_DROPPED;
        let mut end_pos = self.ixs.retained.len();
        for _ in cut_range.min_size..cut_range.max_size {
            // This should not happen since we do `max_size - min_size` iterations,
            // at each iteration gap `end_pos - cut_pos` decreases by one, and at the start
            // `end_pos - cut_pos >= max_size - min_size`.
            debug_assert!(cut_pos < end_pos);
            let ix = self.ixs.retained[cut_pos];
            let pred = self.predictions[ix].as_ref().expect("Unavailable prediction");
            let ln_prob = compare_two_likelihoods(pred, &best);
            if ln_prob < prob_thresh {
                if dropped_counter.next().is_none() { break };
                end_pos -= 1;
                self.ixs.retained.swap(cut_pos, end_pos);
            } else {
                cut_pos += 1;
                if let Some(ref mut v) = out_ln_probs {
                    v.push(ln_prob);
                }
            }
        }

        // Remove dropped assignment counts.
        self.ixs.retained[cut_pos..].iter()
            .filter(|&&ix| Some(ix) != self.ixs.fixed)
            .for_each(|&ix| self.predictions[ix].as_mut().unwrap().assgn_counts = None);
        self.ixs.truncate(cut_pos, &cut_range, &self.predictions, verbose);
        if let Some(ref mut v) = out_ln_probs && self.ixs.retained.len() == cut_pos + 1 {
            // One fixed genotype remains.
            let last_ix = *self.ixs.retained.last().unwrap();
            let fixed_pred = self.predictions[last_ix].as_ref().expect("Unavailable prediction");
            v.push(compare_two_likelihoods(&fixed_pred, &best));
        }
    }

    fn produce_result(mut self, data: &Data) -> Genotyping {
        // ln(1e-5).
        const THRESH: f64 = -11.512925464970229;
        // Output at most this number of genotypes.
        const MIN_GENOTYPES: usize = 10;
        const MAX_GENOTYPES: usize = 100;

        let params = &data.assgn_params;
        let genotypes = &data.genotypes;
        let mut ln_probs = Vec::new();
        self.discard_improbable_genotypes(genotypes, THRESH.min(params.prob_thresh),
            max(MIN_GENOTYPES, params.out_bams), MAX_GENOTYPES, 1, false, Some(&mut ln_probs));
        let n = self.ixs.len();
        assert_eq!(ln_probs.len(), n);

        let quality = if n < 2 {
            log::warn!("Single genotype remaining, quality will be undefined");
            f64::NAN
        } else {
            let norm_fct = Ln::sum(&ln_probs);
            ln_probs.iter_mut().for_each(|v| *v -= norm_fct);
            Phred::from_ln_prob(Ln::sum(&ln_probs[1..])).min(1e9)
        };

        let mut write_n = MIN_GENOTYPES;
        if n > write_n && let Some(j) = self.ixs.fixed
                .and_then(|fixed_ix| self.ixs.retained[write_n..].iter().position(|&ix| ix == fixed_ix)) {
            self.ixs.retained.swap(write_n, j + write_n);
            ln_probs.swap(write_n, j + write_n);
            write_n += 1;
        }
        ln_probs.truncate(write_n);
        self.ixs.retained.truncate(write_n);
        Genotyping {
            genotypes: self.ixs.retained.iter().map(|&ix| genotypes[ix].genotype.clone()).collect(),
            predictions: self.ixs.retained.iter().map(|&ix| self.predictions[ix].take().unwrap()).collect(),
            warnings: Vec::new(),
            weighted_dist: None,
            total_reads: data.all_alns.reads().len() as u32,
            distances: None,
            true_edit_distances: data.true_edit_distances,
            unexpl_reads: None,
            tag: data.contigs.tag().to_owned(),
            ln_probs, quality,
        }
    }
}

#[derive(PartialEq, Clone, Copy)]
pub enum GenotypingWarning {
    /// There are very no reads available for genotyping.
    NoReads,
    /// There are some reads, but fewer than expected (number of reads).
    FewReads(u32),
    /// There are more reads than expected (number of reads).
    #[allow(unused)]
    TooManyReads(u32),
    /// Even the best genotype has very low quality.
    NoProbableGenotype,
}

impl fmt::Display for GenotypingWarning {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::NoReads => write!(f, "NoReads"),
            Self::FewReads(count) => write!(f, "FewReads({})", count),
            Self::TooManyReads(count) => write!(f, "TooManyReads({})", count),
            Self::NoProbableGenotype => write!(f, "NoProbableGenotype"),
        }
    }
}

impl Into<json::JsonValue> for GenotypingWarning {
    fn into(self) -> json::JsonValue {
        json::JsonValue::String(self.to_string())
    }
}

/// Genotyping result.
pub struct Genotyping {
    tag: String,
    /// Filtered list with the best genotypes.
    /// Genotypes are sorted from best to worst.
    genotypes: Vec<Genotype>,
    /// Normalized ln-probabilities for each genotype.
    ln_probs: Vec<f64>,
    /// Predictions (likelihood and read counts) for each of the genotypes,
    predictions: Vec<Prediction>,
    /// Weighted distance from the best to all other genotypes.
    weighted_dist: Option<f64>,
    /// Distance from the secondary to the primary genotype prediction.
    distances: Option<Vec<Option<u32>>>,
    true_edit_distances: bool,
    total_reads: u32,
    /// Fraction of unexplained reads.
    unexpl_reads: Option<u32>,
    /// Quality of the top genotype.
    quality: f64,
    /// Warnings, issued for this sample and this locus.
    warnings: Vec<GenotypingWarning>,
}

impl Genotyping {
    /// Creates empty result.
    pub fn empty_result(tag: String, warnings: Vec<GenotypingWarning>) -> Self {
        Self {
            tag, warnings,
            genotypes: Vec::new(),
            predictions: Vec::new(),
            ln_probs: Vec::new(),
            weighted_dist: None,
            distances: None,
            true_edit_distances: false,
            total_reads: 0,
            unexpl_reads: None,
            quality: 0.0,
        }
    }

    pub fn print_log(&self) {
        let pred0 = &self.predictions[0];
        log::info!("    [{}] GT = {},  lik. = {:.2} ± {:.2}, confidence = {:.4}%, weight.dist = {:.4}",
            self.tag, self.genotypes[0], Ln::to_log10(pred0.lik_mean), Ln::to_log10(pred0.lik_var),
            100.0 * self.ln_probs[0].exp(), self.weighted_dist.unwrap_or(f64::NAN));
        if !self.warnings.is_empty() {
            let warns_str = self.warnings.iter().map(GenotypingWarning::to_string).collect::<Vec<_>>().join(", ");
            log::debug!("        Warnings: {}", warns_str);
        }
    }

    /// Calculate weighted distance to all other genotypes.
    pub fn find_weighted_dist(&mut self, dist_matrix: &TriangleMatrix<Option<u32>>) {
        let Some(gt0) = self.genotypes.first() else { return };
        let mut sum_prob = 0.0;
        let mut sum_dist = Some(0.0);
        let mut distances = Vec::with_capacity(self.genotypes.len());
        for (i, (gt, &ln_prob)) in self.genotypes.iter().zip(&self.ln_probs).enumerate() {
            let prob = ln_prob.exp();
            sum_prob += prob;
            let dist = if i > 0 { genotype_distance(gt0, gt, dist_matrix) } else { Some(0) };
            sum_dist = sum_dist.zip(dist).map(|(s, d)| s + prob * f64::from(d));
            distances.push(dist);
        }
        self.weighted_dist = sum_dist.map(|s| s / sum_prob);
        self.distances = Some(distances);
    }

    pub fn check_first_prob(&mut self) {
        let ln_prob0 = *self.ln_probs.first().expect("Expected at least one genotype");
        // < 0.01
        if ln_prob0.is_nan() || ln_prob0 < -2.0 * Ln::LN10 {
            self.warnings.push(GenotypingWarning::NoProbableGenotype);
            log::warn!("[{}] Best genotype {} is improbable (10^{:.5} = {:.5})",
                self.tag, self.genotypes[0], Ln::to_log10(ln_prob0), ln_prob0.exp());
        }
    }

    /// Tests the total number of reads and issues a warning if the number of reads
    /// is significantly smaller than expected.
    pub fn check_num_of_reads(
        &mut self,
        n_reads: u32,
        #[allow(unused_variables)]
        bg_depth: &ReadDepth,
        #[allow(unused_variables)]
        contig_infos: &ContigInfos,
    ) {
        let gt = &self.genotypes[0];
        let ploidy = gt.ploidy() as u32;
        if n_reads < ploidy {
            log::warn!("[{}] Impossible to find {}-ploid genotype from {} read(s)", self.tag, ploidy, n_reads);
            self.warnings.push(GenotypingWarning::FewReads(n_reads));
            return;
        } else if ploidy > 1 && n_reads < ploidy * 10 {
            let k = f64::from(ploidy);
            let n = f64::from(n_reads);
            // Expected number of zeros in the multinomial distribution with k possibilities, all probabilities = 1/k,
            // and n observations: (k-1)^n / k^(n-1).
            let exp_zeros = ((k - 1.0).ln() * n - k.ln() * (n - 1.0)).exp();
            if exp_zeros > 0.1 {
                log::warn!("[{}] Too few reads ({}) to find {}-ploid genotype: \
                    expected number of unobserved contigs = {:.4}", self.tag, ploidy, n_reads, exp_zeros);
                self.warnings.push(GenotypingWarning::FewReads(n_reads));
                return;
            } else if exp_zeros > 0.01 {
                log::debug!("    [{}] Possibly too few reads ({}) to find {}-ploid genotype: \
                    expected number of unobserved contigs = {:.4}. Continuing", self.tag, ploidy, n_reads, exp_zeros);
            }
        }

        // TODO: Refine this comparison.
        // let mut lbound = 0.0;
        // let mut ubound = 0.0;
        // for contig_id in gt.ids() {
        //     let contig_info = &contig_infos[contig_id.ix()];
        //     for (start, end) in contig_info.default_windows() {
        //         let window_chars = contig_info.window_characteristics(start, end);
        //         // Note: number of reads is actually the number of read pairs.
        //         // However, bg_depth also only counts the number of read pairs.
        //         let m = bg_depth.depth_distribution(window_chars.gc_content).mean();
        //         lbound += window_chars.weight * m;
        //         ubound += m;
        //     }
        // }

        // // Simple heuristic:
        // //    at <= 10 reads:   lbound *= 0;   ubound *= 10,
        // //    at >= 1000 reads: lbound *= 0.5; ubound *=  2.
        // const X1: f64 =   10.0;
        // const X2: f64 = 1000.0;
        // const LY1: f64 = 0.0;
        // const LY2: f64 = 0.5;
        // const UY1: f64 = 10.0;
        // const UY2: f64 =  2.0;
        // let lbound2 = math::interpolate((X1, X2), (LY1, LY2), lbound.clamp(X1, X2)) * lbound;
        // let ubound2 = math::interpolate((X1, X2), (UY1, UY2), ubound.clamp(X1, X2)) * ubound.max(1.5);

        // if f64::from(n_reads) < lbound2 {
        //     log::warn!("[{}] There are fewer reads ({}) than expected (≥ {:.2} for {})",
        //         self.tag, n_reads, lbound, gt);
        //     self.warnings.push(GenotypingWarning::FewReads(n_reads));
        // } else if f64::from(n_reads) > ubound2 {
        //     log::warn!("[{}] There are much more reads ({}) than expected ({:.2} for {})",
        //         self.tag, n_reads, ubound, gt);
        //     self.warnings.push(GenotypingWarning::TooManyReads(n_reads));
        // }
    }

    /// Checks reads that are much less likely at some other contig compared to the predicted genotype.
    pub fn count_unexplained_reads(&mut self, all_alns: &AllAlignments) {
        // Correct for rounding errors.
        const EPS: f64 = 1e-8;
        let mut unexplained = 0_u32;
        let gt = &self.genotypes[0];
        for read in all_alns.reads() {
            let best_lik = IterExt::max(gt.ids().iter().map(|&id| read.best_at_contig(id)));
            unexplained += u32::from(best_lik < read.unmapped_prob() + EPS);
        }
        self.unexpl_reads = Some(unexplained);
    }

    /// Converts genotyping result into JSON format.
    pub fn to_json(&self) -> json::JsonValue {
        let mut res = json::object! {
            total_reads: self.total_reads,
            quality: self.quality,
        };
        if self.distances.is_some() {
            res.insert("dist_type", if self.true_edit_distances { "edit" } else { "minim-div" }).unwrap();
        }
        if let Some(dist) = self.weighted_dist {
            res.insert("weight_dist", dist).unwrap();
        }
        if let Some(n) = self.unexpl_reads {
            res.insert("unexpl_reads", n).unwrap();
        }

        if !self.genotypes.is_empty() {
            res.insert("genotype", self.genotypes[0].name()).unwrap();
            let options: Vec<_> = self.genotypes.iter().enumerate().map(|(i, gt)| {
                let ln_prob = self.ln_probs[i];
                let pred = &self.predictions[i];
                let mut obj = json::object! {
                    genotype: gt.name(),
                    lik_mean: Ln::to_log10(pred.lik_mean),
                    lik_sd: Ln::to_log10(pred.lik_var),
                    prob: ln_prob.exp(),
                    log10_prob: Ln::to_log10(ln_prob),
                };
                if let Some(distances) = &self.distances {
                    match distances[i] {
                        Some(d) => obj.insert("dist_to_primary", d).unwrap(),
                        None => obj.insert("dist_to_primary", "unknown").unwrap(),
                    }
                }
                obj
            }).collect();
            res.insert("options", options).unwrap();
        }
        if !self.warnings.is_empty() {
            res.insert::<&[GenotypingWarning]>("warnings", &self.warnings).unwrap();
        }
        res
    }

    /// Top predicted genotypes (from best to worst).
    #[allow(unused)]
    pub fn genotypes(&self) -> &[Genotype] {
        &self.genotypes
    }

    /// More detailed predictions for top genotypes (same order as .genotypes).
    #[allow(unused)]
    pub fn predictions(&self) -> &[Prediction] {
        &self.predictions
    }
}

/// Returns: vector of likelihood means and variances (same order ).
fn solve_single_thread(
    data: &Data,
    rng: &mut XoshiroRng,
    predictions: &mut Predictions,
    sol_writer: &mut BrFile,
    mut depth_writer: Option<ThreadFile<BrFile>>,
    mut ext_sol_writer: Option<ThreadFile<BrFile>>,
) -> crate::Result<()>
{
    const ONE_THREAD: usize = 1;
    let total_genotypes = data.genotypes.len();
    let mut logger = Logger::new(total_genotypes);

    let distr_cache = data.distr_cache.deref();
    for stage in &data.scheme.stages {
        describe_stage(Some(stage));
        let n_gts = predictions.ixs.len();
        let out_size = data.scheme.stages.get(stage.ix + 1).map(|next_stage| next_stage.in_size);
        if !(data.assgn_params.dont_skip || out_size.map(|s| s < n_gts).unwrap_or(true)) {
            log::info!("    Skipping stage, not enough genotypes");
            continue;
        }

        let mut curr_sol_writer: Option<&mut BrFile> = (data.debug != DebugLvl::None || out_size.is_none())
            .then_some(sol_writer);

        logger.reset(n_gts);
        let mut liks = vec![f64::NAN; usize::from(stage.attempts)];
        for ix in predictions.ixs.iter() {
            let info = &data.genotypes[ix];
            let prefix = format!("{}\t{}\t", stage.ix + 1, info.genotype);
            let mut gt_alns = GenotypeAlignments::new(info.genotype.clone(), &data.contig_infos, &data.all_alns,
                &data.assgn_params);
            let mut counts = gt_alns.create_counts();

            for (attempt, lik) in liks.iter_mut().enumerate() {
                gt_alns.apply_tweak(rng, distr_cache, &data.assgn_params);
                let assgns = stage.solver.solve(&gt_alns, rng)?;
                *lik = info.prior + assgns.likelihood();
                let ext_prefix = format!("{}{}\t", prefix, attempt + 1);
                if let Some(w) = &mut depth_writer {
                    assgns.write_depth(w, &ext_prefix).map_err(add_path!(!))?;
                }
                if let Some(w) = &mut ext_sol_writer {
                    assgns.summarize(w, &ext_prefix).map_err(add_path!(!))?;
                }
                assgns.update_counts(&mut counts);
            }
            let (lik_mean, lik_var) = F64Ext::mean_variance_or_nan(&liks);
            logger.update(&info.genotype, lik_mean);
            predictions.predictions[ix] = Some(Prediction::new(lik_mean, lik_var, stage.attempts, counts));
            if let Some(ref mut w) = curr_sol_writer {
                writeln!(w, "{}{:.4}", prefix, Ln::to_log10(lik_mean)).map_err(add_path!(!))?;
            }
        }
        if let Some(min_size) = out_size {
            predictions.discard_improbable_genotypes(&data.genotypes, data.assgn_params.prob_thresh,
                min_size, min(2 * min_size, SOPHISTICATED_COUNT), ONE_THREAD, true, None);
        }
        logger.finish_stage();
    }
    Ok(())
}

/// Creates debug files for read depth and for extended solution information.
fn create_debug_files(
    locus_dir: &Path,
    threads: usize,
    debug: DebugLvl,
) -> (ThreadFiles, ThreadFiles) {
    let depth_writers = if debug == DebugLvl::Full {
        ThreadFiles::new(locus_dir.join("depth.csv.br"), threads,
            Some(format!("stage\tgenotype\tattempt\t{}", ReadAssignment::DEPTH_CSV_HEADER)))
    } else {
        ThreadFiles::disabled(threads)
    };
    let ext_sol_writers = if debug >= DebugLvl::Some {
        ThreadFiles::new(locus_dir.join("sol_ext.csv.br"), threads,
            Some(format!("stage\tgenotype\tattempt\t{}", ReadAssignment::SUMMARY_HEADER)))
    } else {
        ThreadFiles::disabled(threads)
    };
    (depth_writers, ext_sol_writers)
}

pub fn solve(
    data: &Arc<Data>,
    bg_distr: &BgDistr,
    locus_dir: &Path,
    rng: &mut XoshiroRng,
) -> crate::Result<Genotyping> {
    let n_gts = data.genotypes.len();
    assert!(n_gts > 0);
    log::info!("=== Genotyping {}: {} possible genotypes", data.contigs.tag(), n_gts);

    let mut sol_writer = ext::sys::create_brotli(&locus_dir.join("sol.csv.br"))?;
    writeln!(sol_writer, "stage\tgenotype\tscore").map_err(add_path!(!))?;
    let mut predictions = Predictions::new(n_gts, data.fixed_gt_ix);
    let out_size0 = data.scheme.stages[0].in_size;
    if data.assgn_params.dont_skip || out_size0 < n_gts {
        let curr_sol_writer = (data.debug != DebugLvl::None).then_some(&mut sol_writer);
        run_filter(&data.contigs, &data.genotypes, &mut predictions.ixs, &data.all_alns, curr_sol_writer,
            &data.assgn_params, out_size0, data.threads)?;
    }

    let (mut depth_writers, mut ext_sol_writers) = create_debug_files(locus_dir, data.threads, data.debug);
    if data.threads == 1 {
        solve_single_thread(&data, rng, &mut predictions, &mut sol_writer,
            depth_writers.create_optional_next()?, ext_sol_writers.create_optional_next()?)?;
    } else {
        let main_worker = MainWorker::new(Arc::clone(data), rng, &mut depth_writers, &mut ext_sol_writers)?;
        main_worker.run(&mut sol_writer, rng, &mut predictions)?;
    }
    depth_writers.merge()?;
    ext_sol_writers.merge()?;
    std::mem::drop(sol_writer);

    let data = &data;
    let mut genotyping = predictions.produce_result(&data);
    if data.assgn_params.out_bams > 0 {
        let bam_dir = locus_dir.join(crate::command::paths::ALNS_DIR);
        ext::sys::mkdir(&bam_dir)?;
        let mut contig_to_tid = Vec::new();
        log::info!("    Writing output alignments to {}", ext::fmt::path(&bam_dir));
        let out_gts = data.assgn_params.out_bams.min(genotyping.genotypes.len());
        let width = max(2, math::num_digits(out_gts as f64) as usize);
        for (i, (gt, pred)) in genotyping.genotypes.iter().zip(&genotyping.predictions)
                .take(data.assgn_params.out_bams).enumerate() {
            let bam_path = bam_dir.join(format!("{:0width$}.bam", i));
            model::bam::write_bam(&bam_path, gt, data, &mut contig_to_tid, pred.attempts,
                pred.assgn_counts.as_ref().expect("Assignment counts undefined"))?
        }
    }
    if let Some(dist_matrix) = &data.contig_distances {
        genotyping.find_weighted_dist(dist_matrix);
    }
    genotyping.check_first_prob();
    genotyping.check_num_of_reads(data.all_alns.reads().len() as u32, bg_distr.depth(), &data.contig_infos);
    genotyping.count_unexplained_reads(&data.all_alns);
    genotyping.print_log();
    Ok(genotyping)
}

/// Task, sent to the workers: stage index and a vector of genotype indices.
type Task = (usize, Vec<usize>);
/// Genotype index and calculated likelihood mean and variance.
type Solution = (usize, Prediction);

struct MainWorker {
    data: Arc<Data>,
    senders: Vec<Sender<Task>>,
    receivers: Vec<Receiver<Solution>>,
    handles: Vec<thread::JoinHandle<crate::Result<()>>>,
}

impl MainWorker {
    fn new(
        data: Arc<Data>,
        rng: &mut XoshiroRng,
        depth_writers: &mut ThreadFiles,
        ext_sol_writers: &mut ThreadFiles,
    ) -> crate::Result<Self>
    {
        let n_workers = data.threads;
        let mut senders = Vec::with_capacity(n_workers);
        let mut receivers = Vec::with_capacity(n_workers);
        let mut handles = Vec::with_capacity(n_workers);

        for _ in 0..data.threads {
            let (task_sender, task_receiver) = mpsc::channel();
            let (sol_sender, sol_receiver) = mpsc::channel();
            let worker = Worker {
                data: Arc::clone(&data),
                rng: rng.clone(),
                receiver: task_receiver,
                sender: sol_sender,
                depth_writer: depth_writers.create_optional_next()?,
                ext_sol_writer: ext_sol_writers.create_optional_next()?,
            };
            rng.jump();
            senders.push(task_sender);
            receivers.push(sol_receiver);
            handles.push(thread::spawn(|| worker.run()));
        }
        Ok(MainWorker { data, senders, receivers, handles })
    }

    fn run(self,
        sol_writer: &mut BrFile,
        rng: &mut impl Rng,
        predictions: &mut Predictions,
    ) -> crate::Result<()>
    {
        let n_workers = self.handles.len();
        let data = self.data.deref();
        let scheme = data.scheme.deref();
        let total_genotypes = data.genotypes.len();
        let mut logger = Logger::new(total_genotypes);
        let mut rem_jobs = Vec::with_capacity(n_workers);

        for stage in &scheme.stages {
            describe_stage(Some(stage));
            let n_gts = predictions.ixs.len();
            let out_size = data.scheme.stages.get(stage.ix + 1).map(|next_stage| next_stage.in_size);
            if !(data.assgn_params.dont_skip || out_size.map(|s| s < n_gts).unwrap_or(true)) {
                log::info!("    Skipping stage, not enough genotypes");
                continue;
            }

            logger.reset(n_gts);
            rem_jobs.clear();
            let mut start = 0;
            // Shuffle indices to better distribute them to threads.
            predictions.ixs.retained.shuffle(rng);
            for (i, sender) in self.senders.iter().enumerate() {
                if start == n_gts {
                    break;
                }
                let rem_workers = n_workers - i;
                let curr_jobs = (n_gts - start).fast_ceil_div(rem_workers);
                let task = predictions.ixs.retained[start..start + curr_jobs].to_vec();
                sender.send((stage.ix, task)).expect("Genotyping worker has failed!");
                start += curr_jobs;
                rem_jobs.push(curr_jobs);
            }
            assert_eq!(start, n_gts);

            let mut curr_sol_writer: Option<&mut BrFile> = (data.debug != DebugLvl::None || out_size.is_none())
                .then_some(sol_writer);
            while logger.solved_genotypes < n_gts {
                for (receiver, jobs) in self.receivers.iter().zip(rem_jobs.iter_mut()) {
                    if *jobs > 0 {
                        let (ix, pred) = receiver.recv().expect("Genotyping worker has failed!");
                        let info = &data.genotypes[ix];
                        logger.update(&info.genotype, pred.lik_mean);
                        if let Some(ref mut w) = curr_sol_writer {
                            writeln!(w, "{}\t{}\t{:.4}", stage.ix + 1, info.genotype, Ln::to_log10(pred.lik_mean))
                                .map_err(add_path!(!))?;
                        }
                        predictions.predictions[ix] = Some(pred);
                        *jobs -= 1;
                    }
                }
            }
            if let Some(min_size) = out_size {
                predictions.discard_improbable_genotypes(&data.genotypes, data.assgn_params.prob_thresh,
                    min_size, min(2 * min_size, SOPHISTICATED_COUNT), data.threads, true, None);
            }
            logger.finish_stage();
        }
        std::mem::drop(self.senders);
        self.handles.into_iter()
            .flat_map(|handle| handle.join())
            .collect()
    }
}

struct Worker {
    data: Arc<Data>,
    rng: XoshiroRng,
    receiver: Receiver<Task>,
    sender: Sender<Solution>,
    depth_writer: Option<ThreadFile<BrFile>>,
    ext_sol_writer: Option<ThreadFile<BrFile>>,
}

impl Worker {
    fn run(mut self) -> crate::Result<()> {
        let data = self.data.deref();
        let scheme = data.scheme.deref();
        let distr_cache = data.distr_cache.deref();

        // Block thread and wait for the shipment.
        while let Ok((stage_ix, task)) = self.receiver.recv() {
            let stage = &scheme.stages[stage_ix];
            let n = task.len();
            assert_ne!(n, 0, "Received empty task");

            let mut liks = vec![f64::NAN; usize::from(stage.attempts)];
            for ix in task.into_iter() {
                let info = &data.genotypes[ix];
                let prefix = format!("{}\t{}\t", stage_ix + 1, info.genotype);
                let mut gt_alns = GenotypeAlignments::new(info.genotype.clone(), &data.contig_infos, &data.all_alns,
                    &data.assgn_params);
                let mut counts = gt_alns.create_counts();
                for (attempt, lik) in liks.iter_mut().enumerate() {
                    gt_alns.apply_tweak(&mut self.rng, distr_cache, &data.assgn_params);
                    let assgns = stage.solver.solve(&gt_alns, &mut self.rng)?;
                    *lik = info.prior + assgns.likelihood();
                    let ext_prefix = format!("{}{}\t", prefix, attempt + 1);
                    if let Some(w) = &mut self.depth_writer {
                        assgns.write_depth(w, &ext_prefix).map_err(add_path!(!))?;
                    }
                    if let Some(w) = &mut self.ext_sol_writer {
                        assgns.summarize(w, &ext_prefix).map_err(add_path!(!))?;
                    }
                    assgns.update_counts(&mut counts);
                }
                let (lik_mean, lik_var) = F64Ext::mean_variance_or_nan(&liks);
                let pred = Prediction::new(lik_mean, lik_var, stage.attempts, counts);
                if let Err(_) = self.sender.send((ix, pred)) {
                    log::error!("Read recruitment: main thread stopped before the child thread.");
                    break;
                }
            }
        }
        Ok(())
    }
}

struct Logger {
    solved_genotypes: usize,
    curr_genotypes: usize,
    num_width: usize,

    best_lik: f64,
    best_str: String,

    timer: Instant,
    stage_start: Duration,
    last_msg: Duration,
}

impl Logger {
    fn new(total_genotypes: usize) -> Self {
        Self {
            solved_genotypes: 0,
            curr_genotypes: total_genotypes,
            num_width: math::num_digits(total_genotypes as f64) as usize,

            best_lik: f64::NEG_INFINITY,
            best_str: String::new(),

            timer: Instant::now(),
            stage_start: Duration::default(),
            last_msg: Duration::default(),
        }
    }

    fn reset(&mut self, curr_genotypes: usize) {
        self.stage_start = self.timer.elapsed();
        assert!(curr_genotypes <= self.curr_genotypes);
        if curr_genotypes < self.curr_genotypes {
            self.curr_genotypes = curr_genotypes;
            self.num_width = math::num_digits(self.curr_genotypes as f64) as usize;
        }
        self.solved_genotypes = 0;
        self.last_msg = self.stage_start;
    }

    fn update(&mut self, gt: &Genotype, lik: f64) {
        if lik > self.best_lik {
            self.best_lik = lik;
            self.best_str = gt.name().to_owned();
        }
        self.solved_genotypes += 1;
        let now_dur = self.timer.elapsed();
        // Update frequency in seconds.
        if (now_dur - self.last_msg).as_secs() >= UPDATE_SECS {
            let speed = (now_dur.as_secs_f64() - self.stage_start.as_secs_f64()) / self.solved_genotypes as f64;
            log::debug!("        [{:width$}/{}, {:7.4} s/gt]  Best: {} -> {:8.0}", self.solved_genotypes,
                self.curr_genotypes, speed, self.best_str, Ln::to_log10(self.best_lik), width = self.num_width);
            self.last_msg = now_dur;
        }
    }

    fn finish_stage(&mut self) {
        self.last_msg = self.timer.elapsed();
        let speed = (self.last_msg.as_secs_f64() - self.stage_start.as_secs_f64()) / self.solved_genotypes as f64;
        log::info!("      * Finished in {} ({:.4} s/gt)",
            ext::fmt::Duration(self.last_msg - self.stage_start), speed);
    }
}
