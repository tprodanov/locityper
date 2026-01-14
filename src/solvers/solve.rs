//! Series of solvers.
//! Each next solver is expected to take more time, therefore is called on a smaller subset of haplotypes.

use std::{
    thread, iter, fmt,
    cmp::{min, max},
    io::Write,
    time::{Instant, Duration},
    path::{Path, PathBuf},
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
        sys::GzFile,
        vec::{F64Ext, IterExt},
        rand::XoshiroRng,
        fmt::PrettyUsize,
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

/// Sorts indices by score (largest = best), and truncates the list if needed.
fn truncate_ixs(
    ixs: &mut Vec<usize>,
    scores: &[f64],
    genotypes: &[Genotype],
    filt_diff: f64,
    min_size: usize,
    threads: usize,
) {
    ixs.sort_unstable_by(|&i, &j| scores[j].total_cmp(&scores[i]));
    let n = ixs.len();

    let best = scores[ixs[0]];
    let worst = scores[ixs[n - 1]];
    log::debug!("        Worst: {} ({:.0}), Best: {} ({:.0})",
        genotypes[ixs[n - 1]], Ln::to_log10(worst), genotypes[ixs[0]], Ln::to_log10(best));
    let mut thresh = best - filt_diff;
    if min_size >= n || worst >= thresh {
        log::debug!("        Keep {n}/{n} genotypes (100.0%)");
        return;
    }

    let mut m = ixs.partition_point(|&i| scores[i] >= thresh);
    if m < min_size {
        thresh = scores[ixs[min_size - 1]];
        m = ixs.partition_point(|&i| scores[i] >= thresh);
    }

    // Increase up to `threads`.
    m = m.max(threads).min(n);
    ixs.truncate(m);
    log::debug!("        Keep {}/{} genotypes ({:.1}%), threshold: {:.0} ({:.1}%)",
        m, n, 100.0 * m as f64 / n as f64, Ln::to_log10(thresh), 100.0 * (thresh - worst) / (best - worst));
}

/// Filter genotypes based on read alignments alone (without accounting for read depth).
fn run_filter(
    contigs: &ContigNames,
    genotypes: &[Genotype],
    ixs: &mut Vec<usize>,
    priors: &[f64],
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
    for &i in ixs.iter() {
        let gt = &genotypes[i];
        let prior = priors[i];
        let gt_ids = gt.ids();
        gt_best_probs.copy_from_slice(&best_aln_matrix[gt_ids[0].ix()]);
        for id in gt_ids[1..].iter() {
            gt_best_probs.iter_mut().zip(&best_aln_matrix[id.ix()])
                .for_each(|(best_prob, &curr_prob)| *best_prob = best_prob.max(curr_prob));
        }
        let score = prior + gt_best_probs.iter().sum::<f64>();
        if let Some(ref mut w) = writer {
            writeln!(w, "0\t{}\t{:.3}", gt, Ln::to_log10(score)).map_err(add_path!(!))?;
        }
        scores[i] = score;
    }
    truncate_ixs(ixs, &scores, genotypes, params.filt_diff, out_size, threads);
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

/// Various structures, needed for sample genotyping.
pub struct Data {
    /// Solving scheme.
    pub scheme: Arc<Scheme>,
    pub contigs: Arc<ContigNames>,
    pub contig_distances: Option<TriangleMatrix<u32>>,
    /// All read alignments, groupped by contigs.
    pub all_alns: AllAlignments,
    pub contig_infos: ContigInfos,
    pub distr_cache: Arc<DistrCache>,

    /// Genotypes and their priors.
    pub genotypes: Vec<Genotype>,
    pub priors: Vec<f64>,

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
pub fn genotype_distance(gt1: &Genotype, gt2: &Genotype, distances: &TriangleMatrix<u32>) -> u32 {
    let mut min_dist = u32::MAX;
    ext::vec::gen_permutations(gt1.ids(), |perm_ids1| {
        let dist: u32 = perm_ids1.iter().zip(gt2.ids())
            .map(|(i, j)| distances.get_symmetric(i.ix(), j.ix()).copied().unwrap_or(0))
            .sum();
        min_dist = min(min_dist, dist);
    });
    min_dist
}

/// Prediction for a single genotype.
#[derive(Clone)]
struct Prediction {
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
}

struct Predictions {
    predictions: Vec<Option<Prediction>>,
    /// Remaining indices after the latest iteration.
    ixs: Vec<usize>,
}

impl Predictions {
    fn new(n_gts: usize) -> Self {
        Self {
            predictions: vec![None; n_gts],
            ixs: (0..n_gts).collect(),
        }
    }

    /// Number of genotypes, still in contest.
    fn curr_len(&self) -> usize {
        self.ixs.len()
    }

    fn sort_indices(&mut self) {
        // Decreasing sort by average likelihood.
        self.ixs.sort_unstable_by(|&i, &j|
            self.predictions[j].as_ref().map(|pred| pred.lik_mean).expect("Prediction undefined").total_cmp(
            &self.predictions[i].as_ref().map(|pred| pred.lik_mean).expect("Prediction undefined")));
    }

    fn discard_improbable_genotypes(
        &mut self,
        genotypes: &[Genotype],
        prob_thresh: f64,
        mut out_size: usize,
        threads: usize,
    ) {
        let n = self.ixs.len();
        out_size = out_size.max(threads);
        if prob_thresh == f64::NEG_INFINITY || out_size >= n {
            return;
        }

        self.sort_indices();
        let best_ix = self.ixs[0];
        let worst_ix = self.ixs[n - 1];
        let best = self.predictions[best_ix].as_ref().unwrap().cheap_clone();
        let worst = self.predictions[worst_ix].as_ref().unwrap().cheap_clone();
        log::debug!("        Worst: {} ({:.0}), Best: {} ({:.0})",
            genotypes[worst_ix], Ln::to_log10(worst.lik_mean),
            genotypes[best_ix], Ln::to_log10(best.lik_mean));

        // Only run sophisticated genotype comparison if output size is at most 500.
        const SOPHISTICATED_COUNT: usize = 500;

        let mut new_ixs = self.ixs[..out_size].to_vec();
        let mut rem_ixs = self.ixs[out_size..].iter().copied();
        if out_size <= SOPHISTICATED_COUNT {
            // Stop when encountered 5 unlikely genotypes.
            const STOP_COUNT: u32 = 5;
            let mut dropped = 0;
            while let Some(ix) = rem_ixs.next() {
                let pred = self.predictions[ix].as_ref().unwrap();
                let ln_pval = compare_two_likelihoods(pred, &best);
                if ln_pval >= prob_thresh {
                    new_ixs.push(ix);
                } else {
                    dropped += 1;
                    self.predictions[ix].as_mut().unwrap().assgn_counts = None;
                    if dropped >= STOP_COUNT {
                        break;
                    }
                }
            }
        }
        // Remove dropped assignment counts.
        for ix in rem_ixs {
            self.predictions[ix].as_mut().unwrap().assgn_counts = None;
        }
        let m = new_ixs.len();
        let last = self.predictions[new_ixs[m - 1]].as_ref().unwrap();
        log::debug!("        Keep {}/{} genotypes ({:.1}%), smallest lik: {:.0} ({:.1}%)",
            m, n, 100.0 * m as f64 / n as f64, Ln::to_log10(last.lik_mean),
            100.0 * (last.lik_mean - worst.lik_mean) / (best.lik_mean - worst.lik_mean).max(1.0));
        self.ixs = new_ixs;
    }

    fn produce_result(mut self, data: &Data) -> Genotyping {
        // ln(1e-5).
        const THRESH: f64 = -11.512925464970229;
        // Output at most 100 genotypes.
        const MAX_GENOTYPES: usize = 50;

        let params = &data.assgn_params;
        let min_output = max(4, params.out_bams);
        let thresh_prob = THRESH.min(params.prob_thresh);
        self.sort_indices();
        let mut n = min(self.ixs.len(), MAX_GENOTYPES);
        if n < 2 {
            log::warn!("Only {} genotype(s) remaining, quality will be undefined", n);
        }

        let genotypes = &data.genotypes;
        let mut ln_probs = vec![0.0; n];
        let mut predictions = Vec::with_capacity(n);
        let mut out_genotypes = Vec::with_capacity(n);
        let mut i = 0;
        // Use while instead of for, as upper bound can change.
        while i < n {
            let u = self.ixs[i];
            let pred_i = self.predictions[u].take().unwrap();
            out_genotypes.push(genotypes[u].clone());
            for j in i + 1..n {
                let pred_j = self.predictions[self.ixs[j]].as_ref().unwrap();
                let prob_j = compare_two_likelihoods(pred_j, &pred_i);
                if i == 0 && j >= min_output && prob_j < thresh_prob {
                    n = j;
                    break;
                }
                ln_probs[i] += (-prob_j.exp()).ln_1p();
                ln_probs[j] += prob_j;
            }
            predictions.push(pred_i);
            i += 1;
        }
        ln_probs.truncate(n);
        let norm_fct = Ln::sum(&ln_probs);
        ln_probs.iter_mut().for_each(|v| *v -= norm_fct);
        let quality = Phred::from_ln_prob(Ln::sum(&ln_probs[1..])).min(1e9);
        Genotyping {
            genotypes: out_genotypes,
            warnings: Vec::new(),
            weighted_dist: None,
            total_reads: data.all_alns.reads().len() as u32,
            distances: None,
            unexpl_reads: None,
            tag: data.contigs.tag().to_owned(),
            predictions, ln_probs, quality,
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
    distances: Option<Vec<u32>>,
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
    pub fn find_weighted_dist(&mut self, dist_matrix: &TriangleMatrix<u32>) {
        let Some(gt0) = self.genotypes.first() else { return };
        let mut sum_prob = 0.0;
        let mut sum_dist = 0.0;
        let mut distances = Vec::with_capacity(self.genotypes.len());
        for (i, (gt, &ln_prob)) in self.genotypes.iter().zip(&self.ln_probs).enumerate() {
            let prob = ln_prob.exp();
            sum_prob += prob;
            let dist = if i > 0 { genotype_distance(gt0, gt, dist_matrix) } else { 0 };
            sum_dist += prob * f64::from(dist);
            distances.push(dist);
        }
        self.weighted_dist = Some(sum_dist / sum_prob);
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
                    obj.insert("dist_to_primary", distances[i]).unwrap();
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
}

/// Returns: vector of likelihood means and variances (same order ).
fn solve_single_thread(
    data: &Data,
    rng: &mut XoshiroRng,
    predictions: &mut Predictions,
    sol_writer: &mut GzFile,
    mut files: ThreadDebugFiles,
) -> crate::Result<()>
{
    const ONE_THREAD: usize = 1;
    let total_genotypes = data.genotypes.len();
    let mut logger = Logger::new(total_genotypes);

    let distr_cache = data.distr_cache.deref();
    for stage in &data.scheme.stages {
        describe_stage(Some(stage));
        let n_gts = predictions.curr_len();
        let out_size = data.scheme.stages.get(stage.ix + 1).map(|next_stage| next_stage.in_size);
        if !(data.assgn_params.dont_skip || out_size.map(|s| s < n_gts).unwrap_or(true)) {
            log::info!("    Skipping stage, not enough genotypes");
            continue;
        }

        let mut curr_sol_writer: Option<&mut GzFile> =
            if data.debug == DebugLvl::None && out_size.is_some() { None } else { Some(sol_writer) };

        logger.reset(n_gts);
        let mut liks = vec![f64::NAN; usize::from(stage.attempts)];
        for &ix in predictions.ixs.iter() {
            let gt = &data.genotypes[ix];
            let prefix = format!("{}\t{}\t", stage.ix + 1, gt);
            let prior = data.priors[ix];
            let mut gt_alns = GenotypeAlignments::new(gt.clone(), &data.contig_infos, &data.all_alns,
                &data.assgn_params);
            let mut counts = gt_alns.create_counts();

            for (attempt, lik) in liks.iter_mut().enumerate() {
                gt_alns.apply_tweak(rng, distr_cache, &data.assgn_params);
                let assgns = stage.solver.solve(&gt_alns, rng)?;
                *lik = prior + assgns.likelihood();
                let ext_prefix = format!("{}{}\t", prefix, attempt + 1);
                if let Some(w) = &mut files.depth_writer {
                    assgns.write_depth(w, &ext_prefix).map_err(add_path!(!))?;
                }
                if let Some(w) = &mut files.ext_sol_writer {
                    assgns.summarize(w, &ext_prefix).map_err(add_path!(!))?;
                }
                assgns.update_counts(&mut counts);
            }
            let (lik_mean, lik_var) = F64Ext::mean_variance_or_nan(&liks);
            logger.update(gt, lik_mean);
            predictions.predictions[ix] = Some(Prediction::new(lik_mean, lik_var, stage.attempts, counts));
            if let Some(ref mut w) = curr_sol_writer {
                writeln!(w, "{}{:.4}", prefix, Ln::to_log10(lik_mean)).map_err(add_path!(!))?;
            }
        }
        if let Some(s) = out_size {
            predictions.discard_improbable_genotypes(&data.genotypes, data.assgn_params.prob_thresh, s, ONE_THREAD);
        }
        logger.finish_stage();
    }
    Ok(())
}

/// Creates `threads` filenames, one for each thread.
fn csv_filenames(prefix: &Path, threads: usize) -> Vec<PathBuf> {
    (0..threads).map(|i| if i == 0 {
        ext::sys::path_append(prefix, ".csv.gz")
    } else {
        ext::sys::path_append(prefix, format!(".{}.csv.gz", i))
    }).collect()
}

/// Opens one gzip file for each filename.
fn open_gzips(filenames: &[PathBuf]) -> crate::Result<Vec<GzFile>> {
    filenames.iter().map(|path| ext::sys::create_gzip(path)).collect()
}

struct ThreadDebugFiles {
    depth_writer: Option<GzFile>,
    ext_sol_writer: Option<GzFile>,
}

struct DebugFiles {
    depth_filenames: Option<Vec<PathBuf>>,
    depth_writers: Option<Vec<GzFile>>,

    ext_sol_filenames: Option<Vec<PathBuf>>,
    ext_sol_writers: Option<Vec<GzFile>>,
}

impl DebugFiles {
    fn new(locus_dir: &Path, threads: usize, debug: DebugLvl) -> crate::Result<Self> {
        let mut depth_filenames = None;
        let mut depth_writers = None;
        if debug == DebugLvl::Full {
            depth_filenames = Some(csv_filenames(&locus_dir.join("depth"), threads));
            depth_writers = depth_filenames.as_ref().map(|filenames| open_gzips(filenames)).transpose()?;
            writeln!(depth_writers.as_mut().unwrap()[0],
                "stage\tgenotype\tattempt\t{}", ReadAssignment::DEPTH_CSV_HEADER).map_err(add_path!(!))?;
        }

        let mut ext_sol_filenames = None;
        let mut ext_sol_writers = None;
        if debug >= DebugLvl::Some {
            ext_sol_filenames = Some(csv_filenames(&locus_dir.join("sol_ext"), threads));
            ext_sol_writers = ext_sol_filenames.as_ref().map(|filenames| open_gzips(filenames)).transpose()?;
            writeln!(ext_sol_writers.as_mut().unwrap()[0],
                "stage\tgenotype\tattempt\t{}", ReadAssignment::SUMMARY_HEADER).map_err(add_path!(!))?;
        }

        Ok(Self {
            depth_filenames, depth_writers, ext_sol_filenames, ext_sol_writers,
        })
    }

    /// Takes files and splits them between threads.
    fn into_threads(&mut self) -> impl Iterator<Item = ThreadDebugFiles> + '_ {
        let depth_iter = self.depth_writers.take().into_iter()
            .flatten().map(Some).chain(iter::repeat_with(|| None));
        let ext_sol_iter = self.ext_sol_writers.take().into_iter()
            .flatten().map(Some).chain(iter::repeat_with(|| None));
        depth_iter.zip(ext_sol_iter)
            .map(|(depth_writer, ext_sol_writer)| ThreadDebugFiles { depth_writer, ext_sol_writer })
    }

    fn merge(&self) -> crate::Result<()> {
        debug_assert!(self.depth_writers.is_none() && self.ext_sol_writers.is_none());
        if let Some(filenames) = &self.depth_filenames {
            ext::sys::merge_files(&filenames[0], &filenames[1..])?;
        }
        if let Some(filenames) = &self.ext_sol_filenames {
            ext::sys::merge_files(&filenames[0], &filenames[1..])?;
        }
        Ok(())
    }
}

pub fn solve(
    mut data: Data,
    bg_distr: &BgDistr,
    locus_dir: &Path,
    rng: &mut XoshiroRng,
) -> crate::Result<Genotyping>
{
    let n_gts = data.genotypes.len();
    assert!(n_gts > 0);
    data.threads = min(data.threads, n_gts);
    log::info!("    Genotyping {}: {} possible genotypes", data.contigs.tag(), n_gts);

    let mut sol_writer = ext::sys::create_gzip(&locus_dir.join("sol.csv.gz"))?;
    writeln!(sol_writer, "stage\tgenotype\tscore").map_err(add_path!(!))?;
    let mut predictions = Predictions::new(n_gts);
    let out_size0 = data.scheme.stages[0].in_size;
    if data.assgn_params.dont_skip || out_size0 < n_gts {
        let curr_sol_writer = if data.debug == DebugLvl::None { None } else { Some(&mut sol_writer) };
        run_filter(&data.contigs, &data.genotypes, &mut predictions.ixs,
            &data.priors, &data.all_alns, curr_sol_writer, &data.assgn_params, out_size0, data.threads)?;
    }

    let mut files = DebugFiles::new(locus_dir, data.threads, data.debug)?;
    let data = Arc::new(data);
    if data.threads == 1 {
        solve_single_thread(&data, rng, &mut predictions, &mut sol_writer, files.into_threads().next().unwrap())?;
    } else {
        let main_worker = MainWorker::new(Arc::clone(&data), rng, files.into_threads());
        main_worker.run(&mut sol_writer, rng, &mut predictions)?;
        files.merge()?;
    }
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
        mut files: impl Iterator<Item = ThreadDebugFiles>,
    ) -> Self
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
                files: files.next().expect("Not enough debug files"),
            };
            rng.jump();
            senders.push(task_sender);
            receivers.push(sol_receiver);
            handles.push(thread::spawn(|| worker.run()));
        }
        MainWorker { data, senders, receivers, handles }
    }

    fn run(self,
        sol_writer: &mut GzFile,
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
            let n_gts = predictions.curr_len();
            let out_size = data.scheme.stages.get(stage.ix + 1).map(|next_stage| next_stage.in_size);
            if !(data.assgn_params.dont_skip || out_size.map(|s| s < n_gts).unwrap_or(true)) {
                log::info!("    Skipping stage, not enough genotypes");
                continue;
            }

            logger.reset(n_gts);
            rem_jobs.clear();
            let mut start = 0;
            // Shuffle indices to better distribute them to threads.
            predictions.ixs.shuffle(rng);
            for (i, sender) in self.senders.iter().enumerate() {
                if start == n_gts {
                    break;
                }
                let rem_workers = n_workers - i;
                let curr_jobs = (n_gts - start).fast_ceil_div(rem_workers);
                let task = predictions.ixs[start..start + curr_jobs].to_vec();
                sender.send((stage.ix, task)).expect("Genotyping worker has failed!");
                start += curr_jobs;
                rem_jobs.push(curr_jobs);
            }
            assert_eq!(start, n_gts);

            let mut curr_sol_writer: Option<&mut GzFile> =
                if data.debug == DebugLvl::None && out_size.is_some() { None } else { Some(sol_writer) };
            while logger.solved_genotypes < n_gts {
                for (receiver, jobs) in self.receivers.iter().zip(rem_jobs.iter_mut()) {
                    if *jobs > 0 {
                        let (ix, pred) = receiver.recv().expect("Genotyping worker has failed!");
                        let gt = &data.genotypes[ix];
                        logger.update(gt, pred.lik_mean);
                        if let Some(ref mut w) = curr_sol_writer {
                            writeln!(w, "{}\t{}\t{:.4}", stage.ix + 1, gt, Ln::to_log10(pred.lik_mean))
                                .map_err(add_path!(!))?;
                        }
                        predictions.predictions[ix] = Some(pred);
                        *jobs -= 1;
                    }
                }
            }
            if let Some(s) = out_size {
                predictions.discard_improbable_genotypes(&data.genotypes, data.assgn_params.prob_thresh,
                    s, data.threads);
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
    files: ThreadDebugFiles,
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
                let gt = &data.genotypes[ix];
                let prefix = format!("{}\t{}\t", stage_ix + 1, gt);
                let prior = data.priors[ix];
                let mut gt_alns = GenotypeAlignments::new(gt.clone(), &data.contig_infos, &data.all_alns,
                    &data.assgn_params);
                let mut counts = gt_alns.create_counts();
                for (attempt, lik) in liks.iter_mut().enumerate() {
                    gt_alns.apply_tweak(&mut self.rng, distr_cache, &data.assgn_params);
                    let assgns = stage.solver.solve(&gt_alns, &mut self.rng)?;
                    *lik = prior + assgns.likelihood();
                    let ext_prefix = format!("{}{}\t", prefix, attempt + 1);
                    if let Some(w) = &mut self.files.depth_writer {
                        assgns.write_depth(w, &ext_prefix).map_err(add_path!(!))?;
                    }
                    if let Some(w) = &mut self.files.ext_sol_writer {
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
