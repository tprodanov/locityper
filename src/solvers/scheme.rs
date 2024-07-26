//! Series of solvers.
//! Each next solver is expected to take more time, therefore is called on a smaller subset of haplotypes.

use std::{
    thread, fmt,
    cmp::{min, max},
    io::{self, Write},
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
        windows::ContigInfo,
        distr_cache::DistrCache,
    },
    command::DebugLvl,
};
use super::Solver;

pub struct SchemeParams {
    pub stages: String,
    pub greedy_params: Vec<String>,
    pub anneal_params: Vec<String>,
    pub highs_params: Vec<String>,
    pub gurobi_params: Vec<String>,
}

impl Default for SchemeParams {
    fn default() -> Self {
        Self {
            stages: "filter,anneal".to_owned(),
            greedy_params: Vec::new(),
            anneal_params: Vec::new(),
            highs_params: Vec::new(),
            gurobi_params: Vec::new(),
        }
    }
}

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
fn prefilter_genotypes(
    contigs: &ContigNames,
    genotypes: &[Genotype],
    priors: &[f64],
    all_alns: &AllAlignments,
    writer: &mut impl Write,
    params: &AssgnParams,
    threads: usize,
) -> io::Result<Vec<usize>>
{
    let n = genotypes.len();
    let mut ixs = (0..n).collect();
    log::info!("*** Filtering genotypes based on read alignment likelihood");

    let contig_ids: Vec<_> = contigs.ids().collect();
    let best_aln_matrix = all_alns.best_aln_matrix(&contig_ids);
    // Vector (genotype index, score).
    let mut scores = Vec::with_capacity(n);
    let mut gt_best_probs = vec![0.0_f64; all_alns.reads().len()];
    for (gt, &prior) in genotypes.iter().zip(priors) {
        let gt_ids = gt.ids();
        gt_best_probs.copy_from_slice(&best_aln_matrix[gt_ids[0].ix()]);
        for id in gt_ids[1..].iter() {
            gt_best_probs.iter_mut().zip(&best_aln_matrix[id.ix()])
                .for_each(|(best_prob, &curr_prob)| *best_prob = best_prob.max(curr_prob));
        }
        let score = prior + gt_best_probs.iter().sum::<f64>();
        writeln!(writer, "{}\t{:.3}", gt, Ln::to_log10(score))?;
        scores.push(score);
    }
    truncate_ixs(&mut ixs, &scores, genotypes, params.filt_diff, params.min_gts, threads);
    Ok(ixs)
}

/// Calculates average sum window weights across remaining genotypes.
fn calc_average_sum_weight(genotypes: &[Genotype], ixs: &[usize], contig_infos: &[Arc<ContigInfo>]) -> f64 {
    let mut s = 0.0;
    for &i in ixs {
        s += genotypes[i].ids().iter()
            .map(|id| contig_infos[id.ix()].default_weights().iter().sum::<f64>()).sum::<f64>();
    }
    s / ixs.len() as f64
}

const MAX_STAGES: usize = 10;
const LATIN_NUMS: [&'static str; MAX_STAGES] = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"];

/// Solver scheme.
/// Consists of multiple solvers, each executed on (possibly) smaller subset of genotypes.
pub struct Scheme {
    /// Is pre-filtering enabled?
    filter: bool,
    stages: Vec<Box<dyn Solver>>,
}

impl Scheme {
    pub fn create(params: &SchemeParams) -> crate::Result<Self> {
        let mut filter = false;
        let mut stages = Vec::new();
        for (i, stage) in params.stages.split(',').enumerate() {
            let (mut solver, solver_params): (Box<dyn Solver>, _) = match &stage.trim().to_lowercase() as &str {
                "filter" | "filtering" => {
                    validate_param!(i == 0, "Filtering stage must go first");
                    filter = true;
                    continue;
                }
                "greedy" => (Box::new(super::GreedySolver::default()), &params.greedy_params),
                "anneal" | "simanneal" | "annealing" | "simannealing"
                    => (Box::new(super::SimAnneal::default()), &params.anneal_params),
                "highs" => {
                    #[cfg(feature = "highs")]
                    { (Box::new(super::HighsSolver::default()), &params.highs_params) }
                    #[cfg(not(feature = "highs"))]
                    { return Err(error!(RuntimeError,
                        "HiGHS solver is disabled. Please recompile with `highs` feature.")) }
                }
                "gurobi" => {
                    #[cfg(feature = "gurobi")]
                    { (Box::new(super::GurobiSolver::default()), &params.gurobi_params) }
                    #[cfg(not(feature = "gurobi"))]
                    { return Err(error!(RuntimeError,
                        "Gurobi solver is disabled. Please recompile with `gurobi` feature.")) }
                }
                _ => return Err(error!(InvalidInput, "Unknown solver {:?}", stage)),
            };
            solver.set_params(solver_params)?;
            stages.push(solver);
        }
        validate_param!(stages.len() <= MAX_STAGES, "Too many solution stages ({}), allowed at most {}",
            stages.len(), MAX_STAGES);
        if stages.is_empty() {
            log::warn!("No stages selected, solving disabled!");
        }
        Ok(Self { filter, stages })
    }

    fn iter(&self) -> std::slice::Iter<'_, Box<dyn Solver>> {
        self.stages.iter()
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
    pub all_contig_infos: Vec<Arc<ContigInfo>>,
    pub distr_cache: Arc<DistrCache>,

    /// Genotypes and their priors.
    pub genotypes: Vec<Genotype>,
    pub priors: Vec<f64>,

    pub assgn_params: AssgnParams,
    pub debug: DebugLvl,
    pub threads: usize,
    pub is_paired_end: bool,
}

const ALNS_CSV_HEADER: &'static str = "read_hash\tcontig\tpos1\tposs2\tlik\tselected";

/// Write reads and their assignments to a CSV file in the following format (tab-separated):
/// `prefix  read_hash  aln1  aln2  w1  w2  log10-prob  selected`
fn write_alns(
    f: &mut impl io::Write,
    prefix: &str,
    gt_alns: &GenotypeAlignments,
    all_alns: &AllAlignments,
    assgn_counts: &[u16],
    attempts: u16,
) -> io::Result<()> {
    let attempts = f32::from(attempts);
    let mut counts_iter = assgn_counts.iter();
    for (rp, groupped_alns) in all_alns.reads().iter().enumerate() {
        let hash = groupped_alns.read_data().name_hash();
        for curr_windows in gt_alns.possible_read_alns(rp) {
            write!(f, "{}{}\t", prefix, hash)?;
            match curr_windows.parent() {
                None => write!(f, "*\t*\t*\t")?,
                Some((contig_ix, aln_pair)) => {
                    write!(f, "{}\t", contig_ix + 1)?;
                    match aln_pair.ix1() {
                        Some(i) => write!(f, "{}\t", groupped_alns.ith_aln(i).interval().start() + 1)?,
                        None => write!(f, "*\t")?,
                    };
                    match aln_pair.ix2() {
                        Some(i) => write!(f, "{}\t", groupped_alns.ith_aln(i).interval().start() + 1)?,
                        None => write!(f, "*\t")?,
                    };
                }
            }
            let prob = Ln::to_log10(curr_windows.ln_prob());
            let count = *counts_iter.next().expect("Not enough assignment counts");
            writeln!(f, "{:.3}\t{:.2}", prob, f32::from(count) / attempts)?;
        }
    }
    assert!(counts_iter.next().is_none(), "Too many assignment counts");
    Ok(())
}

/// Returns probability that the first mean is larger than the second.
/// If variances are well defined, returns ln(p-value) of the t-test.
/// Otherwise, returns probability based on the difference between the two means.
fn compare_two_likelihoods(mean1: f64, var1: f64, mean2: f64, var2: f64, attempts: f64) -> f64 {
    // Simple normalization.
    let simple_norm = mean1 - Ln::add(mean1, mean2);
    if var1.is_normal() && var2.is_normal() {
        simple_norm.max(
            math::unpaired_onesided_t_test::<false>(mean1, var2, mean2, var2, attempts).ln())
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

struct Likelihoods {
    /// Mean and variance of various genotypes.
    likelihoods: Vec<(f64, f64)>,
    /// Genotype indices after the latest iteration.
    ixs: Vec<usize>,
    /// Count how many times each read assignment was selected.
    assgn_counts: Vec<Option<Vec<u16>>>,
}

impl Likelihoods {
    fn new(n_gts: usize, ixs: Vec<usize>) -> Self {
        Self {
            likelihoods: vec![(f64::NEG_INFINITY, 0.0); n_gts],
            ixs,
            assgn_counts: vec![None; n_gts],
        }
    }

    /// Number of genotypes, still in contest.
    fn curr_len(&self) -> usize {
        self.ixs.len()
    }

    fn discard_improbable_genotypes(
        &mut self,
        genotypes: &[Genotype],
        params: &AssgnParams,
        threads: usize,
    ) {
        let n = self.ixs.len();
        let min_gts = params.min_gts.max(threads);
        if params.prob_thresh == f64::NEG_INFINITY || min_gts >= n {
            return;
        }

        let attempts = f64::from(params.attempts);
        // Decreasing sort by average likelihood.
        self.ixs.sort_unstable_by(|&i, &j| self.likelihoods[j].0.total_cmp(&self.likelihoods[i].0));
        let best_ix = self.ixs[0];
        let worst_ix = self.ixs[n - 1];
        let (mean1, var1) = self.likelihoods[best_ix];
        let (worst, _) = self.likelihoods[worst_ix];
        log::debug!("        Worst: {} ({:.0}), Best: {} ({:.0})",
            genotypes[worst_ix], Ln::to_log10(worst),
            genotypes[best_ix], Ln::to_log10(mean1));

        // Stop when encountered 5 unlikely genotypes.
        const STOP_COUNT: u32 = 5;
        let mut dropped = 0;
        let mut new_ixs = self.ixs[..min_gts].to_vec();

        let mut ixs_iter = self.ixs[min_gts..].iter().copied();
        while let Some(ix) = ixs_iter.next() {
            let (curr_mean, curr_var) = self.likelihoods[ix];
            let ln_pval = compare_two_likelihoods(curr_mean, curr_var, mean1, var1, attempts);
            if ln_pval >= params.prob_thresh {
                new_ixs.push(ix);
            } else {
                dropped += 1;
                self.assgn_counts[ix] = None;
                if dropped >= STOP_COUNT {
                    break;
                }
            }
        }
        // Remove dropped assignment counts.
        for ix in ixs_iter {
            self.assgn_counts[ix] = None;
        }
        let m = new_ixs.len();
        let (thresh, _) = self.likelihoods[new_ixs[m - 1]];
        log::debug!("        Keep {}/{} genotypes ({:.1}%), smallest lik: {:.0} ({:.1}%)",
            m, n, 100.0 * m as f64 / n as f64, Ln::to_log10(thresh),
            100.0 * (thresh - worst) / (mean1 - worst).max(1.0));
        self.ixs = new_ixs;
    }

    fn produce_result(&mut self, data: &Data) -> Genotyping {
        // ln(1e-5).
        const THRESH: f64 = -11.512925464970229;
        let params = &data.assgn_params;
        let min_output = max(4, params.out_bams);
        let thresh_prob = THRESH.min(params.prob_thresh);
        let attempts = f64::from(params.attempts);
        self.ixs.sort_unstable_by(|&i, &j| self.likelihoods[j].0.total_cmp(&self.likelihoods[i].0));
        let mut n = self.ixs.len();
        if n < 2 {
            log::warn!("Only {} genotype(s) remaining, quality will be undefined", n);
        }

        let genotypes = &data.genotypes;
        let mut ln_probs = vec![0.0; n];
        let mut out_genotypes = Vec::with_capacity(n);
        let mut assgn_counts = Vec::with_capacity(n);
        let mut mean_sds = Vec::with_capacity(n);
        let mut i = 0;
        // Use while instead of for, as upper bound can change.
        while i < n {
            let (mean_i, var_i) = self.likelihoods[self.ixs[i]];
            out_genotypes.push(genotypes[self.ixs[i]].clone());
            assgn_counts.push(self.assgn_counts[self.ixs[i]].take().expect("Assignment counts undefined"));
            mean_sds.push((mean_i, var_i.sqrt()));
            for j in i + 1..n {
                let (mean_j, var_j) = self.likelihoods[self.ixs[j]];
                let prob_j = compare_two_likelihoods(mean_j, var_j, mean_i, var_i, attempts);
                if i == 0 && j >= min_output && prob_j < thresh_prob {
                    n = j;
                    break;
                }
                ln_probs[i] += (-prob_j.exp()).ln_1p();
                ln_probs[j] += prob_j;
            }
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
            assgn_counts, mean_sds, ln_probs, quality,
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
    /// Count how many times each read assignment was selected.
    assgn_counts: Vec<Vec<u16>>,
    /// Mean and standard deviation, divided by sqrt(attempts).
    mean_sds: Vec<(f64, f64)>,
    /// Normalized ln-probabilities for each genotype.
    ln_probs: Vec<f64>,
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
            assgn_counts: Vec::new(),
            mean_sds: Vec::new(),
            ln_probs: Vec::new(),
            weighted_dist: None,
            distances: None,
            total_reads: 0,
            unexpl_reads: None,
            quality: 0.0,
        }
    }

    pub fn print_log(&self) {
        let (mean, sd) = self.mean_sds[0];
        log::info!("    [{}] GT = {},  lik. = {:.2} ± {:.2}, confidence = {:.4}%, weight.dist = {:.4}",
            self.tag, self.genotypes[0], Ln::to_log10(mean), Ln::to_log10(sd),
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
        all_contig_infos: &[Arc<ContigInfo>],
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
        //     let contig_info = &all_contig_infos[contig_id.ix()];
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
    pub fn count_unexplained_reads(&mut self, all_alns: &AllAlignments, unmapped_penalty: f64) {
        // Correct for rounding errors.
        let unmapped_penalty = unmapped_penalty + 1e-8;
        let mut unexplained = 0_u32;
        let gt = &self.genotypes[0];
        for read in all_alns.reads() {
            let best_lik = IterExt::max(gt.ids().iter().map(|&id| read.best_at_contig(id)));
            // Unmapped penalty is negative.
            if best_lik < read.max_lik() + unmapped_penalty {
                unexplained += 1;
            }
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
                let mut obj = json::object! {
                    genotype: gt.name(),
                    lik_mean: Ln::to_log10(self.mean_sds[i].0),
                    lik_sd: Ln::to_log10(self.mean_sds[i].1),
                    prob: self.ln_probs[i].exp(),
                    log10_prob: Ln::to_log10(self.ln_probs[i]),
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
    mean_weight: f64,
    rng: &mut XoshiroRng,
    likelihoods: &mut Likelihoods,
    mut sol_writer: impl Write,
    mut depth_writer: Option<GzFile>,
    mut aln_writer: Option<GzFile>,
) -> crate::Result<()>
{
    let total_genotypes = data.genotypes.len();
    let mut logger = Logger::new(&data.scheme, total_genotypes);

    let n_stages = data.scheme.stages.len();
    let attempts = data.assgn_params.attempts;
    let distr_cache = data.distr_cache.deref();
    for (stage_ix, solver) in data.scheme.iter().enumerate() {
        logger.start_stage(stage_ix, likelihoods.curr_len());
        if stage_ix + 1 < n_stages && likelihoods.curr_len() <= data.assgn_params.min_gts {
            log::warn!("    Skip stage {} (too few genotypes)", LATIN_NUMS[stage_ix]);
            continue;
        }

        let mut liks = vec![f64::NAN; usize::from(attempts)];
        for &ix in likelihoods.ixs.iter() {
            let gt = &data.genotypes[ix];
            let prefix = format!("{}\t{}\t", stage_ix + 1, gt);
            let prior = data.priors[ix];
            let mut gt_alns = GenotypeAlignments::new(gt.clone(), &data.all_contig_infos, &data.all_alns,
                &data.assgn_params);
            let mut counts = gt_alns.create_counts();

            for (attempt, lik) in liks.iter_mut().enumerate() {
                gt_alns.apply_tweak(rng, distr_cache, mean_weight, &data.assgn_params);
                let assgns = solver.solve(&gt_alns, rng)?;
                *lik = prior + assgns.likelihood();
                let ext_prefix = format!("{}{}\t", prefix, attempt + 1);
                if let Some(writer) = depth_writer.as_mut() {
                    assgns.write_depth(writer, &ext_prefix).map_err(add_path!(!))?;
                }
                assgns.update_counts(&mut counts);
                assgns.summarize(&mut sol_writer, &ext_prefix).map_err(add_path!(!))?;
            }
            if let Some(writer) = aln_writer.as_mut() {
                write_alns(writer, &prefix, &gt_alns, &data.all_alns, &counts, attempts).map_err(add_path!(!))?;
            }
            let (lik_mean, lik_var) = F64Ext::mean_variance_or_nan(&liks);
            logger.update(gt.name(), lik_mean);
            likelihoods.likelihoods[ix] = (lik_mean, lik_var);
            likelihoods.assgn_counts[ix] = Some(counts);
        }
        if stage_ix + 1 < n_stages {
            likelihoods.discard_improbable_genotypes(&data.genotypes, &data.assgn_params, 1);
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

/// If there is more than one file, consecutively append all of them to the first file, and keep only the first one.
fn merge_files(filenames: &[PathBuf]) -> crate::Result<()> {
    if filenames.len() > 1 {
        // By this point, all depth_writers should be already dropped.
        let mut file1 = std::fs::OpenOptions::new().append(true).open(&filenames[0]).map_err(add_path!(filenames[0]))?;
        ext::sys::concat_files(filenames[1..].iter(), &mut file1)?;
        filenames[1..].iter().map(|path| std::fs::remove_file(path).map_err(add_path!(path)))
            .collect::<crate::Result<()>>()?;
    }
    Ok(())
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

    let mut depth_filenames = None;
    let mut depth_writers = None;
    let mut aln_filenames = None;
    let mut aln_writers = None;
    if data.debug >= DebugLvl::Some {
        depth_filenames = Some(csv_filenames(&locus_dir.join("depth"), data.threads));
        depth_writers = depth_filenames.as_ref().map(|filenames| open_gzips(filenames)).transpose()?;
        writeln!(depth_writers.as_mut().unwrap()[0], "stage\tgenotype\tattempt\t{}", ReadAssignment::DEPTH_CSV_HEADER)
            .map_err(add_path!(!))?;
    }

    if data.debug >= DebugLvl::Full {
        aln_filenames = Some(csv_filenames(&locus_dir.join("alns"), data.threads));
        aln_writers = aln_filenames.as_ref().map(|filenames| open_gzips(filenames)).transpose()?;
        writeln!(aln_writers.as_mut().unwrap()[0], "stage\tgenotype\t{}", ALNS_CSV_HEADER).map_err(add_path!(!))?;
    }

    let rem_ixs = if data.scheme.filter && (data.debug >= DebugLvl::Some || n_gts > data.assgn_params.min_gts) {
        let filt_filename = locus_dir.join("filter.csv.gz");
        let mut filt_writer = ext::sys::create_gzip(&filt_filename)?;
        writeln!(filt_writer, "genotype\tscore").map_err(add_path!(filt_filename))?;
        prefilter_genotypes(&data.contigs, &data.genotypes, &data.priors, &data.all_alns, &mut filt_writer,
            &data.assgn_params, data.threads).map_err(add_path!(filt_filename))?
    } else {
        (0..n_gts).collect()
    };
    let mean_weight = if data.assgn_params.depth_norm_power != 0.0 {
        calc_average_sum_weight(&data.genotypes, &rem_ixs, &data.all_contig_infos)
    } else {
        f64::NAN
    };

    let sol_filename = locus_dir.join("sol.csv.gz");
    let mut sol_writer = ext::sys::create_gzip(&sol_filename)?;
    writeln!(sol_writer, "stage\tgenotype\tattempt\t{}", ReadAssignment::SUMMARY_HEADER)
        .map_err(add_path!(sol_filename))?;
    let mut likelihoods = Likelihoods::new(n_gts, rem_ixs);
    let data = Arc::new(data);
    if data.threads == 1 {
        solve_single_thread(&data, mean_weight, rng, &mut likelihoods, sol_writer,
            depth_writers.map(|mut writers| writers.pop().unwrap()),
            aln_writers.map(|mut writers| writers.pop().unwrap()))?;
    } else {
        let main_worker = MainWorker::new(Arc::clone(&data), mean_weight, rng, depth_writers, aln_writers);
        main_worker.run(sol_writer, rng, &mut likelihoods)?;
    }
    if let Some(filenames) = depth_filenames {
        merge_files(&filenames)?;
    }
    if let Some(filenames) = aln_filenames {
        merge_files(&filenames)?;
    }

    let data = &data;
    let mut genotyping = likelihoods.produce_result(&data);
    if data.assgn_params.out_bams > 0 {
        let bam_dir = locus_dir.join(crate::command::paths::ALNS_DIR);
        ext::sys::mkdir(&bam_dir)?;
        let mut contig_to_tid = Vec::new();
        log::info!("    Writing output alignments to {}", ext::fmt::path(&bam_dir));
        let out_gts = data.assgn_params.out_bams.min(genotyping.genotypes.len());
        let width = max(2, math::num_digits(out_gts as f64) as usize);
        for (i, (gt, assgn_counts)) in genotyping.genotypes.iter().zip(&genotyping.assgn_counts)
                .take(data.assgn_params.out_bams).enumerate() {
            let bam_path = bam_dir.join(format!("{:0width$}.bam", i));
            model::bam::write_bam(&bam_path, gt, data, &mut contig_to_tid, assgn_counts)?
        }
    }
    if let Some(dist_matrix) = &data.contig_distances {
        genotyping.find_weighted_dist(dist_matrix);
    }
    genotyping.check_first_prob();
    genotyping.check_num_of_reads(data.all_alns.reads().len() as u32, bg_distr.depth(), &data.all_contig_infos);
    genotyping.count_unexplained_reads(&data.all_alns, data.assgn_params.unmapped_penalty);
    genotyping.print_log();
    Ok(genotyping)
}

/// Task, sent to the workers: stage index and a vector of genotype indices.
type Task = (usize, Vec<usize>);
/// Genotype index and calculated likelihood mean and variance.
type Solution = (usize, f64, f64, Vec<u16>);

struct MainWorker {
    data: Arc<Data>,
    senders: Vec<Sender<Task>>,
    receivers: Vec<Receiver<Solution>>,
    handles: Vec<thread::JoinHandle<crate::Result<Vec<u8>>>>,
}

impl MainWorker {
    fn new(
        data: Arc<Data>,
        mean_weight: f64,
        rng: &mut XoshiroRng,
        depth_writers: Option<Vec<GzFile>>,
        aln_writers: Option<Vec<GzFile>>,
    ) -> Self
    {
        let n_workers = data.threads;
        let mut senders = Vec::with_capacity(n_workers);
        let mut receivers = Vec::with_capacity(n_workers);
        let mut handles = Vec::with_capacity(n_workers);

        let mut depth_writers_iter = depth_writers.into_iter().flatten();
        let mut aln_writers_iter = aln_writers.into_iter().flatten();
        for _ in 0..data.threads {
            let (task_sender, task_receiver) = mpsc::channel();
            let (sol_sender, sol_receiver) = mpsc::channel();
            let worker = Worker {
                data: Arc::clone(&data),
                rng: rng.clone(),
                receiver: task_receiver,
                sender: sol_sender,
                sol_buffer: Vec::new(),
                depth_writer: depth_writers_iter.next(),
                aln_writer: aln_writers_iter.next(),
                mean_weight,
            };
            rng.jump();
            senders.push(task_sender);
            receivers.push(sol_receiver);
            handles.push(thread::spawn(|| worker.run()));
        }
        MainWorker { data, senders, receivers, handles }
    }

    fn run(self,
        mut sol_writer: impl Write,
        rng: &mut impl Rng,
        likelihoods: &mut Likelihoods,
    ) -> crate::Result<()>
    {
        let n_workers = self.handles.len();
        let data = self.data.deref();
        let scheme = data.scheme.deref();
        let total_genotypes = data.genotypes.len();
        let mut logger = Logger::new(&scheme, total_genotypes);
        let mut rem_jobs = Vec::with_capacity(n_workers);

        let n_stages = scheme.stages.len();
        for stage_ix in 0..n_stages {
            logger.start_stage(stage_ix, likelihoods.curr_len());
            let m = likelihoods.curr_len();
            if stage_ix + 1 < n_stages && m <= data.assgn_params.min_gts.max(data.threads) {
                log::warn!("    Skip stage {} (too few genotypes)", LATIN_NUMS[stage_ix]);
                continue;
            }

            rem_jobs.clear();
            let mut start = 0;
            // Shuffle indices to better distribute them to threads.
            likelihoods.ixs.shuffle(rng);
            for (i, sender) in self.senders.iter().enumerate() {
                if start == m {
                    break;
                }
                let rem_workers = n_workers - i;
                let curr_jobs = (m - start).fast_ceil_div(rem_workers);
                let task = likelihoods.ixs[start..start + curr_jobs].to_vec();
                sender.send((stage_ix, task)).expect("Genotyping worker has failed!");
                start += curr_jobs;
                rem_jobs.push(curr_jobs);
            }
            assert_eq!(start, m);

            while logger.solved_genotypes < m {
                for (receiver, jobs) in self.receivers.iter().zip(rem_jobs.iter_mut()) {
                    if *jobs > 0 {
                        let (ix, lik_mean, lik_var, counts) = receiver.recv().expect("Genotyping worker has failed!");
                        logger.update(data.genotypes[ix].name(), lik_mean);
                        likelihoods.likelihoods[ix] = (lik_mean, lik_var);
                        likelihoods.assgn_counts[ix] = Some(counts);
                        *jobs -= 1;
                    }
                }
            }
            if stage_ix + 1 < n_stages {
                likelihoods.discard_improbable_genotypes(&data.genotypes, &data.assgn_params, data.threads);
            }
            logger.finish_stage();
        }
        std::mem::drop(self.senders);
        for handle in self.handles.into_iter() {
            let sol_bytes = handle.join().expect("Process failed for unknown reason")?;
            sol_writer.write_all(&sol_bytes).map_err(add_path!(!))?;
        }
        Ok(())
    }
}

struct Worker {
    data: Arc<Data>,
    rng: XoshiroRng,
    receiver: Receiver<Task>,
    sender: Sender<Solution>,
    sol_buffer: Vec<u8>,
    depth_writer: Option<GzFile>,
    aln_writer: Option<GzFile>,
    mean_weight: f64,
}

impl Worker {
    fn run(mut self) -> crate::Result<Vec<u8>> {
        let data = self.data.deref();
        let scheme = data.scheme.deref();
        let attempts = data.assgn_params.attempts;
        let distr_cache = data.distr_cache.deref();

        // Block thread and wait for the shipment.
        while let Ok((stage_ix, task)) = self.receiver.recv() {
            let solver = &scheme.stages[stage_ix];
            let n = task.len();
            assert_ne!(n, 0, "Received empty task");

            let mut liks = vec![f64::NAN; usize::from(attempts)];
            for ix in task.into_iter() {
                let gt = &data.genotypes[ix];
                let prefix = format!("{}\t{}\t", stage_ix + 1, gt);
                let prior = data.priors[ix];
                let mut gt_alns = GenotypeAlignments::new(gt.clone(), &data.all_contig_infos, &data.all_alns,
                    &data.assgn_params);
                let mut counts = gt_alns.create_counts();
                for (attempt, lik) in liks.iter_mut().enumerate() {
                    gt_alns.apply_tweak(&mut self.rng, distr_cache, self.mean_weight, &data.assgn_params);
                    let assgns = solver.solve(&gt_alns, &mut self.rng)?;
                    *lik = prior + assgns.likelihood();
                    let ext_prefix = format!("{}{}\t", prefix, attempt + 1);
                    if let Some(writer) = self.depth_writer.as_mut() {
                        assgns.write_depth(writer, &ext_prefix).map_err(add_path!(!))?;
                    }
                    assgns.update_counts(&mut counts);
                    assgns.summarize(&mut self.sol_buffer, &ext_prefix).map_err(add_path!(!))?;
                }
                if let Some(writer) = self.aln_writer.as_mut() {
                    write_alns(writer, &prefix, &gt_alns, &data.all_alns, &counts, attempts).map_err(add_path!(!))?;
                }
                let (lik_mean, lik_var) = F64Ext::mean_variance_or_nan(&liks);
                if let Err(_) = self.sender.send((ix, lik_mean, lik_var, counts)) {
                    log::error!("Read recruitment: main thread stopped before the child thread.");
                    break;
                }
            }
        }
        Ok(self.sol_buffer)
    }
}

struct Logger<'a> {
    scheme: &'a Scheme,
    stage_ix: usize,
    solved_genotypes: usize,
    curr_genotypes: usize,
    num_width: usize,

    best_lik: f64,
    best_str: String,

    timer: Instant,
    stage_start: Duration,
    last_msg: Duration,
}

impl<'a> Logger<'a> {
    fn new(scheme: &'a Scheme, total_genotypes: usize) -> Self {
        Self {
            scheme,
            stage_ix: 0,
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

    fn start_stage(&mut self, stage_ix: usize, curr_genotypes: usize) {
        self.stage_start = self.timer.elapsed();
        self.stage_ix = stage_ix;
        assert!(curr_genotypes <= self.curr_genotypes);
        if curr_genotypes < self.curr_genotypes {
            self.curr_genotypes = curr_genotypes;
            self.num_width = math::num_digits(self.curr_genotypes as f64) as usize;
        }
        self.solved_genotypes = 0;
        self.last_msg = self.stage_start;
        log::info!("*** Stage {:>3}.  {} on {} genotypes",
            LATIN_NUMS[stage_ix], self.scheme.stages[stage_ix], self.curr_genotypes);
    }

    fn update(&mut self, gt_name: &str, lik: f64) {
        if lik > self.best_lik {
            self.best_lik = lik;
            self.best_str = gt_name.to_owned();
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
