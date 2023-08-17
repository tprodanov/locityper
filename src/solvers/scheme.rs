//! Series of solvers.
//! Each next solver is expected to take more time, therefore is called on a smaller subset of haplotypes.

use std::{
    thread,
    cmp::min,
    io::{self, Write},
    time::{Instant, Duration},
    path::{Path, PathBuf},
    sync::{
        Arc,
        mpsc::{self, Sender, Receiver},
    },
    ops::Deref,
};
use statrs::distribution::Poisson;
use rand::{
    Rng,
    prelude::SliceRandom,
    distributions::Distribution,
};
use crate::{
    ext::{
        self,
        vec::F64Ext,
        rand::XoshiroRng,
    },
    math::{self, Ln, Phred},
    err::{Error, validate_param, add_path},
    seq::contigs::{ContigNames, Genotype},
    model::{
        Params as AssgnParams,
        locs::{AllAlignments, Pair},
        assgn::{GenotypeAlignments, ReadAssignment, SelectedCounter},
        windows::ContigWindows,
    },
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
    score_thresh: f64,
    min_size: usize,
    threads: usize,
) {
    ixs.sort_unstable_by(|&i, &j| scores[j].total_cmp(&scores[i]));
    let n = ixs.len();

    // ~ 10^-5.
    const MIN_DIST: f64 = 11.51;
    let best = scores[ixs[0]];
    let worst = scores[ixs[n - 1]];
    log::debug!("        Worst: {} ({:.0}), Best: {} ({:.0})",
        genotypes[ixs[n - 1]], Ln::to_log10(worst), genotypes[ixs[0]], Ln::to_log10(best));
    if min_size >= n || worst >= best - MIN_DIST {
        log::debug!("        Keep {}/{} genotypes (100.0%)", n, n);
        return;
    }

    let range = best - worst;
    let mut thresh = (worst + range * score_thresh).min(best - MIN_DIST);
    let mut m = ixs.partition_point(|&i| scores[i] >= thresh);
    if m < min_size {
        thresh = scores[ixs[min_size - 1]];
        m = ixs.partition_point(|&i| scores[i] >= thresh);
    }

    // Increase up to `threads`.
    m = m.max(threads).min(n);
    ixs.truncate(m);
    log::debug!("        Keep {}/{} genotypes ({:.1}%), threshold: {:.0} ({:.1}%)",
        m, n, 100.0 * m as f64 / n as f64, Ln::to_log10(thresh), 100.0 * (thresh - worst) / range);
}

/// Filter genotypes based on read alignments alone (without accounting for read depth).
fn prefilter_genotypes(
    contigs: &ContigNames,
    genotypes: &[Genotype],
    priors: &[f64],
    all_alns: &AllAlignments,
    lik_writer: &mut impl Write,
    params: &AssgnParams,
    debug: bool,
    threads: usize,
) -> io::Result<Vec<usize>>
{
    let n = genotypes.len();
    let mut ixs = (0..n).collect();
    // Need to run anyway in case of `debug`, so that we output all values.
    if params.min_gts.0 >= n && !debug {
        return Ok(ixs);
    }
    log::info!("*** Filtering genotypes based on read alignment likelihood");

    let contig_ids: Vec<_> = contigs.ids().collect();
    let best_aln_matrix = all_alns.best_aln_matrix(&contig_ids);
    // Vector (genotype index, score).
    let mut scores = Vec::with_capacity(n);
    let mut gt_best_probs = vec![0.0_f64; all_alns.len()];
    for (gt, &prior) in genotypes.iter().zip(priors) {
        let gt_ids = gt.ids();
        gt_best_probs.copy_from_slice(&best_aln_matrix[gt_ids[0].ix()]);
        for id in gt_ids[1..].iter() {
            gt_best_probs.iter_mut().zip(&best_aln_matrix[id.ix()])
                .for_each(|(best_prob, &curr_prob)| *best_prob = best_prob.max(curr_prob));
        }
        let score = prior + gt_best_probs.iter().sum::<f64>();
        writeln!(lik_writer, "0\t{}\t{:.3}", gt, Ln::to_log10(score))?;
        scores.push(score);
    }

    let (mean, var) = F64Ext::mean_variance(&scores);
    log::debug!("    Aln   mean {:10.3},  sd {:10.3}", Ln::to_log10(mean), Ln::to_log10(var.sqrt()));
    let mut scores_sort = scores.clone();
    scores_sort.sort_unstable_by(|x, y| y.partial_cmp(x).unwrap());
    for &i in &[50, 100, 500, 1000] {
        let (mean, var) = F64Ext::mean_variance(&scores_sort[..min(i, n)]);
        log::debug!("    Top {:4} mean {:10.3},  sd {:10.3}", min(i, n), Ln::to_log10(mean), Ln::to_log10(var.sqrt()));
    }
    for &i in &[2, 3, 5, 10, 50, 100] {
        log::debug!("    Top {:2} diff {:10.3}", i, Ln::to_log10(scores_sort[0] - scores_sort[i - 1]));
    }
    truncate_ixs(&mut ixs, &scores, genotypes, params.score_thresh, params.min_gts.0, threads);
    Ok(ixs)
}

fn calc_depth_variance(
    contig_windows: &[ContigWindows],
    mut rng: &mut impl Rng,
    genotypes: &[Genotype],
    n_reads: usize,
    n_tries: usize,
) {
    assert!(n_reads > 0);
    assert!(n_tries > 2);
    let mut probs = Vec::with_capacity(n_tries);
    for _ in 0..n_tries {
        let gt = &genotypes[rng.gen_range(0..genotypes.len())];
        let usable_windows: u32 = gt.ids().iter().map(|id| contig_windows[id.ix()].usable_windows()).sum();
        let poisson = Poisson::new(n_reads as f64 / f64::from(usable_windows)).expect("Incorrect Poisson distribution");
        let mut poisson_gen = poisson.sample_iter(&mut rng);

        let mut total_prob = 0.0;
        for id in gt.ids() {
            let windows = &contig_windows[id.ix()];
            for (&usew, distr) in windows.use_window().iter().zip(windows.depth_distrs()) {
                if usew {
                    total_prob += distr.ln_pmf(poisson_gen.next().unwrap().round() as u32);
                }
            }
        }
        log::debug!("    {} -> {:10.3}", gt, Ln::to_log10(total_prob));
        probs.push(total_prob);
    }
    let (mean, var) = F64Ext::mean_variance(&probs);
    log::debug!("    Depth mean {:10.3},  sd {:10.3}", Ln::to_log10(mean), Ln::to_log10(var.sqrt()));

    probs.sort_unstable_by(|x, y| y.partial_cmp(x).unwrap());
    for &i in &[2, 3, 5, 10, 50, 100] {
        log::debug!("    Top {:2} diff {:10.3}", i, Ln::to_log10(probs[0] - probs[i - 1]));
    }
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
    pub fn create(params: &SchemeParams) -> Result<Self, Error> {
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
                    { return Err(Error::RuntimeError(
                        "HiGHS solver is disabled. Please recompile with `highs` feature.".to_owned())) }
                }
                "gurobi" => {
                    #[cfg(feature = "gurobi")]
                    { (Box::new(super::GurobiSolver::default()), &params.gurobi_params) }
                    #[cfg(not(feature = "gurobi"))]
                    { return Err(Error::RuntimeError(
                        "Gurobi solver is disabled. Please recompile with `gurobi` feature.".to_owned())) }
                }
                "ilp" => {
                    #[cfg(feature = "highs")]
                    { (Box::new(super::HighsSolver::default()), &params.highs_params) }
                    #[cfg(all(not(feature = "highs"), feature = "gurobi"))]
                    { (Box::new(super::GurobiSolver::default()), &params.gurobi_params) }
                    #[cfg(all(not(feature = "highs"), not(feature = "gurobi")))]
                    { return Err(Error::RuntimeError(
                        "Both ILP solvers are disabled. Please recompile with `gurobi` or `highs` feature."
                        .to_owned())) }
                }
                _ => return Err(Error::InvalidInput(format!("Unknown solver {:?}", stage))),
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
    /// All read alignments, groupped by contigs.
    pub all_alns: AllAlignments,
    pub contig_windows: Vec<ContigWindows>,

    /// Genotypes and their priors.
    pub genotypes: Vec<Genotype>,
    pub priors: Vec<f64>,

    pub assgn_params: AssgnParams,
    pub debug: bool,
    pub threads: usize,
}

const ALNS_CSV_HEADER: &'static str = "read_hash\taln1\taln2\tlik\tselected";

/// Write reads and their assignments to a CSV file in the following format (tab-separated):
/// `prefix  read_hash  aln1  aln2  w1  w2  log10-prob  selected`
fn write_alns(
    f: &mut impl io::Write,
    prefix: &str,
    gt_alns: &GenotypeAlignments,
    all_alns: &AllAlignments,
    mut selected: impl Iterator<Item = f64>,
) -> io::Result<()> {
    for (rp, paired_alns) in all_alns.iter().enumerate() {
        let hash = paired_alns.name_hash();
        for curr_windows in gt_alns.possible_read_alns(rp) {
            write!(f, "{}\t{:X}\t", prefix, hash)?;
            let aln_ix = curr_windows.aln_ix();
            if aln_ix == u32::MAX {
                write!(f, "*\t*\t")?;
            } else {
                match paired_alns.ith_aln(aln_ix as usize).intervals() {
                    Pair::Both(interv1, interv2) => write!(f, "{}\t{}\t", interv1, interv2),
                    Pair::First(interv1) => write!(f, "{}\t*\t", interv1),
                    Pair::Second(interv2) => write!(f, "*\t{}\t", interv2),
                }?;
            }
            let prob = Ln::to_log10(curr_windows.ln_prob());
            writeln!(f, "{:.3}\t{:.2}", prob, selected.next().expect("Not enough `selected` values"))?;
        }
    }
    assert!(selected.next().is_none(), "Too many `selected` values");
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

struct Likelihoods {
    /// Mean and variance of various genotypes.
    likelihoods: Vec<(f64, f64)>,
    /// Genotype indices after the latest iteration.
    ixs: Vec<usize>,
}

impl Likelihoods {
    fn new(n_gts: usize, ixs: Vec<usize>) -> Self {
        Self {
            likelihoods: vec![(f64::NEG_INFINITY, 0.0); n_gts],
            ixs,
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
        let min_gts = params.min_gts.0.max(threads);
        if params.prob_thresh == f64::NEG_INFINITY || min_gts >= n {
            return;
        }

        let attempts = params.attempts as f64;
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
        for &i in self.ixs[min_gts..].iter() {
            let (mean_i, var_i) = self.likelihoods[i];
            let ln_pval = compare_two_likelihoods(mean_i, var_i, mean1, var1, attempts);
            if ln_pval >= params.prob_thresh {
                new_ixs.push(i);
            } else {
                dropped += 1;
                if dropped >= STOP_COUNT {
                    break;
                }
            }
        }
        let m = new_ixs.len();
        let (thresh, _) = self.likelihoods[new_ixs[m - 1]];
        log::debug!("        Keep {}/{} genotypes ({:.1}%), smallest lik: {:.0} ({:.1}%)",
            m, n, 100.0 * m as f64 / n as f64, Ln::to_log10(thresh),
            100.0 * (thresh - worst) / (mean1 - worst).max(1.0));
        self.ixs = new_ixs;
    }

    fn produce_result(
        &mut self,
        tag: String,
        genotypes: &[Genotype],
        params: &AssgnParams,
    ) -> Genotyping {
        // ln(1e-5).
        const THRESH: f64 = -11.512925464970229;
        const MIN_OUTPUT: usize = 4;
        let thresh_prob = THRESH.min(params.prob_thresh);
        let attempts = params.attempts as f64;
        self.ixs.sort_unstable_by(|&i, &j| self.likelihoods[j].0.total_cmp(&self.likelihoods[i].0));
        let mut n = self.ixs.len();
        if n < 2 {
            log::warn!("For some reason, only {} genotypes remaining, quality will be undefined", n);
        }

        let mut ln_probs = vec![0.0; n];
        let mut out_genotypes = Vec::with_capacity(n);
        let mut mean_sds = Vec::with_capacity(n);
        let mut i = 0;
        // Use while instead of for, as upper bound can change.
        while i < n {
            let (mean_i, var_i) = self.likelihoods[self.ixs[i]];
            out_genotypes.push(genotypes[self.ixs[i]].clone());
            mean_sds.push((mean_i, var_i.sqrt()));
            for j in i + 1..n {
                let (mean_j, var_j) = self.likelihoods[self.ixs[j]];
                let prob_j = compare_two_likelihoods(mean_j, var_j, mean_i, var_i, attempts);
                if i == 0 && j >= MIN_OUTPUT && prob_j < thresh_prob {
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
        let quality = Phred::from_ln_prob(Ln::sum(&ln_probs[1..]));
        Genotyping {
            genotypes: out_genotypes,
            tag, mean_sds, ln_probs, quality,
        }
    }
}

/// Genotyping result.
pub struct Genotyping {
    tag: String,
    /// Filtered list with the best genotypes.
    /// Genotypes are sorted from best to worst.
    genotypes: Vec<Genotype>,
    /// Mean and standard deviation, divided by sqrt(attempts).
    mean_sds: Vec<(f64, f64)>,
    /// Normalized ln-probabilities for each genotype.
    ln_probs: Vec<f64>,
    /// Quality of the top genotype.
    quality: f64,
}

impl Genotyping {
    pub fn print_log(&self) {
        let (mean, sd) = self.mean_sds[0];
        log::info!("    Best genotype for {}: {},  likelihood = {:.2} ± {:.2}, quality = {:.1}, confidence = {:.4}%",
            self.tag, self.genotypes[0], Ln::to_log10(mean), Ln::to_log10(sd),
            self.quality, 100.0 * self.ln_probs[0].exp());
    }

    pub fn to_json(&self) -> json::JsonValue {
        let options: Vec<_> = self.genotypes.iter().enumerate().map(|(i, gt)| {
            json::object! {
                genotype: gt.name(),
                lik_mean: Ln::to_log10(self.mean_sds[i].0),
                lik_sd: Ln::to_log10(self.mean_sds[i].1),
                prob: self.ln_probs[i].exp(),
                log10_prob: Ln::to_log10(self.ln_probs[i]),
            }
        }).collect();
        json::object! {
            locus: &self.tag as &str,
            genotype: self.genotypes[0].name(),
            quality: self.quality,
            options: options,
        }
    }
}

/// Returns: vector of likelihood means and variances (same order ).
fn solve_single_thread(
    data: &Data,
    likelihoods: &mut Likelihoods,
    mut lik_writer: impl Write,
    mut depth_writer: impl Write,
    mut aln_writer: impl Write,
    rng: &mut XoshiroRng,
) -> Result<(), Error>
{
    let tweak = data.assgn_params.tweak.unwrap();
    let total_genotypes = data.genotypes.len();
    let mut logger = Logger::new(&data.scheme, total_genotypes);

    let n_stages = data.scheme.stages.len();
    let mut selected = SelectedCounter::default();
    for (stage_ix, solver) in data.scheme.iter().enumerate() {
        logger.start_stage(stage_ix, likelihoods.curr_len());
        if stage_ix + 1 < n_stages && likelihoods.curr_len() <= data.assgn_params.min_gts.0 {
            log::warn!("    Skip stage {} (too few genotypes)", LATIN_NUMS[stage_ix]);
            continue;
        }

        let mut liks = vec![f64::NAN; data.assgn_params.attempts];
        for &ix in likelihoods.ixs.iter() {
            let gt = &data.genotypes[ix];
            let prior = data.priors[ix];
            let mut gt_alns = GenotypeAlignments::new(gt.clone(), &data.contig_windows, &data.all_alns,
                &data.assgn_params);
            if data.debug {
                selected.reset(&gt_alns);
            }

            for (attempt, lik) in liks.iter_mut().enumerate() {
                gt_alns.define_read_windows(tweak, rng);
                let assgns = solver.solve(&gt_alns, rng)?;
                *lik = prior + assgns.likelihood();
                if data.debug {
                    let ext_prefix = format!("{}\t{}\t{}", stage_ix + 1, gt, attempt + 1);
                    assgns.write_depth(&mut depth_writer, &ext_prefix).map_err(add_path!(!))?;
                    if data.debug {
                        selected.update(&assgns);
                    }
                }
            }
            let prefix = format!("{}\t{}", stage_ix + 1, gt);
            if data.debug {
                write_alns(&mut aln_writer, &prefix, &gt_alns, &data.all_alns, selected.fractions())
                    .map_err(add_path!(!))?;
            }
            let (lik_mean, lik_var) = F64Ext::mean_variance_or_nan(&liks);
            logger.update(gt.name(), lik_mean);
            writeln!(lik_writer, "{}\t{:.3}\t{:.3}", prefix, Ln::to_log10(lik_mean), Ln::to_log10(lik_var.sqrt()))
                .map_err(add_path!(!))?;
            likelihoods.likelihoods[ix] = (lik_mean, lik_var);
        }
        if stage_ix + 1 < n_stages {
            likelihoods.discard_improbable_genotypes(&data.genotypes, &data.assgn_params, 1);
        }
        logger.finish_stage();
    }
    Ok(())
}

/// Creates `threads` debug files for writing read assignments, and returns their file names.
/// If `debug` is false, returns a vector of `threads` sink objects + an empty vector of names.
fn create_debug_files(
    prefix: &Path,
    threads: usize,
    debug: bool,
) -> Result<(Vec<Box<dyn Write + Send>>, Vec<PathBuf>), Error>
{
    if !debug {
        return Ok(((0..threads).map(|_| Box::new(io::sink()) as Box<dyn Write + Send>).collect(), vec![]));
    }
    let filenames: Vec<_> = (0..threads).map(|i| if i == 0 {
        ext::sys::path_append(prefix, "csv.gz")
    } else {
        ext::sys::path_append(prefix, format!("{}.csv.gz", i))
    }).collect();

    let mut writers = Vec::with_capacity(threads);
    for filename in filenames.iter() {
        // Using `for` instead of `collect`, as it becomes hard to cast Box and process Errors at the same time.
        writers.push(Box::new(ext::sys::create_gzip(filename)?) as Box<dyn Write + Send>);
    }
    Ok((writers, filenames))
}

/// Merge output files if there is debug output and more than one thread.
fn merge_dbg_files(filenames: &[PathBuf]) -> Result<(), Error> {
    if filenames.len() > 1 {
        // By this point, all depth_writers should be already dropped.
        let mut file1 = std::fs::OpenOptions::new().append(true).open(&filenames[0]).map_err(add_path!(filenames[0]))?;
        ext::sys::concat_files(filenames[1..].iter(), &mut file1)?;
        filenames[1..].iter().map(|path| std::fs::remove_file(path).map_err(add_path!(path)))
            .collect::<Result<(), Error>>()?;
    }
    Ok(())
}

pub fn solve(
    mut data: Data,
    locus_dir: &Path,
    rng: &mut XoshiroRng,
) -> Result<(), Error>
{
    let n_gts = data.genotypes.len();
    assert!(n_gts > 0);
    data.threads = min(data.threads, n_gts);
    log::info!("    Genotyping {}: {} possible genotypes", data.contigs.tag(), n_gts);
    let (mut depth_writers, depth_filenames) = create_debug_files(
        &locus_dir.join("depth."), data.threads, data.debug)?;
    writeln!(depth_writers[0], "stage\tgenotype\tattempt\t{}", ReadAssignment::DEPTH_CSV_HEADER).map_err(add_path!(!))?;

    let (mut aln_writers, aln_filenames) = create_debug_files(
        &locus_dir.join("alns."), data.threads, data.debug)?;
    writeln!(aln_writers[0], "stage\tgenotype\t{}", ALNS_CSV_HEADER).map_err(add_path!(!))?;

    let lik_filename = locus_dir.join("lik.csv.gz");
    let mut lik_writer = ext::sys::create_gzip(&lik_filename)?;
    writeln!(lik_writer, "stage\tgenotype\tlik\tlik_std").map_err(add_path!(lik_filename))?;
    let rem_ixs = if data.scheme.filter {
        prefilter_genotypes(&data.contigs, &data.genotypes, &data.priors, &data.all_alns, &mut lik_writer,
            &data.assgn_params, data.debug, data.threads).map_err(add_path!(!))?
    } else {
        (0..n_gts).collect()
    };

    calc_depth_variance(&data.contig_windows, rng, &data.genotypes, data.all_alns.len(), 1000);

    let mut likelihoods = Likelihoods::new(n_gts, rem_ixs);
    let data = Arc::new(data);
    if data.threads == 1 {
        solve_single_thread(&data, &mut likelihoods, lik_writer,
            depth_writers.pop().unwrap(), aln_writers.pop().unwrap(), rng)?;
    } else {
        let main_worker = MainWorker::new(Arc::clone(&data), rng, depth_writers, aln_writers);
        main_worker.run(rng, &mut likelihoods, lik_writer)?;
    }
    merge_dbg_files(&depth_filenames)?;
    merge_dbg_files(&aln_filenames)?;

    let data = &data;
    let genotyping = likelihoods.produce_result(data.contigs.tag().to_owned(), &data.genotypes, &data.assgn_params);
    genotyping.print_log();
    let json_filename = locus_dir.join("res.json.gz");
    let mut json_writer = ext::sys::create_gzip(&json_filename)?;
    genotyping.to_json().write_pretty(&mut json_writer, 4).map_err(add_path!(json_filename))?;
    Ok(())
}

/// Task, sent to the workers: stage index and a vector of genotype indices.
type Task = (usize, Vec<usize>);
/// Genotype index and calculated likelihood mean and variance.
type Solution = (usize, f64, f64);

struct MainWorker {
    data: Arc<Data>,
    senders: Vec<Sender<Task>>,
    receivers: Vec<Receiver<Solution>>,
    handles: Vec<thread::JoinHandle<Result<(), Error>>>,
}

impl MainWorker {
    fn new(
        data: Arc<Data>,
        rng: &mut XoshiroRng,
        depth_writers: Vec<impl Write + Send + 'static>,
        aln_writers: Vec<impl Write + Send + 'static>,
    ) -> Self
    {
        let n_workers = data.threads;
        let mut senders = Vec::with_capacity(n_workers);
        let mut receivers = Vec::with_capacity(n_workers);
        let mut handles = Vec::with_capacity(n_workers);
        debug_assert!(n_workers == depth_writers.len() && n_workers == aln_writers.len());

        for (depth_writer, aln_writer) in depth_writers.into_iter().zip(aln_writers.into_iter()) {
            let (task_sender, task_receiver) = mpsc::channel();
            let (sol_sender, sol_receiver) = mpsc::channel();
            let worker = Worker {
                data: Arc::clone(&data),
                rng: rng.clone(),
                receiver: task_receiver,
                sender: sol_sender,
                depth_writer, aln_writer,
            };
            rng.jump();
            senders.push(task_sender);
            receivers.push(sol_receiver);
            handles.push(thread::spawn(|| worker.run()));
        }
        MainWorker { data, senders, receivers, handles }
    }

    fn run(self,
        rng: &mut impl Rng,
        likelihoods: &mut Likelihoods,
        mut lik_writer: impl Write,
    ) -> Result<(), Error>
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
            if stage_ix + 1 < n_stages && m <= data.assgn_params.min_gts.0.max(data.threads) {
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
                // Ceiling division.
                let curr_jobs = ((m - start) + rem_workers - 1) / rem_workers;
                let task = likelihoods.ixs[start..start + curr_jobs].to_vec();
                sender.send((stage_ix, task)).expect("Genotyping worker has failed!");
                start += curr_jobs;
                rem_jobs.push(curr_jobs);
            }
            assert_eq!(start, m);

            while logger.solved_genotypes < m {
                for (receiver, jobs) in self.receivers.iter().zip(rem_jobs.iter_mut()) {
                    if *jobs > 0 {
                        let (ix, lik_mean, lik_var) = receiver.recv().expect("Genotyping worker has failed!");
                        logger.update(data.genotypes[ix].name(), lik_mean);
                        writeln!(lik_writer, "{}\t{}\t{:.3}\t{:.3}", stage_ix + 1, data.genotypes[ix],
                            Ln::to_log10(lik_mean), Ln::to_log10(lik_var.sqrt())).map_err(add_path!(!))?;
                        likelihoods.likelihoods[ix] = (lik_mean, lik_var);
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
            handle.join().expect("Process failed for unknown reason")?;
        }
        Ok(())
    }
}

struct Worker<W, U> {
    data: Arc<Data>,
    rng: XoshiroRng,
    receiver: Receiver<Task>,
    sender: Sender<Solution>,
    depth_writer: W,
    aln_writer: U,
}

impl<W: Write, U: Write> Worker<W, U> {
    fn run(mut self) -> Result<(), Error> {
        let data = self.data.deref();
        let tweak = data.assgn_params.tweak.unwrap();
        let scheme = data.scheme.deref();

        // Block thread and wait for the shipment.
        while let Ok((stage_ix, task)) = self.receiver.recv() {
            let solver = &scheme.stages[stage_ix];
            let n = task.len();
            assert_ne!(n, 0, "Received empty task");

            let mut liks = vec![f64::NAN; data.assgn_params.attempts];
            let mut selected = SelectedCounter::default();
            for ix in task.into_iter() {
                let gt = &data.genotypes[ix];
                let prior = data.priors[ix];
                let mut gt_alns = GenotypeAlignments::new(gt.clone(), &data.contig_windows, &data.all_alns,
                    &data.assgn_params);
                if data.debug {
                    selected.reset(&gt_alns);
                }
                for (attempt, lik) in liks.iter_mut().enumerate() {
                    gt_alns.define_read_windows(tweak, &mut self.rng);
                    let assgns = solver.solve(&gt_alns, &mut self.rng)?;
                    *lik = prior + assgns.likelihood();
                    if data.debug {
                        let prefix = format!("{}\t{}\t{}", stage_ix + 1, gt, attempt + 1);
                        assgns.write_depth(&mut self.depth_writer, &prefix).map_err(add_path!(!))?;
                        selected.update(&assgns);
                    }
                }
                if data.debug {
                    write_alns(&mut self.aln_writer, &format!("{}\t{}", stage_ix + 1, gt),
                        &gt_alns, &data.all_alns, selected.fractions()).map_err(add_path!(!))?;
                }
                let (lik_mean, lik_var) = F64Ext::mean_variance_or_nan(&liks);
                if let Err(_) = self.sender.send((ix, lik_mean, lik_var)) {
                    log::error!("Read recruitment: main thread stopped before the child thread.");
                    break;
                }
            }
        }
        Ok(())
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
        const UPDATE_SECS: u64 = 10;
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
