//! Series of solvers.
//! Each next solver is expected to take more time, therefore is called on a smaller subset of haplotypes.

use std::{
    thread,
    io::{self, Write},
    time::{Instant, Duration},
    path::{Path, PathBuf},
    sync::{
        Arc,
        mpsc::{self, Sender, Receiver},
    },
};
use json::JsonValue;
use crate::{
    ext::{
        self,
        vec::{F64Ext, Tuples},
        rand::XoshiroRng,
    },
    math::{self, Ln, Phred},
    err::{Error, validate_param, add_path},
    bg::ser::json_get,
    seq::{ContigId, ContigNames},
    model::{
        Params,
        locs::{AllAlignments, Pair},
        windows::{ContigWindows, GenotypeWindows},
        assgn::{GenotypeAlignments, ReadAssignment, SelectedCounter},
        dp_cache::CachedDepthDistrs,
    },
};
use super::Solver;

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

fn parse_averaging_mode(val: &JsonValue) -> Result<f64, Error> {
    if val.is_null() {
        Ok(0.0)
    } else if let Some(s) = val.as_str() {
        match &s.to_lowercase() as &str {
            "min" | "minimum" | "-inf" | "-infinity" => Ok(f64::NEG_INFINITY),
            "max" | "maximum" | "inf" | "infinity" => Ok(f64::INFINITY),
            "mean" | "average" =>
                Err(Error::InvalidInput(format!("Averaging mode {} is ambiguous, please use 0 or 1", s))),
            _ => Err(Error::InvalidInput(format!("Unknown averaging mode {}", s))),
        }
    } else {
        val.as_f64().ok_or_else(||
            Error::InvalidInput(format!("Unknown averaging mode {:?}", val)))
    }
}

/// One stage of the scheme: a solver and fraction of haplotypes, on which it is executed.
#[derive(Clone)]
struct Stage {
    solver: Box<dyn Solver>,
    /// Fraction of best genotypes, used for analysis.
    fraction: f64,
    /// Write debug information about read assignments?
    write: bool,
    /// Number of attempts.
    attempts: u16,
    /// Averaging function: generalized mean with power `p`.
    aver_power: f64,
}

impl Stage {
    fn new<S: Solver + 'static>(solver: S) -> Self {
        Self {
            solver: Box::new(solver),
            fraction: 1.0,
            write: false,
            attempts: 5,
            aver_power: 1.0,
        }
    }

    fn set_fraction(&mut self, fraction: f64) -> Result<&mut Self, Error> {
        validate_param!(fraction >= 0.0 && fraction <= 1.0, "Ratio ({}) must be within [0, 1]", fraction);
        self.fraction = fraction;
        Ok(self)
    }

    fn set_write(&mut self, write: bool) -> &mut Self {
        self.write = write;
        self
    }

    fn set_tries(&mut self, attempts: u16) -> Result<&mut Self, Error> {
        validate_param!(attempts > 0, "Number of attempts must be over 1");
        self.attempts = attempts;
        Ok(self)
    }

    fn set_aver_power(&mut self, aver_power: f64) -> &mut Self {
        self.aver_power = aver_power;
        self
    }

    fn from_json(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> solver (as_str));
        let mut stage = match &solver.to_lowercase() as &str {
            "greedy" => Stage::new(super::GreedySolver::default()),
            "anneal" | "simanneal" | "simulatedannealing" => Stage::new(super::SimAnneal::default()),
            "highs" => {
                #[cfg(feature = "highs")]
                { Stage::new(super::HighsSolver::default()) }
                #[cfg(not(feature = "highs"))]
                panic!("HiGHS feature is disabled. Please recompile with `highs` feature.")
            }
            "gurobi" => {
                #[cfg(feature = "gurobi")]
                { Stage::new(super::GurobiSolver::default()) }
                #[cfg(not(feature = "gurobi"))]
                panic!("Gurobi feature is disabled. Please recompile with `gurobi` feature.")
            }
            _ => panic!("Unknown solver '{}'", solver),
        };

        json_get!(obj -> fraction? (as_f64), write? (as_bool), attempts? (as_u16));
        stage.solver.set_params(obj)?;
        if let Some(frac) = fraction {
            stage.set_fraction(frac)?;
        }
        if let Some(write) = write {
            stage.set_write(write);
        }
        if let Some(attempts) = attempts {
            stage.set_tries(attempts)?;
        }
        stage.set_aver_power(parse_averaging_mode(&obj["aver"])?);
        Ok(stage)
    }
}

/// Solver scheme.
/// Consists of multiple solvers, each executed on (possibly) smaller subset of haplotypes.
/// First stage must have fraction = 1.
#[derive(Clone)]
pub struct Scheme(Vec<Stage>);

const DEFAULT_FRAC2: f64 = 0.05;

impl Default for Scheme {
    #[cfg(feature = "highs")]
    fn default() -> Self {
        Self(vec![
            // First, run Greedy solver on all haplotypes.
            Stage {
                solver: Box::new(super::SimAnneal::default()),
                fraction: 1.0,
                write: false,
                attempts: 5,
                aver_power: f64::INFINITY, // Take maximum across all random attempts.
            },
            // Then, run HiGHS solver on the best 3%.
            Stage {
                solver: Box::new(super::HighsSolver::default()),
                fraction: DEFAULT_FRAC2,
                write: false,
                attempts: 5,
                aver_power: 1.0, // Take arithmetic mean of all likelihoods.
            },
        ])
    }

    #[cfg(all(not(feature = "highs"), feature = "gurobi"))]
    fn default() -> Self {
        Self(vec![
            // First, run Greedy solver on all haplotypes.
            Stage {
                solver: Box::new(super::SimAnneal::default()),
                fraction: 1.0,
                write: false,
                attempts: 5,
                aver_power: f64::INFINITY, // Take maximum across all random attempts.
            },
            // Then, run HiGHS solver on the best 3%.
            Stage {
                solver: Box::new(super::GurobiSolver::default()),
                fraction: DEFAULT_FRAC2,
                write: false,
                attempts: 5,
                aver_power: 1.0, // Take arithmetic mean of all likelihoods.
            },
        ])
    }

    #[cfg(all(not(feature = "highs"), not(feature = "gurobi")))]
    fn default() -> Self {
        Self(vec![
            // First, run Greedy solver on all haplotypes.
            Stage {
                solver: Box::new(super::GreedySolver::default()),
                fraction: 1.0,
                write: false,
                attempts: 5,
                aver_power: f64::INFINITY, // Take maximum across all random attempts.
            },
            // Then, run HiGHS solver on the best 3%.
            Stage {
                solver: Box::new(super::SimAnneal::default()),
                fraction: DEFAULT_FRAC2,
                write: false,
                attempts: 10,
                aver_power: 1.0, // Take arithmetic mean of all likelihoods.
            },
        ])
    }
}

const MAX_STAGES: usize = 10;
const LATIN_NUMS: [&'static str; MAX_STAGES] = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"];

impl Scheme {
    pub fn from_json(obj: &json::JsonValue) -> Result<Self, Error> {
        let json_arr = match obj {
            json::JsonValue::Array(arr) => arr,
            _ => return Err(Error::JsonLoad(format!("Failed to parse '{}': must be an array", obj))),
        };
        if json_arr.len() > MAX_STAGES {
            return Err(Error::JsonLoad(
                format!("Cannot create solver scheme: number of stages ({}) is bigger than allowed ({})",
                json_arr.len(), MAX_STAGES)));
        }
        let stages: Vec<_> = json_arr.iter().map(|v| Stage::from_json(v)).collect::<Result<_, _>>()?;
        if stages.is_empty() {
            return Err(Error::JsonLoad(format!("Failed to parse '{}': empty array is not allowed", obj)));
        }
        for (i, el) in stages.iter().enumerate() {
            if i == 0 && el.fraction < 1.0 {
                return Err(Error::JsonLoad(
                    "Cannot create solver scheme: first fraction is smaller than one".to_owned()));
            } else if i > 0 && el.fraction > stages[i - 1].fraction {
                return Err(Error::JsonLoad(
                    "Cannot create solver scheme: fractions must be non-increasing".to_owned()));
            }
        }
        Ok(Self(stages))
    }

    /// Returns the number of stages in the scheme.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> std::slice::Iter<'_, Stage> {
        self.0.iter()
    }

    /// Returns true if debug information is written at any of the stages.
    pub fn has_dbg_output(&self) -> bool {
        self.iter().any(|stage| stage.write)
    }
}

/// Various structures, needed for sample genotyping.
pub struct Data {
    /// Solving scheme.
    pub scheme: Scheme,
    pub contigs: Arc<ContigNames>,
    /// Genotypes, for which likelihoods need to be calculated.
    pub genotypes: Tuples<ContigId>,
    pub all_alns: AllAlignments,
    pub tweak: u32,
    pub debug: bool,
}

fn solve_single_thread(
    data: Data,
    mut all_gt_alns: Vec<Option<GenotypeAlignments>>,
    lik_writer: impl Write,
    mut depth_writer: impl Write,
    mut aln_writer: impl Write,
    rng: &mut XoshiroRng,
) -> Result<(), Error>
{
    let total_genotypes = data.genotypes.len();
    let mut helper = Helper::new(&data.scheme, data.contigs.tag(), lik_writer, total_genotypes)?;
    let mut rem_ixs = (0..total_genotypes).collect();

    let mut selected = SelectedCounter::default();
    for (stage_ix, stage) in data.scheme.iter().enumerate() {
        helper.start_stage(stage_ix, &mut rem_ixs, &mut all_gt_alns);
        let solver = &stage.solver;
        let mean = F64Ext::generalized_ln_mean(stage.aver_power);
        let mut liks = vec![f64::NAN; usize::from(stage.attempts)];

        for &ix in rem_ixs.iter() {
            let gt_alns = all_gt_alns[ix].as_mut().expect("Genotype alignments must be present");
            if stage.write && data.debug {
                selected.reset(gt_alns);
            }

            for (attempt, lik) in liks.iter_mut().enumerate() {
                gt_alns.define_read_windows(data.tweak, rng);
                let assgns = solver.solve(&gt_alns, rng)?;
                *lik = assgns.likelihood();
                if stage.write {
                    let prefix = format!("{}\t{}\t{}", stage_ix + 1, assgns.gt_name(), attempt + 1);
                    assgns.write_depth(&mut depth_writer, &prefix).map_err(add_path!(!))?;
                    if data.debug {
                        selected.update(&assgns);
                    }
                }
            }
            if stage.write && data.debug {
                let prefix = format!("{}\t{}", stage_ix + 1, gt_alns.gt_name());
                write_alns(&mut aln_writer, &prefix, &gt_alns, &data.all_alns, selected.fractions())
                    .map_err(add_path!(!))?;
            }
            helper.update(ix, gt_alns.gt_name(), mean(&liks))?;
        }
        helper.finish_stage();
    }
    helper.finish_overall();
    Ok(())
}

/// Creates `threads` debug files for writing read assignments, and returns their file names.
/// If `has_dbg_output` is false, returns a vector of `threads` sink objects + an empty vector of names.
fn create_debug_files(
    prefix: &Path,
    threads: u16,
    has_dbg_output: bool,
) -> Result<(Vec<Box<dyn Write + Send>>, Vec<PathBuf>), Error>
{
    if !has_dbg_output {
        return Ok(((0..threads).map(|_| Box::new(io::sink()) as Box<dyn Write + Send>).collect(), vec![]));
    }
    let filenames: Vec<_> = (0..threads).map(|i| if i == 0 {
        ext::sys::path_append(prefix, "csv.gz")
    } else {
        ext::sys::path_append(prefix, format!("{}.csv.gz", i))
    }).collect();

    let mut writers = Vec::with_capacity(usize::from(threads));
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

// `gt_alns`: Read alignments to all genotypes, in the same order as `data,genotypes`.
pub fn solve(
    data: Data,
    gt_alns: Vec<Option<GenotypeAlignments>>,
    lik_writer: impl Write,
    locus_dir: &Path,
    rng: &mut XoshiroRng,
    mut threads: u16,
) -> Result<(), Error>
{
    if usize::from(threads) > data.genotypes.len() {
        threads = data.genotypes.len() as u16;
    }
    log::info!("    [{}] Genotyping complex locus in {} stages and {} threads across {} possible genotypes",
        data.contigs.tag(), data.scheme.len(), threads, data.genotypes.len());
    let has_dbg_output = data.scheme.has_dbg_output();
    let (mut depth_writers, depth_filenames) = create_debug_files(
        &locus_dir.join("depth."), threads, has_dbg_output)?;
    writeln!(depth_writers[0], "stage\tgenotype\tattempt\t{}", ReadAssignment::DEPTH_CSV_HEADER).map_err(add_path!(!))?;

    let (mut aln_writers, aln_filenames) = create_debug_files(
        &locus_dir.join("alns."), threads, has_dbg_output && data.debug)?;
    writeln!(aln_writers[0], "stage\tgenotype\t{}", ALNS_CSV_HEADER).map_err(add_path!(!))?;

    if threads == 1 {
        solve_single_thread(data, gt_alns, lik_writer, depth_writers.pop().unwrap(), aln_writers.pop().unwrap(), rng)?;
    } else {
        let main_worker = MainWorker::new(Arc::new(data), gt_alns, rng, threads, depth_writers, aln_writers);
        main_worker.run(lik_writer)?;
    }
    merge_dbg_files(&depth_filenames)?;
    merge_dbg_files(&aln_filenames)?;
    Ok(())
}

/// Task, sent to the workers: stage index and a vector of (genotype index, genotype alignments).
type Task = (usize, Vec<(usize, GenotypeAlignments)>);
/// Tuple index, returned genotype alignments and calculated likelihood.
type Solution = (usize, GenotypeAlignments, f64);

struct MainWorker {
    data: Arc<Data>,
    gt_alns: Vec<Option<GenotypeAlignments>>,
    senders: Vec<Sender<Task>>,
    receivers: Vec<Receiver<Solution>>,
    handles: Vec<thread::JoinHandle<Result<(), Error>>>,
}

impl MainWorker {
    fn new(
        data: Arc<Data>,
        gt_alns: Vec<Option<GenotypeAlignments>>,
        rng: &mut XoshiroRng,
        threads: u16,
        depth_writers: Vec<impl Write + Send + 'static>,
        aln_writers: Vec<impl Write + Send + 'static>,
    ) -> Self
    {
        let n_workers = usize::from(threads);
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
        MainWorker { data, gt_alns, senders, receivers, handles }
    }

    fn run(mut self, lik_writer: impl Write) -> Result<(), Error> {
        let n_workers = self.handles.len();
        let scheme = &self.data.scheme;
        let contigs = &self.data.contigs;
        let genotypes = &self.data.genotypes;
        let total_genotypes = genotypes.len();
        let mut helper = Helper::new(&scheme, contigs.tag(), lik_writer, total_genotypes)?;
        let mut rem_jobs = Vec::with_capacity(n_workers);
        let mut rem_ixs = (0..total_genotypes).collect();

        for stage_ix in 0..scheme.len() {
            helper.start_stage(stage_ix, &mut rem_ixs, &mut self.gt_alns);
            rem_jobs.clear();
            let m = rem_ixs.len();
            let mut start = 0;
            for (i, sender) in self.senders.iter().enumerate() {
                if start == m {
                    break;
                }
                let rem_workers = n_workers - i;
                // Ceiling division.
                let curr_jobs = ((m - start) + rem_workers - 1) / rem_workers;
                let task = rem_ixs[start..start + curr_jobs].iter()
                    .map(|&i| (i, self.gt_alns[i].take().expect("Genotype alignments were dropped")))
                    .collect();
                sender.send((stage_ix, task)).expect("Genotyping worker has failed!");
                start += curr_jobs;
                rem_jobs.push(curr_jobs);
            }
            assert_eq!(start, m);

            while helper.solved_genotypes < m {
                for (receiver, jobs) in self.receivers.iter().zip(rem_jobs.iter_mut()) {
                    if *jobs > 0 {
                        let (ix, gt_alns, lik) = receiver.recv().expect("Genotyping worker has failed!");
                        helper.update(ix, gt_alns.gt_name(), lik)?;
                        // Do not remove assert!
                        assert!(self.gt_alns[ix].replace(gt_alns).is_none(), "Solved genotype twice");
                        *jobs -= 1;
                    }
                }
            }
            helper.finish_stage();
        }
        std::mem::drop(self.senders);
        for handle in self.handles.into_iter() {
            handle.join().expect("Process failed for unknown reason")?;
        }
        helper.finish_overall();
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
        let all_alns = &self.data.all_alns;
        let tweak = self.data.tweak;

        // Block thread and wait for the shipment.
        while let Ok((stage_ix, task)) = self.receiver.recv() {
            let stage = &self.data.scheme.0[stage_ix];
            let solver = &stage.solver;
            let debug = self.data.debug && stage.write;
            let n = task.len();
            assert_ne!(n, 0, "Received empty task");

            let mut liks = vec![f64::NAN; usize::from(stage.attempts)];
            let mut selected = SelectedCounter::default();
            let mean = F64Ext::generalized_ln_mean(stage.aver_power);
            for (ix, mut gt_alns) in task.into_iter() {
                if debug {
                    selected.reset(&gt_alns);
                }
                for (attempt, lik) in liks.iter_mut().enumerate() {
                    gt_alns.define_read_windows(tweak, &mut self.rng);
                    let assgns = solver.solve(&gt_alns, &mut self.rng)?;
                    *lik = assgns.likelihood();
                    if stage.write {
                        let prefix = format!("{}\t{}\t{}", stage_ix + 1, assgns.gt_name(), attempt + 1);
                        assgns.write_depth(&mut self.depth_writer, &prefix).map_err(add_path!(!))?;
                        if debug {
                            selected.update(&assgns);
                        }
                    }
                }
                if debug {
                    write_alns(&mut self.aln_writer, &format!("{}\t{}", stage_ix + 1, gt_alns.gt_name()),
                        &gt_alns, all_alns, selected.fractions()).map_err(add_path!(!))?;
                }
                if let Err(_) = self.sender.send((ix, gt_alns, mean(&liks))) {
                    log::error!("Read recruitment: main thread stopped before the child thread.");
                    break;
                }
            }
        }
        Ok(())
    }
}

struct Helper<'a, W> {
    scheme: &'a Scheme,
    tag: &'a str,
    lik_writer: W,

    stage_ix: usize,
    solved_genotypes: usize,
    curr_genotypes: usize,
    num_width: usize,

    likelihoods: Vec<f64>,
    best_lik: f64,
    best_ix: usize,
    best_str: String,

    timer: Instant,
    stage_start: Duration,
    last_msg: Duration,
}

impl<'a, W: Write> Helper<'a, W> {
    fn new(scheme: &'a Scheme, tag: &'a str, mut lik_writer: W, total_genotypes: usize) -> Result<Self, Error> {
        writeln!(lik_writer, "stage\tgenotype\tlik").map_err(add_path!(!))?;
        Ok(Self {
            scheme, tag, lik_writer,

            stage_ix: 0,
            solved_genotypes: 0,
            curr_genotypes: total_genotypes,
            num_width: math::num_digits(total_genotypes as f64) as usize,

            likelihoods: vec![f64::NEG_INFINITY; total_genotypes],
            best_lik: f64::NEG_INFINITY,
            best_ix: 0,
            best_str: String::new(),

            timer: Instant::now(),
            stage_start: Duration::default(),
            last_msg: Duration::default(),
        })
    }

    /// Does several things: logs the start of the stage, truncates the list of indices,
    /// and drops discarded genotype alignments.
    fn start_stage(
        &mut self,
        stage_ix: usize,
        rem_ixs: &mut Vec<usize>,
        gt_alns: &mut [Option<GenotypeAlignments>],
    ) {
        self.stage_start = self.timer.elapsed();
        self.stage_ix = stage_ix;
        let stage = &self.scheme.0[stage_ix];
        let total_genotypes = self.likelihoods.len();
        // Consider at least one tuple irrespective of the fraction.
        self.curr_genotypes = ((total_genotypes as f64 * stage.fraction).ceil() as usize)
            .clamp(1, self.curr_genotypes);
        if self.curr_genotypes < rem_ixs.len() {
            rem_ixs.sort_unstable_by(|&i, &j| self.likelihoods[j].total_cmp(&self.likelihoods[i]));
            for &i in rem_ixs[self.curr_genotypes..].iter() {
                // Do not remove this assert!
                assert!(gt_alns[i].take().is_some(), "Genotype alignments are already dropped");
            }
            rem_ixs.truncate(self.curr_genotypes);
            self.num_width = math::num_digits(self.curr_genotypes as f64) as usize;
        }

        self.solved_genotypes = 0;
        self.last_msg = self.stage_start;
        log::info!("    [{}] Stage {:>3}.  {}.  {} genotypes, {} attempts, averaging power: {:.0}",
            self.tag, LATIN_NUMS[stage_ix], stage.solver,
            self.curr_genotypes, stage.attempts, stage.aver_power);
    }

    fn update(&mut self, ix: usize, gt_name: &str, lik: f64) -> Result<(), Error> {
        writeln!(self.lik_writer, "{}\t{}\t{:.3}", self.stage_ix + 1, gt_name, Ln::to_log10(lik))
            .map_err(add_path!(!))?;
        let stored_lik = &mut self.likelihoods[ix];
        *stored_lik = stored_lik.max(lik);
        if lik > self.best_lik {
            self.best_lik = lik;
            self.best_ix = ix;
            self.best_str = gt_name.to_owned();
        }
        self.solved_genotypes += 1;
        let now_dur = self.timer.elapsed();
        // Update frequency in seconds.
        const UPDATE_SECS: u64 = 10;
        if (now_dur - self.last_msg).as_secs() >= UPDATE_SECS {
            let speed = (now_dur.as_secs_f64() - self.stage_start.as_secs_f64()) / self.solved_genotypes as f64;
            log::debug!("        [{:width$}/{},  {:.4} s/tuple]  Best: {} -> {:11.2}", self.solved_genotypes,
                self.curr_genotypes, speed, self.best_str, Ln::to_log10(self.best_lik), width = self.num_width);
            self.last_msg = now_dur;
        }
        Ok(())
    }

    fn finish_stage(&mut self) {
        self.last_msg = self.timer.elapsed();
        log::info!("    [{}] Stage {:>3} finished in {}.  Best: {} -> {:11.2}", self.tag, LATIN_NUMS[self.stage_ix],
            ext::fmt::Duration(self.last_msg - self.stage_start), self.best_str, Ln::to_log10(self.best_lik));
    }

    /// Returns index of the best tuple, as well as its quality.
    fn finish_overall(mut self) -> (usize, f64) {
        let norm_fct = Ln::sum(&self.likelihoods);
        self.likelihoods.iter_mut().for_each(|v| *v -= norm_fct);
        let quality = Phred::from_likelihoods(&mut self.likelihoods, self.best_ix);
        log::info!("    [{}] Overall best: {}  (Quality = {:.1}, Confidence = {:.4}%)", self.tag, self.best_str,
            quality, self.likelihoods[self.best_ix].exp() * 100.0);
        (self.best_ix, quality)
    }
}
