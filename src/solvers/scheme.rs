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
use crate::{
    ext::{
        self,
        vec::Tuples,
        rand::XoshiroRng,
    },
    math::{self, Ln, Phred},
    err::{Error, validate_param},
    bg::ser::json_get,
    seq::{ContigId, ContigNames},
    model::{
        Params,
        locs::AllPairAlignments,
        windows::{ContigWindows, MultiContigWindows},
        assgn::ReadAssignment,
        dp_cache::CachedDepthDistrs,
    },
};
use super::Solver;

/// One stage of the scheme: a solver and fraction of haplotypes, on which it is executed.
#[derive(Clone)]
struct Stage {
    solver: Box<dyn Solver>,
    /// Fraction of best tuples, used for analysis.
    fraction: f64,
    /// Write debug information about read assignments?
    write: bool,
}

impl Stage {
    fn with_default<S: Solver + Default + 'static>(fraction: f64, write: bool) -> Result<Self, Error> {
        Self::new(Box::new(S::default()), fraction, write)
    }

    fn new(solver: Box<dyn Solver>, fraction: f64, write: bool) -> Result<Self, Error> {
        validate_param!(fraction >= 0.0 && fraction <= 1.0, "Ratio ({}) must be within [0, 1].", fraction);
        Ok(Self { solver, fraction, write })
    }

    fn from_json(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> solver (as_str), fraction? (as_f64), write? (as_bool));
        let mut solver: Box<dyn Solver> = match &solver.to_lowercase() as &str {
            "greedy" => Box::new(super::GreedySolver::default()),
            "simanneal" | "simulatedannealing" => Box::new(super::SimAnneal::default()),
            "highs" =>
                if cfg!(feature = "highs") {
                    Box::new(super::HighsSolver::default())
                } else {
                    panic!("HiGHS feature is disabled. Consider recompiling with `highs` feature enabled.");
                },
            "gurobi" =>
                if cfg!(feature = "gurobi") {
                    Box::new(super::GurobiSolver::default())
                } else {
                    panic!("Gurobi feature is disabled. Consider recompiling with `gurobi` feature enabled.");
                },
            _ => panic!("Unknown solver '{}'", solver),
        };

        if obj.has_key("params") {
            solver.set_params(&obj["params"])?;
        }
        Self::new(solver, fraction.unwrap_or(1.0), write.unwrap_or(false))
    }
}

/// Solver scheme.
/// Consists of multiple solvers, each executed on (possibly) smaller subset of haplotypes.
/// First stage must have fraction = 1.
#[derive(Clone)]
pub struct Scheme(Vec<Stage>);

impl Default for Scheme {
    fn default() -> Self {
        Self(vec![
            // First, run Greedy solver on all haplotypes.
            Stage::with_default::<super::GreedySolver>(1.0, false).unwrap(),
            // Then, run HiGHS solver on the best 3%.
            Stage::with_default::<super::HighsSolver>(0.03, false).unwrap(),
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
    pub scheme: Scheme,
    pub all_alns: AllPairAlignments,
    pub contig_windows: Vec<ContigWindows>,
    pub contigs: Arc<ContigNames>,
    pub cached_distrs: Arc<CachedDepthDistrs>,
    pub tuples: Tuples<ContigId>,
    pub params: Params,
    pub debug: bool,
}

fn solve_single_thread(
    data: Data,
    lik_writer: impl Write,
    mut depth_writer: impl Write,
    mut aln_writer: impl Write,
    rng: &mut XoshiroRng,
) -> Result<(), Error>
{
    let total_tuples = data.tuples.len();
    let mut helper = Helper::new(&data.scheme, data.contigs.tag(), lik_writer, total_tuples)?;
    let mut rem_ixs = (0..total_tuples).collect();

    for (stage_ix, stage) in data.scheme.iter().enumerate() {
        helper.start_stage(stage_ix, &mut rem_ixs);
        let solver = &stage.solver;
        for &ix in rem_ixs.iter() {
            let mcontig_windows = MultiContigWindows::new(&data.tuples[ix], &data.contig_windows);
            let mut assgn = ReadAssignment::new(mcontig_windows, &data.all_alns, &data.cached_distrs, &data.params);
            assgn.define_read_windows(data.params.tweak, rng);
            let lik = solver.solve(&mut assgn, rng)?;
            let contigs_str = data.contigs.get_names(data.tuples[ix].iter().copied());
            if stage.write {
                let prefix = format!("{}\t{}", stage_ix + 1, contigs_str);
                assgn.write_depth(&mut depth_writer, &prefix)?;
                if data.debug {
                    assgn.write_alns(&mut aln_writer, &prefix, &data.all_alns)?;
                }
            }
            helper.update(ix, contigs_str, lik)?;
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
) -> io::Result<(Vec<Box<dyn Write + Send>>, Vec<PathBuf>)>
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
fn merge_dbg_files(filenames: &[PathBuf]) -> io::Result<()> {
    if filenames.len() > 1 {
        // By this point, all depth_writers should be already dropped.
        let mut file1 = std::fs::OpenOptions::new().append(true).open(&filenames[0])?;
        ext::sys::concat_files(filenames[1..].iter(), &mut file1)?;
        filenames[1..].iter().map(std::fs::remove_file).collect::<io::Result<()>>()?;
    }
    Ok(())
}

pub fn solve(
    data: Data,
    lik_writer: impl Write,
    locus_dir: &Path,
    rng: &mut XoshiroRng,
    threads: u16,
) -> Result<(), Error>
{
    log::info!("    [{}] Genotyping complex locus in {} stages and {} threads across {} possible tuples",
        data.contigs.tag(), data.scheme.len(), threads, data.tuples.len());
    let has_dbg_output = data.scheme.has_dbg_output();
    let (mut depth_writers, depth_filenames) = create_debug_files(
        &locus_dir.join("depth."), threads, has_dbg_output)?;
    writeln!(depth_writers[0], "stage\tgenotype\t{}", ReadAssignment::DEPTH_CSV_HEADER)?;

    let (mut aln_writers, aln_filenames) = create_debug_files(
        &locus_dir.join("alns."), threads, has_dbg_output && data.debug)?;
    writeln!(aln_writers[0], "stage\tgenotype\t{}", ReadAssignment::ALNS_CSV_HEADER)?;

    if threads == 1 {
        solve_single_thread(data, lik_writer, depth_writers.pop().unwrap(), aln_writers.pop().unwrap(), rng)?;
    } else {
        let main_worker = MainWorker::new(Arc::new(data), rng, threads, depth_writers, aln_writers);
        main_worker.run(lik_writer)?;
    }
    merge_dbg_files(&depth_filenames)?;
    merge_dbg_files(&aln_filenames)?;
    Ok(())
}

/// Task, sent to the workers: stage index + tuple indices to solve.
type Task = (usize, Vec<usize>);
/// Tuple index and likelihood.
type Solution = (usize, f64);

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
        MainWorker { data, senders, receivers, handles }
    }

    fn run(self, lik_writer: impl Write) -> Result<(), Error> {
        let n_workers = self.handles.len();
        let total_tuples = self.data.tuples.len();
        let scheme = &self.data.scheme;
        let tuples = &self.data.tuples;
        let contigs = &self.data.contigs;
        let mut helper = Helper::new(&scheme, self.data.contigs.tag(), lik_writer, total_tuples)?;
        let mut rem_jobs = vec![0; n_workers];
        let mut rem_ixs = (0..total_tuples).collect();

        for stage_ix in 0..scheme.len() {
            helper.start_stage(stage_ix, &mut rem_ixs);
            rem_jobs.fill(0);
            let m = rem_ixs.len();
            let mut start = 0;
            for (i, (sender, jobs)) in self.senders.iter().zip(rem_jobs.iter_mut()).enumerate() {
                if start == m {
                    break;
                }
                let rem_workers = n_workers - i;
                // Ceiling division.
                *jobs = ((m - start) + rem_workers - 1) / rem_workers;
                sender.send((stage_ix, rem_ixs[start..start + *jobs].to_vec()))
                    .expect("Genotyping worker has failed!");
                start += *jobs;
            }
            assert_eq!(start, m);

            while helper.solved_tuples < m {
                for (receiver, jobs) in self.receivers.iter().zip(rem_jobs.iter_mut()) {
                    if *jobs > 0 {
                        let (ix, lik) = receiver.recv().expect("Genotyping worker has failed!");
                        let contigs_str = contigs.get_names(tuples[ix].iter().copied());
                        helper.update(ix, contigs_str, lik)?;
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
        // Block thread and wait for the shipment.
        while let Ok((stage_ix, task)) = self.receiver.recv() {
            let n = task.len();
            assert_ne!(n, 0);
            let all_alns = &self.data.all_alns;
            let stage = &self.data.scheme.0[stage_ix];
            let solver = &stage.solver;
            let tuples = &self.data.tuples;
            let contigs = &self.data.contigs;
            let contig_windows = &self.data.contig_windows;
            let cached_distrs = &self.data.cached_distrs;
            let params = &self.data.params;

            for ix in task.into_iter() {
                let mcontig_windows = MultiContigWindows::new(&tuples[ix], contig_windows);
                let mut assgn = ReadAssignment::new(mcontig_windows, all_alns, cached_distrs, params);
                assgn.define_read_windows(params.tweak, &mut self.rng);
                let lik = solver.solve(&mut assgn, &mut self.rng)?;
                if stage.write {
                    let prefix = format!("{}\t{}", stage_ix + 1, contigs.get_names(tuples[ix].iter().copied()));
                    assgn.write_depth(&mut self.depth_writer, &prefix)?;
                    if self.data.debug {
                        assgn.write_alns(&mut self.aln_writer, &prefix, all_alns)?;
                    }
                }

                if let Err(_) = self.sender.send((ix, lik)) {
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
    solved_tuples: usize,
    curr_tuples: usize,
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
    fn new(scheme: &'a Scheme, tag: &'a str, mut lik_writer: W, total_tuples: usize) -> io::Result<Self> {
        writeln!(lik_writer, "stage\tgenotype\tlik")?;
        Ok(Self {
            scheme, tag, lik_writer,

            stage_ix: 0,
            solved_tuples: 0,
            curr_tuples: total_tuples,
            num_width: math::num_digits(total_tuples as u64),

            likelihoods: vec![f64::NEG_INFINITY; total_tuples],
            best_lik: f64::NEG_INFINITY,
            best_ix: 0,
            best_str: String::new(),

            timer: Instant::now(),
            stage_start: Duration::default(),
            last_msg: Duration::default(),
        })
    }

    fn start_stage(&mut self, stage_ix: usize, rem_ixs: &mut Vec<usize>) {
        self.stage_start = self.timer.elapsed();
        self.stage_ix = stage_ix;
        let stage = &self.scheme.0[stage_ix];
        let total_tuples = self.likelihoods.len();
        // Consider at least one tuple irrespective of the fraction.
        self.curr_tuples = ((total_tuples as f64 * stage.fraction).ceil() as usize).clamp(1, self.curr_tuples);
        if self.curr_tuples < rem_ixs.len() {
            rem_ixs.sort_unstable_by(|&i, &j| self.likelihoods[j].total_cmp(&self.likelihoods[i]));
            rem_ixs.truncate(self.curr_tuples);
            self.num_width = math::num_digits(self.curr_tuples as u64);
        }

        self.solved_tuples = 0;
        self.last_msg = self.stage_start;
        log::info!("    [{}] Stage {:>3}.  {}  ({} tuples)", self.tag, LATIN_NUMS[stage_ix], stage.solver,
            self.curr_tuples);
    }

    fn update(&mut self, ix: usize, tuple_str: String, lik: f64) -> io::Result<()> {
        writeln!(self.lik_writer, "{}\t{}\t{:.3}", self.stage_ix + 1, &tuple_str, Ln::to_log10(lik))?;
        let stored_lik = &mut self.likelihoods[ix];
        *stored_lik = stored_lik.max(lik);
        if lik > self.best_lik {
            self.best_lik = lik;
            self.best_ix = ix;
            self.best_str = tuple_str;
        }
        self.solved_tuples += 1;
        let now_dur = self.timer.elapsed();
        const UPDATE_FREQ: u64 = 2;
        if (now_dur - self.last_msg).as_secs() >= UPDATE_FREQ {
            let speed = (now_dur.as_secs_f64() - self.stage_start.as_secs_f64()) / self.solved_tuples as f64;
            log::debug!("        [{:width$}/{},  {:.4} s/tuple]  Best: {} -> {:11.2}", self.solved_tuples,
                self.curr_tuples, speed, self.best_str, Ln::to_log10(self.best_lik), width = self.num_width);
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
