//! Series of solvers.
//! Each next solver is expected to take more time, therefore is called on a smaller subset of haplotypes.

use std::{
    thread,
    io::Write,
    time::Instant,
    cmp::min,
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

/// One stage of the scheme: a solver and ratio of haplotypes, on which it is executed.
#[derive(Clone)]
struct Stage {
    solver: Box<dyn Solver>,
    ratio: f64,
}

impl Stage {
    fn with_default<S: Solver + Default + 'static>(ratio: f64) -> Result<Self, Error> {
        Self::new(Box::new(S::default()), ratio)
    }

    fn new(solver: Box<dyn Solver>, ratio: f64) -> Result<Self, Error> {
        validate_param!(ratio >= 0.0 && ratio <= 1.0, "Ratio ({}) must be within [0, 1].", ratio);
        Ok(Self { solver, ratio })
    }

    fn from_json(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> solver (as_str));
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

        let ratio = if obj.has_key("ratio") {
            json_get!(obj -> ratio (as_f64));
            ratio
        } else {
            1.0
        };
        if obj.has_key("params") {
            solver.set_params(&obj["params"])?;
        }
        Self::new(solver, ratio)
    }
}

/// Solver scheme.
/// Consists of multiple solvers, each executed on (possibly) smaller subset of haplotypes.
/// First stage must have ratio = 1.
#[derive(Clone)]
pub struct Scheme(Vec<Stage>);

impl Default for Scheme {
    fn default() -> Self {
        Self(vec![
            // First, run Greedy solver on all haplotypes.
            Stage::with_default::<super::GreedySolver>(1.0).unwrap(),
            // Then, run HiGHS solver on the best 3%.
            Stage::with_default::<super::HighsSolver>(0.03).unwrap(),
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
            if i == 0 && el.ratio < 1.0 {
                return Err(Error::JsonLoad(
                    "Cannot create solver scheme: first ratio is smaller than one".to_owned()));
            } else if i > 0 && el.ratio > stages[i - 1].ratio {
                return Err(Error::JsonLoad(
                    "Cannot create solver scheme: ratios must be non-increasing".to_owned()));
            }
        }
        Ok(Self(stages))
    }

    fn solve_single_thread(
        &self,
        all_alns: AllPairAlignments,
        contig_windows: Vec<ContigWindows>,
        contigs: &Arc<ContigNames>,
        cached_distrs: &CachedDepthDistrs,
        tuples: Tuples<ContigId>,
        mut writer: impl Write,
        params: &Params,
        rng: &mut XoshiroRng,
    ) -> Result<(), Error>
    {
        let tag = contigs.tag();
        let n = tuples.len();
        let width2 = 12 * tuples.tup_len();
        let mut rem_ixs: Vec<_> = (0..n).collect();
        // Best likelihoods and read assignments for each tuple.
        // Do not store actual best assignments.
        // Instead, if needed, we can store Rng states, in order to reconstruct solutions later.
        let mut likelihoods = vec![f64::NEG_INFINITY; n];
        let mut best_lik = f64::NEG_INFINITY;
        let mut best_ix = 0;
        let mut best_str = String::new();

        for (stage_ix, stage) in self.0.iter().enumerate() {
            let m = filter_ixs(&mut rem_ixs, &likelihoods, stage.ratio);
            let width1 = math::num_digits(m as u64);
            let solver = &stage.solver;
            let timer = Instant::now();
            log::info!("    [{}] Stage {}. {}  ({} tuples, 1 thread)", tag, LATIN_NUMS[stage_ix], solver, m);
            for (i, &ix) in rem_ixs.iter().enumerate() {
                let mcontig_windows = MultiContigWindows::new(&tuples[ix], &contig_windows);
                let contigs_str = mcontig_windows.ids_str(contigs);
                let mut assgn = ReadAssignment::new(mcontig_windows, &all_alns, cached_distrs, params);
                let lik = solver.solve(&mut assgn, rng)?;
                writeln!(writer, "{}\t{}\t{:.3}", stage_ix + 1, contigs_str, Ln::to_log10(lik))?;

                let stored_lik = &mut likelihoods[ix];
                *stored_lik = stored_lik.max(lik);
                log::debug!("        [{:width1$} / {}] {:width2$}  -> {:11.2}", i + 1, m, contigs_str,
                    Ln::to_log10(*stored_lik));

                if lik > best_lik {
                    best_lik = lik;
                    best_ix = ix;
                    best_str = contigs_str;
                }
            }
            log::info!("    [{}] Stage {} finished in {}. Best: {} ({:.2})", tag, LATIN_NUMS[stage_ix],
                ext::fmt::Duration(timer.elapsed()), best_str, Ln::to_log10(best_lik))
        }
        let norm_fct = Ln::sum(&likelihoods);
        likelihoods.iter_mut().for_each(|v| *v -= norm_fct);
        log::info!("    [{}] Overall best: {}  (Quality = {:.1}, Norm.Likelihood = {})",
            tag, best_str, Phred::from_likelihoods(&mut likelihoods, best_ix), Ln::to_log10(likelihoods[best_ix]));
        Ok(())
    }

    fn solve_multi_thread(
        &self,
        all_alns: AllPairAlignments,
        contig_windows: Vec<ContigWindows>,
        contigs: &Arc<ContigNames>,
        cached_distrs: &Arc<CachedDepthDistrs>,
        tuples: Tuples<ContigId>,
        writer: impl Write,
        params: &Params,
        rng: &mut XoshiroRng,
        threads: u16,
    ) -> Result<(), Error>
    {
        let main_worker = MainWorker::new(self.clone(), Arc::new(all_alns), Arc::new(contig_windows),
            Arc::clone(&contigs), cached_distrs, Arc::new(tuples), writer, params, rng, threads);
        main_worker.run()
    }

    pub fn solve(
        &self,
        all_alns: AllPairAlignments,
        contig_windows: Vec<ContigWindows>,
        contigs: &Arc<ContigNames>,
        cached_distrs: &Arc<CachedDepthDistrs>,
        tuples: Tuples<ContigId>,
        writer: impl Write,
        params: &Params,
        rng: &mut XoshiroRng,
        threads: u16,
    ) -> Result<(), Error>
    {
        if threads == 1 {
            self.solve_single_thread(all_alns, contig_windows, contigs, &cached_distrs, tuples, writer, params, rng)
        } else {
            self.solve_multi_thread(all_alns, contig_windows, contigs, cached_distrs, tuples, writer,
                params, rng, threads)
        }
    }
}

///Keep only the best `ratio` indices (out of the total number).
fn filter_ixs(rem_ixs: &mut Vec<usize>, likelihoods: &[f64], ratio: f64) -> usize {
    let n = likelihoods.len();
    // Consider at least one tuple irrespective of the ratio.
    let m = ((n as f64 * ratio).ceil() as usize).clamp(1, n);
    if m < rem_ixs.len() {
        rem_ixs.sort_unstable_by(|&i, &j| likelihoods[j].total_cmp(&likelihoods[i]));
        rem_ixs.truncate(m);
    }
    m
}

/// Task, sent to the workers: stage index + tuple indices to solve.
type Task = (usize, Vec<usize>);
type Solution = Vec<f64>;

struct MainWorker<W> {
    scheme: Scheme,
    contigs: Arc<ContigNames>,
    tuples: Arc<Tuples<ContigId>>,
    senders: Vec<Sender<Task>>,
    receivers: Vec<Receiver<Solution>>,
    handles: Vec<thread::JoinHandle<Result<(), Error>>>,
    writer: W,
}

impl<W: Write> MainWorker<W> {
    fn new(
        scheme: Scheme,
        all_alns: Arc<AllPairAlignments>,
        contig_windows: Arc<Vec<ContigWindows>>,
        contigs: Arc<ContigNames>,
        cached_distrs: &Arc<CachedDepthDistrs>,
        tuples: Arc<Tuples<ContigId>>,
        writer: W,
        params: &Params,
        rng: &mut XoshiroRng,
        threads: u16,
    ) -> Self {
        let n_workers = usize::from(threads);
        let mut senders = Vec::with_capacity(n_workers);
        let mut receivers = Vec::with_capacity(n_workers);
        let mut handles = Vec::with_capacity(n_workers);

        for _ in 0..n_workers {
            let (task_sender, task_receiver) = mpsc::channel();
            let (sol_sender, sol_receiver) = mpsc::channel();
            let worker = Worker {
                scheme: scheme.clone(),
                rng: rng.clone(),
                all_alns: Arc::clone(&all_alns),
                contig_windows: Arc::clone(&contig_windows),
                tuples: Arc::clone(&tuples),
                cached_distrs: Arc::clone(&cached_distrs),
                params: params.clone(),
                receiver: task_receiver,
                sender: sol_sender,
            };
            rng.jump();
            senders.push(task_sender);
            receivers.push(sol_receiver);
            handles.push(thread::spawn(|| worker.run()));
        }
        MainWorker { scheme, contigs, tuples, senders, receivers, handles, writer }
    }

    fn run(mut self) -> Result<(), Error> {
        let n_workers = self.handles.len();
        let tag = self.contigs.tag();
        let n = self.tuples.len();
        let mut rem_ixs: Vec<_> = (0..n).collect();
        // Best likelihoods and read assignments for each tuple.
        // Do not store actual best assignments.
        // Instead, if needed, we can store Rng states, in order to reconstruct solutions later.
        let mut likelihoods = vec![f64::NEG_INFINITY; n];
        let mut best_lik = f64::NEG_INFINITY;
        let mut best_ix = 0;
        let mut best_str = String::new();

        for (stage_ix, stage) in self.scheme.0.iter().enumerate() {
            let timer = Instant::now();
            let m = filter_ixs(&mut rem_ixs, &likelihoods, stage.ratio);
            log::info!("    [{}] Stage {}. {}  ({} tuples, {} threads)", tag, LATIN_NUMS[stage_ix], stage.solver,
                m, n_workers);

            // Ceiling division.
            let per_worker = (m + n_workers - 1) / n_workers;
            for (i, sender) in self.senders.iter().enumerate() {
                let start = per_worker * i;
                if start < m {
                    let end = min(start + per_worker, m);
                    sender.send((stage_ix, rem_ixs[start..end].to_vec())).expect("Recruitment worker has failed!");
                }
            }

            for (i, receiver) in self.receivers.iter().enumerate() {
                let start = per_worker * i;
                if start < m {
                    let end = min(start + per_worker, m);
                    let worker_likelihoods = receiver.recv().expect("Recruitment worker has failed!");
                    assert_eq!(worker_likelihoods.len(), end - start);
                    for (&ix, &lik) in rem_ixs[start..end].iter().zip(&worker_likelihoods) {
                        let contigs_str = self.contigs.get_names(self.tuples[ix].iter().copied());
                        writeln!(self.writer, "{}\t{}\t{:.3}", stage_ix + 1, contigs_str, Ln::to_log10(lik))?;
                        let stored_lik = &mut likelihoods[ix];
                        *stored_lik = stored_lik.max(lik);
                        if lik > best_lik {
                            best_lik = lik;
                            best_ix = ix;
                            best_str = contigs_str;
                        }
                    }
                }
            }

            log::info!("    [{}] Stage {} finished in {}. Best: {} ({:.2})", tag, LATIN_NUMS[stage_ix],
                ext::fmt::Duration(timer.elapsed()), best_str, Ln::to_log10(best_lik))
        }
        let norm_fct = Ln::sum(&likelihoods);
        likelihoods.iter_mut().for_each(|v| *v -= norm_fct);
        log::info!("    [{}] Overall best: {}  (Quality = {:.1}, Norm.Likelihood = {})",
            tag, best_str, Phred::from_likelihoods(&mut likelihoods, best_ix), Ln::to_log10(likelihoods[best_ix]));

        std::mem::drop(self.senders);
        for handle in self.handles.into_iter() {
            handle.join().expect("Process failed for unknown reason")?;
        }
        Ok(())
    }
}

struct Worker {
    scheme: Scheme,
    rng: XoshiroRng,
    all_alns: Arc<AllPairAlignments>,
    contig_windows: Arc<Vec<ContigWindows>>,
    tuples: Arc<Tuples<ContigId>>,
    cached_distrs: Arc<CachedDepthDistrs>,
    params: Params,
    receiver: Receiver<Task>,
    sender: Sender<Solution>,
}

impl Worker {
    fn run(mut self) -> Result<(), Error> {
        // Block thread and wait for the shipment.
        while let Ok((stage_ix, task)) = self.receiver.recv() {
            let n = task.len();
            assert_ne!(n, 0);
            let all_alns = &self.all_alns;
            let solver = &self.scheme.0[stage_ix].solver;
            let tuples = &self.tuples;
            let contig_windows = &self.contig_windows;
            let mut likelihoods = Vec::with_capacity(n);

            for ix in task.into_iter() {
                let mcontig_windows = MultiContigWindows::new(&tuples[ix], contig_windows);
                let mut assgn = ReadAssignment::new(mcontig_windows, all_alns, &self.cached_distrs, &self.params);
                likelihoods.push(solver.solve(&mut assgn, &mut self.rng)?);
            }
            if let Err(_) = self.sender.send(likelihoods) {
                log::error!("Read recruitment: main thread stopped before the child thread.");
                break;
            }
        }
        Ok(())
    }
}
