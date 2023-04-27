//! Series of solvers.
//! Each next solver is expected to take more time, therefore is called on a smaller subset of haplotypes.

use std::{
    rc::Rc,
    time::Instant,
};
use crate::{
    ext,
    math::{self, Ln},
    err::{Error, validate_param},
    bg::ser::json_get,
    seq::{ContigId, ContigNames},
    model::{
        Params,
        locs::AllPairAlignments,
        windows::{ContigWindows, MultiContigWindows},
        assgn::{CachedDepthDistrs, ReadAssignment},
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
        json_get!(obj -> solver_name (as_str));
        let mut solver: Box<dyn Solver> = match &solver_name.to_lowercase() as &str {
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
            _ => panic!("Unknown solver '{}'", solver_name),
        };

        let ratio = if obj.has_key("ratio") {
            json_get!(obj -> val (as_f64));
            val
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
            // Then, run HiGHS solver on the best 1%.
            Stage::with_default::<super::HighsSolver>(0.01).unwrap(),
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
        all_alns: &AllPairAlignments,
        contig_windows: &[ContigWindows],
        contigs: &Rc<ContigNames>,
        cached_distrs: &CachedDepthDistrs,
        tuples: &ext::vec::Tuples<ContigId>,
        params: &Params,
        rng: &mut ext::rand::XoshiroRng,
    ) -> Result<(), Error>
    {
        let tag = contigs.tag();
        let n = tuples.len();
        let mut width1 = math::num_digits(n as u64);
        let width2 = 12 * tuples.tup_len();
        let mut rem_ixs: Vec<_> = (0..n).collect();
        // Best likelihoods and read assignments for each tuple.
        // Do not store actual best assignments.
        // Instead, if needed, we can store Rng states, in order to reconstruct solutions later.
        let mut likelihoods = vec![f64::NAN; n];
        let mut best_lik = f64::NEG_INFINITY;
        let mut best_str = "???".to_owned();

        for (stage_ix, stage) in self.0.iter().enumerate() {
            // Consider at least one tuple irrespective of the ratio.
            let m = ((n as f64 * stage.ratio).ceil() as usize).clamp(1, n);
            if m < rem_ixs.len() {
                rem_ixs.sort_unstable_by(|&i, &j| likelihoods[j].total_cmp(&likelihoods[i]));
                rem_ixs.truncate(m);
                width1 = math::num_digits(m as u64);
            }

            let solver = &stage.solver;
            let timer = Instant::now();
            log::info!("    [{}] Stage {}. {}  ({} tuples, 1 thread)", tag, LATIN_NUMS[stage_ix], solver, m);
            for (i, &ix) in rem_ixs.iter().enumerate() {
                let mcontig_windows = MultiContigWindows::new(&tuples[ix], contig_windows);
                let contigs_str = mcontig_windows.ids_str(contigs);
                let mut assgn = ReadAssignment::new(mcontig_windows, all_alns, cached_distrs, params);
                let lik = solver.solve(&mut assgn, rng)?;
                let stored_lik = &mut likelihoods[ix];
                *stored_lik = stored_lik.max(lik);
                log::debug!("        [{:width1$} / {}] {:width2$}  -> {:11.2}", i + 1, m, contigs_str,
                    Ln::to_log10(*stored_lik));

                if lik > best_lik {
                    best_lik = lik;
                    best_str = contigs_str;
                }
            }
            log::info!("    [{}] Stage {} finished in {}. Best: {} ({:.2})", tag, LATIN_NUMS[stage_ix],
                ext::fmt::Duration(timer.elapsed()), best_str, Ln::to_log10(best_lik))
        }
        Ok(())
    }

    pub fn solve(
        &self,
        all_alns: &AllPairAlignments,
        contig_windows: &[ContigWindows],
        contigs: &Rc<ContigNames>,
        cached_distrs: &CachedDepthDistrs,
        tuples: &ext::vec::Tuples<ContigId>,
        params: &Params,
        rng: &mut ext::rand::XoshiroRng,
        threads: u16,
    ) -> Result<(), Error>
    {
        if threads == 1 {
            self.solve_single_thread(all_alns, contig_windows, contigs, cached_distrs, tuples, params, rng)
        } else {
            unimplemented!("Multi-thread solution is not implemented")
        }
    }
}