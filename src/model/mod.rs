pub mod windows;
pub mod locs;
pub mod assgn;
pub mod dp_cache;

use crate::{
    math::Ln,
    err::{Error, validate_param},
    bg::err_prof,
    ext::vec::UsizeOrInf,
};

/// Read depth model parameters.
#[derive(Clone, Debug)]
pub struct Params {
    /// Boundary size: ignore left- and right-most `boundary_size` bp.
    pub boundary_size: u32,
    /// Read depth likelihood contribution, relative to read alignment likelihood.
    pub depth_contrib: f64,
    /// For each read pair, all alignments to a specific genotype,
    /// less probable than `best_prob - prob_diff` are discarded.
    pub prob_diff: f64,
    /// Unmapped reads receive minimum alignment probability for this read PLUS this penalty.
    pub unmapped_penalty: f64,
    /// Edit distance p-value thresholds, first for good alignments, second for passable.
    pub edit_pvals: (f64, f64),

    /// See `windows::WeightCalculator` for more details.
    pub weight_breakpoint: f64,
    pub weight_power: f64,
    /// Ignore reads and windows with weight under this value.
    pub min_weight: f64,

    /// Randomly move read middle by `tweak` bp into one of the directions.
    /// None: half window size.
    pub tweak: Option<u32>,
    /// Alternative hypotheses copy number values (main hypothesis is 1).
    pub alt_cn: (f64, f64),
    /// Use paired-end reads that do not have a mapped mate within any contig.
    pub use_unpaired: bool,
    /// Use at most this number of alignments (subsample if necessary).
    pub max_alns: usize,

    /// Minimal number of genotypes after each step of solving.
    pub min_gts: UsizeOrInf,
    /// Discard low likelihood genotypes based on this threshold (`min + score_thresh * (max - min)`).
    pub score_thresh: f64,
    /// Discard genotypes with low probability to be best (after each step).
    pub prob_thresh: f64,
    /// Number of attempts with different tweak sizes.
    pub attempts: usize,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            boundary_size: 200,
            depth_contrib: 50.0,
            prob_diff: Ln::from_log10(5.0),
            unmapped_penalty: Ln::from_log10(-5.0),
            edit_pvals: err_prof::DEF_EDIT_PVAL,

            weight_breakpoint: 0.2,
            weight_power: 2.0,
            min_weight: 0.001,

            tweak: None,
            alt_cn: (0.5, 1.5),
            use_unpaired: false,
            max_alns: 500_000,

            min_gts: UsizeOrInf(10),
            score_thresh: 0.9,
            prob_thresh: Ln::from_log10(-4.0),
            attempts: 5,
        }
    }
}

impl Params {
    pub fn validate(&mut self) -> Result<(), Error> {
        validate_param!(self.boundary_size > 0, "Boundary size ({}) cannot be zero.", self.boundary_size);
        validate_param!(!self.prob_diff.is_nan() && self.prob_diff.is_finite(),
            "Unexpected probability difference ({:.4}) value", self.prob_diff);
        validate_param!(self.depth_contrib > 0.0,
            "Depth contribution ({}) must be positive", self.depth_contrib);
        self.prob_diff = self.prob_diff.abs();
        if self.prob_diff < 1.0 {
            log::warn!("Note that probability difference ({}) is in log-10 space", Ln::to_log10(self.prob_diff));
        }
        validate_param!(self.unmapped_penalty < 0.0 && self.unmapped_penalty.is_normal(),
            "Unmapped penalty ({:.4}) must be negative", Ln::to_log10(self.unmapped_penalty));
        for &pval in &[self.edit_pvals.0, self.edit_pvals.1] {
            validate_param!(0.0 < pval && pval < 0.5, "p-value threshold ({}) must be within (0, 0.5)", pval);
        }
        validate_param!(self.edit_pvals.1 <= self.edit_pvals.0,
            "Second p-value threshold ({}) must not be greater than the first threshold ({})",
            self.edit_pvals.1, self.edit_pvals.0);

        validate_param!(0.0 < self.weight_breakpoint && self.weight_breakpoint < 1.0,
            "Weight breakpoint ({}) must be within (0, 1)", self.weight_breakpoint);
        validate_param!(0.5 <= self.weight_power && self.weight_power <= 50.0,
            "Weight power ({}) must be within [0.5, 50]", self.weight_power);
        validate_param!(0.0 <= self.min_weight && self.min_weight <= 0.5,
            "Minimal weight ({}) must be within [0, 0.5].", self.min_weight);

        validate_param!(self.alt_cn.0 > 0.0 && self.alt_cn.0 < 1.0,
            "Alternative copy number #1 ({}) must be in (0, 1).", self.alt_cn.0);
        validate_param!(self.alt_cn.1 > 1.0,
            "Alternative copy number #2 ({}) must be over 1.", self.alt_cn.1);

        validate_param!(self.max_alns > 0, "Number of alignments must be positive");
        if self.max_alns < 1000 {
            log::warn!("Number of alignments ({}) is too small", self.max_alns);
        }
        validate_param!(0.0 <= self.score_thresh && self.score_thresh <= 1.0,
            "Score threshold ({}) must be within [0, 1]", self.score_thresh);
        validate_param!(self.prob_thresh.is_normal() && self.prob_thresh < 0.0,
            "Probability threshold ({}) must be negative (log10-space)", Ln::to_log10(self.prob_thresh));
        validate_param!(self.attempts > 0, "Number of attempts must be positive");
        if self.attempts < 3 {
            log::warn!("At least 3 attempts are advised, otherwise qualities are inaccurate!");
        }
        validate_param!(self.min_gts.0 > 1,
            "Minimal number of genotypes ({}) must be at least 2", self.min_gts);
        Ok(())
    }

    pub fn set_tweak_size(&mut self, window_size: u32) -> Result<(), Error> {
        const MAX_THEORETICAL_TWEAK: u32 = u16::MAX as u32 / 2 - 1;
        const MAX_TWEAK: u32 = 100;
        if self.tweak.is_none() {
            self.tweak = Some((window_size / 3).min(MAX_TWEAK).min(self.boundary_size.saturating_sub(1)));
        }

        let tweak = self.tweak.unwrap();
        validate_param!(tweak < self.boundary_size, "Boundary size ({}) must be greater than tweak size ({}).",
            self.boundary_size, tweak);
        validate_param!(tweak < MAX_THEORETICAL_TWEAK, "Tweaking size ({}) is too large (max = {})",
            tweak, MAX_THEORETICAL_TWEAK);
        if tweak > MAX_TWEAK {
            log::warn!("Tweaking size ({}) may be too large", tweak);
        }
        Ok(())
    }
}
