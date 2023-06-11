pub mod windows;
pub mod locs;
pub mod assgn;
pub mod dp_cache;

use crate::{
    math::Ln,
    err::{Error, validate_param},
};

/// Read depth model parameters.
#[derive(Clone, Debug)]
pub struct Params {
    /// Boundary size: ignore left- and right-most `boundary_size` bp.
    pub boundary_size: u32,
    /// Read depth likelihood contribution, relative to read alignment likelihood.
    pub depth_contrib: f64,
    /// For each read pair, all alignments less probable than `best_prob - prob_diff` are discarded.
    pub prob_diff: f64,
    /// Unmapped reads receive this penalty.
    pub unmapped_penalty: f64,

    /// Average k-mer frequency is calculated for a window in question.
    /// If the value does not exceed `rare_kmer`, the window received a weight = 1.
    /// If the value equals to `semicommon_kmer`, weight would be 0.5.
    pub rare_kmer: f64,
    pub semicommon_kmer: f64,

    /// Randomly move read middle by `tweak` bp into one of the directions.
    pub tweak: u32,
    /// Alternative hypotheses copy number values (main hypothesis is 1).
    pub alt_cn: (f64, f64),
}

impl Default for Params {
    fn default() -> Self {
        Self {
            boundary_size: 200,
            depth_contrib: 2.0,
            prob_diff: Ln::from_log10(5.0),
            unmapped_penalty: Ln::from_log10(-10.0),
            rare_kmer: 3.0,
            semicommon_kmer: 5.0,

            tweak: 70,
            alt_cn: (0.01, 2.0),
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
        validate_param!(self.rare_kmer > 1.0, "First k-mer frequency threshold ({:.4}) must be over 1",
            self.rare_kmer);
        validate_param!(self.rare_kmer < self.semicommon_kmer,
            "k-mer frequency thresholds ({:.4}, {:.4}) are non-increasing", self.rare_kmer, self.semicommon_kmer);

        validate_param!(self.tweak < self.boundary_size, "Boundary size ({}) must be greater than tweak size ({}).",
            self.boundary_size, self.tweak);
        validate_param!(self.alt_cn.0 > 0.0 && self.alt_cn.0 < 1.0,
            "Alternative copy number #1 ({}) must be in (0, 1).", self.alt_cn.0);
        validate_param!(self.alt_cn.1 > 1.0,
            "Alternative copy number #2 ({}) must be over 1.", self.alt_cn.1);
        Ok(())
    }
}
