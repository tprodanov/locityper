pub mod windows;
pub mod locs;
pub mod assgn;
pub mod distr_cache;
pub mod bam;

use std::{
    fmt,
    str::FromStr,
};
use crate::{
    math::{
        Ln,
        distr::bayes,
    },
    err::{validate_param},
    seq::kmers::Kmer,
    bg::Technology,
};
use windows::WeightCalculator;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Polarity {
    Best,
    Worst,
}

impl Polarity {
    pub fn to_str(self) -> &'static str {
        match self {
            Self::Best => "best",
            Self::Worst => "worst",
        }
    }
}

impl FromStr for Polarity {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match &s.to_lowercase() as &str {
            "best" => Ok(Self::Best),
            "worst" => Ok(Self::Worst),
            _ => Err(format!("Unknown polarity {:?}", s)),
        }
    }
}

impl fmt::Display for Polarity {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(self.to_str())
    }
}

pub fn default_unmapped_penalty(tech: Technology) -> (Polarity, f64) {
    match tech {
        Technology::Illumina => (Polarity::Best, Ln::from_log10(-10.0)),
        Technology::HiFi | Technology::PacBio | Technology::Nanopore => (Polarity::Worst, Ln::from_log10(-20.0)),
    }
}

/// Read depth model parameters.
#[derive(Clone, Debug)]
pub struct Params {
    /// Boundary size: ignore left- and right-most `boundary_size` bp.
    pub boundary_size: u32,
    /// Read depth likelihood contribution, relative to read alignment likelihood, must be in [0, 2].
    pub lik_skew: f64,
    /// For each read pair, all alignments to a specific genotype,
    /// less probable than `best_prob - prob_diff` are discarded.
    pub prob_diff: f64,
    /// Unmapped reads are penalized by float value either relative to the best or to the worst retained alignment.
    pub unmapped_penalty: (Polarity, f64),
    /// Calculate linguistic complexity as fraction of non-repetitive k-mers.
    pub complexity_k: u8,
    /// Poor complexity threshold.
    pub poor_compl: f64,
    /// Edit distance threshold as a fraction of read length in regions of low complexity.
    pub poor_compl_edit: f64,
    /// Two window weight calculators: one for linguistic complexity of the window
    /// and another for the fraction of region-specific k-mers.
    pub compl_weight_calc: Option<WeightCalculator>,
    pub kmers_weight_calc: Option<WeightCalculator>,
    /// Soft and hard thresholds for the number of unique k-mers per read/read pair.
    pub kmer_soft_thresh: u16,
    pub kmer_hard_thresh: u16,
    /// Ignore reads and windows with weight under this value.
    pub min_weight: f64,

    /// Randomly move read middle by `tweak` bp into one of the directions.
    /// None: half window size.
    pub tweak: Option<u32>,
    /// Alternative hypotheses copy number values (main hypothesis is 1).
    pub default_alt_cn: &'static str,
    pub alt_cn: Vec<f64>,

    /// Discard low likelihood genotypes with alignment likelihood smaller than `best - filt_diff`.
    pub filt_diff: f64,
    /// Discard genotypes with low probability to be best (after each step).
    pub prob_thresh: f64,

    /// Do not skip intermediate solvers.
    pub dont_skip: bool,
    /// Generate BAM files for this many best genotypes.
    pub out_bams: usize,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            boundary_size: 200,
            lik_skew: 0.85,
            prob_diff: Ln::from_log10(10.0),
            unmapped_penalty: (Polarity::Best, f64::NAN),
            complexity_k: 5,
            compl_weight_calc: Some(WeightCalculator::new(0.5, 4.0).unwrap()),
            kmers_weight_calc: Some(WeightCalculator::new(0.2, 4.0).unwrap()),
            kmer_soft_thresh: 5,
            kmer_hard_thresh: 1,
            poor_compl: 0.5,
            poor_compl_edit: 0.7,
            min_weight: 0.001,

            tweak: None,
            default_alt_cn: "0.3,2,3,4,5",
            alt_cn: vec![0.3, 2.0, 3.0, 4.0, 5.0],

            filt_diff: Ln::from_log10(100.0),
            prob_thresh: Ln::from_log10(-4.0),
            dont_skip: false,
            out_bams: 0,
        }
    }
}

/// If tweak size is bigger than this, there may appear problems with overflowing.
const MAX_ALLOWED_TWEAK: u32 = u16::MAX as u32 / 2 - 1;
/// Automatically select tweak size as 1/2 of the window size.
const AUTO_TWEAK_MULT: f64 = 0.5;
/// Automatically crop tweak size to 200 bp.
const AUTO_TWEAK_MAX: u32 = 200;

impl Params {
    pub fn validate(&mut self) -> crate::Result<()> {
        validate_param!(self.boundary_size > 0, "Boundary size ({}) cannot be zero.", self.boundary_size);
        validate_param!(!self.prob_diff.is_nan() && self.prob_diff.is_finite(),
            "Unexpected probability difference ({:.4}) value", self.prob_diff);
        validate_param!(self.lik_skew >= -1.0 + 1e-10 && self.lik_skew <= 1.0 - 1e-10,
            "Likelihood skew ({}) must be within (-1, 1)", self.lik_skew);
        self.prob_diff = self.prob_diff.abs();
        if self.prob_diff < 1.0 {
            log::warn!("Note that probability difference ({}) is in log-10 space", Ln::to_log10(self.prob_diff));
        }
        validate_param!(0.0 <= self.min_weight && self.min_weight <= 0.5,
            "Minimal weight ({}) must be within [0, 0.5].", self.min_weight);
        validate_param!(0 < self.complexity_k && self.complexity_k <= u32::MAX_KMER_SIZE,
            "Complexity k size ({}) must be between 1 and {}", self.complexity_k, u32::MAX_KMER_SIZE);
        validate_param!(0.0 <= self.poor_compl && self.poor_compl <= 1.0,
            "Complexity threshold ({}) must be between 0 and 1.", self.poor_compl);
        validate_param!(0.0 <= self.poor_compl_edit && self.poor_compl_edit <= 1.0,
            "Edit distance threshold ({}) in low complexity regions must be between 0 and 1.", self.poor_compl_edit);

        for &cn in &self.alt_cn {
            validate_param!(cn > 0.0 && cn != 1.0,
                "Alternative copy number ({}) must be positive and not equal to one.", cn);
        }
        validate_param!(self.alt_cn.len() <= bayes::N_ALTS, "Too many alternative CNs ({}, max {})",
            self.alt_cn.len(), bayes::N_ALTS);

        validate_param!(self.kmer_hard_thresh <= self.kmer_soft_thresh,
            "Hard k-mer threshold ({}) should not be greater than the soft k-mer threshold ({})",
            self.kmer_hard_thresh, self.kmer_soft_thresh);

        validate_param!(self.filt_diff >= 0.0,
            "Filtering likelihood difference ({}) must be non-negative (log10-space)", Ln::to_log10(self.filt_diff));
        validate_param!(self.prob_thresh < 0.0,
            "Probability threshold ({}) must be negative (log10-space)", Ln::to_log10(self.prob_thresh));
        Ok(())
    }

    pub fn set_tweak_size(&mut self, window_size: u32) -> crate::Result<()> {
        if self.tweak.is_none() {
            self.tweak = Some((
                (window_size as f64 * AUTO_TWEAK_MULT).round() as u32)
                    .min(AUTO_TWEAK_MAX)
                    .min(self.boundary_size.saturating_sub(1))
            );
        }

        let tweak = self.tweak.unwrap();
        validate_param!(tweak < self.boundary_size, "Boundary size ({}) must be greater than tweak size ({}).",
            self.boundary_size, tweak);
        validate_param!(tweak <= MAX_ALLOWED_TWEAK, "Tweaking size ({}) is too large (max = {})",
            tweak, MAX_ALLOWED_TWEAK);
        if tweak > AUTO_TWEAK_MAX {
            log::warn!("Tweaking size ({}) may be too large", tweak);
        }
        Ok(())
    }
}
