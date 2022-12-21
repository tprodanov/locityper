//! Common math functions.
//! Everywhere in comments `log` represents natural log.

/// Constant log(10).
pub const LOG10: f64 = 2.302585092994045684_f64;
/// Constant 1 / log(10).
pub const INV_LOG10: f64 = 0.4342944819032518277_f64;

/// Convert log10 value to natural log value.
#[inline]
pub fn log10_to_ln(v: f64) -> f64 {
    v * LOG10
}

/// Convert natural log value to log10 value.
#[inline]
pub fn ln_to_log10(v: f64) -> f64 {
    v * INV_LOG10
}

/// Static Phred class with useful conversion methods.
pub struct Phred;

impl Phred {
    /// Convert regular probability to a Phred quality.
    pub fn from_prob(prob: f64) -> f64 {
        debug_assert!(0.0 <= prob && prob <= 1.0);
        -10.0 * prob.log10()
    }

    /// Convert log probability to a Phred quality.
    pub fn from_lprob(lprob: f64) -> f64 {
        debug_assert!(lprob <= 0.0);
        -10.0 * INV_LOG10 * lprob
    }

    /// Convert phred quality into regular probability (0 - 1).
    pub fn to_prob(phred: f64) -> f64 {
        debug_assert!(phred >= 0.0);
        10.0_f64.powf(-0.1 * phred)
    }

    /// Convert phred quality into regular probability (0 - 1).
    pub fn to_lprob(phred: f64) -> f64 {
        debug_assert!(phred >= 0.0);
        -0.1 * phred * LOG10
    }
}