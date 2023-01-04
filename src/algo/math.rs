//! Common math functions.
//! Everywhere in comments `log` represents natural log.

use crate::algo::vec_ext::*;

/// Constant log(10).
pub const LOG10: f64 = 2.302585092994045684_f64;
/// Constant 1 / log(10).
pub const INV_LOG10: f64 = 0.4342944819032518277_f64;

pub struct Ln;

impl Ln {
    /// Converts log10 value into natural log value.
    #[inline]
    pub fn from_log10(l10: f64) -> f64 {
        l10 * LOG10
    }

    /// Converts natural log value into log10 value.
    #[inline]
    pub fn to_log10(ln: f64) -> f64 {
        ln * INV_LOG10
    }

    /// Calculates *log(exp(a) + exp(b))*.
    pub fn add(a: f64, b: f64) -> f64 {
        if a >= b {
            b + (a - b).exp().ln_1p()
        } else {
            a + (b - a).exp().ln_1p()
        }
    }

    /// Calculates *log(exp(a) - exp(b))*.
    pub fn sub(a: f64, b: f64) -> f64 {
        let c = a - b;
        debug_assert!(c >= 0.0);
        b + (c.exp() - 1.0).ln()
    }

    /// Calculates *log(sum(exp(values)))*.
    pub fn sum(values: &[f64]) -> f64 {
        let m = values.min();
        let s = values.iter().fold(0.0_f64, |acc, v| acc + (v - m).exp());
        m + s.ln()
    }
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