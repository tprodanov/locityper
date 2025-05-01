pub mod distr;
pub mod frac;

use statrs::distribution::{StudentsT, ContinuousCDF};

pub struct Ln;
pub use frac::Fraction;

impl Ln {
    /// log(10).
    pub const LN10: f64 = 2.302585092994045684_f64;
    /// 1 / log(10).
    pub const INV_LN10: f64 = 0.4342944819032518277_f64;

    /// Converts log10 value into natural log value.
    #[inline]
    pub fn from_log10(l10: f64) -> f64 {
        l10 * Self::LN10
    }

    /// Converts natural log value into log10 value.
    #[inline]
    pub fn to_log10(ln: f64) -> f64 {
        ln * Self::INV_LN10
    }

    /// Calculates *log(exp(a) + exp(b))*.
    pub fn add(a: f64, b: f64) -> f64 {
        if a >= b {
            if b == f64::NEG_INFINITY { a } else { b + (a - b).exp().ln_1p() }
        } else {
            if a == f64::NEG_INFINITY { b } else { a + (b - a).exp().ln_1p() }
        }
    }

    /// Calculates *log(exp(a) - exp(b))*.
    #[allow(unused)]
    pub fn sub(a: f64, b: f64) -> f64 {
        if b == f64::NEG_INFINITY {
            a
        } else {
            let c = a - b;
            assert!(c >= 0.0, "Ln::sub({}, {}) is impossible!", a, b);
            b + (c.exp() - 1.0).ln()
        }
    }

    /// Calculates logsumexp: *log(sum(exp(values)))*.
    pub fn sum(values: &[f64]) -> f64 {
        Ln::map_sum(values, f64::clone)
    }

    /// Similar to `Ln::sum`, but starts summation with an `init` value (ln-space).
    /// Calculates *log( exp(init) + sum(exp(values)) )*.
    pub fn sum_init(values: &[f64], init: f64) -> f64 {
        Ln::map_sum_init(values, f64::clone, init)
    }

    /// Calculates logsumexp, but applies function `f` to each element of `a` before calculating sum.
    /// WARNING: Function `f` is applied to each element twice!
    pub fn map_sum<T, F: FnMut(&T) -> f64>(a: &[T], mut f: F) -> f64 {
        match a.len() {
            0 => f64::NEG_INFINITY,
            1 => f(&a[0]),
            _ => {
                let m = a.iter().map(|v| f(v)).fold(f64::NEG_INFINITY, f64::max);
                if m.is_infinite() {
                    m
                } else {
                    let s = a.iter().fold(0.0_f64, |acc, el| acc + (f(&el) - m).exp());
                    m + s.ln()
                }
            },
        }
    }

    /// Similar to `Ln::map_sum`, but start summation with an `init` value (ln-space).
    /// Calculates *log( exp(init) + sum(exp(f(a))) )*
    pub fn map_sum_init<T, F: FnMut(&T) -> f64>(a: &[T], mut f: F, init: f64) -> f64 {
        match a.len() {
            0 => init,
            1 => Ln::add(init, f(&a[0])),
            _ => {
                let m = a.iter().map(|v| f(v)).fold(init, f64::max);
                if m.is_infinite() {
                    m
                } else {
                    let s = a.iter().fold((init - m).exp(), |acc, el| acc + (f(&el) - m).exp());
                    m + s.ln()
                }
            },
        }
    }
}

/// Static Phred class with useful conversion methods.
pub struct Phred;

impl Phred {
    /// Convert regular probability to a Phred quality.
    #[allow(unused)]
    pub fn from_prob(prob: f64) -> f64 {
        debug_assert!(0.0 <= prob && prob <= 1.0);
        -10.0 * prob.log10()
    }

    /// Convert log probability to a Phred quality.
    pub fn from_ln_prob(lprob: f64) -> f64 {
        debug_assert!(lprob <= 0.0);
        -10.0 * Ln::to_log10(lprob)
    }

    /// Convert phred quality into regular probability (0 - 1).
    #[allow(unused)]
    pub fn to_prob(phred: f64) -> f64 {
        debug_assert!(phred >= 0.0);
        10.0_f64.powf(-0.1 * phred)
    }

    /// Convert phred quality into regular probability (0 - 1).
    #[allow(unused)]
    pub fn to_ln_prob(phred: f64) -> f64 {
        debug_assert!(phred >= 0.0);
        -0.1 * Ln::from_log10(phred)
    }

    /// Returns Phred score ln-likelihood with index `ix` across all `likelihoods`.
    #[allow(unused)]
    pub fn from_likelihoods(likelihoods: &mut [f64], ix: usize) -> f64 {
        let stored = likelihoods[ix];
        likelihoods[ix] = f64::NEG_INFINITY;
        let error_prob = Ln::sum(likelihoods);
        likelihoods[ix] = stored;
        Self::from_ln_prob(error_prob)
    }
}

/// Index of the leading in the 10-base representation of x.
pub fn num_digits(x: f64) -> i32 {
    x.abs().log10().floor() as i32 + 1
}

/// Round the number to the corresponding number of significant digits.
pub fn round_signif(x: impl Into<f64>, digits: u8) -> f64 {
    assert!(digits > 0, "Number of significant digits cannot be zero");
    let x = x.into();
    if x == 0.0 {
        0.0
    } else {
        let shift = num_digits(x) - i32::from(digits);
        let fct = 10_f64.powi(shift);
        (x / fct).round() * fct
    }
}

/// Converts number into string with the corresponding number of significant digits.
/// Is not very fast.
pub fn fmt_signif(x: impl Into<f64>, digits: u8) -> String {
    assert!(digits > 0, "Number of significant digits cannot be zero");
    let x = x.into();
    if x == 0.0 {
        "0".to_owned()
    } else {
        let shift = num_digits(x) - i32::from(digits);
        let fct = 10_f64.powi(shift);
        if shift < 0 {
            format!("{:.prec$}", x, prec = (-shift) as usize)
                .trim_end_matches('0').trim_end_matches('.').to_owned()
        } else {
            (((x / fct).round() * fct).round() as i64).to_string()
        }
    }
}

/// Computes unpaired t-test p-value.
/// Returns probability of observing t-statistic under the null hypothesis: `mean1 >= mean2`.
/// Null hipothesis: first mean is equal or larger than the second mean.
/// Sample sizes should be the same.
pub fn unpaired_onesided_t_test<const EQ_VAR: bool>(
    mean1: f64,
    var1: f64,
    mean2: f64,
    var2: f64,
    n: f64,
) -> f64
{
    let var_sum = var1 + var2;
    let t_stat = (mean1 - mean2) * (n / var_sum).sqrt();
    let freedom = if EQ_VAR {
        2.0 * n - 2.0
    } else {
        (n - 1.0) * var_sum * var_sum / (var1 * var1 + var2 * var2)
    };
    let t_distr = StudentsT::new(0.0, 1.0, freedom).expect("Could not create Student's t distribution");
    t_distr.cdf(t_stat)
}

/// Same as `unpaired_onesided_t_test`, but allowing for different sample sizes.
pub fn unpaired_onesided_t_test_diffsizes<const EQ_VAR: bool>(
    mean1: f64,
    var1: f64,
    mean2: f64,
    var2: f64,
    n1: f64,
    n2: f64,
) -> f64
{
    let nvar1 = var1 / n1;
    let nvar2 = var2 / n2;
    let sum_nvar = nvar1 + nvar2;
    let t_stat = (mean1 - mean2) / sum_nvar.sqrt();
    let freedom = if EQ_VAR {
        n1 + n2 - 2.0
    } else {
        sum_nvar * sum_nvar / (nvar1 * nvar1 / (n1 - 1.0) + nvar2 * nvar2 / (n2 - 1.0))
    };
    let t_distr = StudentsT::new(0.0, 1.0, freedom).expect("Could not create Student's t distribution");
    t_distr.cdf(t_stat)
}

/// Rounding divisions.
pub trait RoundDiv {
    /// Correct division, without overflows.
    #[allow(unused)]
    fn correct_ceil_div(self, rhs: Self) -> Self;

    /// Fast ceiling division. Can give incorrect results due to overflowing.
    fn fast_ceil_div(self, rhs: Self) -> Self;

    /// Fast rounding division. Can give incorrect results due to overflowing.
    fn fast_round_div(self, rhs: Self) -> Self;
}

impl<T> RoundDiv for T
where T: num_traits::PrimInt + num_traits::Unsigned + num_traits::ConstZero + num_traits::ConstOne
{
    #[inline]
    fn correct_ceil_div(self, rhs: Self) -> Self {
        let d = self / rhs;
        if self % rhs == Self::ZERO { d } else { d + Self::ONE }
    }

    #[inline]
    fn fast_ceil_div(self, rhs: Self) -> Self {
        (self + rhs - Self::ONE) / rhs
    }

    #[inline]
    fn fast_round_div(self, rhs: Self) -> Self {
        (self + (rhs >> 1)) / rhs
    }
}

/// Linear extrapolation. `x1` must not be equal to `x2`.
#[allow(unused)]
pub fn interpolate((x1, x2): (f64, f64), (y1, y2): (f64, f64), x: f64) -> f64 {
    x * (y2 - y1) / (x2 - x1) + y1
}

/// Logical operator `a => b` = `!a || b`.
#[inline(always)]
pub fn implies(a: bool, b: bool) -> bool {
    !a || b
}
