pub mod nbinom;

/// Constant log(10).
pub const LOG10: f64 = 2.302585092994045684_f64;
/// Constant 1 / log(10).
pub const INV_LOG10: f64 = 0.4342944819032518277_f64;

/// Returns the smaller value between `a` and `b`.
/// However, if only one of them equals -INF, returns the other one.
fn non_inf_min(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        b
    } else if b == f64::NEG_INFINITY {
        a
    } else {
        a.min(b)
    }
}

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
            if b == f64::NEG_INFINITY { a } else { b + (a - b).exp().ln_1p() }
        } else {
            if a == f64::NEG_INFINITY { b } else { a + (b - a).exp().ln_1p() }
        }
    }

    /// Calculates *log(exp(a) - exp(b))*.
    pub fn sub(a: f64, b: f64) -> f64 {
        if b == f64::NEG_INFINITY {
            a
        } else {
            let c = a - b;
            assert!(c >= 0.0, "Ln::sub({}, {}) impossible!", a, b);
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
                let m = a.iter().map(|v| f(v)).fold(f64::INFINITY, non_inf_min);
                let s = a.iter().fold(0.0_f64, |acc, el| acc + (f(&el) - m).exp());
                m + s.ln()
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
                let m = a.iter().map(|v| f(v)).fold(init, non_inf_min);
                let s = a.iter().fold((init - m).exp(), |acc, el| acc + (f(&el) - m).exp());
                m + s.ln()
            },
        }
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