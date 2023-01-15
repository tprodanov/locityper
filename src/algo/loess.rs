//! Module to calculate local regression (LOESS / LOWESS).

use std::cmp::min;
use itertools::izip;
use nalgebra::{DVector, DMatrix, SVD};
use crate::algo::{
    vec_ext::*,
    bisect,
};

/// Calculates local regression (LOESS / LOWESS).
/// Given input arrays `x`, `y` and, optionally weights `w`, finds best `y_out` values for each `xout` value.
/// If the `xout` array is not set, finds best `y_out` value for each input `x` value.
///
/// The structure allows to set custom `frac` option (2/3 by default):
/// For each `xout` point, looks at `frac * x.len()` input points.
///
/// Additionally, the structure allows to set the degree of the local polynomials (1 by default).
///
/// Usage:
/// ```rust
/// Loess::new().set_degree(2).calculate(&x, &y);
/// Loess::new().set_frac(0.3).set_xout(xout).calculate_weighted(&x, &y, &w);
/// ```
#[derive(Debug, Clone)]
pub struct Loess {
    frac: f64,
    degree: usize,
    xout: Option<Vec<f64>>,
}

impl Loess {
    /// Creates a new Loess class with default parameters.
    pub fn new(frac: f64, degree: usize) -> Loess {
        assert!(frac > 0.0 && frac <= 1.0, "Fraction must be in (0, 1].");
        Loess {
            frac,
            degree,
            xout: None,
        }
    }

    pub fn frac(&self) -> f64 {
        self.frac
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn xout(&self) -> Option<&[f64]> {
        match &self.xout {
            Some(xout) => Some(&xout),
            None => None,
        }
    }

    pub fn set_xout(&mut self, xout: Vec<f64>) -> &mut Loess {
        assert!(!xout.is_empty());
        self.xout = Some(xout);
        self
    }

    pub fn calculate(&self, x: &[f64], y: &[f64]) -> Vec<f64> {
        loess(x, y, None, self.xout().unwrap_or(x), self.frac, self.degree)
    }

    pub fn calculate_weighted(&self, x: &[f64], y: &[f64], w: &[f64]) -> Vec<f64> {
        loess(x, y, Some(w), self.xout().unwrap_or(x), self.frac, self.degree)
    }
}

impl Default for Loess {
    #[inline]
    fn default() -> Self {
        Self::new(2.0 / 3.0, 1)
    }
}

/// Calculates local regression (LOESS / LOWESS).
/// Given input arrays `x`, `y` and, optionally weights `w`, finds best `y` values for each `xout` value.
/// For each `xout` point, looks at `frac * x.len()` input points.
/// Uses local polynomials of degree `deg`.
fn loess(x: &[f64], y: &[f64], w: Option<&[f64]>, xout: &[f64], frac: f64, deg: usize) -> Vec<f64> {
    let n = x.len();
    assert!(n == y.len(), "Cannot calculate LOESS on vectors of different length ({} and {})", n, y.len());
    let ixs = x.argsort();
    let x = x.reorder(&ixs);
    let y = y.reorder(&ixs);
    let w = w.map(|values| values.reorder(&ixs));

    let n_frac = (n as f64 * frac).round().max(1.0) as usize;
    let range = x[n - 1] - x[0];
    assert!(range > 0.0, "Cannot calculate LOESS: x contains a single value {}", x[0]);

    let mut y_out = Vec::with_capacity(xout.len());
    for &xval in xout.iter() {
        let mut a = bisect::left(&x, &xval);
        let mut b = bisect::right_at(&x, &xval, a, n);
        let curr_n = b - a;
        if curr_n >= n_frac {
            y_out.push(y[a..b].iter().sum::<f64>() / curr_n as f64);
            continue
        }

        let rem = n_frac - curr_n;
        let left;
        let right;
        if a < n - b {
            left = min(a, rem / 2);
            right = min(n - b, rem - left);
        } else {
            right = min(n - b, rem / 2);
            left = min(a, rem - right);
        }
        a -= left;
        b += right;

        let sub_x = &x[a..b];
        let sub_y = &y[a..b];
        // Calculates (sub_x - xval) / range.
        let norm_sub_x = sub_x.sc_add_mul(-xval, 1.0 / range);
        let mut weight = tricube_kernel(&norm_sub_x);

        if let Some(in_weight) = &w {
            weight.ew_mul_assign(&in_weight[a..b]);
        }
        let coefs = polyfit(sub_x, sub_y, &weight, deg).unwrap();
        y_out.push(polyval(&coefs, xval));
    }
    y_out
}

/// Fit a polynomial of degree `deg` to points `x`, `y`. Points are weighted with the slice `w`.
/// Returns a vector of coefficients, starting with x^0 coefficient and ending with x^deg.
pub fn polyfit(x: &[f64], y: &[f64], w: &[f64], deg: usize) ->
        Result<Vec<f64>, &'static str> {
    let nrow = x.len();
    let ncol = deg + 1;
    let mut a = DMatrix::zeros(nrow, ncol);
    let mut b = DVector::zeros(nrow);

    for (i, (&xi, &yi, &wi)) in izip!(x, y, w).enumerate() {
        a[(i, 0)] = wi;
        b[i] = yi * wi;

        let mut xpow = 1.0;
        for j in 1..ncol {
            xpow *= xi;
            a[(i, j)] = xpow * wi;
        }
    }

    let decomp = SVD::new(a, true, true);
    decomp.solve(&b, 1e-18).map(|v| v.data.into())
}

/// Given coefficients a0, a1, a2, ..., an calculates
/// a0 + a1 * x + a2 * x^2 + ... + an * x^n.
pub fn polyval(coefs: &[f64], x: f64) -> f64 {
    let mut curr_x = 1.0;
    let mut y = 0.0;
    for coef in coefs.iter() {
        y += coef * curr_x;
        curr_x *= x;
    }
    y
}

/// Calculates tricube kernel: **70/81 (1 - |v|^3)^3**.
/// Returns 0 for values out of range `[-1, 1]`.
pub fn tricube_kernel(values: &[f64]) -> Vec<f64> {
    const MULT: f64 = 70.0 / 81.0;
    let mut y = Vec::with_capacity(values.len());
    for &v in values {
        y.push(MULT * (1.0 - v.abs().min(1.0).powi(3)).powi(3))
    }
    y
}
