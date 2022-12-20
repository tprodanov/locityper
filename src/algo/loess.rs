use std::cmp::{Ordering, min};

use nalgebra::{DVector, DMatrix, SVD};

// pub struct Builder {
//     x: DVector<f64>,
//     y: DVector<f64>,
//     w: Option<DVector<f64>>,
//     frac: f64,
//     degree: u8,
// }


// def loess(x, y, xout, frac=2/3, deg=1, w=None):
//     ixs = np.argsort(x)
//     x = x[ixs]
//     y = y[ixs]
//     in_weight = w[ixs] if w is not None else None

//     n = len(x)
//     n_frac = int(round(n * frac))
//     size = x[-1] - x[0]
//     assert size > 0.0

//     res = np.full(len(xout), np.nan)
//     for i, xval in enumerate(xout):
//         a = x.searchsorted(xval, 'left')
//         b = x.searchsorted(xval, 'right')
//         if b - a >= n_frac:
//             res[i] = np.mean(y[a:b])
//             continue

//         if b - a < n_frac:
//             rem = n_frac - b + a
//             left = min(a, (rem // 2))
//             right = min(n - b, rem - left)
//             a -= left
//             b += right
//             assert a >= 0 and b <= n

//         sub_x = x[a:b]
//         sub_y = y[a:b]
//         weight = common.tricube_kernel((sub_x - xval) / size)
//         if in_weight is not None:
//             weight *= in_weight[a:b]
//         coef = np.polyfit(sub_x, sub_y, deg=deg, w=weight)
//         res[i] = np.polyval(coef, xval)
//     return res

/// Reorders vector by taking elements in the order of `ixs`.
pub fn reorder(x: &DVector<f64>, ixs: &[usize]) -> DVector<f64> {
    let n = ixs.len();
    let mut x2 = DVector::zeros(n);
    for i in 0..n {
        x2[i] = x[ixs[i]];
    }
    x2
}

/// Performs binary search and finds the index `i` such that `v[i-1] < target <= v[i]`.
pub fn binary_search_left(v: &[f64], target: f64) -> usize {
    let mut low = 0;
    let mut high = v.len();
    while low < high {
        let mid = (low + high) / 2;
        match v[mid].total_cmp(&target) {
            Ordering::Equal if mid == 0 || v[mid - 1] < target => return mid,
            Ordering::Less => low = mid + 1,
            _ => high = mid,
        }
    }
    low
}

/// Performs binary search and finds the index `i` such that `v[i-1] <= target < v[i]`.
pub fn binary_search_right(v: &[f64], target: f64) -> usize {
    let n = v.len();
    let mut low = 0;
    let mut high = n;
    while low < high {
        let mid = (low + high) / 2;
        if v[mid] > target {
            if mid == 0 || v[mid - 1] == target {
                return mid;
            }
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    low
}

/// Calculates tricube kernel: **70/81 (1 - |v|^3)^3**.
/// Returns 0 for values out of range `[-1, 1]`.
pub fn tricube_kernel(values: &[f64]) -> DVector<f64> {
    const MULT: f64 = 70.0 / 81.0;
    let n = values.len();
    let mut res = DVector::zeros(n);
    for (i, &v) in values.iter().enumerate() {
        res[i] = MULT * (1.0 - v.abs().min(1.0).powi(3)).powi(3);
    }
    res
}

pub fn loess(x: &DVector<f64>, y: &DVector<f64>, w: Option<&DVector<f64>>,
        xout: &DVector<f64>, frac: f64, degree: usize) -> DVector<f64> {
    let n = x.len();
    assert!(n == y.len(), "Cannot calculate LOESS on vectors of different length ({} and {})", n, y.len());
    let mut ixs: Vec<usize> = (0..n).collect();
    ixs.sort_by(|&i, &j| x[i].total_cmp(&x[j]));
    let x = reorder(x, &ixs);
    let y = reorder(y, &ixs);
    let w = match w {
        Some(values) => Some(reorder(values, &ixs)),
        None => None,
    };

    let n_frac = (n as f64 * frac).round() as usize;
    let range = x[n - 1] - x[0];
    assert!(range > 0.0, "Cannot calculate LOESS: x contains a single value {}", x[0]);

    let mut res = DVector::zeros(xout.len());
    for (i, &xval) in xout.iter().enumerate() {
        let mut a = binary_search_left((&x).into(), xval);
        let mut b = binary_search_right((&x).into(), xval);
        let curr_n = b - a;
        if curr_n >= n_frac {
            res[i] = x.index((a..b, 0)).sum() / curr_n as f64;
            continue
        }

        let rem = n - curr_n;
        let left = min(a, rem / 2);
        let right = min(n - b, rem - left);
        a -= left;
        b += right;

        let sub_x = x.index((a..b, 0));
        let sub_y = y.index((a..b, 0));
        let mut norm_sub_x = sub_x.add_scalar(-xval);
        norm_sub_x.unscale_mut(range); // Divides vector by range.
        let mut weight = tricube_kernel((&norm_sub_x).into());
        if let Some(in_weight) = &w {
            weight.component_mul_assign(in_weight);
        }
    }
    res
    // for i, xval in enumerate(xout):
    //     a = x.searchsorted(xval, 'left')
    //     b = x.searchsorted(xval, 'right')
    //     if b - a >= n_frac:
    //         res[i] = np.mean(y[a:b])
    //         continue

    //     if b - a < n_frac:
    //         rem = n_frac - b + a
    //         left = min(a, (rem // 2))
    //         right = min(n - b, rem - left)
    //         a -= left
    //         b += right
    //         assert a >= 0 and b <= n

    //     sub_x = x[a:b]
    //     sub_y = y[a:b]
    //     weight = common.tricube_kernel((sub_x - xval) / size)
    //     if in_weight is not None:
    //         weight *= in_weight[a:b]
    //     coef = np.polyfit(sub_x, sub_y, deg=deg, w=weight)
    //     res[i] = np.polyval(coef, xval)
    // return res
}

pub fn polyfit(x : &DVector<f64>, y: &DVector<f64>, w: &DVector<f64>, degree: usize) ->
        Result<DVector<f64>, &'static str> {
    let nrow = x.len();
    let ncol = degree + 1;

    let mut a = DMatrix::zeros(nrow, ncol);
    let mut col = DVector::from_element(nrow, 1.0);
    a.set_column(0, &col);
    for i in 1..ncol {
        col *= x;
        a.set_column(i, &col);
    }

    let decomp = SVD::new(a, true, true);
    decomp.solve(&y, 1e-18)
}
