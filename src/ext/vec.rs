use std::{
    fmt::{self, Write},
    ops::{AddAssign},
};

/// Static methods extending slices.
pub struct VecExt;

impl VecExt {
    /// Return indices, sorted according to the vector values.
    pub fn argsort<T: PartialOrd>(a: &[T]) -> Vec<usize> {
        let mut ixs: Vec<usize> = (0..a.len()).collect();
        ixs.sort_unstable_by(|&i, &j| a[i].partial_cmp(&a[j])
            .expect("Error in `argsort`: elements are not comparable"));
        ixs
    }

    /// Reorders array by taking elements with indices `ixs` and cloning them. Returns a new vector.
    pub fn reorder<T: Clone>(a: &[T], ixs: &[usize]) -> Vec<T> {
        let mut res = Vec::with_capacity(ixs.len());
        for &i in ixs {
            res.push(a[i].clone());
        }
        res
    }

    /// Sort vector/slice using `partial_ord`.
    pub fn sort<T: PartialOrd>(a: &mut [T]) {
        a.sort_unstable_by(|x, y| x.partial_cmp(y).expect("Error in sort: elements are not comparable"));
    }

    /// Finds `q`-th quantile in a sorted array.
    /// Does not use linear interpolation, just returns closest of the two values.
    fn quantile_sorted<T>(a: &[T], q: f64) -> T
    where T: Copy,
    {
        let n = a.len();
        assert!(0.0 <= q && q <= 1.0, "Quantile must be within [0, 1]!");
        assert!(n > 0, "Cannot find quantile on an empty array!");
        let i = (((a.len() - 1) as f64 * q).round() as usize).min(n - 1);
        a[i]
    }

    pub fn quantile<T>(a: &mut [T], q: f64) -> T
    where T: Copy + PartialOrd,
    {
        Self::sort(a);
        Self::quantile_sorted(a, q)
    }
}

/// Static methods related to vectors over f64.
pub struct F64Ext;

impl F64Ext {
    /// Calculate sample mean.
    pub fn mean<T>(a: &[T]) -> f64
    where T: Into<f64> + Copy,
    {
        a.iter().copied().map(T::into).sum::<f64>() / a.len() as f64
    }

    /// Calculates geometric mean.
    pub fn geom_mean(a: &[f64]) -> f64 {
        (a.iter().copied().map(f64::ln).sum::<f64>() / a.len() as f64).exp()
    }

    /// Returns variance for already calculated mean.
    pub fn fast_variance(a: &[f64], mean: f64) -> f64 {
        let n = a.len();
        assert!(n > 1, "Cannot calculate variance from less than 2 elements!");
        a.iter().fold(0.0, |acc, x| {
            let diff = x - mean;
            acc + diff * diff
        }) / (n - 1) as f64
    }

    /// Calculate sample variance.
    pub fn variance(a: &[f64]) -> f64 {
        Self::fast_variance(a, Self::mean(a))
    }

    /// Returns mean and variance together.
    pub fn mean_variance(a: &[f64]) -> (f64, f64) {
        let mean = Self::mean(a);
        (mean, Self::fast_variance(a, mean))
    }

    /// Returns mean and variance, but, if the number of elements is too low, returns NaN.
    pub fn mean_variance_or_nan(a: &[f64]) -> (f64, f64) {
        let n = a.len();
        match n {
            0 => (f64::NAN, f64::NAN),
            1 => (Self::mean(a), f64::NAN),
            _ => Self::mean_variance(a),
        }
    }

    /// Returns minimal vector value.
    pub fn min(a: &[f64]) -> f64 {
        a.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Returns maximal vector value.
    pub fn max(a: &[f64]) -> f64 {
        IterExt::max(a.iter().copied())
    }

    /// Returns the index of the minimal value and the value itself.
    pub fn argmin(a: &[f64]) -> (usize, f64) {
        IterExt::argmin(a.iter().copied())
    }

    /// Returns the index of the maximal value and the value itself.
    pub fn argmax(a: &[f64]) -> (usize, f64) {
        IterExt::argmax(a.iter().copied())
    }

    /// Finds `q`-th quantile in a sorted array.
    /// Uses linear interpolation, if the quantile is between two elements.
    fn interpol_quantile_sorted(a: &[f64], q: f64) -> f64 {
        assert!(0.0 <= q && q <= 1.0, "Quantile must be within [0, 1]!");
        assert!(!a.is_empty(), "Cannot find quantile on an empty array!");
        let f = (a.len() - 1) as f64 * q;
        let i = f as usize;
        let r = f.fract();
        debug_assert_eq!(i as f64 + r, f);

        let x = a[i];
        if r < 1e-6 {
            x
        } else {
            let y = a[i + 1];
            debug_assert!(y >= x);
            x + (y - x) * r
        }
    }

    /// Finds `q`-th quantile in an unsorted array.
    /// Uses linear interpolation, if the quantile is between two elements.
    pub fn interpol_quantile(a: &mut [f64], q: f64) -> f64 {
        VecExt::sort(a);
        Self::interpol_quantile_sorted(a, q)
    }

    /// Converts vector to string with given precision and width.
    pub fn to_str(a: &[f64], width: usize, precision: usize) -> String {
        let mut buffer = String::new();
        buffer.write_char('[').unwrap();
        for (i, v) in a.iter().enumerate() {
            if i > 0 {
                buffer.write_str(", ").unwrap();
            }
            write!(buffer, "{:width$.precision$}", v).unwrap();
        }
        buffer.write_char(']').unwrap();
        buffer
    }

    /// Converts vector to string with precision = 5 and width = 10.
    #[inline]
    pub fn to_str5(a: &[f64]) -> String {
        F64Ext::to_str(a, 10, 5)
    }
}

/// Static methods related to iterators.
pub struct IterExt;

impl IterExt {
    /// Finds an index of an optimal value `opt`, according `is_better`.
    /// Unless stochastic, `is_better` will return `is_better(opt, e)` for all other `e` in the iterator.
    pub fn arg_optimal<T, I, F>(mut it: I, is_better: F) -> (usize, T)
    where T: Copy,
          I: Iterator<Item = T>,
          F: Fn(T, T) -> bool,
    {
        let mut k = 0;
        let mut opt = it.next().expect("arg_best: empty iterator!");
        for (i, e) in it.enumerate() {
            if is_better(e, opt) {
                k = i + 1;
                opt = e;
            }
        }
        (k, opt)
    }

    /// Returns maximal f64 value across the iterator.
    #[inline]
    pub fn max(it: impl Iterator<Item = f64>) -> f64 {
        it.fold(f64::NEG_INFINITY, f64::max)
    }

    /// Finds an index of the minimal value in the iterator and the value itself.
    /// If the minimal value appears several times, returns the index of the first value.
    /// Panics on an empty iterator.
    pub fn argmin(it: impl Iterator<Item = f64>) -> (usize, f64) {
        Self::arg_optimal(it, |opt, e| opt < e)
    }

    /// Finds an index of the maximal value in the iterator and the value itself.
    /// If the maximal value appears several times, returns the index of the first value.
    /// Panics on an empty iterator.
    pub fn argmax(it: impl Iterator<Item = f64>) -> (usize, f64) {
        Self::arg_optimal(it, |opt, e| opt > e)
    }

    /// Calculates cumulative sums over iterator and returns the vector of length `len(iter) + 1`.
    /// First element = 0.
    pub fn cumul_sums<T, U>(it: impl Iterator<Item = T>) -> Vec<U>
    where U: From<T> + num_traits::Zero + AddAssign + Copy,
    {
        let mut c = U::zero();
        let mut res = vec![c];
        for el in it {
            c += U::from(el);
            res.push(c);
        }
        res
    }
}

/// Trait that encodes either `&mut Vec<T>` or `()`, depending on if the information needs to be saved or not.
pub trait VecOrNone<T> : fmt::Debug {
    /// True for `()` - meaning that the object does not actually do any work.
    const IS_SINK: bool;

    /// Push a new element.
    fn push(&mut self, val: T);

    /// Returns length, if available.
    fn try_len(&self) -> Option<usize>;
}

impl<T> VecOrNone<T> for () {
    const IS_SINK: bool = true;

    fn push(&mut self, _val: T) {}

    fn try_len(&self) -> Option<usize> { None }
}

impl<T> VecOrNone<T> for &mut Vec<T>
where T: fmt::Debug
{
    const IS_SINK: bool = false;

    fn push(&mut self, val: T) {
        Vec::push(*self, val);
    }

    fn try_len(&self) -> Option<usize> {
        Some(self.len())
    }
}

/// Count the number of combinations.
/// Taken from https://stackoverflow.com/questions/65561566/number-of-combinations-permutations.
pub fn count_combinations(n: usize, r: usize) -> usize {
    if r > n {
        0
    } else {
        (1..=r.min(n - r)).fold(1, |acc, val| acc * (n - val + 1) / val)
    }
}

/// Count the number of combinations with replacement.
pub fn count_combinations_with_repl(n: usize, r: usize) -> usize {
    count_combinations(n + r - 1, r)
}

fn recursive_combinations_with_repl<T, F>(
    v: &[T],
    buffer: &mut [T],
    start: usize,
    depth: usize,
    size: usize,
    action: &mut F,
)
where T: Copy,
      F: FnMut(&[T]),
{
    if depth + 1 == size {
        for &el in &v[start..] {
            buffer[depth] = el;
            action(buffer);
        }
    } else {
        for (i, &el) in v[start..].iter().enumerate() {
            buffer[depth] = el;
            recursive_combinations_with_repl(v, buffer, start + i, depth + 1, size, action);
        }
    }
}

/// Generates combinations with replacement, and calls `action` for each generated slice.
pub fn gen_combinations_with_repl<T, F>(v: &[T], size: usize, mut action: F)
where T: Copy,
      F: FnMut(&[T]),
{
    if v.is_empty() || size == 0 {
        // Do nothing.
    } else if size == 1 {
        let mut buffer = [v[0]];
        for &el in v {
            buffer[0] = el;
            action(&buffer);
        }
    } else {
        let mut buffer = vec![v[0]; size];
        recursive_combinations_with_repl(v, &mut buffer, 0, 0, size, &mut action);
    }
}

/// Use Heap's algorithm to generate all permutations. For each permutation, `action` is called.
pub fn gen_permutations<T, F>(v: &[T], mut action: F)
where T: Copy,
      F: FnMut(&[T]),
{
    let n = v.len();
    if n == 0 {
        // Do nothing.
    } else if n == 1 {
        action(v);
    } else if n == 2 {
        action(v);
        action(unsafe { &[*v.get_unchecked(1), *v.get_unchecked(0)] });
    } else {
        let mut buffer = v.to_vec();
        let mut c = vec![0; n];
        let mut i = 1;
        while i < n {
            let ci = unsafe { c.get_unchecked_mut(i) };
            if *ci < i {
                // 0 if i is even, *ci if i is odd.
                buffer.swap(i, *ci * (i % 2));
                action(&buffer);
                *ci += 1;
                i = 1;
            } else {
                *ci = 0;
                i += 1;
            }
        }
    }
}
