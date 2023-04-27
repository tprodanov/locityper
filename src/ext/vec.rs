use std::{
    fmt::{self, Write},
    ops::{Add, Sub},
};

/// Static methods extending slices.
pub struct VecExt;

impl VecExt {
    /// Return indices, sorted according to the vector values.
    pub fn argsort<T: PartialOrd>(a: &[T]) -> Vec<usize> {
        let mut ixs: Vec<usize> = (0..a.len()).collect();
        ixs.sort_unstable_by(|&i, &j| a[i].partial_cmp(&a[j]).expect("Error in `argsort`: elements are not comparable"));
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

    /// For each moving window of size `w`, calculates sum over elements in `a`.
    /// Returns `a.len() - w + 1` sums.
    pub fn moving_window_sums<T>(a: &[T], w: usize) -> Vec<T>
    where T: Add<Output = T> + Sub<Output = T> + Copy
    {
        assert!(w >= 1, "moving_window_sums: window size ({}) must be at least 1.", w);
        let n = a.len();
        if n < w {
            return Vec::new();
        }
        let mut sums = Vec::with_capacity(n - w + 1);
        let mut acc = a[..w].iter().copied().reduce(|acc, e| acc + e).unwrap();
        sums.push(acc);

        for i in w..n {
            acc = acc + a[i] - a[i - w];
            sums.push(acc);
        }
        debug_assert_eq!(sums.len(), n - w + 1);
        sums
    }
}

/// Static methods related to vectors over f64.
pub struct F64Ext;

impl F64Ext {
    /// Calculate sample mean.
    pub fn mean<T>(a: &[T]) -> f64
    where T: Into<f64> + Copy,
    {
        a.iter().fold(0.0_f64, |acc, &v| acc + v.into()) / a.len() as f64
    }

    /// Calculate sample variance.
    pub fn variance(a: &[f64], mean: Option<f64>) -> f64 {
        let mean = match mean {
            Some(val) => val,
            None => F64Ext::mean(a),
        };
        let n = a.len();
        assert!(n > 1, "Cannot calculate variance from less than 2 elements!");
        a.iter().fold(0.0, |acc, x| {
            let diff = x - mean;
            acc + diff * diff
        }) / (n - 1) as f64
    }

    /// Returns minimal vector value.
    pub fn min(a: &[f64]) -> f64 {
        a.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Returns maximal vector value.
    pub fn max(a: &[f64]) -> f64 {
        a.iter().copied().fold(f64::NEG_INFINITY, f64::max)
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
    pub fn quantile_sorted(a: &[f64], q: f64) -> f64 {
        assert!(0.0 <= q && q <= 1.0, "Quantile must be within [0, 1]!");
        assert!(a.len() > 0, "Cannot find quantile on an empty array!");
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
    /// Consume iterator and find `q`-th quantile in it. Takes *O(n log n)*.
    pub fn quantile(it: impl Iterator<Item = f64>, q: f64) -> f64 {
        let mut v: Vec<f64> = it.collect();
        VecExt::sort(&mut v);
        F64Ext::quantile_sorted(&v, q)
    }

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
fn count_combinations(n: u64, r: u64) -> u64 {
    if r > n {
        0
    } else {
        (1..=r.min(n - r)).fold(1, |acc, val| acc * (n - val + 1) / val)
    }
}

/// Count the number of combinations with replacement.
fn count_repl_combinations(n: u64, r: u64) -> u64 {
    count_combinations(n + r - 1, r)
}

/// Structure that stores tuples of fixed size.
#[derive(Clone)]
pub struct Tuples<T> {
    /// Linear storage of all tuples.
    data: Vec<T>,
    /// Size of a single tuple.
    tup_len: usize,
}

impl<T> Tuples<T> {
    pub fn new(tup_len: usize) -> Self {
        assert!(tup_len > 0, "Cannot construct tuples with empty length.");
        Self {
            data: Vec::new(),
            tup_len,
        }
    }

    pub fn with_capacity(tup_len: usize, capacity: usize) -> Self {
        assert!(tup_len > 0, "Cannot construct tuples with empty length.");
        Self {
            data: Vec::with_capacity(capacity * tup_len),
            tup_len,
        }
    }

    /// Returns the number of tuples in the set.
    pub fn len(&self) -> usize {
        self.data.len() / self.tup_len
    }

    /// Returns tuple length.
    pub fn tup_len(&self) -> usize {
        self.tup_len
    }

    /// Iterates over all tuples in the set.
    pub fn iter(&self) -> std::slice::ChunksExact<'_, T> {
        self.data.chunks_exact(self.tup_len)
    }
}

impl<T: Copy> Tuples<T> {
    fn init_repl_combinations(&mut self, v: &[T], subtuple: &mut [T], depth: usize, start_ix: usize) {
        if depth + 1 == self.tup_len {
            for &el in &v[start_ix..] {
                self.data.extend_from_slice(subtuple);
                self.data.push(el);
            }
        } else {
            for (i, &el) in v[start_ix..].iter().enumerate() {
                subtuple[depth] = el;
                self.init_repl_combinations(v, subtuple, depth + 1, start_ix + i);
            }
        }
    }

    /// Stores all tuples of size `tup_len`, constructed from vector `v` through all combinations with replacement.
    pub fn repl_combinations(v: &[T], tup_len: usize) -> Self {
        assert!(!v.is_empty(), "Cannot construct combinations on an empty vector or empty tuple size");
        let count = usize::try_from(count_repl_combinations(v.len() as u64, tup_len as u64)).unwrap();
        assert!(count < 10_000_000, "Number of possible tuples ({}) is too great", count);
        let mut res = Tuples::with_capacity(tup_len, count);

        let mut subtuple = vec![v[0]; tup_len - 1];
        if tup_len == 1 {
            res.data.extend_from_slice(v);
        } else {
            res.init_repl_combinations(v, &mut subtuple, 0, 0);
        }
        assert_eq!(res.len(), count);
        res
    }
}

impl<T> std::ops::Index<usize> for Tuples<T> {
    type Output = [T];

    fn index(&self, index: usize) -> &[T] {
        let i = index * self.tup_len;
        let j = i + self.tup_len;
        &self.data[i..j]
    }
}
