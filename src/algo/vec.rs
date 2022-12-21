//! Simple vector algorithms and traits.

use std::cmp::Ordering;
use std::ops::{Add, AddAssign, Mul, MulAssign};

use itertools::izip;

/// Reorders `v` by taking elements with indices `ixs`. Returns `Vec`.
pub fn reorder<T>(v: &[T], ixs: &[usize]) -> Vec<T>
where T: Copy,
{
    let mut res = Vec::with_capacity(ixs.len());
    for &i in ixs {
        res.push(v[i]);
    }
    res
}

/// Performs binary search and finds the index `i` such that `v[i-1] < target <= v[i]`.
#[inline]
pub fn binary_search_left<T: PartialOrd>(v: &[T], target: T) -> usize {
    binary_search_left_at(v, target, 0, v.len())
}

/// Performs binary search between indices `low` and `high`
/// and finds the index `i` such that `v[i-1] < target <= v[i]`.
pub fn binary_search_left_at<T: PartialOrd>(v: &[T], target: T, mut low: usize, mut high: usize) -> usize {
    while low < high {
        let mid = (low + high) / 2;
        match v[mid].partial_cmp(&target).expect("Incomparable values in binary_search_left") {
            Ordering::Equal if mid == 0 || v[mid - 1] < target => return mid,
            Ordering::Less => low = mid + 1,
            _ => high = mid,
        }
    }
    low
}

/// Performs binary search and finds the index `i` such that `v[i-1] <= target < v[i]`.
#[inline]
pub fn binary_search_right<T: PartialOrd>(v: &[T], target: T) -> usize {
    binary_search_right_at(v, target, 0, v.len())
}

/// Performs binary search between indices `low` and `high`
/// and finds the index `i` such that `v[i-1] <= target < v[i]`.
pub fn binary_search_right_at<T: PartialOrd>(v: &[T], target: T, mut low: usize, mut high: usize) -> usize {
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

pub trait F64Vec {
    fn min(&self) -> f64;
    fn max(&self) -> f64;
    fn argsort(&self) -> Vec<usize>;
}

pub trait F64MutVec {
    fn sort(&mut self);
    fn sort_rev(&mut self);
}

impl F64Vec for &[f64] {
    fn min(&self) -> f64 {
        self.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    fn max(&self) -> f64 {
        self.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }

    fn sort(&mut self) {
        self.sort_by(|a, b| a.total_cmp(&b));
    }

    fn sort_rev(&mut self) {
        self.sort_by(|a, b| b.total_cmp(&a));
    }

    fn argsort(&self) -> Vec<usize> {
        let mut ixs: Vec<usize> = (0..self.len()).collect();
        ixs.sort_by(|&i, &j| self[i].total_cmp(&self[j]));
        ixs
    }
}

/// Element-wise expansion over vectors.
/// All operators modify the vector itself and return a mutable reference to itself.
pub trait ElementWise<T, Rhs = T> {
    /// Adds a scalar or another vector to self.
    fn add(&mut self, oth: Rhs) -> &mut Self;

    /// Multiplies the vector by another vector or a scalar.
    fn mul(&mut self, oth: Rhs) -> &mut Self;

    /// First, adds a vector (or scalar), and then multiplies by another vector (or scalar).
    fn add_mul(&mut self, add: Rhs, mul: Rhs) -> &mut Self;
}

impl<T> ElementWise<T, T> for Vec<T>
where T: Copy + Add<Output = T> + AddAssign + Mul<Output = T> + MulAssign
{
    fn add(&mut self, scalar: T) -> &mut Self {
        for v in self.iter_mut() {
            *v += scalar;
        }
        self
    }

    fn mul(&mut self, scalar: T) -> &mut Self {
        for v in self.iter_mut() {
            *v += scalar;
        }
        self
    }

    fn add_mul(&mut self, add: T, mul: T) -> &mut Self {
        for v in self.iter_mut() {
            *v = (*v + add) * mul;
        }
        self
    }
}

impl<T> ElementWise<T, &[T]> for Vec<T>
where T: Copy + Add<Output = T> + AddAssign + Mul<Output = T> + MulAssign
{
    fn add(&mut self, oth: &[T]) -> &mut Self {
        assert_eq!(self.len(), oth.len(), "ElementWise::add - vector sizes are different!");
        for (a, &b) in self.iter_mut().zip(oth) {
            *a += b;
        }
        self
    }

    fn mul(&mut self, oth: &[T]) -> &mut Self {
        assert_eq!(self.len(), oth.len(), "ElementWise::mul - vector sizes are different!");
        for (a, &b) in self.iter_mut().zip(oth) {
            *a *= b;
        }
        self
    }

    fn add_mul(&mut self, add: &[T], mul: &[T]) -> &mut Self {
        assert_eq!(self.len(), add.len(), "ElementWise::add_mul - unexpected `add` vector size!");
        assert_eq!(self.len(), mul.len(), "ElementWise::add_mul - unexpected `mul` vector size!");
        for (v, &a, &m) in izip!(self.iter_mut(), add, mul) {
            *v = (*v + a) * m;
        }
        self
    }
}

impl<T> ElementWise<T, &Vec<T>> for Vec<T>
where T: Copy + Add<Output = T> + AddAssign + Mul<Output = T> + MulAssign
{
    #[inline]
    fn add(&mut self, oth: &Vec<T>) -> &mut Self {
        self.add(oth as &[T])
    }

    #[inline]
    fn mul(&mut self, oth: &Vec<T>) -> &mut Self {
        self.mul(oth as &[T])
    }

    #[inline]
    fn add_mul(&mut self, add: &Vec<T>, mul: &Vec<T>) -> &mut Self {
        self.add_mul(add as &[T], mul as &[T])
    }
}