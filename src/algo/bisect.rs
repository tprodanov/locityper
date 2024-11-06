//! Various functions related by binary search.

use std::cmp::Ordering;

/// Performs binary search and finds index `i` such that `a[i-1] < target <= a[i].
#[inline]
pub fn left<T: PartialOrd>(a: &[T], target: &T) -> usize {
    left_by(a, |v| v.partial_cmp(target).expect("Bisect failed: elements are not comparable"))
}

/// Performs binary search between indices `lo` and `hi`
/// and finds index `i` such that `a[i-1] < target <= a[i].
#[inline]
#[allow(unused)]
pub fn left_at<T: PartialOrd>(a: &[T], target: &T, lo: usize, hi: usize) -> usize {
    left_by_at(a, |v| v.partial_cmp(target).expect("Bisect failed: elements are not comparable"), lo, hi)
}

/// Performs binary search and finds index `i` such that `a[i-1] <= target < a[i].
#[inline]
pub fn right<T: PartialOrd>(a: &[T], target: &T) -> usize {
    right_by(a, |v| v.partial_cmp(target).expect("Bisect failed: elements are not comparable"))
}

/// Performs binary search between indices `lo` and `hi`
/// and finds index `i` such that `a[i-1] <= target < a[i].
#[inline]
pub fn right_at<T: PartialOrd>(a: &[T], target: &T, lo: usize, hi: usize) -> usize {
    right_by_at(a, |v| v.partial_cmp(target).expect("Bisect failed: elements are not comparable"), lo, hi)
}

/// Performs binary search
/// and finds the index `i` such that `f(a[i-1]) -> Less` and `f(a[i]) -> Equal | Greater`.
#[inline]
pub fn left_by<T, F: FnMut(&T) -> Ordering>(a: &[T], f: F) -> usize {
    left_by_at(a, f, 0, a.len())
}

/// Performs binary search between indices `lo` and `hi`
/// and finds the index `i` such that `f(a[i-1]) -> Less` and `f(a[i]) -> Equal | Greater`.
pub fn left_by_at<T, F: FnMut(&T) -> Ordering>(a: &[T], mut f: F, mut lo: usize, mut hi: usize) -> usize {
    assert!(lo <= hi && hi <= a.len(), "Cannot perform binary search on indices {}, {} (len: {})", lo, hi, a.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if f(unsafe { a.get_unchecked(mid) }) == Ordering::Less {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    // debug_assert!(lo == 0 || f(&a[lo - 1]) == Ordering::Less,
    //     "Binary search left: Incorrect output index {}: f({}) = {:?}", lo, lo - 1, f(&a[lo - 1]));
    // debug_assert!(lo == a.len() || f(&a[lo]) != Ordering::Less,
    //     "Binary search left: Incorrect output index {}: f({}) = {:?}", lo, lo, f(&a[lo]));
    lo
}

/// Performs binary search
/// and finds the index `i` such that `f(a[i-1]) -> Less | Equal` and `f(a[i]) -> Greater`.
#[inline]
pub fn right_by<T, F: FnMut(&T) -> Ordering>(a: &[T], f: F) -> usize {
    right_by_at(a, f, 0, a.len())
}

/// Performs binary search between indices `lo` and `hi`
/// and finds the index `i` such that `f(a[i-1]) -> Less | Equal` and `f(a[i]) -> Greater`.
pub fn right_by_at<T, F: FnMut(&T) -> Ordering>(a: &[T], mut f: F, mut lo: usize, mut hi: usize) -> usize {
    assert!(lo <= hi && hi <= a.len(), "Cannot perform binary search on indices {}, {} (len: {})", lo, hi, a.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if f(unsafe { a.get_unchecked(mid) }) == Ordering::Greater {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

/// Finds the first index `i` (between `lo` and `hi`) where `! f(a[i])`.
/// Returns `hi` if all `f(a[i])` for all `i`.
///
/// This function takes $O(hi - lo)$, but can still be faster than binary search if expected $i - lo$ is small.
pub fn right_boundary<T, F>(a: &[T], mut f: F, lo: usize, hi: usize) -> usize
where F: FnMut(&T) -> bool,
{
    assert!(hi <= a.len(), "Cannot perform linear search on indices {}, {} (len: {})", lo, hi, a.len());
    for i in lo..hi {
        if !f(unsafe { a.get_unchecked(i) }) {
            return i;
        }
    }
    hi
}
