use std::cmp::Ordering;
use std::fmt::{Display, Debug};
use std::ops::Index;

use num_traits::Zero;
use nalgebra::DVector;

/// Reorders `storage` by taking elements in the order of `ixs`.
/// Returns reordered iterator.
pub fn reorder<T, V, I>(storage: V, ixs: I) -> Reorder<T, V, I>
where T: Copy,
      V: Index<usize, Output = T>,
      I: Iterator<Item = usize>
{
    Reorder { storage, ixs }
}

pub struct Reorder<T, V, I>
where T: Copy,
      V: Index<usize, Output = T>,
      I: Iterator<Item = usize>
{
    storage: V,
    ixs: I,
}

impl<T, V, I> Iterator for Reorder<T, V, I>
where T: Copy,
      V: Index<usize, Output = T>,
      I: Iterator<Item = usize>
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        match self.ixs.next() {
            Some(ix) => Some(*self.storage.index(ix)),
            None => None,
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.ixs.size_hint()
    }
}

impl<T, V, I> ExactSizeIterator for Reorder<T, V, I>
where T: Copy,
      V: Index<usize, Output = T>,
      I: Iterator<Item = usize> + ExactSizeIterator
{
    #[inline]
    fn len(&self) -> usize {
        self.ixs.len()
    }
}

impl<T, V, I> std::iter::FusedIterator for Reorder<T, V, I>
where T: Copy,
      V: Index<usize, Output = T>,
      I: Iterator<Item = usize> + std::iter::FusedIterator
{ }

/// Performs binary search and finds the index `i` such that `v[i-1] < target <= v[i]`.
pub fn binary_search_left<T: PartialOrd + Display>(v: &[T], target: T) -> usize {
    let mut low = 0;
    let mut high = v.len();
    while low < high {
        let mid = (low + high) / 2;
        match v[mid].partial_cmp(&target) {
            None => panic!("Error in binary_search: values v[{}] = {} and target = {} cannot be compared",
                mid, v[mid], target),
            Some(Ordering::Equal) if mid == 0 || v[mid - 1] < target => return mid,
            Some(Ordering::Less) => low = mid + 1,
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