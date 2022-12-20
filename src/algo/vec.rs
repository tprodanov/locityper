use std::cmp::Ordering;
use std::fmt::{Display, Debug};
use std::ops::{Index, IndexMut};

use num_traits::Zero;
use nalgebra::DVector;

/// Reorders `storage` by taking elements in the order of `ixs`.
/// Returns `Vec`.
pub fn reorder_vec<T, V, I>(storage: V, ixs: I) -> Vec<T>
where T: Copy,
      V: Index<usize, Output = T>,
      I: Iterator<Item = usize> + ExactSizeIterator
{
    let n = ixs.len();
    let mut res = Vec::with_capacity(n);
    for i in ixs {
        res.push(storage[i]);
    }
    res
}

/// Reorders `storage` by taking elements in the order of `ixs`.
/// Returns `DVector`.
pub fn reorder_dvec<T, V, I>(storage: V, ixs: I) -> DVector<T>
where T: Copy + Debug + Eq + Zero + 'static,
      V: Index<usize, Output = T>,
      I: Iterator<Item = usize> + ExactSizeIterator
{
    let n = ixs.len();
    let mut res = DVector::zeros(n);
    for (i, j) in ixs.enumerate() {
        res[i] = storage[j]
    }
    res
}

pub trait Scalar : nalgebra::Scalar + Zero + Copy { }

pub trait GeneralVec<T: Scalar> : Index<usize, Output = T> + IndexMut<usize, Output = T> {
    fn create_vec(n: usize) -> Self;
}

impl<T: Scalar> GeneralVec<T> for Vec<T> {
    fn create_vec(n: usize) -> Self {
        vec![T::zero(); n]
    }
}

impl<T: Scalar> GeneralVec<T> for DVector<T> {
    fn create_vec(n: usize) -> Self {
        DVector::zeros(n)
    }
}

/// Reorders `storage` by taking elements in the order of `ixs`.
/// Returns `DVector`.
pub fn reorder<T, V, I, W>(storage: V, ixs: I) -> W
where T: Scalar,
      V: Index<usize, Output = T>,
      I: Iterator<Item = usize> + ExactSizeIterator,
      W: GeneralVec<T>
{
    let n = ixs.len();
    let mut res = W::create_vec(n);
    for (i, j) in ixs.enumerate() {
        res[i] = storage[j]
    }
    res
}

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

// /// Performs binary search and finds the index `i` such that `v[i-1] <= target < v[i]`.
// pub fn binary_search_right<T: PartialOrd>(v: &[T], target: T) -> usize {
//     let n = v.len();
//     let mut low = 0;
//     let mut high = n;
//     while low < high {
//         let mid = (low + high) / 2;
//         if v[mid] > target {
//             if mid == 0 || v[mid - 1] == target {
//                 return mid;
//             }
//             high = mid;
//         } else {
//             low = mid + 1;
//         }
//     }
//     low
// }