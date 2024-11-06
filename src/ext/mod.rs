//! Extending standard functions and structures.

pub mod vec;
pub mod fmt;
pub mod sys;
pub mod rand;

use std::{
    cmp::Ordering,
    ops::{Index, IndexMut},
};

/// Upper triangle matrix, excluding the diagonal (i < j).
pub struct TriangleMatrix<T> {
    side: usize,
    data: Vec<T>,
}

impl TriangleMatrix<()> {
    /// Returns iterator over all pairs `(i, j)` such that `0 <= i < j < side`.
    pub fn indices(side: usize) -> impl Iterator<Item = (usize, usize)> {
        (0..side - 1).flat_map(move |i| (i + 1..side).map(move |j| (i, j)))
    }

    /// Returns iterator over all pairs `(i, j)` such that `0 <= i < j < side`.
    pub fn indices_u32(side: usize) -> impl Iterator<Item = (u32, u32)> {
        let side = side as u32;
        (0..side - 1).flat_map(move |i| (i + 1..side).map(move |j| (i, j)))
    }
}

impl<T> TriangleMatrix<T> {
    /// Returns 0-based index (i,j) in the regular matrix n×n based on linear index k.
    #[allow(unused)]
    pub fn from_linear_index(&self, k: usize) -> (usize, usize) {
        assert!(k < self.data.len());
        let under_root = (8 * self.data.len()).checked_sub(8 * k + 7).unwrap();
        let i = self.side.checked_sub(2 + (0.5 * (under_root as f64).sqrt() - 0.5).floor() as usize).unwrap();
        let j = (k + i * (i + 3) / 2 + 1).checked_sub(self.side * i).unwrap();
        (i, j)
    }

    #[inline]
    pub fn expected_len(side: usize) -> usize {
        side.saturating_sub(1) * side / 2
    }

    #[inline]
    fn to_linear_index(&self, i: usize, j: usize) -> usize {
        assert!(i < j && j < self.side, "Incorrect indices ({}, {}) to triangle matrix", i, j);
        (2 * self.side - 3 - i) * i / 2 + j - 1
    }

    /// Creates triangle matrix from linear storage (must have correct order: sorted first by row, then by column).
    #[allow(unused)]
    pub fn from_linear(side: usize, data: Vec<T>) -> Self {
        assert_eq!(data.len(), Self::expected_len(side), "Incorrect triangle matrix size");
        Self { side, data }
    }

    /// Creates triangle matrix by running `f(i, j)` for corresponding indices.
    #[allow(unused)]
    pub fn create(side: usize, f: impl FnMut((usize, usize)) -> T) -> Self {
        Self {
            side,
            data: TriangleMatrix::indices(side).map(f).collect(),
        }
    }

    /// Total number of elements in the triangle matrix.
    #[allow(unused)]
    pub fn linear_len(&self) -> usize {
        self.data.len()
    }

    #[allow(unused)]
    pub fn linear_data(&self) -> &[T] {
        &self.data
    }

    #[allow(unused)]
    pub fn take_linear(self) -> Vec<T> {
        self.data
    }

    /// Size of the matrix side.
    pub fn side(&self) -> usize {
        self.side
    }

    /// Linear iterator over elements (first by row, second by column).
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.data.iter()
    }

    /// Returns `self[(i, j)]` if `i < j`, `self[(j, i)] if `j < i`, and `None` when `i == j`.
    pub fn get_symmetric(&self, i: usize, j: usize) -> Option<&T> {
        match i.cmp(&j) {
            Ordering::Less => Some(self.index((i, j))),
            Ordering::Equal => None,
            Ordering::Greater => Some(self.index((j, i))),
        }
    }
}

impl<T: Clone> TriangleMatrix<T> {
    #[allow(unused)]
    pub fn new(side: usize, val: T) -> Self {
        Self {
            side,
            data: vec![val; Self::expected_len(side)],
        }
    }
}

impl<T> Index<(usize, usize)> for TriangleMatrix<T> {
    type Output = T;

    #[inline]
    fn index(&self, (i, j): (usize, usize)) -> &T {
        self.data.index(self.to_linear_index(i, j))
    }
}

impl<T> IndexMut<(usize, usize)> for TriangleMatrix<T> {
    #[inline]
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut T {
        self.data.index_mut(self.to_linear_index(i, j))
    }
}
