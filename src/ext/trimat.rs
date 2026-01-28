use std::{
    cmp::Ordering,
    ops::{Index, IndexMut},
};

/// Upper triangle matrix (i < j), diagonal is excluded.
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

    #[inline]
    pub fn calc_linear_len(side: usize) -> usize {
        side.saturating_sub(1) * side / 2
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
    fn to_linear_index(&self, i: usize, j: usize) -> usize {
        assert!(i < j && j < self.side, "Incorrect indices ({}, {}) to triangle matrix", i, j);
        (2 * self.side - 3 - i) * i / 2 + j - 1
    }

    /// Creates triangle matrix from linear storage (must have correct order: sorted first by row, then by column).
    pub fn from_linear_data(side: usize, data: Vec<T>) -> Self {
        assert_eq!(data.len(), TriangleMatrix::calc_linear_len(side), "Incorrect triangle matrix size");
        Self { side, data }
    }

    /// Creates triangle matrix by running `f(i, j)` for corresponding indices.
    #[allow(unused)]
    pub fn new_with(side: usize, f: impl FnMut((usize, usize)) -> T) -> Self {
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

    /// Linear iterator over elements (first by row, second by column).
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, T> {
        self.data.iter_mut()
    }

    /// Returns `self[(i, j)]` if `i < j`, `self[(j, i)] if `j < i`, and `None` when `i == j`.
    pub fn get_symmetric(&self, i: usize, j: usize) -> Option<&T> {
        match i.cmp(&j) {
            Ordering::Less => Some(self.index((i, j))),
            Ordering::Equal => None,
            Ordering::Greater => Some(self.index((j, i))),
        }
    }

    /// Returns `self[(i, j)]` if `i < j`, `self[(j, i)] if `j < i`, and `None` when `i == j`.
    pub fn get_mut_symmetric(&mut self, i: usize, j: usize) -> Option<&mut T> {
        match i.cmp(&j) {
            Ordering::Less => Some(self.index_mut((i, j))),
            Ordering::Equal => None,
            Ordering::Greater => Some(self.index_mut((j, i))),
        }
    }
}

impl<T: Clone> TriangleMatrix<T> {
    pub fn new(side: usize, val: T) -> Self {
        Self {
            side,
            data: vec![val; TriangleMatrix::calc_linear_len(side)],
        }
    }

    /// Returns submatrix only with the given indices, in the given order.
    pub fn thin_out<G>(&self, ixs: &[G]) -> TriangleMatrix<T>
    where G: TryInto<usize> + Copy,
          G::Error: std::fmt::Debug,
    {
        Self {
            side: ixs.len(),
            data: TriangleMatrix::indices(ixs.len())
                .map(|(i, j)| {
                    let k = ixs[i].try_into().expect("Matrix index out of usize range");
                    let l = ixs[j].try_into().expect("Matrix index out of usize range");
                    self.get_symmetric(k, l).unwrap().clone()
                }).collect(),
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
