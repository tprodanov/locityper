use std::ops::{Index, Add, AddAssign, Mul, MulAssign, Sub, SubAssign, Div, DivAssign};
use std::fmt::Write;
use itertools::{izip, Itertools};

/// General vector extension.
pub trait VecExt<T> : Index<usize, Output = T>
where T: Clone + PartialOrd
{
    /// Returns vector length.
    fn len(&self) -> usize;

    /// Return indices, sorted according to the vector values.
    fn argsort(&self) -> Vec<usize> {
        let mut ixs: Vec<usize> = (0..self.len()).collect();
        ixs.sort_by(|&i, &j| self[i].partial_cmp(&self[j]).expect("Error in `argsort`: elements re not comparable"));
        ixs
    }

    /// Reorders `self` by taking elements with indices `ixs`. Returns a new vector.
    fn reorder(&self, ixs: &[usize]) -> Vec<T> {
        let mut res = Vec::with_capacity(ixs.len());
        for &i in ixs {
            res.push(self[i].clone());
        }
        res
    }
}

impl<T: Clone + PartialOrd> VecExt<T> for Vec<T> {
    #[inline]
    fn len(&self) -> usize {
        self.len()
    }
}

impl<T: Clone + PartialOrd> VecExt<T> for [T] {
    #[inline]
    fn len(&self) -> usize {
        self.len()
    }
}

/// `f64` vector extension.
pub trait F64VecExt : VecExt<f64> {
    /// Sort f64 vector/slice.
    fn sort(&mut self);

    /// Calculate sample mean.
    fn mean(&self) -> f64;

    /// Calculate sample variance.
    fn variance(&self, mean: Option<f64>) -> f64;

    /// Returns minimal vector value.
    fn min(&self) -> f64;

    /// Returns maximal vector value.
    fn max(&self) -> f64;

    /// Finds `q`-th quantile in a sorted array.
    /// Uses linear interpolation, if the quantile is between two elements.
    fn quantile_sorted(&self, q: f64) -> f64;

    /// Converts vector to string with given precision and width.
    fn to_str(&self, width: usize, precision: usize) -> String;

    /// Converts vector to string with precision = 5 and width = 10.
    #[inline]
    fn to_str5(&self) -> String {
        self.to_str(10, 5)
    }

    /// Element-wise applies function to every element and returns a new vector.
    fn ew_apply1<F: FnMut(f64) -> f64>(&self, f: F) -> Vec<f64>;

    /// Element-wise applies function to `self` and `oth`, and returns a new vector.
    fn ew_apply2<F: FnMut(f64, f64) -> f64>(&self, oth: &[f64], f: F) -> Vec<f64>;

    /// Element-wise applies function to `self`, `a` and `b`, and returns a new vector.
    fn ew_apply3<F: FnMut(f64, f64, f64) -> f64>(&self, a: &[f64], b: &[f64], f: F) -> Vec<f64>;

    /// Element-wise modifies the vector by applying function `f` to every element.
    fn ew_assign1<F: FnMut(&mut f64)>(&mut self, f: F);

    /// Element-wise modifies the vector by applying function `f` to `self` and `oth`.
    fn ew_assign2<F: FnMut(&mut f64, f64)>(&mut self, oth: &[f64], f: F);

    /// Element-wise modifies the vector by applying function `f` to `self`, `a` and `b`.
    fn ew_assign3<F: FnMut(&mut f64, f64, f64)>(&mut self, a: &[f64], b: &[f64], f: F);

    /// Adds two vectors element-wise and returns a new vector.
    #[inline]
    fn ew_add(&self, oth: &[f64]) -> Vec<f64> {
        self.ew_apply2(oth, f64::add)
    }

    /// Element-wise adds `oth` to `self`, updates `self`.
    #[inline]
    fn ew_add_assign(&mut self, oth: &[f64]) {
        self.ew_assign2(oth, f64::add_assign)
    }

    /// Multiplies two vectors element-wise and returns a new vector.
    #[inline]
    fn ew_mul(&self, oth: &[f64]) -> Vec<f64> {
        self.ew_apply2(oth, f64::mul)
    }

    /// Element-wise multiplies `oth` by `self`, updates `self`.
    #[inline]
    fn ew_mul_assign(&mut self, oth: &[f64]) {
        self.ew_assign2(oth, f64::mul_assign)
    }

    #[inline]
    fn ew_sub(&self, oth: &[f64]) -> Vec<f64> {
        self.ew_apply2(oth, f64::sub)
    }

    #[inline]
    fn ew_sub_assign(&mut self, oth: &[f64]) {
        self.ew_assign2(oth, f64::sub_assign)
    }

    #[inline]
    fn ew_div(&self, oth: &[f64]) -> Vec<f64> {
        self.ew_apply2(oth, f64::div)
    }

    #[inline]
    fn ew_div_assign(&mut self, oth: &[f64]) {
        self.ew_assign2(oth, f64::div_assign)
    }

    #[inline]
    fn ew_add_mul(&self, a: &[f64], b: &[f64]) -> Vec<f64> {
        self.ew_apply3(a, b, |x, y, z| (x + y) * z)
    }

    #[inline]
    fn ew_add_mul_assign(&mut self, a: &[f64], b: &[f64]) {
        self.ew_assign3(a, b, |x, y, z| *x = (*x + y) * z)
    }

    #[inline]
    fn sc_add(&self, scalar: f64) -> Vec<f64> {
        self.ew_apply1(|x| x + scalar)
    }

    #[inline]
    fn sc_add_assign(&mut self, scalar: f64) {
        self.ew_assign1(|x| *x += scalar)
    }

    #[inline]
    fn sc_sub(&self, scalar: f64) -> Vec<f64> {
        self.ew_apply1(|x| x - scalar)
    }

    #[inline]
    fn sc_sub_assign(&mut self, scalar: f64) {
        self.ew_assign1(|x| *x -= scalar)
    }

    #[inline]
    fn sc_mul(&self, scalar: f64) -> Vec<f64> {
        self.ew_apply1(|x| x * scalar)
    }

    #[inline]
    fn sc_mul_assign(&mut self, scalar: f64) {
        self.ew_assign1(|x| *x *= scalar)
    }

    #[inline]
    fn sc_div(&self, scalar: f64) -> Vec<f64> {
        self.sc_mul(1.0 / scalar)
    }

    #[inline]
    fn sc_div_assign(&mut self, scalar: f64) {
        self.sc_mul_assign(1.0 / scalar)
    }

    #[inline]
    fn sc_add_mul(&self, a: f64, b: f64) -> Vec<f64> {
        self.ew_apply1(|x| (x + a) * b)
    }

    #[inline]
    fn sc_add_mul_assign(&mut self, a: f64, b: f64) {
        self.ew_assign1(|x| *x = (*x + a) * b)
    }
}

macro_rules! f64_vec_ext_impl {
    (
        $coll:ty
    ) => {
        impl F64VecExt for $coll {
            fn sort(&mut self) {
                self.sort_by(|a, b| a.total_cmp(b));
            }

            fn mean(&self) -> f64 {
                self.iter().sum::<f64>() / self.len() as f64
            }

            fn variance(&self, mean: Option<f64>) -> f64 {
                let mean = match mean {
                    Some(val) => val,
                    None => self.mean(),
                };
                let n = self.len();
                assert!(n > 1, "Cannot calculate variance from less than 2 elements!");
                self.iter().fold(0.0, |acc, x| {
                    let diff = x - mean;
                    acc + diff * diff
                }) / (n - 1) as f64
            }

            #[inline]
            fn min(&self) -> f64 {
                self.iter().cloned().fold(f64::INFINITY, f64::min)
            }

            #[inline]
            fn max(&self) -> f64 {
                self.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            }

            fn quantile_sorted(&self, q: f64) -> f64 {
                assert!(0.0 <= q && q <= 1.0, "Quantile must be within [0, 1]!");
                assert!(self.len() > 0, "Cannot find quantile on an empty array!");
                let f = (self.len() - 1) as f64 * q;
                let i = f as usize;
                let r = f.fract();
                debug_assert_eq!(i as f64 + r, f);

                let x = self[i];
                if r < 1e-6 {
                    x
                } else {
                    let y = self[i + 1];
                    debug_assert!(y >= x);
                    x + (y - x) * r
                }
            }

            fn to_str(&self, width: usize, precision: usize) -> String {
                let mut buffer = String::new();
                buffer.write_char('[').unwrap();
                for (i, v) in self.iter().enumerate() {
                    if i > 0 {
                        buffer.write_str(", ").unwrap();
                    }
                    write!(buffer, "{:width$.precision$}", v).unwrap();
                }
                buffer.write_char(']').unwrap();
                buffer
            }

            #[inline]
            fn ew_apply1<F: FnMut(f64) -> f64>(&self, f: F) -> Vec<f64> {
                self.iter().cloned().map(f).collect()
            }

            fn ew_apply2<F: FnMut(f64, f64) -> f64>(&self, oth: &[f64], mut f: F) -> Vec<f64> {
                let n = self.len();
                assert_eq!(n, oth.len(), "Element-wise operation is not permitted: different array lengths!");
                let mut res = Vec::with_capacity(n);
                for (&a, &b) in self.iter().zip(oth) {
                    res.push(f(a, b));
                }
                res
            }

            fn ew_apply3<F: FnMut(f64, f64, f64) -> f64>(&self, a: &[f64], b: &[f64], mut f: F) -> Vec<f64> {
                let n = self.len();
                assert_eq!(n, a.len(), "Element-wise operation is not permitted: different array lengths!");
                assert_eq!(n, b.len(), "Element-wise operation is not permitted: different array lengths!");
                let mut res = Vec::with_capacity(n);
                for (&v, &a, &b) in izip!(self.iter(), a, b) {
                    res.push(f(v, a, b));
                }
                res
            }

            fn ew_assign1<F: FnMut(&mut f64)>(&mut self, mut f: F) {
                for v in self.iter_mut() {
                    f(v)
                }
            }

            fn ew_assign2<F: FnMut(&mut f64, f64)>(&mut self, oth: &[f64], mut f: F) {
                assert_eq!(self.len(), oth.len(),
                    "Element-wise operation is not permitted: different array lengths!");
                for (a, &b) in self.iter_mut().zip(oth) {
                    f(a, b);
                }
            }

            fn ew_assign3<F: FnMut(&mut f64, f64, f64)>(&mut self, a: &[f64], b: &[f64], mut f: F) {
                let n = self.len();
                assert_eq!(n, a.len(), "Element-wise operation is not permitted: different array lengths!");
                assert_eq!(n, b.len(), "Element-wise operation is not permitted: different array lengths!");
                for (v, &a, &b) in izip!(self.iter_mut(), a, b) {
                    f(v, a, b)
                }
            }
        }
    }
}

f64_vec_ext_impl!{Vec<f64>}
f64_vec_ext_impl!{[f64]}

/// Static methods related to iterators.
pub struct IterExt;

impl IterExt {
    /// Sort iterator over f64.
    pub fn sorted(it: impl Iterator<Item = f64>) -> Vec<f64> {
        it.sorted_by(f64::total_cmp).collect()
    }

    /// Consume iterator and find `q`-th quantile in it. Takes *O(n log n)*.
    pub fn quantile(it: impl Iterator<Item = f64>, q: f64) -> f64 {
        let v = IterExt::sorted(it);
        v.quantile_sorted(q)
    }
}
