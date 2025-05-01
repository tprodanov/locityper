use std::{
    fmt::{self, Display, Debug},
    cmp::{PartialOrd, Ordering},
};
use num_traits::{
    ConstZero, ConstOne, CheckedAdd, CheckedMul,
    cast::{AsPrimitive, FromPrimitive},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Fraction<T> {
    numerator: T,
    denominator: T,
}

impl<T> Fraction<T> {
    #[inline]
    pub fn new(numerator: T, denominator: T) -> Self {
        Self { numerator, denominator }
    }
}

impl<T: Copy> Fraction<T> {
    #[inline]
    pub fn numerator(self) -> T {
        self.numerator
    }

    #[inline]
    pub fn denominator(self) -> T {
        self.denominator
    }
}

impl<T: AsPrimitive<f64>> Fraction<T> {
    #[allow(dead_code)]
    pub fn eval(self) -> f64 {
        self.numerator.as_() / self.denominator.as_()
    }
}

impl<T> Fraction<T>
where T: Copy + PartialOrd + FromPrimitive + AsPrimitive<f64>
          + ConstZero + ConstOne + CheckedAdd + CheckedMul,
{
    /// Approximate x with a rational fraction where both numerator and denominator fit into `T`.
    /// Uses continued fraction algorithm, generally very fast.
    pub fn approximate(x: f64) -> Self {
        assert!(x >= 0.0);

        // a_{k - 2}, a_{k - 1}, b_{k - 2} and b_{k - 1}.
        let mut a2 = T::ONE;
        let mut a1 = T::from_f64(x.floor()).unwrap();
        let mut b2 = T::ZERO;
        let mut b1 = T::ONE;

        let mut xk = x;
        const ITERATIONS: u32 = 20;
        for _ in 0..ITERATIONS {
            let numer = xk - xk.floor();
            if numer <= f64::EPSILON { break; }
            xk = 1.0 / numer;
            let Some(floor) = T::from_f64(xk.floor()) else { break };
            let Some(a0) = floor.checked_mul(&a1).and_then(|a| a.checked_add(&a2)) else { break };
            let Some(b0) = floor.checked_mul(&b1).and_then(|b| b.checked_add(&b2)) else { break };
            a2 = a1;
            a1 = a0;
            b2 = b1;
            b1 = b0;
            if (a1.as_() / b1.as_() - x).abs() <= f64::EPSILON { break };
        }
        Fraction::new(a1, b1)
    }
}

// impl<T> PartialOrd for Fraction<T>
// where T: PartialOrd + num_traits::CheckedMul + num_traits::sign::Unsigned,
// {
//     fn partial_cmp(&self, oth: &Self) -> Ordering {
//         self.numerator.checked_mul(oth.denominator).unwrap()
//             .partial_cmp(&oth.numerator.checked_mul(self.denominator))
//     }
// }

impl PartialOrd for Fraction<u16> {
    fn partial_cmp(&self, oth: &Self) -> Option<Ordering> {
        // Go to higher size to have no overflow
        (self.numerator as u32 * oth.denominator as u32)
            .partial_cmp(&(oth.numerator as u32 * self.denominator as u32))
    }
}

impl PartialOrd for Fraction<u32> {
    fn partial_cmp(&self, oth: &Self) -> Option<Ordering> {
        // Go to higher size to have no overflow
        (self.numerator as u64 * oth.denominator as u64)
            .partial_cmp(&(oth.numerator as u64 * self.denominator as u64))
    }
}

impl<T: Display> Display for Fraction<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} / {}", self.numerator, self.denominator)
    }
}
