use std::{
    cmp::Ord,
    ops::{Shl, Shr, BitOr, BitAnd},
};
use ruint::aliases::U256;
use crate::algo::IntSet;

/// Store k-mers in integers of different length (15-mers in u32, 31-mers in u64, 63-mers in u128).
pub trait Kmer: Copy + Ord + Eq
    + std::fmt::Display + std::fmt::Binary + std::fmt::LowerHex + std::fmt::UpperHex
    + Shl<u8, Output = Self>
    + Shr<u8, Output = Self>
    + BitOr<Self, Output = Self>
    + BitAnd<Self, Output = Self>
{
    const ZERO: Self;

    const ONE: Self;

    /// Maximum k-mer size `(BITS / 2 - 1)`.
    /// It is actually possible to store `BITS / 2` k-mers in case of minimizers and canonical k-mers, but not in case
    /// of non-canonical k-mers.
    /// To be on the safe side, we forbid this k-mer size.
    const MAX_KMER_SIZE: u8;

    /// Undefined k-mer/minimizer.
    /// In case of minimizers, it is possible that hash function would produce `UNDEF` value inadvertently.
    /// Nevertheless, such minimizers would have the highest hash function anyway, and should not go to the output.
    const UNDEF: Self;

    fn create_mask(k: u8) -> Self;

    fn from_u8(val: u8) -> Self;
}

macro_rules! impl_kmer {
    ($prim:ty, $zero:expr, $one:expr) => {
        impl Kmer for $prim {
            const ZERO: Self = $zero;

            const ONE: Self = $one;

            const MAX_KMER_SIZE: u8 = (Self::BITS / 2 - 1) as u8;

            const UNDEF: Self = Self::MAX;

            #[inline]
            fn create_mask(k: u8) -> Self {
                // if k == 16 { -1_i32 as u32 } else { (1_u32 << 2 * k) - 1 }
                // Already know that k is not `Self::BITS / 2`.
                (Self::ONE << 2 * k) - Self::ONE
            }

            #[inline]
            fn from_u8(val: u8) -> Self {
                Self::from(val)
            }
        }
    };

    ($prim:ty) => { impl_kmer!($prim, 0, 1); }
}

impl_kmer!(u8);
impl_kmer!(u16);
impl_kmer!(u32);
impl_kmer!(u64);
impl_kmer!(u128);
impl_kmer!(U256, U256::ZERO, U256::ONE);

pub trait Minimizer: Kmer {
    /// Quickly calculate hash function.
    /// For now, use Murmur3 function, this can change later.
    fn fast_hash(self) -> Self;
}

impl Minimizer for u32 {
    /// Murmur3 for uint32.
    #[inline]
    fn fast_hash(mut self) -> u32 {
        // Use bitwise not to convert 0 to u32::MAX.
        // This is needed because otherwise 0 would produce 0.
        self = !self;
        self ^= self >> 16;
        self = self.wrapping_mul(0x85ebca6b);
        self ^= self >> 13;
        // self = self.wrapping_mul(0xc2b2ae35);
        // self ^= self >> 16;
        self
    }
}

impl Minimizer for u64 {
    /// fasthash (https://github.com/rurban/smhasher/blob/master/fasthash.cpp) mix function.
    #[inline]
    fn fast_hash(mut self) -> u64 {
        self = !self;
        self ^= self >> 23;
        self = self.wrapping_mul(0x2127599bf4325c37);
        self ^= self >> 47;
        self
    }
}

impl Minimizer for u128 {
    #[inline]
    fn fast_hash(self) -> Self {
        let [h1, h2] = unsafe { std::mem::transmute::<u128, [u64; 2]>(self) };
        let h1 = h1.fast_hash();
        let h2 = h2.fast_hash();
        unsafe { std::mem::transmute::<[u64; 2], u128>([h1, h2]) }
    }
}

/// Output data structure for k-mers or minimizers.
/// Either Vec<(position, kmer)>, Vec<kmer> or IntSet<kmer>.
pub trait Output<K: Kmer> {
    fn push(&mut self, pos: u32, kmer: K, forward: bool);
}

impl<K: Kmer> Output<K> for Vec<(u32, K)> {
    #[inline(always)]
    fn push(&mut self, pos: u32, kmer: K, _forward: bool) {
        self.push((pos, kmer));
    }
}

impl<K: Kmer> Output<K> for Vec<K> {
    #[inline(always)]
    fn push(&mut self, _pos: u32, kmer: K, _forward: bool) {
        self.push(kmer);
    }
}

impl<K: Kmer> Output<K> for Vec<(u32, K, bool)> {
    #[inline(always)]
    fn push(&mut self, pos: u32, kmer: K, forward: bool) {
        self.push((pos, kmer, forward));
    }
}

impl<K: Kmer> Output<K> for Vec<(K, bool)> {
    #[inline(always)]
    fn push(&mut self, _pos: u32, kmer: K, forward: bool) {
        self.push((kmer, forward));
    }
}

impl<K: Kmer + std::hash::Hash + nohash::IsEnabled> Output<K> for IntSet<K> {
    #[inline(always)]
    fn push(&mut self, _pos: u32, kmer: K, _forward: bool) {
        self.insert(kmer);
    }
}

pub const CANONICAL: bool = true;
pub const NON_CANONICAL: bool = false;

/// Returns all k-mers in the sequence.
/// Returns `K::UNDEF` for k-mers containing N.
///
/// If `CANONICAL` is true, returns canonical k-mers (min between forward and reverse-complement k-mer).
pub fn kmers<K, O, const CANON: bool>(seq: &[u8], k: u8, output: &mut O)
where K: Kmer,
      O: Output<K>,
{
    debug_assert!(k <= K::MAX_KMER_SIZE, "k-mer size ({}) can be at most {}", k, K::MAX_KMER_SIZE);
    let mask = K::create_mask(k);
    let rv_shift = if CANON { 2 * k - 2 } else { 0 };
    let mut fw_kmer = K::ZERO;
    let mut rv_kmer = K::ZERO;

    let k = u32::from(k);
    let k_1 = k - 1;
    let mut reset = k_1;

    for (i, &nt) in seq.iter().enumerate() {
        let i = i as u32;
        let fw_enc: u8 = match nt {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                reset = i + k;
                if i + 1 >= k {
                    output.push(i - k_1, K::UNDEF, true);
                }
                continue;
            },
        };
        fw_kmer = ((fw_kmer << 2) | K::from_u8(fw_enc)) & mask;
        if CANON { rv_kmer = (rv_kmer >> 2) | (K::from_u8(3 - fw_enc) << rv_shift); }

        if i >= reset {
            let (kmer, forward) = if CANON && rv_kmer < fw_kmer { (rv_kmer, false) } else { (fw_kmer, true) };
            output.push(i - k_1, kmer, forward);
        } else if i + 1 >= k {
            output.push(i - k_1, K::UNDEF, true);
        }
    }
}

/// Biggest allowed window size (number of consecutive k-mer). Must be the power of two.
pub const MAX_MINIMIZER_W: u8 = 64;
const MOD_MAXW: u32 = MAX_MINIMIZER_W as u32 - 1;

/// Circular array of fixed size.
struct CircArray<T> {
    arr: [T; MAX_MINIMIZER_W as usize],
}

impl<T: Copy> CircArray<T> {
    #[inline(always)]
    fn new(val: T) -> Self {
        const _: () = assert!(MAX_MINIMIZER_W.count_ones() == 1);
        Self {
            arr: [val; MAX_MINIMIZER_W as usize],
        }
    }
}

impl<T> std::ops::Index<u32> for CircArray<T> {
    type Output = T;

    #[inline(always)]
    fn index(&self, i: u32) -> &T {
        unsafe { self.arr.get_unchecked((i & MOD_MAXW) as usize) }
    }
}

impl<T> std::ops::IndexMut<u32> for CircArray<T> {
    #[inline(always)]
    fn index_mut(&mut self, i: u32) -> &mut T {
        unsafe { self.arr.get_unchecked_mut((i & MOD_MAXW) as usize) }
    }
}

/// Function that goes over non-empty slice `start..end`, and returns minimal position and value.
#[inline(always)]
fn find_min<T: Copy + Ord>(
    arr: &CircArray<T>,
    start: u32,
    end: u32,
) -> (u32, T)
{
    debug_assert!(start < end);
    let mut pos = start;
    let mut min = arr[start];
    for j in start + 1..end {
        let v = arr[j];
        if v < min {
            pos = j;
            min = v;
        }
    }
    (pos, min)
}

/// Finds sequence minimizers.
/// Minimizer is a k-mer with the smallest hash value across `w` consecutive k-mers.
/// Output vector should have type `Vec<K>` or `Vec<(u32, K)>`, then minimizers are saved together with their positions.
///
/// NOTE: minimizers should be cleared in advance.
pub fn minimizers<K, O, const CANON: bool>(seq: &[u8], k: u8, w: u8, output: &mut O)
where K: Minimizer,
      O: Output<K>,
{
    debug_assert!(0 < k && k <= K::MAX_KMER_SIZE, "k-mer size ({}) can be at most {}", k, K::MAX_KMER_SIZE);
    debug_assert!(0 < w && w < MAX_MINIMIZER_W, "Minimizer window ({}) must be smaller than {}", w, MAX_MINIMIZER_W);
    let mask = K::create_mask(k);
    let rv_shift = if CANON { 2 * k - 2 } else { 0 };
    let mut fw_kmer = K::ZERO;
    let mut rv_kmer = K::ZERO;

    let k = u32::from(k);
    let k_1 = k - 1;
    let w = u32::from(w);
    let w_1 = w - 1;

    // Hashes in a window, stored in a cycling array.
    let mut hashes = CircArray::new(K::UNDEF);
    let mut forward = CircArray::new(true);
    let mut last_pos = -1_i32;
    let mut best_pos = 0;
    let mut best_hash = K::UNDEF;

    let mut first_kmer = k_1;
    let mut first_window = k_1 + w_1;

    for (i, &nt) in seq.iter().enumerate() {
        let i = i as u32;
        let (fw_enc, rv_enc): (u8, u8) = match nt {
            b'A' => (0, 3),
            b'C' => (1, 2),
            b'G' => (2, 1),
            b'T' => (3, 0),
            _ => {
                first_kmer = i + k;
                (0, 0)
            },
        };
        fw_kmer = ((fw_kmer << 2) | K::from_u8(fw_enc)) & mask;
        if CANON { rv_kmer = (rv_kmer >> 2) | (K::from_u8(rv_enc) << rv_shift); }
        let (kmer, fw) = if CANON && rv_kmer < fw_kmer { (rv_kmer, false) } else { (fw_kmer, true) };
        let h = if i < first_kmer { K::UNDEF } else { kmer.fast_hash() };
        hashes[i] = h;
        forward[i] = fw;

        if h < best_hash {
            best_hash = h;
            best_pos = i;
        }
        if i < first_window {
            continue;
        }

        let start = i - w_1;
        if best_pos < start {
            (best_pos, best_hash) = find_min(&hashes, start, i + 1);
            if best_hash == K::UNDEF {
                first_window = first_window + w_1;
                continue;
            }
        }
        if best_pos as i32 > last_pos {
            last_pos = best_pos as i32;
            output.push(best_pos - k_1, best_hash, forward[best_pos]);
        }
    }
}

/// Find canonical minimizers (wrapper around `minimizers`).
#[inline]
pub fn canon_minimizers<K, O>(seq: &[u8], k: u8, w: u8, output: &mut O)
where K: Minimizer,
      O: Output<K>,
{
    minimizers::<K, O, { CANONICAL }>(seq, k, w, output)
}
