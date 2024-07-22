use std::{
    cmp::{min, Ord},
    ops::{Shl, Shr, BitOr, BitAnd},
};

/// Store k-mers in integers of different length (15-mers in u32, 31-mers in u64, 63-mers in u128).
pub trait Kmer: From<u8>
    + Copy + Ord + Eq
    + std::fmt::Display + std::fmt::Binary + std::fmt::LowerHex + std::fmt::UpperHex
    + Shl<u8, Output = Self>
    + Shr<u8, Output = Self>
    + BitOr<Self, Output = Self>
    + BitAnd<Self, Output = Self>
{
    /// Zero value.
    const ZERO: Self;

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
}

macro_rules! impl_kmer {
    ($prim:ty) => {
        impl Kmer for $prim {
            const ZERO: Self = 0;

            const MAX_KMER_SIZE: u8 = (Self::BITS / 2 - 1) as u8;

            const UNDEF: Self = Self::MAX;

            #[inline]
            fn create_mask(k: u8) -> Self {
                // if k == 16 { -1_i32 as u32 } else { (1_u32 << 2 * k) - 1 }
                // Already know that k is not `Self::BITS / 2`.
                (1 << 2 * k) - 1
            }
        }
    }
}

impl_kmer!(u8);
impl_kmer!(u16);
impl_kmer!(u32);
impl_kmer!(u64);
impl_kmer!(u128);

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

pub const CANONICAL: bool = true;
pub const NON_CANONICAL: bool = false;

/// Returns all k-mers in the sequence.
/// Returns `K::UNDEF` for k-mers containing N.
///
/// If `CANONICAL` is true, returns canonical k-mers (min between forward and reverse-complement k-mer).
pub fn kmers<K: Kmer, const CANON: bool>(seq: &[u8], k: u8, output: &mut Vec<K>) {
    debug_assert!(k <= K::MAX_KMER_SIZE, "k-mer size ({}) can be at most {}", k, K::MAX_KMER_SIZE);
    let mask = K::create_mask(k);
    let rv_shift = if CANON { 2 * k - 2 } else { 0 };
    let mut fw_kmer = K::ZERO;
    let mut rv_kmer = K::ZERO;

    let k = u32::from(k);
    let mut reset = k - 1;

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
                    output.push(K::UNDEF);
                }
                continue;
            },
        };
        fw_kmer = ((fw_kmer << 2) | K::from(fw_enc)) & mask;
        if CANON { rv_kmer = (rv_kmer >> 2) | (K::from(3 - fw_enc) << rv_shift); }

        if i >= reset {
            output.push(if CANON { min(fw_kmer, rv_kmer) } else { fw_kmer });
        } else if i + 1 >= k {
            output.push(K::UNDEF);
        } // else { }
    }
}

/// Trait that summarizes tuple `(u32, kmer)` and simply `kmer`.
/// In the first case, minimizers are saved together with their positions,
/// in the other, only kmers are saved.
pub trait PosKmer<K> {
    fn transform(pos: u32, kmer: K) -> Self;

    /// Replace position with reverse complement position.
    /// n = seq_len - k.
    fn rc_pos(&mut self, n: u32);
}

impl<K: Kmer> PosKmer<K> for K {
    #[inline(always)]
    fn transform(_: u32, kmer: K) -> K {
        kmer
    }

    #[inline(always)]
    fn rc_pos(&mut self, _: u32) {}
}

impl<K: Kmer> PosKmer<K> for (u32, K) {
    #[inline(always)]
    fn transform(pos: u32, kmer: K) -> (u32, K) {
        (pos, kmer)
    }

    #[inline(always)]
    fn rc_pos(&mut self, n: u32) {
        self.0 = n - self.0;
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

/// Function that goes over indices `start..end`, and returns new `start`.
/// Additionally, the function pushes the new minimizer to the output, if it is not `UNDEF`.
#[inline]
fn select_minimizer<K: Kmer, P: PosKmer<K>>(
    start: u32,
    end: u32,
    k: u32,
    hashes: &mut [K; MAX_MINIMIZER_W as usize],
    output: &mut Vec<P>,
) -> u32 {
    let mut minimizer = K::UNDEF;
    let mut kmer_end = end;
    for j in start..end {
        let h = hashes[(j & MOD_MAXW) as usize];
        if h < minimizer {
            minimizer = h;
            kmer_end = j + 1;
        }
    }
    if minimizer != K::UNDEF {
        output.push(P::transform(kmer_end - k, minimizer));
    }
    kmer_end
}

/// Finds sequence minimizers.
/// Minimizer is a k-mer with the smallest hash value across `w` consecutive k-mers.
/// Output vector should have type `Vec<K>` or `Vec<(u32, K)>`, then minimizers are saved together with their positions.
///
/// NOTE: minimizers should be cleared in advance.
pub fn minimizers0<K, P, const CANON: bool>(seq: &[u8], k: u8, w: u8, output: &mut Vec<P>)
where K: Minimizer,
      P: PosKmer<K>,
{
    debug_assert!(k <= K::MAX_KMER_SIZE, "k-mer size ({}) can be at most {}", k, K::MAX_KMER_SIZE);
    debug_assert!(w <= MAX_MINIMIZER_W, "Minimizer window ({}) can be at most {}", w, MAX_MINIMIZER_W);
    let mask = K::create_mask(k);
    let rv_shift = 2 * k - 2;
    // Hashes in a window, stored in a cycling array.
    let mut hashes = [K::UNDEF; MAX_MINIMIZER_W as usize];

    let mut fw_kmer = K::ZERO;
    let mut rv_kmer = K::ZERO;

    let k = u32::from(k);
    let w = u32::from(w);
    // At what index will the first k-mer be available.
    let mut reset = k - 1;
    // Start of the window with consecutive k-mers.
    let mut start = reset;

    for (i, &nt) in seq.iter().enumerate() {
        let i = i as u32;
        let (fw_enc, rv_enc): (u8, u8) = match nt {
            b'A' => (0, 3),
            b'C' => (1, 2),
            b'G' => (2, 1),
            b'T' => (3, 0),
            _ => {
                reset = i + k;
                if reset >= start + w {
                    if i > start {
                        select_minimizer(start, i, k, &mut hashes, output);
                    }
                    start = reset;
                }
                hashes[(i & MOD_MAXW) as usize] = K::UNDEF;
                continue;
            },
        };
        fw_kmer = ((fw_kmer << 2) | K::from(fw_enc)) & mask;
        if CANON {
            rv_kmer = (K::from(rv_enc) << rv_shift) | (rv_kmer >> 2);
        }
        if i < reset {
            hashes[(i & MOD_MAXW) as usize] = K::UNDEF;
            continue;
        }

        let kmer = if CANON {
            min(fw_kmer, rv_kmer)
        } else {
            fw_kmer
        };
        hashes[(i & MOD_MAXW) as usize] = kmer.fast_hash();
        if i == start + w - 1 {
            start = select_minimizer(start, i + 1, k, &mut hashes, output);
        }
    }
    let l = seq.len() as u32;
    if l >= start {
        debug_assert!(l <= start + w - 1);
        select_minimizer(start, l, k, &mut hashes, output);
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
pub fn minimizers<K, P, const CANON: bool>(seq: &[u8], k: u8, w: u8, output: &mut Vec<P>)
where K: Minimizer,
      P: PosKmer<K>,
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

    // Hashes in a window, stored in a cycling array.
    let mut hashes = CircArray::new(K::UNDEF);
    // Position of the best k-mer.
    // In addition, it indicates the position, where first k-mer becomes available.
    let mut hash_pos = k_1;
    let mut min_hash = K::UNDEF;

    // New window starts at least at this position.
    let mut min_start = k_1 + w;
    // and at most at this position.
    let mut max_start = min_start;

    for (i, &nt) in seq.iter().enumerate() {
        let i = i as u32;
        let (fw_enc, rv_enc): (u8, u8) = match nt {
            b'A' => (0, 3),
            b'C' => (1, 2),
            b'G' => (2, 1),
            b'T' => (3, 0),
            _ => {
                if min_hash != K::UNDEF && min_start <= i {
                    output.push(P::transform(hash_pos - k_1, min_hash));
                }
                hash_pos = i + k;
                min_hash = K::UNDEF;
                min_start = hash_pos + w;
                max_start = min_start;
                continue;
            },
        };
        fw_kmer = ((fw_kmer << 2) | K::from(fw_enc)) & mask;
        if CANON { rv_kmer = (rv_kmer >> 2) | (K::from(rv_enc) << rv_shift); }
        if i < hash_pos {
            continue;
        }
        let hash = if CANON { min(fw_kmer, rv_kmer) } else { fw_kmer }.fast_hash();
        hashes[i] = hash;

        if min_hash == K::UNDEF {
            min_hash = hash;
            hash_pos = i;
        } else if max_start <= i {
            output.push(P::transform(hash_pos - k_1, min_hash));
            min_start = hash_pos + w + 1;
            (hash_pos, min_hash) = find_min(&hashes, hash_pos + 1, i + 1);
            max_start = hash_pos + w;
        } else if hash < min_hash {
            if min_start <= i {
                output.push(P::transform(hash_pos - k_1, min_hash));
                min_start = i + 1;
            }
            hash_pos = i;
            min_hash = hash;
        }
    }
    if min_hash != K::UNDEF && min_start <= seq.len() as u32 {
        output.push(P::transform(hash_pos - k_1, min_hash));
    }
}

/// Find canonical minimizers (wrapper around `minimizers`).
#[inline]
pub fn canon_minimizers<K, P>(seq: &[u8], k: u8, w: u8, output: &mut Vec<P>)
where K: Minimizer,
      P: PosKmer<K>,
{
    minimizers::<K, P, { CANONICAL }>(seq, k, w, output)
}

/// Density of mod-minimizers.
fn mod_density(k: u8, w: u8, t: u8) -> f64 {
    let kf = f64::from(k);
    let wf = f64::from(w);
    let tf = f64::from(t);
    let l = wf + kf - 1.0;
    let x = if t == k { 0.0 } else { 1.0 / (l - tf + 1.0) };
    let density = (2.0 + ((l - tf) / wf).ceil() * (1.0 - x)) / (l - tf + 1.0);
    density
}

/// Finds such `t`, that would produce the smallest minimizer density.
/// If `very_small` is true, consider very small `t`, under theoretical limit (practically, should still work).
/// If `fw_only` is set, only allow forward-only minimizer schemes.
pub fn find_best_t(k: u8, w: u8, very_small: bool, fw_only: bool) -> u8 {
    let mut best_t = 0;
    let mut best_density = f64::INFINITY;
    let start_t = if very_small {
        min(k, 3)
    } else {
        (1.5 * (f64::from(w) + f64::from(k) - 1.0).log2()).ceil().min(f64::from(k)) as u8
    };
    for t in start_t..=k {
        if fw_only && t % w != k % w && t % w != (k + 1) % w {
            continue;
        }
        let density = mod_density(k, w, t);
        if density < best_density {
            best_t = t;
            best_density = density;
        }
    }
    best_t
}

/// Finds sequence mod-minimizers (https://doi.org/10.1101/2024.05.25.595898)
/// for the forward sequence.
/// Output vector should have type `Vec<K>` or `Vec<(u32, K)>`, then minimizers are saved together with their positions.
///
/// NOTE: minimizers should be cleared in advance.
pub fn mod_minimizers<K, P, const CANON: bool>(seq: &[u8], k: u8, w: u8, t: u8, output: &mut Vec<P>)
where K: Minimizer,
      P: PosKmer<K>,
{
    println!("mod-minimizers({k}, {w}, {t}) for {}", String::from_utf8_lossy(seq));
    debug_assert!(k <= K::MAX_KMER_SIZE, "k-mer size ({}) can be at most {}", k, K::MAX_KMER_SIZE);
    debug_assert!(t <= k, "t-mer size ({}) must not be greater than k-mer size ({})", t, k);
    debug_assert!(w.checked_add(k - t).unwrap() <= MAX_MINIMIZER_W,
        "Minimizer window w ({}) + k ({}) - t ({}) can be at most {}", w, k, t, MAX_MINIMIZER_W);
    let rv_shift = 2 * k - 2;
    let t_shift = 2 * k - 2 * t;
    let k_mask = K::create_mask(k);
    let t_mask = K::create_mask(t);

    // canonical t-mer hashes and forward t-mers,
    // needed to select the same direction when going on the forward and reverse-complement.
    let mut t_keys = CircArray::new((K::UNDEF, K::UNDEF));
    let mut t_strands = CircArray::new(true);
    let mut kmers = CircArray::new(K::UNDEF);
    let mut fw_kmer = K::ZERO;
    let mut rv_kmer = K::ZERO;

    let k = u32::from(k);
    let k_1 = k - 1;
    let t = u32::from(t);
    let t_1 = t - 1;
    let w = u32::from(w);
    let w_1 = w - 1;
    let kwt_1 = k + w - t - 1;

    // let mut first_kmer = k_1;
    let mut first_tmer = t_1;
    let mut first_window = w_1 + k_1;

    let mut best_key = (K::UNDEF, K::UNDEF);
    let mut best_t_pos = 0;

    for (i, &nt) in seq.iter().enumerate() {
        let i = i as u32;
        let (fw_enc, rv_enc): (u8, u8) = match nt {
            b'A' => (0, 3),
            b'C' => (1, 2),
            b'G' => (2, 1),
            b'T' => (3, 0),
            _ => {
                // First available position will be i + 1.
                // first_kmer = i + k;
                first_tmer = i + t;
                first_window = i + w + k_1;
                best_key = (K::UNDEF, K::UNDEF);
                continue;
            },
        };
        fw_kmer = ((fw_kmer << 2) | K::from(fw_enc)) & k_mask;
        if CANON { rv_kmer = (rv_kmer >> 2) | (K::from(rv_enc) << rv_shift); }
        if i < first_tmer { continue };

        // This may be incorrect for first k-mers, but we should not access them.
        kmers[i] = if CANON { min(fw_kmer, rv_kmer) } else { fw_kmer };
        println!("   {:3} '{}'{:16x}  {:030b}", i - t_1, nt as char, kmers[i], fw_kmer);
        let fw_tmer = fw_kmer & t_mask;
        let fw_t_hash = fw_tmer.fast_hash();
        let (t_hash, strand) = if CANON {
            let rv_t_hash = (rv_kmer >> t_shift).fast_hash();
            println!("       F  {:016x}  {:010b}", fw_t_hash, fw_tmer);
            println!("       B  {:016x}  {:010b}", rv_t_hash, rv_kmer >> t_shift);
            if fw_t_hash <= rv_t_hash { (fw_t_hash, true) } else { (rv_t_hash, false) }
        } else { (fw_t_hash, true) };
        let key = (t_hash, fw_tmer);
        println!("       K  {:016x}  {:010b}  {}", key.0, key.1, if strand { "+" } else { "-" });

        let t_pos = i - t_1;
        t_keys[t_pos] = key;
        t_strands[t_pos] = strand;
        if key < best_key {
            best_t_pos = t_pos;
            best_key = key;
            println!("       Update best");
        }
        if i < first_window { continue };

        let start = i - w_1 - k_1;
        if best_t_pos < start {
            (best_t_pos, best_key) = find_min(&t_keys, start, t_pos + 1);
            println!("       Find min -> {:2}  {:016x}  {:010b}", best_t_pos, best_key.0, best_key.1);
        }

        let j = best_t_pos - start;
        let new_pos = if CANON && !t_strands[best_t_pos] {
            start + w_1 - (kwt_1 - j) % w
        } else {
            start + j % w
        };
        // Replace previous value with UNDEF, so that we don't output the same k-mer twice.
        let sel_kmer = std::mem::replace(&mut kmers[new_pos + k_1], K::UNDEF);
        println!("       Save {:2} -> {:2} ({})", best_t_pos, new_pos, if sel_kmer == K::UNDEF { "dupl" } else { "new" });
        if sel_kmer != K::UNDEF {
            output.push(P::transform(new_pos, sel_kmer));
        }
    }
}

// /// Canonical mod-minimizers.
// /// In contrast to `canon_minimizers`, this function simply concatenates minimizers from forward
// /// and reverse complement strands.
// pub fn canon_mod_minimizers<K, P>(seq: &[u8], k: u8, w: u8, t: u8, output: &mut Vec<P>)
// where K: Minimizer,
//       P: PosKmer<K>,
// {
//     if seq.len() < usize::from(k) {
//         return;
//     }

//     mod_minimizers(seq, k, w, t, output);
//     let s = output.len();
//     mod_minimizers(&super::reverse_complement(seq), k, w, t, output);
//     let n = seq.len() as u32 - u32::from(k);
//     output[s..].iter_mut().for_each(|pos_kmer| pos_kmer.rc_pos(n));
// }

/// Naïve implementation of mod-minimizers.
pub fn naive_mod_minimizers<K, P, const CANON: bool>(seq: &[u8], k: u8, w: u8, t: u8, output: &mut Vec<P>)
where K: Minimizer,
      P: PosKmer<K>,
{
    let k_mask = K::create_mask(k);
    let t_mask = K::create_mask(t);
    let rv_shift = 2 * k - 2;
    let t_shift = 2 * k - 2 * t;
    let mut fw_kmer = K::ZERO;
    let mut rv_kmer = K::ZERO;
    let mut kmers = Vec::new();
    let mut used = Vec::new();
    let mut tmers = Vec::new();
    let mut reset = 0;

    let k = u32::from(k);
    let t = u32::from(t);
    let w = u32::from(w);
    'outer: for (i, &nt) in seq.iter().enumerate() {
        let i = i as u32;
        let (fw_enc, rv_enc): (u8, u8) = match nt {
            b'A' => (0, 3),
            b'C' => (1, 2),
            b'G' => (2, 1),
            b'T' => (3, 0),
            _ => {
                reset = i + 1;
                (0, 0)
            },
        };

        fw_kmer = ((fw_kmer << 2) | K::from(fw_enc)) & k_mask;
        rv_kmer = (rv_kmer >> 2) | (K::from(rv_enc) << rv_shift);
        if i >= k - 1 {
            kmers.push(if i < reset + k - 1 { K::UNDEF } else
                if CANON { min(fw_kmer, rv_kmer) } else { fw_kmer });
            used.push(false);
        }
        if i >= t - 1 {
            let fw_hash = (fw_kmer & t_mask).fast_hash();
            let rv_hash = (rv_kmer >> t_shift).fast_hash();
            let hash = if i < reset + t - 1 {
                (K::UNDEF, true)
            } else if CANON && rv_hash < fw_hash {
                (rv_hash, false)
            } else {
                (fw_hash, true)
            };
            tmers.push(hash);
        }

        let Some(start) = (i + 2).checked_sub(w + k) else { continue };
        let mut best_tmer = K::UNDEF;
        let mut best_pos = 0;
        for j in 0..w + k - t {
            let (tmer, fw) = tmers[(start + j) as usize];
            if tmer == K::UNDEF {
                continue 'outer;
            }
            if tmer < best_tmer {
                best_tmer = tmer;
                best_pos = if fw {
                    start + (j % w)
                } else {
                    let j2 = k + w - j - t - 1;
                    start + w - 1 - (j2 % w)
                };
            }
        }
        if !std::mem::replace(&mut used[best_pos as usize], true) {
            output.push(P::transform(best_pos, kmers[best_pos as usize]));
        }
    }
}
