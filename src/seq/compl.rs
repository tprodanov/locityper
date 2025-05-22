use std::cmp::min;
use crate::{
    seq::kmers::{self, NON_CANONICAL},
    algo::{IntMap},
};

trait Counts<K>
where usize: From<K>,
{
    fn get_mut(&mut self, kmer: K) -> &mut u16;
}

impl Counts<u8> for [u16] {
    #[inline(always)]
    fn get_mut(&mut self, kmer: u8) -> &mut u16 {
        &mut self[usize::from(kmer)]
    }
}

impl<const N: usize> Counts<u8> for [u16; N] {
    #[inline(always)]
    fn get_mut(&mut self, kmer: u8) -> &mut u16 {
        &mut self[usize::from(kmer)]
    }
}

impl Counts<u16> for IntMap<u16, u16> {
    #[inline(always)]
    fn get_mut(&mut self, kmer: u16) -> &mut u16 {
        self.entry(kmer).or_default()
    }
}

/// Adds a new k-mer and removes one k-mer from k-mer counts.
fn remove_add_kmer<K, T>(counts: &mut T, remove: K, add: K, unique: &mut u16)
where usize: From<K>,
      K: PartialEq,
      T: Counts<K>,
{
    if remove != add {
        let c1 = counts.get_mut(add);
        *unique += u16::from(*c1 == 0);
        *c1 += 1;
        let c2 = counts.get_mut(remove);
        *unique -= u16::from(*c2 == 1);
        *c2 -= 1;
    }
}

#[inline]
fn encode_nt(nt: &u8) -> u8 {
    match nt {
        b'A' => 0b00,
        b'C' => 0b01,
        b'G' => 0b10,
        b'T' => 0b11,
        _ => panic!("Cannot calculate linguistic complexity in the presence of unknown nucleotides"),
    }
}

/// Linguistic complexity of a sequence = U1 × U2 × U3, where Uk is the fraction of unique k-mers out of max available.
/// Sequence must not contain Ns.
pub fn linguistic_complexity_123(seq: &[u8], w: usize) -> Vec<f64> {
    const MASK1: u8 = 0b000011;
    const MASK2: u8 = 0b001111;
    const MASK3: u8 = 0b111111;

    let n = seq.len();
    assert!(3 <= w, "Window size ({}) must be over 3", w);
    assert!(w <= n, "Window size ({}) must not be greater than the sequence length ({})", w, n);
    let w_u16 = u16::try_from(w).expect("Window size must fit in two bytes");
    // Maximum achievable count multiplication.
    let divisor = (min(4, w) * min(16, w - 1) * min(64, w - 2)) as f64;

    let mut counts1 = [0; 4];
    let mut counts2 = [0; 16];
    let mut counts3 = [0; 64];
    // Fill counts to account for corresponding numbers of lagging 0 values.
    counts1[0] = w_u16;
    counts2[0] = w_u16 - 1;
    counts3[0] = w_u16 - 2;
    let mut uniq1 = 1;
    let mut uniq2 = 1;
    let mut uniq3 = 1;
    let mut complexities = Vec::with_capacity(n - w + 1);

    let mut lagging_kmer = 0;
    let mut current_kmer = 0;
    let lagging_iter = std::iter::repeat(0).take(w - 2).chain(seq.iter().map(encode_nt));
    let current_iter = seq.iter().map(encode_nt);
    for (i, (lagging_enc, current_enc)) in lagging_iter.zip(current_iter).enumerate() {
        lagging_kmer = (lagging_kmer << 2) | lagging_enc;
        current_kmer = (current_kmer << 2) | current_enc;
        remove_add_kmer(&mut counts1, (lagging_kmer >> 4) & MASK1, current_kmer & MASK1, &mut uniq1);
        remove_add_kmer(&mut counts2, (lagging_kmer >> 2) & MASK2, current_kmer & MASK2, &mut uniq2);
        remove_add_kmer(&mut counts3, lagging_kmer & MASK3, current_kmer & MASK3, &mut uniq3);
        if i >= w - 1 {
            complexities.push(f64::from(uniq1 * uniq2 * uniq3) / divisor);
        }
    }
    assert_eq!(complexities.len(), n - w + 1);
    complexities
}

// /// Counts the fraction of unique k-mers in a sequence out of the total possible.
// /// k should be small enough that k-mer fits into u16 (<= 15), but large enough
// /// that sequence length is significantly larger than 4 ** k (for speed reasons).
// pub fn frac_unique_kmers(seq: &[u8], k: u8, set: &mut IntSet<u16>) -> f64 {
//     debug_assert!(k <= <u16 as kmers::Kmer>::MAX_KMER_SIZE);
//     set.clear();
//     kmers::kmers::<u16, _, NON_CANONICAL>(seq, k, set);
//     set.len() as f64 / (seq.len() + 1 - k as usize) as f64
// }

/// Calculate simple complexity (fraction of unique k-mers) across each moving window of size w.
/// kmers and counts are buffers, and will be completely rewritten.
pub fn simple_complexity(seq: &[u8], k: u8, w: usize) -> Vec<f64> {
    let n = seq.len();
    debug_assert!(k <= <u16 as kmers::Kmer>::MAX_KMER_SIZE);
    debug_assert!(w < n && w < usize::from(u16::MAX) && w > 3 * usize::from(k));

    let mut kmers = Vec::with_capacity(n - k as usize + 1);
    kmers::kmers::<u16, _, NON_CANONICAL>(seq, k, &mut kmers);

    let k = k as usize;
    let mult = 1.0 / (w + 1 - k) as f64;
    let mut counts = IntMap::default();
    let mut unique = 0;
    for &kmer in &kmers[..(w - k + 1)] {
        let c = counts.entry(kmer).or_default();
        unique += u16::from(*c == 0);
        *c += 1;
    }

    let mut res = Vec::with_capacity(n - w + 1);
    res.push(unique as f64 * mult);
    for (&lag_kmer, &kmer) in kmers.iter().zip(&kmers[(w - k + 1)..]) {
        remove_add_kmer(&mut counts, lag_kmer, kmer, &mut unique);
        res.push(unique as f64 * mult);
    }
    res
}
