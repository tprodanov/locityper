use std::cmp::min;

/// Adds a new k-mer to k-mer counts.
#[inline]
fn add_counts(counts: &mut [u8], unique: &mut u8, kmer: u64) {
    if counts[kmer as usize] == 0 {
        *unique += 1;
    }
    counts[kmer as usize] += 1;
}

/// Adds a new k-mer and removes one k-mer from k-mer counts.
#[inline]
fn update_counts(counts: &mut [u8], unique: &mut u8, remove_kmer: u64, add_kmer: u64) {
    if remove_kmer != add_kmer {
        add_counts(counts, unique, add_kmer);
        if counts[remove_kmer as usize] == 1 {
            *unique -= 1;
        }
        counts[remove_kmer as usize] -= 1;
    }
}

const N_CODE: u8 = 4;

/// Encodes nucleotide sequence into two-bits. Stores Ns as `N_CODE`.
#[inline]
fn seq_to_2bit(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&nt| match nt {
        b'A' => 0b00,
        b'C' => 0b01,
        b'G' => 0b10,
        b'T' => 0b11,
        b'N' => N_CODE,
        _ => panic!("Unexpected nucleotide {}", nt as char),
    }).collect()
}

/// Limit k to <= 6, because we use straight vector to count the number of k-mer occurances,
/// and the vector grows exponentially.
const MAX_K_VAL: usize = 6;

/// Calculates the number of unique k-mers in a moving window of size `w`.
/// The resulting numbers for window `i..i+w` are written to `out[i + w/2]` by multiplying by the old value.
fn abundances(seq_2b: &[u8], out: &mut [u16], w: usize, k: usize) {
    let n = seq_2b.len();
    assert!(w >= k && w <= n && w <= u8::MAX as usize + 1 - k && k <= MAX_K_VAL,
        "Cannot calculate k-mer abundance for k = {}, window = {}, and seq. length = {}", k, w, n);
    let halfw = w / 2;
    let old_shift = 2 * (w - k + 1);
    // NOTE: once `generic_const_exprs` becomes stable, this can become a const generic function on k.
    let mut counts = vec![0_u8; 4_usize.pow(k as u32)];
    let mut unique = 0;

    let mask = 2_u64.pow(2 * k as u32) - 1;
    let mut kmer = 0_u64;
    let mut last_reset = 0;
    for (i, &code) in seq_2b.iter().enumerate() {
        if code == N_CODE {
            last_reset = i + 1;
            if i != last_reset {
                counts.fill(0);
                unique = 0;
            } else {
                debug_assert_eq!(unique, 0);
            }
            continue;
        }

        kmer = (kmer << 2) | (code as u64);
        if i >= last_reset + w {
            update_counts(&mut counts, &mut unique, (kmer >> old_shift) & mask, kmer & mask);
        } else if i >= last_reset + k - 1 {
            add_counts(&mut counts, &mut unique, kmer & mask);
        }
        if i >= last_reset + w - 1 {
            out[i - halfw] *= unique as u16;
        }
    }
}

/// Calculates the linguistic complexity of a sequence.
/// Complexity is calculated across moving windows of size `w` (must be odd).
///
/// For each window of size `w`, the function calculates the product of *U_k* of all `k in k1..=k2`,
/// where *U_k* is the number of unique k-mers in the window divided by the maximal possible number of k-mers
/// `min(4^k, w-k+1)`.
///
/// Output complexity is written in the middle of each window (`i -> i + w/2`).
/// Returns NaN for windows close to the edge and windows that contain Ns.
pub fn linguistic_complexity(seq: &[u8], w: usize, k1: usize, k2: usize
) -> impl Iterator<Item = f32> + std::iter::ExactSizeIterator {
    let n = seq.len();
    assert!(w % 2 == 1, "Window size ({}) must be odd", w);
    assert!(k1 <= k2 && k2 <= MAX_K_VAL, "Impossible k values {} and {}", k1, k2);
    assert!(w <= n, "Window size ({}) must not be greater than the sequence length ({})", w, n);
    let seq_2b = seq_to_2bit(seq);

    let mut abund = vec![1_u16; n];
    let halfw = w / 2;
    let mut n_stretch_start = usize::MAX;
    for (i, &code) in seq_2b.iter().enumerate() {
        if code == N_CODE {
            if n_stretch_start == usize::MAX {
                n_stretch_start = i;
            }
        } else if n_stretch_start != usize::MAX {
            abund[n_stretch_start.saturating_sub(halfw)..min(i + halfw, n)].fill(0);
            n_stretch_start = usize::MAX;
        }
    }
    if n_stretch_start != usize::MAX {
        abund[n_stretch_start.saturating_sub(halfw)..n].fill(0);
    }

    let mut divisor = 1.0;
    for k in k1..=k2 {
        abundances(&seq_2b, &mut abund, w, k);
        divisor *= min(4_u32.pow(k as u32), (w + 1 - k) as u32) as f32;
    }
    abund[..halfw].fill(0);
    abund[n - halfw..].fill(0);
    let coeff = 1.0 / divisor;
    abund.into_iter().map(move |a| if a == 0 { f32::NAN } else { coeff * a as f32 })
}

/// Linguistic complexity of a sequence, with k-mer size between 1 and 3 (see `linguistic_complexity`).
#[inline]
pub fn linguistic_complexity_k3(seq: &[u8], w: usize) -> impl Iterator<Item = f32> + std::iter::ExactSizeIterator {
    linguistic_complexity(seq, w, 1, 3)
}