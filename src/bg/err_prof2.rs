use std::cmp::min;

const ABUNDANCE_ERROR: u32 = 0;

/// Calculates: <monomer abundance> * <dimer abundance> * <trimer abundance>.
/// where k-mer abundance is the number of unique k-mers in the sequence.
/// Returns ABUNDANCE_ERROR (= 0) if the sequence contains Ns.
fn kmer_abundance_123(seq: &[u8]) -> u32 {
    let n = seq.len();
    // Resulting value may be too high for sequences of length > 1626.
    // Max number of values is (n-2) * (n-1) * n.
    debug_assert!(n >= 3 && n <= 1626,
        "Cannot calculate k-mer abundance for sequence of length {} (must be between 3 and 1626)", n);

    let mut abundance1 = 0_u8;
    let mut abundance2 = 0_u16;
    let mut abundance3 = 0_u64;
    let mut kmer3 = 0_u8;
    const FOUR_BITS: u8 = 0b00001111;
    for i in 0..n {
        let nt: u8 = match unsafe { *seq.get_unchecked(i) } {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => return ABUNDANCE_ERROR,
            ch => panic!("Unexpected nucleotide {}", ch as char),
        };
        abundance1 |= 1_u8 << nt;
        kmer3 = ((kmer3 & FOUR_BITS) << 2) | nt;
        if i > 1 {
            abundance2 |= 1_u16 << (kmer3 & FOUR_BITS);
            abundance3 |= 1_u64 << kmer3;
        } else if i > 0 {
            abundance2 |= 1_u16 << (kmer3 & FOUR_BITS);
        }
    }
    abundance1.count_ones() * abundance2.count_ones() * abundance3.count_ones()
}

/// Sequence complexity: product of k-mer abundances for k=1,2,3 divided by the total possible number
/// of the corresponding k-mers.
/// Returns NAN if the sequence contains Ns.
fn complexity_123(seq: &[u8]) -> f32 {
    let ab = kmer_abundance_123(seq);
    if ab == ABUNDANCE_ERROR {
        f32::NAN
    } else {
        let n = seq.len() as u16;
        ab as f32 / min(n, 4) as f32 / min(n - 1, 16) as f32 / min(n - 2, 64) as f32
    }
}

/// Calculates sequence complexity for all 7-mers in the sequence.
/// `result[i]` is calculated according to `complexity_123(seq[i - 3..=i + 3])`.
///
/// k-mer abundance in [1, 30] returns 0, and is considered low complexity,
/// abundance in [31, 60] returns 1, and is considered medium complexity,
/// abundance in [61, 120] returns 2, and is considered high complexity.
/// k-mer abundance over 120 is impossible for 7-mers.
///
/// Tails are filled with the last available values.
/// k-mers with Ns are filled with Complexity::Low.
///
/// Output length = interval length.
pub fn sevenmer_complexities(seq: &[u8]) -> Vec<u8> {
    const SIZE: usize = 7;
    assert!(seq.len() >= SIZE, "Cannot calculate seven-mer complexities for sequence of length {}", seq.len());
    // Seven-mer strings have padding = 3 around its center.
    const PADD: usize = SIZE / 2;
    // println!("{}", String::from_utf8_lossy(seq));

    // complexity = min(abundance / DIVISOR, 2).
    // Value 61 would also map to 1 instead of 2, but abundance = 61 cannot be achieved.
    const DIVISOR: u8 = 31;

    let n = seq.len();
    let mut res = vec![0; n];
    // Last calculated index.
    let r = n - PADD - 1;
    for i in PADD..=r {
        res[i] = min(kmer_abundance_123(&seq[i-PADD..=i+PADD]) as u8 / DIVISOR, 2);
        // println!("    [{:2}] {}  ->  {:3}  ~ {}", i, String::from_utf8_lossy(&seq[i-PADD..=i+PADD]),
        //     kmer_abundance_123(&seq[i-PADD..=i+PADD]), res[i]);
    }
    for i in 0..PADD {
        res[i] = res[PADD];
    }
    for i in r + 1..n {
        res[i] = res[r];
    }
    res
}
