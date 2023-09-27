use std::cmp::min;

type Kmer = u8;
type Count = u16;

/// Adds a new k-mer and removes one k-mer from k-mer counts.
fn remove_add_kmer(counts: &mut [Count], remove: Kmer, add: Kmer, unique: &mut u16) {
    if remove != add {
        if counts[add as usize] == 0 {
            *unique += 1;
        }
        counts[add as usize] += 1;
        if counts[remove as usize] == 1 {
            *unique -= 1;
        }
        counts[remove as usize] -= 1;
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
    const MASK1: Kmer = 0b000011;
    const MASK2: Kmer = 0b001111;
    const MASK3: Kmer = 0b111111;

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
