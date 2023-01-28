use std::{
    cmp::{min, max},
    rc::Rc,
    fmt,
};
use htslib::bam::Record;
use intmap::{IntMap, Entry};
use crate::{
    seq::{
        contigs::ContigNames,
        interv::Interval,
        cigar::Operation,
        aln::{ReadEnd, Alignment},
    },
};

const ABUNDANCE_ERROR: u32 = 0;

/// Calculates: <unique 1-mers> * <unique 2-mers> * <unique 3-mers>.
/// Returns ABUNDANCE_ERROR (= 0) if the sequence contains Ns.
///
/// Maximal number if min(4^3, n-2) * min(4^2, n-1) * min(4, n), where n is the length of the sequence.
fn kmer_abundance_123(seq: &[u8]) -> u32 {
    let n = seq.len();
    assert!(n >= 3, "Cannot calculate k-mer abundance for sequence of length {} (must be >= 3)", n);

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

macro_rules! kmer_counts {
    ($name:ident, $k:expr, $mask:expr) => {
        struct $name {
            counts: [u8; 4_usize.pow($k)],
            unique: u8,
            add_countdown: u8,
            remove_countdown: u8,
        }

        impl $name {
            fn new(w: usize) -> Self {
                Self {
                    counts: [0_u8; 4_usize.pow($k)],
                    unique: 0,
                    add_countdown: $k,
                    remove_countdown: $k + u8::try_from(w).unwrap(),
                }
            }

            fn reset(&mut self, w: usize) {
                // Compiler will convert to memset.
                self.counts.iter_mut().for_each(|x| *x = 0);
                self.unique = 0;
                self.countdown = $k;
            }

            fn update(&mut self, remove: u8, add: u8) {
                debug_assert_eq!(self.countdown, 0);
                let remove = remove & $mask;
                let add = add & $mask;
                println!("    k{} {:?}", $k, self.counts);
                print!("    k{}: Remove {:06b}, add {:06b}, unique was {}",
                    $k, remove, add, self.unique);
                if add != remove {
                    if self.counts[add as usize] == 0 {
                        self.unique += 1;
                    }
                    if self.counts[remove as usize] == 1 {
                        self.unique -= 1;
                    }
                    self.counts[add as usize] += 1;
                    self.counts[remove as usize] -= 1;
                }
                println!(" -> {}", self.unique);
            }

            fn add(&mut self, kmer: u64) {
                let kmer = kmer & $mask;
                self.countdown = self.countdown.saturating_sub(1);
                print!("    k{}:   Add  {:06b}, countdown {}, unique was {}",
                    $k, kmer, self.countdown, self.unique);
                if self.countdown == 0 {
                    if self.counts[kmer as usize] == 0 {
                        self.unique += 1;
                    }
                    self.counts[kmer as usize] += 1;
                }
                println!(" -> {}", self.unique);
            }

            fn remove(&mut self, kmer: u64) {
                debug_assert_eq!(self.countdown, 0);
                let kmer = kmer & $mask;
                print!("    k{}: Remove {:06b}, countdown {}, unique was {}",
                    $k, kmer, self.countdown, self.unique);
                if self.countdown == 0 {
                    self.counts[kmer as usize] -= 1;
                    if self.counts[kmer as usize] == 0 {
                        self.unique -= 1;
                    }
                }
                println!(" -> {}", self.unique);
            }
        }
    }
}

kmer_counts!(KmerCounts1, 1, 0b000011);
kmer_counts!(KmerCounts2, 2, 0b001111);
kmer_counts!(KmerCounts3, 3, 0b111111);

pub fn linguistic_complexity_k3(seq: &[u8], w: usize) -> Vec<f32> {
    const KMER_SIZE: usize = 3;
    let n = seq.len();
    assert!(n >= w && w >= KMER_SIZE && w <= u8::MAX as usize,
        "Cannot calculate sequence complexity for sequence length {} and window {}", n, w);
    assert!(w % 2 == 1, "Window size ({}) must be odd!", w);
    let halfw = w / 2;
    let wshift = w * 2;
    let mut res = vec![f32::NAN; n];
    let coef = 1.0_f32 / (min(4, w) * min(16, w - 1) * min(64, w - 2)) as f32;

    let mut counts1 = KmerCounts1::new();
    let mut counts2 = KmerCounts2::new();
    let mut counts3 = KmerCounts3::new();
    let mut kmer = 0_u64;
    let mut i = 0;
    let mut save_ix = w - 1;
    while i < (n - w + 1) {
        println!("i = {},  nt {}", i, seq[i] as char);
        let nt: u64 = match unsafe { *seq.get_unchecked(i) } {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => {
                i += KMER_SIZE;
                save_ix = i + w - 1;
                // Compiler will convert to memset.
                counts1.reset();
                counts2.reset();
                counts3.reset();
                continue;
            },
            ch => panic!("Unexpected nucleotide {}", ch as char),
        };

        if i > save_ix {
            let old_kmer = kmer >> wshift;
            counts1.remove(kmer);
            counts2.remove(kmer);
            counts3.remove(kmer);
        }
        kmer = (kmer << 2) | nt;
        counts1.add(kmer);
        counts2.add(kmer);
        counts3.add(kmer);

        if i >= save_ix {
            res[i + halfw] = f32::from(counts1.unique) * f32::from(counts2.unique) * f32::from(counts3.unique) * coef;
            println!("    save {} ({})   unique {}, {}, {} -> {:.4}", i,
                String::from_utf8_lossy(&seq[i + 1 - w..=i]),
                counts1.unique, counts2.unique, counts3.unique, res[i + halfw]);
        }
        i += 1;
    }
    res
}

// fn add_kmer(counts: &mut [u8], unique: &mut u8, add_kmer: u8) {
//     print!("    Update: add {:06b}, unique was {}", add_kmer, unique);
//     counts[add_kmer as usize] += 1;
//     if counts[add_kmer as usize] == 1 {
//         *unique += 1;
//     }
//     println!(" -> {}", unique);
// }

// fn update_kmer_counts(counts: &mut [u8], unique: &mut u8, remove_kmer: u8, add_kmer: u8) {
//     print!("    Update: remove {:06b}, add {:06b}, unique was {}", remove_kmer, add_kmer, unique);
//     if remove_kmer != add_kmer {
//         counts[remove_kmer as usize] -= 1;
//         if counts[remove_kmer as usize] == 0 {
//             *unique -= 1;
//         }
//         counts[add_kmer as usize] += 1;
//         if counts[add_kmer as usize] == 1 {
//             *unique += 1;
//         }
//     }
//     println!(" -> {}", unique);
// }

// pub fn linguistic_complexity_k3(seq: &[u8], w: usize) -> Vec<f32> {
//     let n = seq.len();
//     assert!(n >= w && w >= 3 && w <= 255,
//         "Cannot calculate sequence complexity for sequence length {} and window {}", n, w);
//     assert!(w % 2 == 1, "Window size ({}) must be odd!", w);
//     let halfw = w / 2;
//     let mut res = vec![f32::NAN; n];
//     let coef = 1.0_f32 / (min(4, w) * min(16, w - 1) * min(64, w - 2)) as f32;

//     let mut counts1 = [0_u8; 4];
//     let mut counts2 = [0_u8; 16];
//     let mut counts3 = [0_u8; 64];
//     let mut unique1 = 0_u8;
//     let mut unique2 = 0_u8;
//     let mut unique3 = 0_u8;

//     const KMER_SIZE: usize = 3;
//     const MASK3: u8 = 0b111111;
//     const MASK2: u8 = 0b001111;
//     const MASK1: u8 = 0b000011;
//     let mut kmer = 0_u8;
//     let mut i = 0;
//     let mut last_reset = 0;

//     while i < (n - w + 1) {
//         println!("i = {},  nt {}", i, seq[i] as char);
//         let nt: u8 = match unsafe { *seq.get_unchecked(i) } {
//             b'A' => 0,
//             b'C' => 1,
//             b'G' => 2,
//             b'T' => 3,
//             b'N' => {
//                 i += KMER_SIZE;
//                 last_reset = i;
//                 // Compiler will convert to memset.
//                 counts1.iter_mut().for_each(|x| *x = 0);
//                 counts2.iter_mut().for_each(|x| *x = 0);
//                 counts3.iter_mut().for_each(|x| *x = 0);
//                 unique1 = 0;
//                 unique2 = 0;
//                 unique3 = 0;
//                 continue;
//             },
//             ch => panic!("Unexpected nucleotide {}", ch as char),
//         };
//         let old_kmer = kmer >> 2;
//         kmer = (kmer << 2) | nt;

//         if i > last_reset {}

//         if i >= update_ix {
//             update_kmer_counts(&mut counts1, &mut unique1, old_kmer & MASK1, kmer & MASK1);
//             update_kmer_counts(&mut counts2, &mut unique2, old_kmer & MASK2, kmer & MASK2);
//             update_kmer_counts(&mut counts3, &mut unique3, old_kmer & MASK3, kmer & MASK3);
//             if i >= save_ix {
//                 res[i + halfw] = f32::from(unique1) * f32::from(unique2) * f32::from(unique3) * coef;
//                 println!("    save {} ({})   unique {}, {}, {} -> {:.4}", i,
//                     String::from_utf8_lossy(&seq[i + 1 - w..=i]),
//                     unique1, unique2, unique3, res[i + halfw]);
//             }
//         } else {
//             add_kmer(&mut counts1, &mut unique1, kmer & MASK1);
//             add_kmer(&mut counts2, &mut unique2, kmer & MASK2);
//             add_kmer(&mut counts3, &mut unique3, kmer & MASK3);
//         }
//         i += 1;
//     }
//     res
// }

const MAX_COMPLEXITY: u8 = 2;
const N_COMPLEXITIES: usize = (MAX_COMPLEXITY + 1) as usize;

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
    // complexity = min(abundance / DIVISOR, MAX_COMPLEXITY).
    // Value 61 would also map to 1 instead of 2, but abundance = 61 cannot be achieved.
    const DIVISOR: u8 = 31;

    let n = seq.len();
    let mut res = vec![0; n];
    // Last calculated index.
    let r = n - PADD - 1;
    for i in PADD..=r {
        let ab = kmer_abundance_123(&seq[i-PADD..=i+PADD]);
        res[i] = if ab <= 15 {
            0
        } else if ab <= 60 {
            1
        } else {
            2
        };
        // res[i] = min(kmer_abundance_123(&seq[i-PADD..=i+PADD]) as u8 / DIVISOR, MAX_COMPLEXITY);
    }
    for i in 0..PADD {
        res[i] = res[PADD];
    }
    for i in r + 1..n {
        res[i] = res[r];
    }
    res
}

/// For a specific (i) read mate, (ii) sequence complexity, and (iii) read position bin,
/// store the total number of matches, mismatches, deletions, and sum insertion size.
#[derive(Default, Clone)]
struct ErrorCounts {
    matches: u32,
    mismatches: u32,
    deletions: u32,
    sum_insertion: u32,
}

impl ErrorCounts {
    fn observations(&self) -> u32 {
        // Insertions always come after another operation, therefore it is not counted here.
        self.matches + self.mismatches + self.deletions
    }
}

impl fmt::Debug for ErrorCounts {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let obs = self.observations();
        write!(f, "{:5} obs.  M: {:5}, X: {:5}, D: {:5}, I: {:5} ({:.10})",
            obs, self.matches, self.mismatches, self.deletions,
            self.sum_insertion, self.sum_insertion as f32 / obs as f32)
    }
}

/// Place read position into one of the bins as
/// `bin <- floor(POSITION_BINS * read_pos / (length - 1))`.
const POSITION_BINS: u32 = 500;

/// For a specific (i) read mate and (ii) sequence complexity,
/// store error counts for each read position bin,
/// as well as different soft clipping sizes.
#[derive(Clone)]
struct SubProfileBuilder {
    /// Vector of error counts. Length = POSITION_BINS + 1.
    err_counts: Vec<ErrorCounts>,
    /// Number of observed soft clippings. Key: `(clipping as f64 / length as f64).to_bits()`.
    clip_counts: IntMap<u32>,
}

impl Default for SubProfileBuilder {
    fn default() -> Self {
        Self {
            err_counts: vec![ErrorCounts::default(); POSITION_BINS as usize + 1],
            clip_counts: IntMap::new(),
        }
    }
}

impl fmt::Debug for SubProfileBuilder {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut clipping: Vec<(f64, u32)> = self.clip_counts.iter()
            .map(|(&key, &count)| (f64::from_bits(key), count)).collect();
        clipping.sort_by(|a, b| a.0.total_cmp(&b.0));
        writeln!(f, "    Clipping: {:?}", clipping)?;
        for i in 0..=POSITION_BINS as usize {
            writeln!(f, "        {:3} {:?}", i, self.err_counts[i])?;
        }
        Ok(())
    }
}

use crate::seq::aln::READ_ENDS;

/// Construct error profiles.
struct ProfileBuilder {
    /// Reference sequence interval.
    interval: Interval,
    /// Sequence complexity around each nucleotide of the input interval.
    /// Stores values in 0..N_COMPLEXITIES.
    complexities: Vec<u8>,
    /// Stratify errors by read mate and sequence complexity.
    profiles: [[SubProfileBuilder; N_COMPLEXITIES]; READ_ENDS],
}

impl ProfileBuilder {
    fn new(interval: &Interval, interval_seq: &[u8]) -> Self {
        assert_eq!(interval.len(), interval_seq.len() as u32);
        Self {
            profiles: Default::default(),
            interval: interval.clone(),
            complexities: sevenmer_complexities(interval_seq),
        }
    }

    fn update(&mut self, aln: &Alignment, read_end: ReadEnd) {
        let aln_interval = aln.ref_interval();
        debug_assert!(aln_interval.overlaps(&self.interval), "Read alignment is out of the reference interval");
        let self_start = self.interval.start();
        let self_end = self.interval.end();
        let ref_start = aln_interval.start();
        let ref_end = aln_interval.end();
        let cigar = aln.cigar();
        let read_len = cigar.query_len();

        let (left_clip, right_clip) = cigar.soft_clipping();
        let profiles = &mut self.profiles[read_end.ix()];

        // Left clipping.
        if ref_start >= self_start {
            let left_clip_key = (left_clip as f64 / read_len as f64).to_bits();
            match profiles[self.complexities[(ref_start - self_start) as usize] as usize]
                    .clip_counts.entry(left_clip_key) {
                Entry::Occupied(mut entry) => *entry.get_mut() += 1,
                Entry::Vacant(entry) => { entry.insert(1); },
            }
        }
        // Right clipping.
        if ref_end <= self_end {
            let right_clip_key = (right_clip as f64 / read_len as f64).to_bits();
            match profiles[self.complexities[(ref_end - self_start - 1) as usize] as usize]
                    .clip_counts.entry(right_clip_key) {
                Entry::Occupied(mut entry) => *entry.get_mut() += 1,
                Entry::Vacant(entry) => { entry.insert(1); },
            }
        }

        let mut ref_pos = ref_start;
        let mut read_pos = 0;
        for item in cigar.iter() {
            let curr_len = item.len();
            match item.operation() {
                Operation::Soft => read_pos += curr_len,
                Operation::Match => panic!("Expected an extended CIGAR, but found M operation."),
                Operation::Equal => {
                    for i in max(ref_pos, self_start)..min(ref_pos + curr_len, self_end) {
                        profiles[self.complexities[(i - self_start) as usize] as usize]
                            .err_counts[(POSITION_BINS * (i - ref_pos + read_pos) / (read_len - 1)) as usize]
                            .matches += 1;
                    }
                    ref_pos += curr_len;
                    read_pos += curr_len;
                },
                Operation::Diff => {
                    for i in max(ref_pos, self_start)..min(ref_pos + curr_len, self_end) {
                        profiles[self.complexities[(i - self_start) as usize] as usize]
                            .err_counts[(POSITION_BINS * (i - ref_pos + read_pos) / (read_len - 1)) as usize]
                            .mismatches += 1;
                    }
                    ref_pos += curr_len;
                    read_pos += curr_len;
                },
                Operation::Ins => {
                    if self_start <= ref_pos && ref_pos < self_end {
                        profiles[self.complexities[(ref_pos - self_start) as usize] as usize]
                            .err_counts[(POSITION_BINS * (read_pos) / (read_len - 1)) as usize]
                            .sum_insertion += curr_len;
                    }
                    read_pos += curr_len;
                },
                Operation::Del => {
                    for i in max(ref_pos, self_start)..min(ref_pos + curr_len, self_end) {
                        profiles[self.complexities[(i - self_start) as usize] as usize]
                            .err_counts[(POSITION_BINS * read_pos / (read_len - 1)) as usize]
                            .deletions += 1;
                    }
                    ref_pos += curr_len;
                },
            }
        }
    }
}

impl fmt::Debug for ProfileBuilder {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Profile builder for {}", self.interval)?;
        for read_end in 0..READ_ENDS as usize {
            for complexity in 0..N_COMPLEXITIES as usize {
                writeln!(f, "Mate {}, complexity {}", read_end + 1, complexity)?;
                self.profiles[read_end][complexity].fmt(f)?;
            }
        }
        Ok(())
    }
}

pub struct ErrorProfile;

impl ErrorProfile {
    /// Create error profile from the iterator over records.
    pub fn estimate<'a, I>(records: I, interval: &Interval, interval_seq: &[u8]) -> ErrorProfile
    where I: Iterator<Item = &'a Record>,
    {
        log::info!("    Estimating read error profiles");
        let mut prof_builder = ProfileBuilder::new(interval, interval_seq);
        let contigs = interval.contigs();
        for record in records {
            // Unmapped or secondary.
            if (record.flags() & 3844) != 0 {
                continue;
            }
            let aln = Alignment::from_record(&record, Rc::clone(contigs));
            prof_builder.update(&aln, ReadEnd::from_record(&record));
        }
        println!("{:?}", prof_builder);
        ErrorProfile
    }
}
