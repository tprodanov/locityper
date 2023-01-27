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
        res[i] = min(kmer_abundance_123(&seq[i-PADD..=i+PADD]) as u8 / DIVISOR, MAX_COMPLEXITY);
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

use crate::seq::aln::READ_ENDS;

/// Construct error profiles.
struct ProfileBuilder {
    /// Stratify errors by read mate and sequence complexity.
    profiles: [[SubProfileBuilder; N_COMPLEXITIES]; READ_ENDS],
    /// Reference sequence interval.
    interval: Interval,
    /// Sequence complexity around each nucleotide of the input interval.
    /// Stores values in 0..N_COMPLEXITIES.
    complexities: Vec<u8>,
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

pub struct ErrorProfile;

impl ErrorProfile {
    /// Create error profile from the iterator over records.
    pub fn estimate<'a, I>(
        records: I,
        contigs: &Rc<ContigNames>,
        interval: &Interval,
        interval_seq: &[u8],
    ) -> ErrorProfile
    where I: Iterator<Item = &'a Record>,
    {
        log::info!("    Estimating read error profiles");
        let mut prof_builder = ProfileBuilder::new(interval, interval_seq);
        for record in records {
            // Unmapped or secondary.
            if (record.flags() & 3844) != 0 {
                continue;
            }
            let aln = Alignment::from_record(&record, Rc::clone(contigs));
            prof_builder.update(&aln, ReadEnd::from_record(&record));
        }
        ErrorProfile
    }
}
