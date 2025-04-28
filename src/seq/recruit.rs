//! Recruit reads from FASTQ files according to their k-mers.

use std::{
    io, thread,
    fmt::{self, Display},
    cmp::{min, max, PartialOrd},
    time::{Instant, Duration},
    sync::mpsc::{self, Sender, Receiver, TryRecvError},
    ops::{Add, Sub, AddAssign},
};
use smallvec::SmallVec;
use rand::Rng;
use crate::{
    err::{validate_param, add_path},
    seq::{
        ContigSet,
        kmers::{self, Kmer},
        counts::KmerCount,
        fastx::{self, FastxRead},
    },
    math::RoundDiv,
    algo::{IntMap, hash_map},
    ext::rand::XoshiroRng,
};

pub type Minimizer = u64;

/// Default minimizer size.
pub const DEFAULT_MINIM_KW: (u8, u8) = (15, 10);
/// Default length, over which a long read must match the target
/// (shorter matches with more minimizer matches will pass as well).
pub const DEFAULT_MATCH_LEN: u32 = 2000;

/// Reads shorter than 500 bp are considered short, larger - long.
/// Note that there must be < 256 minimizers per read end.
const READ_LENGTH_THRESH: u32 = 500;

/// Store this many minimizer thresholds for short reads.
const STORE_THRESHOLDS: usize = READ_LENGTH_THRESH as usize + Minimizer::MAX_KMER_SIZE as usize;

/// Reward matching k-mers by +3, penalize mismatching k-mers by -1.
/// Bigger bonus requires bigger score to achieve `(BONUS + 1) * match_frac - PENALTY`,
/// and therefore promotes longer stretches over short stretches without any mismatches.
const SUBSUM_BONUS: u32 = 3;
const SUBSUM_PENALTY: u32 = 1;

#[derive(Clone)]
pub struct Params {
    /// Recruit reads using k-mers of size `minimizer_k` that have the minimial hash across `minimizer_w`
    /// consecutive k-mers.
    minimizer_k: u8,
    minimizer_w: u8,
    /// Recruit the read if it has at has at least `match_frac` matching k-mers
    /// on a stretch of `min(match_length, read_length)`.
    match_frac: f64,
    match_length: u32,

    /// Use k-mers with multiplicity LESS than this number.
    thresh_kmer_count: KmerCount,

    /// Approximate number of minimizers in a stretch of `match_length`.
    stretch_minims: u32,
    /// If each matching k-mer has bonus +1 and each non-matching is penalized by -1,
    /// stretch of `stretch_minims` k-mers would have a sum of `strech_score`.
    stretch_score: u32,
    /// Threshold number of matches of short reads.
    thresholds: [u16; STORE_THRESHOLDS],
}

impl Params {
    pub fn new(
        (minimizer_k, minimizer_w): (u8, u8),
        match_frac: f64,
        match_length: u32,
        thresh_kmer_count: KmerCount,
    ) -> crate::Result<Self>
    {
        validate_param!(0 < minimizer_k && minimizer_k <= Minimizer::MAX_KMER_SIZE,
            "Minimizer kmer-size must be within [1, {}]", Minimizer::MAX_KMER_SIZE);
        validate_param!(1 < minimizer_w && minimizer_w <= kmers::MAX_MINIMIZER_W,
            "Minimizer window-size must be within [2, {}]", kmers::MAX_MINIMIZER_W);

        let min_frac = f64::from(SUBSUM_PENALTY) / f64::from(SUBSUM_BONUS + 1);
        validate_param!(match_frac >= min_frac && match_frac <= 1.0,
            "Minimizer match fraction ({}) must be in [{:.5}, 1]", match_frac, min_frac);
        if match_frac < 0.3 {
            log::warn!("Small minimizer match fraction ({}) is allowed but not recommended", match_frac);
        }
        validate_param!(match_length >= 200 && match_length <= 100_000,
            "Matching stretch length ({}) should be between 200 and 100,000", match_length);
        if match_length < READ_LENGTH_THRESH {
            log::warn!("Matching stretch length ({}) is too small (note that it is only used for long reads)",
                match_length);
        }
        validate_param!(thresh_kmer_count > 0, "k-mer threshold must be positive");

        let mut thresholds = [0; STORE_THRESHOLDS];
        for i in 0..STORE_THRESHOLDS {
            thresholds[i] = max(1, (i as f64 * match_frac).ceil() as u16);
        }

        // As per https://doi.org/10.1093/bioinformatics/btaa472,
        // there are 2L/(w + 1) minimizers per sequence of length of L.
        let stretch_minims = (2 * match_length).fast_ceil_div(u32::from(minimizer_w) + 1);
        let stretch_score = f64::from(stretch_minims)
            * (f64::from(SUBSUM_BONUS + SUBSUM_PENALTY) * match_frac - f64::from(SUBSUM_PENALTY));
        let stretch_score = stretch_score.max(f64::from(SUBSUM_BONUS)).ceil() as u32;
        Ok(Self {
            minimizer_k, minimizer_w, match_frac, match_length, thresh_kmer_count,
            stretch_minims, stretch_score, thresholds,
        })
    }

    #[inline(always)]
    fn long_read_threshold(&self, n_minims: u32) -> u32 {
        max(1, (f64::from(min(self.stretch_minims, n_minims)) * self.match_frac).ceil() as u32)
    }

    /// k-mer size within the minimizers.
    pub fn minimizer_k(&self) -> u8 {
        self.minimizer_k
    }

    /// Window size for the minimizers.
    #[allow(unused)]
    pub fn minimizer_w(&self) -> u8 {
        self.minimizer_w
    }
}

// Write log messages at most every ten seconds.
pub const UPDATE_SECS: u64 = 10;

// Check timer every `UPDATE_COUNT` processed reads.
const UPDATE_COUNT: u64 = 10_000;

/// Recruitment statistics: how long did it take, how many reads processed and how many recruited.
pub struct Progress {
    timer: Instant,
    recruited: u64,
    processed: u64,
    /// Last log message was at this duration since start.
    last_msg: Duration,
    /// Show recruited reads, not only processed.
    show_recruited: bool,
}

impl Progress {
    pub fn new(show_recruited: bool) -> Self {
        Self {
            timer: Instant::now(),
            recruited: 0,
            processed: 0,
            last_msg: Duration::default(),
            show_recruited,
        }
    }

    /// Create new progress logger, where the number of recruited reads is shown.
    #[allow(unused)]
    pub fn new_recruitment() -> Self {
        Self::new(true)
    }

    /// Create new progress logger, where the number of recruited reads is not shown.
    pub fn new_simple() -> Self {
        Self::new(false)
    }

    /// Add the number of recruited reads.
    /// Must be run before processed reads are added.
    #[inline(always)]
    pub fn add_recruited(&mut self, count: impl Into<u64>) {
        self.recruited += count.into();
    }

    /// Increment the number of processed reads by 1.
    /// Should be run after the number of recruited reads is updated.
    #[inline]
    pub fn inc_processed(&mut self) {
        self.processed += 1;
        if self.processed % UPDATE_COUNT == 0 {
            self.timed_print_message();
        }
    }

    /// Add the number of processed reads.
    /// Should be run after the number of recruited reads is updated.
    #[inline]
    pub fn add_processed(&mut self, count: u64) {
        let old_val = self.processed;
        self.processed += count;
        if old_val / UPDATE_COUNT < self.processed / UPDATE_COUNT {
            self.timed_print_message();
        }
    }

    /// Prints log message if enough time has passed.
    fn timed_print_message(&mut self) {
        let elapsed = self.timer.elapsed();
        if (elapsed - self.last_msg).as_secs() >= UPDATE_SECS {
            let processed = self.processed as f64;
            let speed = 1e-3 * processed / elapsed.as_secs_f64();
            if self.show_recruited {
                log::debug!("    Recruited {:11} /{:8.0}k reads, {:5.1}k reads/s",
                    self.recruited, 1e-3 * processed, speed);
            } else {
                log::debug!("    Processed {:8.0}k reads, {:5.1}k reads/s", 1e-3 * processed, speed);
            }
            self.last_msg = elapsed;
        }
    }

    /// Print final message.
    pub fn final_message(&mut self) {
        let elapsed = self.timer.elapsed();
        let processed = self.processed as f64;
        let speed = 1e-3 * processed / elapsed.as_secs_f64();
        let total_time = crate::ext::fmt::Duration(elapsed);
        if self.show_recruited {
            log::debug!("    Recruited {} / {} reads ({:.4}%) in {} ({:5.1}k reads/s)",
                self.recruited, self.processed, 100.0 * self.recruited as f64 / processed,
                total_time, speed);
        } else {
            log::debug!("    Processed {} reads in {} ({:5.1}k reads/s)", self.processed, total_time, speed);
        }
    }

    /// Total number of recruited reads.
    #[inline]
    #[allow(unused)]
    pub fn recruited(&self) -> u64 {
        self.recruited
    }

    /// Total number of processed reads.
    #[inline]
    pub fn processed(&self) -> u64 {
        self.processed
    }
}

#[repr(C)]
#[derive(Default, Clone, Copy, Debug)]
struct BaseMatchCount<T> {
    /// Four values are: common-backward, common-forward, rare-backward rare-forward.
    arr: [T; 4],
}

impl<T: num_traits::ConstZero + Copy> BaseMatchCount<T> {
    const ZERO: Self = Self {
        arr: [T::ZERO; 4],
    };
}

impl<T: AddAssign + From<bool>> BaseMatchCount<T> {
    #[inline(always)]
    fn inc(&mut self, forward: bool, info: MinimInfo) {
        let i = u8::from(info.rare) << 1;
        // Updating common/rare - backward.
        self.arr[usize::from(i)].add_assign(T::from(info.is_directed_to(!forward)));
        self.arr[usize::from(i | 1)].add_assign(T::from(info.is_directed_to(forward)));
    }
}

impl<T: Copy + Add<Output = T> + Sub<Output = T> + PartialOrd> BaseMatchCount<T> {
    /// Given the total number of minimizer matches, returns
    /// - direction of the match (forward if true),
    /// - divisor and numerators of the match fraction.
    #[inline(always)]
    fn characterize(self, total: T) -> (bool, T, T) {
        let [bw_c, fw_c, bw_r, fw_r] = self.arr;
        if fw_c + fw_r >= bw_c + bw_r {
            (true, fw_r, total - fw_c)
        } else {
            (false, bw_r, total - bw_c)
        }
    }
}

impl<T: Copy + Display> Display for BaseMatchCount<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let [bw_c, fw_c, bw_r, fw_r] = self.arr;
        write!(f, "rare: {}, {}; common: {}, {}", fw_r, bw_r, fw_c, bw_c)
    }
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct ShortMatchCount {
    first: BaseMatchCount<u16>,
    second: BaseMatchCount<u16>,
}

// #[repr(C)]
// #[derive(Clone, Copy, Default)]
// struct LongMatchCount(GenericMatchCount<u32>);

/// During read recruitment, we count minimizer matches within 8 bytes,
/// but use them differently for short and long reads.
#[repr(C)]
#[derive(Clone, Copy)]
pub(crate) union MatchCount {
    short: ShortMatchCount,
    long: BaseMatchCount<u32>,
}

const SHORT_ZERO_MATCHES: MatchCount = MatchCount {
    short: ShortMatchCount {
        first: BaseMatchCount::<u16>::ZERO,
        second: BaseMatchCount::<u16>::ZERO,
    },
};

const LONG_ZERO_MATCHES: MatchCount = MatchCount {
    long: BaseMatchCount::<u32>::ZERO,
};

/// To which loci is the read recruited?
pub(crate) trait Answer: Default + Clone + Sync + Send + 'static {
    type Iter<'a>: Iterator<Item = u16>;

    fn clear(&mut self);

    fn push(&mut self, locus_ix: u16);

    fn not_empty(&self) -> bool;

    /// Iterate over recruited loci.
    fn iter(&self) -> Self::Iter<'_>;
}

impl Answer for Vec<u16> {
   type Iter<'a> = std::iter::Copied<std::slice::Iter<'a, u16>>;

    #[inline(always)]
    fn clear(&mut self) {
        self.clear();
    }

    #[inline(always)]
    fn push(&mut self, locus_ix: u16) {
        self.push(locus_ix);
    }

    #[inline(always)]
    fn not_empty(&self) -> bool {
        !self.is_empty()
    }

    #[inline(always)]
    fn iter(&self) -> Self::Iter<'_> {
        (self as &[u16]).iter().copied()
    }
}

/// Store answer for a single locus.
impl Answer for bool {
    type Iter<'a> = std::option::IntoIter<u16>;

    #[inline(always)]
    fn clear(&mut self) {
        *self = false;
    }

    #[inline(always)]
    fn push(&mut self, locus_ix: u16) {
        debug_assert!(locus_ix == 0);
        *self = true;
    }

    #[inline(always)]
    fn not_empty(&self) -> bool {
        *self
    }

    #[inline(always)]
    fn iter(&self) -> Self::Iter<'_> {
        (if *self { Some(0) } else { None }).into_iter()
    }
}

/// Based on the number of targets, use different buffer structures.
pub(crate) trait MatchesBuffer: Default + Sync + Send + 'static {
    type Iter<'a>: Iterator<Item = (u16, MatchCount)>;

    type Answer: Answer;

    fn clear(&mut self);

    fn get_or_insert(&mut self, locus_ix: u16, val: MatchCount) -> &mut MatchCount;

    fn get_mut(&mut self, locus_ix: u16) -> Option<&mut MatchCount>;

    fn is_empty(&self) -> bool;

    fn iter(&self) -> Self::Iter<'_>;
}

/// When there exist more than one target, use IntMap.
pub(crate) type MatchesMap = IntMap<u16, MatchCount>;

pub(crate) struct MatchesMapIter<'a>(hash_map::Iter<'a, u16, MatchCount>);

impl<'a> Iterator for MatchesMapIter<'a> {
    type Item = (u16, MatchCount);

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|(&locus_ix, &count)| (locus_ix, count))
    }
}

impl MatchesBuffer for MatchesMap {
    type Iter<'a> = MatchesMapIter<'a>;

    type Answer = Vec<u16>;

    #[inline(always)]
    fn clear(&mut self) {
        IntMap::clear(self);
    }

    #[inline(always)]
    fn get_or_insert(&mut self, locus_ix: u16, default: MatchCount) -> &mut MatchCount {
        self.entry(locus_ix).or_insert(default)
    }

    #[inline(always)]
    fn get_mut(&mut self, locus_ix: u16) -> Option<&mut MatchCount> {
        self.get_mut(&locus_ix)
    }

    #[inline(always)]
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    #[inline(always)]
    fn iter(&self) -> MatchesMapIter<'_> {
        MatchesMapIter(self.iter())
    }
}

/// When there exist exactly one target, use simple `Option`.
pub(crate) type SingleMatch = Option<MatchCount>;

impl MatchesBuffer for SingleMatch {
    type Iter<'a> = std::option::IntoIter<(u16, MatchCount)>;

    type Answer = bool;

    #[inline(always)]
    fn clear(&mut self) {
        *self = None;
    }

    #[inline(always)]
    fn get_or_insert(&mut self, locus_ix: u16, default: MatchCount) -> &mut MatchCount {
        debug_assert!(locus_ix == 0);
        self.get_or_insert(default)
    }

    #[inline(always)]
    fn get_mut(&mut self, locus_ix: u16) -> Option<&mut MatchCount> {
        debug_assert!(locus_ix == 0);
        self.as_mut()
    }

    #[inline(always)]
    fn is_empty(&self) -> bool {
        self.is_none()
    }

    #[inline(always)]
    fn iter(&self) -> Self::Iter<'_> {
        self.map(|count| (0, count)).into_iter()
    }
}

/// Trait-extension over single/paired reads.
pub(crate) trait RecruitableRecord : fastx::WritableRecord + Send + 'static {
    /// Recruit a short single-end/paired-end read, or a long read.
    /// `minimizers` and `matches` are buffer structures, that can be freely cleaned and reused.
    fn recruit<M: MatchesBuffer>(&self,
        targets: &Targets,
        answer: &mut M::Answer,
        minimizers: &mut Vec<(Minimizer, bool)>,
        matches: &mut M,
    );
}

impl<T: fastx::SingleRecord + fastx::WritableRecord + Send + 'static> RecruitableRecord for T {
    #[inline]
    fn recruit<M: MatchesBuffer>(&self,
        targets: &Targets,
        answer: &mut M::Answer,
        minimizers: &mut Vec<(Minimizer, bool)>,
        matches: &mut M,
    ) {
        let seq = self.seq();
        if seq.len() as u32 <= READ_LENGTH_THRESH {
            targets.recruit_short_read(seq, answer, minimizers, matches);
        } else {
            targets.recruit_long_read(seq, answer, minimizers, matches);
        }
    }
}

impl<T: fastx::SingleRecord + fastx::WritableRecord + Send + 'static> RecruitableRecord for [T; 2] {
    #[inline]
    fn recruit<M: MatchesBuffer>(&self,
        targets: &Targets,
        answer: &mut M::Answer,
        minimizers: &mut Vec<(Minimizer, bool)>,
        matches: &mut M,
    ) {
        println!("{} {}", std::str::from_utf8(self[0].name()).unwrap(),
            crate::model::locs::NameHash::new(self[0].name()));
        targets.recruit_read_pair(self[0].seq(), self[1].seq(), answer, minimizers, matches);
    }
}

#[derive(Clone, Copy, Debug)]
struct MinimInfo {
    // 0 = 00 - none, 1 = 01 - backward, 2 = 10 - forward, 3 = 11 - both.
    direction: u8,
    rare: bool,
}

impl Default for MinimInfo {
    #[inline(always)]
    fn default() -> Self {
        Self {
            direction: 0,
            rare: true,
        }
    }
}

impl MinimInfo {
    #[inline]
    fn new(forward: bool, rare: bool) -> Self {
        Self {
            direction: 1 + u8::from(forward),
            rare,
        }
    }

    #[inline]
    fn update(&mut self, forward: bool, rare: bool) {
        self.direction |= 1 + u8::from(forward);
        self.rare &= rare;
    }

    #[inline]
    fn is_directed_to(self, forward: bool) -> bool {
        self.direction & (1 + u8::from(forward)) != 0
    }
}

impl Display for MinimInfo {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}{}", if self.rare { 'r' } else { 'c' },
            if self.is_directed_to(false) { "<" } else { "" },
            if self.is_directed_to(true) { ">" } else { "" },
        )
    }
}

const CAPACITY: usize = 4;
/// Key: minimizer,
/// value: vector of (locus index, minimizer information).
type MinimToLoci = IntMap<Minimizer, SmallVec<[(u16, MinimInfo); CAPACITY]>>;
type LocusMinims = Vec<IntMap<Minimizer, MinimInfo>>;

/// Target builder. Can be converted to targets using `finalize()`.
pub struct TargetBuilder {
    params: Params,
    total_seqs: u32,
    buffer: Vec<(u32, Minimizer, bool)>,

    /// Map minimizer -> all locus indices, where the minimizer occurs.
    minim_to_loci: MinimToLoci,
    /// For each locus, map minimizer -> minimizer info without locus index.
    locus_minimizers: LocusMinims,
}

impl TargetBuilder {
    /// Creates a new targets builder.
    pub fn new(params: Params) -> Self {
        Self {
            params,
            total_seqs: 0,
            buffer: Vec::new(),
            minim_to_loci: Default::default(),
            locus_minimizers: Vec::new(),
        }
    }

    /// Add set of locus alleles.
    pub fn add(&mut self, contig_set: &ContigSet, mean_read_len: f64) {
        let locus_ix = u16::try_from(self.locus_minimizers.len())
            .expect(const_format::formatcp!("Too many contig sets (allowed at most {})", u16::MAX));
        let kmer_counts = contig_set.kmer_counts();
        let base_k = kmer_counts.k();
        let shift = if u32::from(self.params.minimizer_k) <= base_k {
            (base_k as usize - usize::from(self.params.minimizer_k)) / 2
        } else {
            usize::from(self.params.minimizer_k) - base_k as usize
        };

        let mut too_short_alleles = 0;
        let mut locus_minimizers = IntMap::<Minimizer, MinimInfo>::default();
        for (seq, counts) in contig_set.seqs().iter().zip(kmer_counts.iter()) {
            assert_eq!((seq.len() + 1).saturating_sub(base_k as usize), counts.len(),
                "Sequence and k-mer lengths do not match");
            too_short_alleles += u32::from(seq.len() < self.params.match_length as usize);
            let n_counts = counts.len();
            self.buffer.clear();
            kmers::canon_minimizers(seq, self.params.minimizer_k, self.params.minimizer_w, &mut self.buffer);
            for &(pos, minimizer, forward) in self.buffer.iter() {
                let pos = pos as usize;
                let rare = if u32::from(self.params.minimizer_k) <= base_k {
                    // Check Jellyfish k-mer that is centered around k-mer at `pos`.
                    counts[min(pos.saturating_sub(shift), n_counts - 1)] < self.params.thresh_kmer_count
                } else {
                    // Compare first and last Jellyfish k-mers contained in the k-mer at `pos`.
                    counts[pos] < self.params.thresh_kmer_count && counts[pos + shift] < self.params.thresh_kmer_count
                };

                let v = self.minim_to_loci.entry(minimizer).or_default();
                match v.last_mut() {
                    Some((last_locus, last_info)) if *last_locus == locus_ix => last_info.update(forward, rare),
                    _ => v.push((locus_ix, MinimInfo::new(forward, rare))),
                }
                locus_minimizers.entry(minimizer).or_default().update(forward, rare);
            }
        }

        self.total_seqs += contig_set.len() as u32;
        self.locus_minimizers.push(locus_minimizers);
        if mean_read_len >= f64::from(self.params.match_length) && too_short_alleles > 0 {
            log::warn!("{}: {}/{} alleles are shorter than match length ({}), not all reads will be recruited",
                contig_set.tag(), too_short_alleles, contig_set.len(), self.params.match_length);
        }
    }

    /// Finalize targets construction.
    /// Remove top `discard_minim` fraction of minimizers.
    pub fn finalize(self) -> Targets {
        let total_minims = self.minim_to_loci.len();
        log::info!("Collected {} minimizers across {} loci and {} sequences",
            total_minims, self.locus_minimizers.len(), self.total_seqs);
        assert!(total_minims > 0, "No minimizers for recruitment");
        println!("Targets:");
        for (minimizer, entries) in self.minim_to_loci.iter() {
            println!("    {:016X} {}", minimizer, entries[0].1);
        }
        println!("Recruit:");
        Targets {
            params: self.params,
            minim_to_loci: self.minim_to_loci,
            locus_minimizers: self.locus_minimizers,
            show_recruited: true,
        }
    }
}

/// Trait for either subsampling, or using all reads.
pub trait Sampling {
    fn next(&mut self) -> bool;
}

pub struct All;

impl Sampling for All {
    #[inline(always)]
    fn next(&mut self) -> bool { true }
}

pub struct Subsampling {
    /// Probability of success, relative to the maximal integer.
    /// Taken from `rand::distributions::Bernoulli`.
    p_int: u64,
    rng: XoshiroRng,
}

impl Subsampling {
    pub fn new(p: f64, rng: XoshiroRng) -> Self {
        const SCALE: f64 = 2.0 * (1u64 << 63) as f64;
        assert!(0.0 < p && p < 1.0, "Subsampling rate ({}) must be in (0, 1)", p);
        Self {
            p_int: (p * SCALE) as u64,
            rng,
        }
    }
}

impl From<(f64, XoshiroRng)> for Subsampling {
    fn from((p, rng): (f64, XoshiroRng)) -> Self {
        Self::new(p, rng)
    }
}

impl Sampling for Subsampling {
    #[inline(always)]
    fn next(&mut self) -> bool {
        self.rng.random::<u64>() < self.p_int
    }
}

// /// Additional recruitment parameters.
// pub struct SupplParams {
//     pub threads: u8,
//     pub chunk_size: usize,
//     /// Maximal number of recruited reads.
//     pub max_recruited: u64,
// }

/// Recruitment targets.
#[derive(Clone)]
pub struct Targets {
    params: Params,
    /// Minimizers appearing across the targets.
    minim_to_loci: MinimToLoci,
    locus_minimizers: LocusMinims,

    show_recruited: bool,
}

impl Targets {
    /// During progress messages, do not show recruited reads, only processed (relevant for preprocessing).
    pub fn hide_recruited_reads(&mut self) {
        self.show_recruited = false;
    }

    /// Record short single-end read to one or more loci.
    fn recruit_short_read<M: MatchesBuffer>(
        &self,
        seq: &[u8],
        answer: &mut M::Answer,
        buffer: &mut Vec<(Minimizer, bool)>,
        matches: &mut M,
    ) {
        buffer.clear();
        kmers::canon_minimizers(seq, self.params.minimizer_k, self.params.minimizer_w, buffer);
        let total = u16::try_from(buffer.len()).expect("Short read has too many minimizers");

        matches.clear();
        for &(minimizer, forward) in buffer.iter() {
            let Some(entries) = self.minim_to_loci.get(&minimizer) else { continue };
            for &(locus_ix, info) in entries {
                let counts = matches.get_or_insert(locus_ix, SHORT_ZERO_MATCHES);
                let counts = unsafe { &mut counts.short };
                // Both backwards means that read and haplotype face the same direction.
                counts.first.inc(forward, info);
            }
        }

        answer.clear();
        for (locus_ix, counts) in matches.iter() {
            let counts = unsafe { counts.short };
            let (_, divisor, numerator) = counts.first.characterize(total);
            if divisor >= self.params.thresholds[usize::from(numerator)] {
                answer.push(locus_ix);
            }
        }
    }

    /// Record one read pair to one or more loci.
    /// The read is recruited when both read mates satisfy the recruitment thresholods.
    fn recruit_read_pair<M: MatchesBuffer>(
        &self,
        seq1: &[u8],
        seq2: &[u8],
        answer: &mut M::Answer,
        buffer: &mut Vec<(Minimizer, bool)>,
        matches: &mut M,
    ) {
        use std::fmt::Write;
        let mut s = String::new();

        answer.clear();
        matches.clear();
        // First mate.
        buffer.clear();
        kmers::canon_minimizers(seq1, self.params.minimizer_k, self.params.minimizer_w, buffer);
        let total1 = u16::try_from(buffer.len()).expect("Paired end read has too many minimizers");
        for &(minimizer, forward) in buffer.iter() {
            if !s.is_empty() {
                write!(s, ", ").unwrap();
            }
            write!(s, "{:016X}{}", minimizer, if forward { '>' } else { '<' }).unwrap();
            let Some(entries) = self.minim_to_loci.get(&minimizer) else { continue };
            write!(s, "!").unwrap();
            for &(locus_ix, info) in entries {
                let counts = matches.get_or_insert(locus_ix, SHORT_ZERO_MATCHES);
                let counts = unsafe { &mut counts.short };
                counts.first.inc(forward, info);
            }
        }
        if matches.is_empty() { return }
        println!("    M {}", s);
        s.clear();

        // Second mate.
        buffer.clear();
        kmers::canon_minimizers(seq2, self.params.minimizer_k, self.params.minimizer_w, buffer);
        let total2 = u16::try_from(buffer.len()).expect("Paired end read has too many minimizers");
        for &(minimizer, forward) in buffer.iter() {
            if !s.is_empty() {
                write!(s, ", ").unwrap();
            }
            write!(s, "{:016X}{}", minimizer, if forward { '>' } else { '<' }).unwrap();
            let Some(entries) = self.minim_to_loci.get(&minimizer) else { continue };
            write!(s, "!").unwrap();
            for &(locus_ix, info) in entries {
                // No reason to insert new loci if they did not match the first read end.
                let Some(counts) = matches.get_mut(locus_ix) else { continue };
                let counts = unsafe { &mut counts.short };
                counts.second.inc(forward, info);
            }
        }

        for (locus_ix, counts) in matches.iter() {
            let counts = unsafe { counts.short };
            println!("    {} / {}", counts.first, total1);
            println!("    M {}", s);
            println!("    {} / {}", counts.second, total2);
            let (fw1, divisor1, numerator1) = counts.first.characterize(total1);
            let (fw2, divisor2, numerator2) = counts.second.characterize(total2);
            if fw1 != fw2
                && divisor1 >= self.params.thresholds[usize::from(numerator1)]
                && divisor2 >= self.params.thresholds[usize::from(numerator2)]
            {
                println!("    OK");
                answer.push(locus_ix);
            }
        }
    }

    /// Compare long read minimizers against locus minimizers.
    /// Unusable matches `locus_minims[minim] = false` are ignored,
    /// Otherwise each mismatch is penalized by -1, each matches rewarded by +1.
    /// Read is recruited if there is a subarray with sum >= req_sum.
    fn has_matching_stretch(
        &self,
        minimizers: &[(Minimizer, bool)],
        locus_ix: u16,
        direction: bool,
    ) -> bool
    {
        let locus_minimizers = &self.locus_minimizers[usize::from(locus_ix)];
        let mut s = 0_u32;
        for &(minimizer, forward) in minimizers {
            // Optimized Kadane's algorithm for finding max subsum.
            match locus_minimizers.get(&minimizer) {
                None => s = s.saturating_sub(SUBSUM_PENALTY),
                Some(info) if info.rare => {
                    if info.is_directed_to(forward == direction) {
                        s += SUBSUM_BONUS;
                        if s >= self.params.stretch_score {
                            return true;
                        }
                    } else {
                        s = s.saturating_sub(SUBSUM_PENALTY);
                    }
                }
                _ => {} // minimizer present, but common.
            }
        }
        false
    }

    /// Record long single-end read to one or more loci.
    fn recruit_long_read<M: MatchesBuffer>(
        &self,
        seq: &[u8],
        answer: &mut M::Answer,
        buffer: &mut Vec<(Minimizer, bool)>,
        matches: &mut M,
    ) {
        buffer.clear();
        kmers::canon_minimizers(seq, self.params.minimizer_k, self.params.minimizer_w, buffer);
        let total = u32::try_from(buffer.len()).expect("Long read has too many minimizers");

        matches.clear();
        for &(minimizer, forward) in buffer.iter() {
            let Some(entries) = self.minim_to_loci.get(&minimizer) else { continue };
            for &(locus_ix, info) in entries {
                let counts = matches.get_or_insert(locus_ix, LONG_ZERO_MATCHES);
                let counts = unsafe { &mut counts.long };
                counts.inc(forward, info);
            }
        }

        answer.clear();
        for (locus_ix, counts) in matches.iter() {
            let counts = unsafe { counts.long };
            let (direction, divisor, numerator) = counts.characterize(total);
            if divisor >= self.params.long_read_threshold(numerator) && (numerator < self.params.stretch_minims
                || self.has_matching_stretch(buffer, locus_ix, direction))
            {
                answer.push(locus_ix);
            }
        }
    }

    fn recruit_single_thread<T, M>(
        &self,
        mut reader: impl FastxRead<Record = T>,
        mut writers: Vec<impl io::Write>,
        mut sampling: impl Sampling,
    ) -> crate::Result<Progress>
    where T: RecruitableRecord,
          M: MatchesBuffer,
    {
        let mut record = T::default();
        let mut answer = M::Answer::default();
        let mut buffer1 = Default::default();
        let mut buffer2 = M::default();

        let mut progress = Progress::new(self.show_recruited);
        while reader.read_next(&mut record)? {
            if sampling.next() {
                record.recruit(self, &mut answer, &mut buffer1, &mut buffer2);
                for locus_ix in answer.iter() {
                    record.write_to(&mut writers[usize::from(locus_ix)]).map_err(add_path!(!))?;
                }
                progress.add_recruited(answer.not_empty());
            }
            progress.inc_processed();
        }
        progress.final_message();
        Ok(progress)
    }

    fn recruit_multi_thread<T, M>(
        &self,
        reader: impl FastxRead<Record = T>,
        writers: Vec<impl io::Write>,
        threads: u16,
        chunk_size: usize,
        sampling: impl Sampling,
    ) -> crate::Result<Progress>
    where T: RecruitableRecord,
          M: MatchesBuffer,
    {
        let n_workers = usize::from(threads - 1);
        log::info!("Starting read recruitment with 1 read/write thread and {} recruitment threads", n_workers);
        let mut main_worker = MainWorker::<T, _, _, M>::new(self, reader, writers, n_workers, chunk_size);
        main_worker.run(sampling)?;
        main_worker.finish()
    }

    /// Recruit reads to the targets, possibly in multiple threads.
    /// Returns the number of recruited reads and the total number of reads.
    pub(crate) fn recruit(
        &self,
        reader: impl FastxRead<Record = impl RecruitableRecord>,
        writers: Vec<impl io::Write>,
        threads: u16,
        chunk_size: usize,
        subsampling: Option<(f64, XoshiroRng)>,
    ) -> crate::Result<Progress>
    {
        assert_eq!(writers.len(), self.locus_minimizers.len(), "Unexpected number of writers");
        assert!(threads >= 1 && self.n_targets() >= 1);

        // As read recruitment is very expensive, split into 8 implementations
        // based on the input parameters.
        match (threads, self.n_targets(), subsampling) {
            (1, 1, None) =>
                self.recruit_single_thread::<_, SingleMatch>(reader, writers, All),
            (1, 1, Some(t)) =>
                self.recruit_single_thread::<_, SingleMatch>(reader, writers, Subsampling::from(t)),
            (1, _, None) =>
                self.recruit_single_thread::<_, MatchesMap>(reader, writers, All),
            (1, _, Some(t)) =>
                self.recruit_single_thread::<_, MatchesMap>(reader, writers, Subsampling::from(t)),
            (_, 1, None) =>
                self.recruit_multi_thread::<_, SingleMatch>(reader, writers, threads, chunk_size, All),
            (_, 1, Some(t)) =>
                self.recruit_multi_thread::<_, SingleMatch>(reader, writers, threads, chunk_size, Subsampling::from(t)),
            (_, _, None) =>
                self.recruit_multi_thread::<_, MatchesMap>(reader, writers, threads, chunk_size, All),
            (_, _, Some(t)) =>
                self.recruit_multi_thread::<_, MatchesMap>(reader, writers, threads, chunk_size, Subsampling::from(t)),
        }
    }

    /// How many targets are there in the set?
    pub fn n_targets(&self) -> usize {
        self.locus_minimizers.len()
    }
}

/// Vector of records and corresponding answers.
/// Is send between threads: main thread reads records and sends the vector to workers.
/// In the meantime, workers receive records and fills corresponding answers (recruitment targets for the record),
/// and then send shipments back to the main thread.
#[derive(Clone)]
struct Shipment<T, A> {
    data: Vec<(T, A)>,
    total: u32,
}

struct Worker<T, M: MatchesBuffer> {
    targets: Targets,
    /// Receives records that need to be recruited.
    receiver: Receiver<Shipment<T, M::Answer>>,
    /// Sends already recruited reads back to the main thread.
    sender: Sender<Shipment<T, M::Answer>>,
    buffer1: Vec<(Minimizer, bool)>,
    buffer2: M,
}

impl<T, M: MatchesBuffer> Worker<T, M> {
    fn new(
        targets: Targets,
        receiver: Receiver<Shipment<T, M::Answer>>,
        sender: Sender<Shipment<T, M::Answer>>,
    ) -> Self
    {
        Self {
            targets, receiver, sender,
            buffer1: Default::default(),
            buffer2: Default::default(),
        }
    }
}

impl<T: RecruitableRecord, M: MatchesBuffer> Worker<T, M> {
    fn run(mut self) {
        // Block thread and wait for the shipment.
        while let Ok(mut shipment) = self.receiver.recv() {
            debug_assert!(!shipment.data.is_empty());
            for (record, answer) in shipment.data.iter_mut() {
                record.recruit(&self.targets, answer, &mut self.buffer1, &mut self.buffer2);
            }
            if let Err(_) = self.sender.send(shipment) {
                log::error!("Read recruitment: main thread stopped before the child thread.");
                return;
            }
        }
    }
}

/// Worker in the main thread, that organizes other workers, as well as reads/writes reads.
struct MainWorker<T, R, W, M>
where R: FastxRead<Record = T>,
      M: MatchesBuffer,
{
    /// Fasta/q reader.
    /// Becomes `None`, once the stream has ended.
    reader: Option<R>,
    /// Fasta/q writers for each of the loci.
    writers: Vec<W>,
    /// Senders from the main thread to the workers. Sends reads to be analyzed.
    senders: Vec<Sender<Shipment<T, M::Answer>>>,
    /// Receivers from workers to the main thread. Receives
    /// - analyzed reads with possible recruited loci,
    /// - bool: true if any of the reads were recruited.
    receivers: Vec<Receiver<Shipment<T, M::Answer>>>,
    /// Thread handles.
    handles: Vec<thread::JoinHandle<()>>,

    /// Send reads between threads in vectors of this size.
    chunk_size: usize,
    /// Does the worker currently recruit reads?
    is_busy: Vec<bool>,
    /// Recruitment statistics.
    progress: Progress,

    /// Chunks of reads that were read from the reader and are ready to be analyzed.
    to_send: Vec<Shipment<T, M::Answer>>,
    /// Chunks of reads that were analyzed and need to be writter to the writers.
    to_write: Vec<Shipment<T, M::Answer>>,
}

impl<T, R, W, M> MainWorker<T, R, W, M>
where T: RecruitableRecord,
      R: FastxRead<Record = T>,
      W: io::Write,
      M: MatchesBuffer
{
    fn new(
        targets: &Targets,
        reader: R,
        writers: Vec<W>,
        n_workers: usize,
        chunk_size: usize,
    ) -> Self
    {
        let mut senders = Vec::with_capacity(n_workers);
        let mut receivers = Vec::with_capacity(n_workers);
        let mut handles = Vec::with_capacity(n_workers);

        for _ in 0..n_workers {
            let (sender1, receiver1) = mpsc::channel();
            let (sender2, receiver2) = mpsc::channel();
            let worker = Worker::<T, M>::new(targets.clone(), receiver1, sender2);
            senders.push(sender1);
            receivers.push(receiver2);
            handles.push(thread::spawn(|| worker.run()));
        }
        Self {
            writers, senders, receivers, handles, chunk_size,
            progress: Progress::new(targets.show_recruited),
            reader: Some(reader),
            is_busy: vec![false; n_workers],
            to_send: Vec::new(),
            to_write: Vec::new(),
        }
    }

    /// Starts the process: provides the first task to each worker.
    fn start(&mut self, sampling: &mut impl Sampling) -> crate::Result<()> {
        for (is_busy, sender) in self.is_busy.iter_mut().zip(&self.senders) {
            let shipment = read_new_shipment(&mut self.reader, self.chunk_size, sampling)?;
            if !shipment.data.is_empty() {
                *is_busy = true;
                sender.send(shipment).expect("Recruitment worker has failed!");
            }
            if self.reader.is_none() {
                break;
            }
        }
        Ok(())
    }

    fn recv_send_iteration(&mut self) -> bool {
        let mut any_action = false;
        for ((receiver, sender), is_busy) in self.receivers.iter().zip(&self.senders)
                .zip(self.is_busy.iter_mut()) {
            if *is_busy {
                match receiver.try_recv() {
                    Ok(recv_shipment) => {
                        any_action = true;
                        self.to_write.push(recv_shipment);
                        if let Some(send_shipment) = self.to_send.pop() {
                            sender.send(send_shipment).expect("Recruitment worker has failed!");
                        } else {
                            *is_busy = false;
                        }
                    }
                    Err(TryRecvError::Empty) => { continue; }
                    Err(TryRecvError::Disconnected) => panic!("Recruitment worker has failed!"),
                }
            } else if let Some(send_shipment) = self.to_send.pop() {
                any_action = true;
                sender.send(send_shipment).expect("Recruitment worker has failed!");
                *is_busy = true;
            }
        }
        any_action
    }

    fn write_read_iteration(&mut self, sampling: &mut impl Sampling) -> crate::Result<()> {
        if self.reader.is_none() {
            return Ok(())
        }
        while let Some(mut shipment) = self.to_write.pop() {
            write_shipment(&mut self.writers, &shipment, &mut self.progress)?;
            fill_shipment(&mut self.reader, &mut shipment, sampling)?;
            if !shipment.data.is_empty() {
                self.to_send.push(shipment);
            }
            if self.reader.is_none() {
                return Ok(());
            }
        }
        if self.to_send.is_empty() {
            let shipment = read_new_shipment(&mut self.reader, self.chunk_size, sampling)?;
            if !shipment.data.is_empty() {
                self.to_send.push(shipment);
            }
        }
        Ok(())
    }

    /// Main part of the multi-thread recruitment.
    /// Iterates until there are any shipments left to read from the input files.
    fn run(&mut self, mut sampling: impl Sampling) -> crate::Result<()> {
        self.start(&mut sampling)?;
        while self.reader.is_some() || !self.to_send.is_empty() {
            self.write_read_iteration(&mut sampling)?;
            // There were no updates, and there are shipments ready to be sent.
            if !self.recv_send_iteration() && !self.to_send.is_empty() {
                const SLEEP: Duration = Duration::from_micros(100);
                thread::sleep(SLEEP);
            }
        }
        Ok(())
    }

    /// Finish the main thread: write all remaining shipments to the output files, and stop worker threads.
    fn finish(mut self) -> crate::Result<Progress> {
        assert!(self.reader.is_none() && self.to_send.is_empty());
        for shipment in self.to_write.into_iter() {
            write_shipment(&mut self.writers, &shipment, &mut self.progress)?;
        }
        for (&is_busy, receiver) in self.is_busy.iter().zip(&self.receivers) {
            if is_busy {
                // Block thread and wait for the task completion.
                let shipment = receiver.recv().expect("Recruitment worker has failed!");
                write_shipment(&mut self.writers, &shipment, &mut self.progress)?;
            }
        }
        std::mem::drop(self.senders);
        for handle in self.handles.into_iter() {
            handle.join().expect("Process failed for unknown reason");
        }
        self.progress.final_message();
        Ok(self.progress)
    }
}

/// Fills `shipment` from the reader.
/// Output shipment may be empty, if the stream has ended.
fn fill_shipment<T, R, A>(
    opt_reader: &mut Option<R>,
    shipment: &mut Shipment<T, A>,
    sampling: &mut impl Sampling,
) -> crate::Result<()>
where R: FastxRead<Record = T>,
{
    let reader = opt_reader.as_mut().expect("fill_shipment: reader must not be None");
    let mut total = 0;
    let mut new_len = 0;

    'outer: for (record, _) in shipment.data.iter_mut() {
        while reader.read_next(record)? {
            total += 1;
            if sampling.next() {
                new_len += 1;
                continue 'outer;
            }
        }
        // Reader has ended.
        shipment.data.truncate(new_len);
        *opt_reader = None;
        break;
    }
    shipment.total = total;
    Ok(())
}

#[inline]
fn read_new_shipment<T, R, A>(
    opt_reader: &mut Option<R>,
    chunk_size: usize,
    sampling: &mut impl Sampling,
) -> crate::Result<Shipment<T, A>>
where T: Clone + Default,
      A: Clone + Default,
      R: FastxRead<Record = T>,
{
    let mut shipment = Shipment {
        data: vec![Default::default(); chunk_size],
        total: 0,
    };
    fill_shipment(opt_reader, &mut shipment, sampling)?;
    Ok(shipment)
}

/// Writes recruited records to the output files.
fn write_shipment<T, A>(
    writers: &mut [impl io::Write],
    shipment: &Shipment<T, A>,
    progress: &mut Progress,
) -> crate::Result<()>
where T: fastx::WritableRecord,
      A: Answer,
{
    for (record, answer) in shipment.data.iter() {
        progress.add_recruited(answer.not_empty());
        for locus_ix in answer.iter() {
            record.write_to(&mut writers[usize::from(locus_ix)]).map_err(add_path!(!))?;
        }
    }
    progress.add_processed(u64::from(shipment.total));
    Ok(())
}
