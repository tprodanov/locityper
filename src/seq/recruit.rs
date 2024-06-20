//! Recruit reads from FASTQ files according to their k-mers.

use std::{
    io, thread,
    cmp::{min, max},
    collections::hash_map::Entry,
    time::{Instant, Duration},
    sync::mpsc::{self, Sender, Receiver, TryRecvError},
};
use nohash::IntMap;
use smallvec::{smallvec, SmallVec};
use crate::{
    err::{Error, validate_param, add_path},
    seq::{
        ContigSet,
        kmers::{self, Kmer, KmerCount},
        fastx::{self, FastxRead, UPDATE_SECS},
    },
    math::RoundDiv,
};

pub type Minimizer = u64;

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
    ) -> Result<Self, Error>
    {
        validate_param!(0 < minimizer_k && minimizer_k <= Minimizer::MAX_KMER_SIZE,
            "Minimizer kmer-size must be within [1, {}]", Minimizer::MAX_KMER_SIZE);
        validate_param!(1 < minimizer_w && minimizer_w <= kmers::MAX_MINIMIZER_W,
            "Minimizer window-size must be within [2, {}]", kmers::MAX_MINIMIZER_W);

        let min_frac = f64::from(SUBSUM_PENALTY) / f64::from(SUBSUM_BONUS + 1);
        validate_param!(match_frac >= min_frac && match_frac <= 1.0,
            "Minimizer match fraction ({}) must be in [{:.5}, 1]", match_frac, min_frac);
        if match_frac < 0.5 {
            log::warn!("Minimizer match fraction ({}) under 0.5 is allowed but not recommended", match_frac);
        }
        validate_param!(match_length >= 200 && match_length <= 100_000,
            "Matching stretch length ({}) should be between 200 and 100,000", match_length);
        if match_length < READ_LENGTH_THRESH {
            log::warn!("Matching stretch length ({}) is too small, it is not be used for short reads in any case",
            match_length);
        }

        let mut thresholds = [0; STORE_THRESHOLDS];
        for i in 0..STORE_THRESHOLDS {
            thresholds[i] = max(1, (i as f64 * match_frac).ceil() as u16);
        }

        // As per https://doi.org/10.1093/bioinformatics/btaa472,
        // there are 2L/(w + 1) minimizers per sequence of length of L.
        let stretch_minims = (2 * match_length).fast_ceil_div(u32::from(minimizer_w) + 1);
        let stretch_score = f64::from(stretch_minims)
            * (f64::from(SUBSUM_BONUS + SUBSUM_PENALTY) * match_frac - f64::from(SUBSUM_PENALTY));
        assert!(stretch_score > 0.0);
        let stretch_score = stretch_score.ceil() as u32;
        Ok(Self { minimizer_k, minimizer_w, match_frac, match_length, stretch_minims, stretch_score, thresholds })
    }

    /// Returns true if the read has enough matching minimizers to be recruited.
    #[inline]
    fn short_read_passes(&self, usable_matches: u16, unusable_matches: u16, total_minims: u16) -> bool {
        usable_matches >= self.thresholds[usize::from(total_minims - unusable_matches)]
    }

    #[inline]
    fn long_read_threshold(&self, n_minims: u32) -> u32 {
        max(1, (f64::from(n_minims) * self.match_frac).ceil() as u32)
    }
}

/// Recruitment statistics: how long did it take, how many reads processed and how many recruited.
struct Stats {
    timer: Instant,
    recruited: u64,
    processed: u64,
    /// Last log message was at this duration since start.
    last_msg: Duration,
}

impl Stats {
    fn new() -> Self {
        Self {
            timer: Instant::now(),
            recruited: 0,
            processed: 0,
            last_msg: Duration::default(),
        }
    }

    /// Prints log message if enough time has passed.
    fn timed_print_log(&mut self) {
        let elapsed = self.timer.elapsed();
        if (elapsed - self.last_msg).as_secs() >= UPDATE_SECS {
            self.print_log_always(elapsed);
        }
    }

    fn print_log_always(&mut self, elapsed: Duration) {
        let processed = self.processed as f64;
        let speed = 1e-3 * processed / elapsed.as_secs_f64();
        log::debug!("    Recruited {:11} /{:8.0}k reads, {:5.1}k reads/s", self.recruited, 1e-3 * processed, speed);
        self.last_msg = elapsed;
    }

    fn finish(&mut self) {
        let elapsed = self.timer.elapsed();
        self.print_log_always(elapsed);
        log::info!("Finished recruitment in {}", crate::ext::fmt::Duration(elapsed));
    }
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct ShortMatchCount {
    usable1: u16,
    unusable1: u16,
    usable2: u16,
    unusable2: u16,
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct LongMatchCount {
    usable: u32,
    unusable: u32,
}

/// During read recruitment, we count minimizer matches within 8 bytes,
/// but use them differently for short and long reads.
#[repr(C)]
#[derive(Clone, Copy)]
pub(crate) union MatchCount {
    short: ShortMatchCount,
    long: LongMatchCount,
}

const SHORT_ZERO_MATCHES: MatchCount = MatchCount {
    short: ShortMatchCount {
        usable1: 0,
        unusable1: 0,
        usable2: 0,
        unusable2: 0,
    },
};

const LONG_ZERO_MATCHES: MatchCount = MatchCount {
    long: LongMatchCount {
        usable: 0,
        unusable: 0,
    },
};

pub(crate) trait MatchesBuffer: Default + Sync + Send + 'static {
    /// How many loci can be processed with this buffer?
    const MAX_LOCI: usize;

    type Iter<'a>: Iterator<Item = (u16, MatchCount)>;

    fn clear(&mut self);

    fn get_or_insert(&mut self, locus_ix: u16, val: MatchCount) -> &mut MatchCount;

    fn get_mut(&mut self, locus_ix: u16) -> Option<&mut MatchCount>;

    fn iter(&self) -> Self::Iter<'_>;
}

pub(crate) type MatchesMap = IntMap<u16, MatchCount>;

pub(crate) struct MatchesMapIter<'a>(std::collections::hash_map::Iter<'a, u16, MatchCount>);

impl<'a> Iterator for MatchesMapIter<'a> {
    type Item = (u16, MatchCount);

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|(&locus_ix, &count)| (locus_ix, count))
    }
}

impl MatchesBuffer for MatchesMap {
    const MAX_LOCI: usize = u16::MAX as usize;

    type Iter<'a> = MatchesMapIter<'a>;

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
    fn iter(&self) -> MatchesMapIter<'_> {
        MatchesMapIter(self.iter())
    }
}

pub(crate) type SingleMatch = Option<MatchCount>;

pub(crate) struct SingleMatchIter<'a>(std::option::Iter<'a, MatchCount>);

impl<'a> Iterator for SingleMatchIter<'a> {
    type Item = (u16, MatchCount);

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|&count| (0, count))
    }
}

impl MatchesBuffer for SingleMatch {
    const MAX_LOCI: usize = 1;

    type Iter<'a> = SingleMatchIter<'a>;

    #[inline(always)]
    fn clear(&mut self) {
        *self = None;
    }

    #[inline(always)]
    fn get_or_insert(&mut self, locus_ix: u16, default: MatchCount) -> &mut MatchCount {
        assert!(locus_ix == 0);
        self.get_or_insert(default)
    }

    #[inline(always)]
    fn get_mut(&mut self, locus_ix: u16) -> Option<&mut MatchCount> {
        assert!(locus_ix == 0);
        self.as_mut()
    }

    #[inline(always)]
    fn iter(&self) -> Self::Iter<'_> {
        SingleMatchIter(self.iter())
    }
}

/// Trait-extension over single/paired reads.
pub(crate) trait RecruitableRecord : fastx::WritableRecord + Send + 'static {
    /// Recruit a short single-end/paired-end read, or a long read.
    /// `minimizers` and `matches` are buffer structures, that can be freely cleaned and reused.
    fn recruit<M: MatchesBuffer>(&self,
        targets: &Targets,
        answer: &mut Answer,
        minimizers: &mut Vec<Minimizer>,
        matches: &mut M,
    );
}

impl<T: fastx::SingleRecord + fastx::WritableRecord + Send + 'static> RecruitableRecord for T {
    #[inline]
    fn recruit<M: MatchesBuffer>(&self,
        targets: &Targets,
        answer: &mut Answer,
        minimizers: &mut Vec<Minimizer>,
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
        answer: &mut Answer,
        minimizers: &mut Vec<Minimizer>,
        matches: &mut M,
    ) {
        targets.recruit_read_pair(self[0].seq(), self[1].seq(), answer, minimizers, matches);
    }
}

const CAPACITY: usize = 4;
/// Key: minimizer,
/// value: vector of loci indices, where the minimizer appears + is k-mer usable (not common).
type MinimToLoci = IntMap<Minimizer, SmallVec<[(u16, bool); CAPACITY]>>;
/// Vector of loci indices, to which the read was recruited.
type Answer = Vec<u16>;

/// Target builder. Can be converted to targets using `finalize()`.
pub struct TargetBuilder {
    params: Params,
    thresh_kmer_count: KmerCount,
    total_seqs: u32,
    buffer: Vec<(u32, Minimizer)>,

    /// Map minimizer -> all locus indices, where the minimizer occurs.
    minim_to_loci: MinimToLoci,
    /// For each locus, map minimizer -> usable/unusable k-mer (entry is present only if k-mer appears in the locus).
    locus_minimizers: Vec<IntMap<Minimizer, bool>>,
}

impl TargetBuilder {
    /// Creates a new targets builder.
    /// Discard k-mers that appear off target more than `thresh_kmer_count` times.
    pub fn new(params: Params, thresh_kmer_count: KmerCount) -> Self {
        Self {
            params, thresh_kmer_count,
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
        let mut locus_minimizers = IntMap::default();
        for (seq, counts) in contig_set.seqs().iter().zip(kmer_counts.iter()) {
            too_short_alleles += u32::from(seq.len() < self.params.match_length as usize);
            let n_counts = counts.len();
            self.buffer.clear();
            kmers::canon_minimizers(seq, self.params.minimizer_k, self.params.minimizer_w, &mut self.buffer);
            for &(pos, minimizer) in self.buffer.iter() {
                let pos = pos as usize;
                let is_usable = if u32::from(self.params.minimizer_k) <= base_k {
                    // Check Jellyfish k-mer that is centered around k-mer at `pos`.
                    counts[min(pos.saturating_sub(shift), n_counts - 1)] < self.thresh_kmer_count
                } else {
                    // Compare first and last Jellyfish k-mers contained in the k-mer at `pos`.
                    counts[pos] < self.thresh_kmer_count && counts[pos + shift] < self.thresh_kmer_count
                };

                match self.minim_to_loci.entry(minimizer) {
                    Entry::Occupied(entry) => {
                        let loci_ixs = entry.into_mut();
                        let last = loci_ixs.last_mut().unwrap();
                        if last.0 == locus_ix {
                            last.1 &= is_usable; // All counts must be usable.
                        } else {
                            loci_ixs.push((locus_ix, is_usable));
                        }
                    }
                    Entry::Vacant(entry) => { entry.insert(smallvec![(locus_ix, is_usable)]); }
                }
                locus_minimizers.entry(minimizer)
                    .and_modify(|was_usable| *was_usable &= is_usable).or_insert(is_usable);
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
        Targets {
            params: self.params,
            minim_to_loci: self.minim_to_loci,
            locus_minimizers: self.locus_minimizers,
        }
    }
}

/// Recruitment targets.
#[derive(Clone)]
pub struct Targets {
    params: Params,
    /// Minimizers appearing across the targets.
    minim_to_loci: MinimToLoci,
    locus_minimizers: Vec<IntMap<Minimizer, bool>>,
}

impl Targets {
    /// Record short single-end read to one or more loci.
    fn recruit_short_read<M: MatchesBuffer>(
        &self,
        seq: &[u8],
        answer: &mut Answer,
        minimizers: &mut Vec<Minimizer>,
        matches: &mut M,
    ) {
        minimizers.clear();
        kmers::canon_minimizers(seq, self.params.minimizer_k, self.params.minimizer_w, minimizers);
        let total = u16::try_from(minimizers.len()).expect("Short read has too many minimizers");

        matches.clear();
        for minimizer in minimizers.iter() {
            let Some(loci_ixs) = self.minim_to_loci.get(minimizer) else { continue };
            for &(locus_ix, is_usable) in loci_ixs.iter() {
                let counts = matches.get_or_insert(locus_ix, SHORT_ZERO_MATCHES);
                let counts = unsafe { &mut counts.short };
                if is_usable { counts.usable1 += 1 } else { counts.unusable1 += 1 };
            }
        }

        answer.clear();
        for (locus_ix, counts) in matches.iter() {
            let counts = unsafe { counts.short };
            if self.params.short_read_passes(counts.usable1, counts.unusable1, total) {
                answer.push(locus_ix);
            }
        }
    }

    /// Compare long read minimizers against locus minimizers.
    /// Unusable matches `locus_minims[minim] = false` are ignored,
    /// Otherwise each mismatch is penalized by -1, each matches rewarded by +1.
    /// Read is recruited if there is a subarray with sum >= req_sum.
    fn has_matching_stretch(&self, minimizers: &[Minimizer], locus_ix: u16) -> bool {
        let locus_minimizers = &self.locus_minimizers[usize::from(locus_ix)];
        let mut s = 0;
        for minim in minimizers {
            // Optimized Kadane's algorithm for finding max subsum.
            match locus_minimizers.get(minim) {
                Some(false) => {}
                Some(true) => {
                    s += SUBSUM_BONUS;
                    if s >= self.params.stretch_score {
                        return true;
                    }
                }
                None => s = s.saturating_sub(SUBSUM_PENALTY),
            }
        }
        false
    }

    /// Record long single-end read to one or more loci.
    fn recruit_long_read<M: MatchesBuffer>(
        &self,
        seq: &[u8],
        answer: &mut Answer,
        minimizers: &mut Vec<Minimizer>,
        matches: &mut M,
    ) {
        minimizers.clear();
        kmers::canon_minimizers(seq, self.params.minimizer_k, self.params.minimizer_w, minimizers);
        let total_minims = u32::try_from(minimizers.len()).expect("Long read has too many minimizers");

        matches.clear();
        for minimizer in minimizers.iter() {
            let Some(loci_ixs) = self.minim_to_loci.get(minimizer) else { continue };
            for &(locus_ix, is_usable) in loci_ixs.iter() {
                let counts = matches.get_or_insert(locus_ix, LONG_ZERO_MATCHES);
                let counts = unsafe { &mut counts.long };
                if is_usable { counts.usable += 1 } else { counts.unusable += 1 };
            }
        }

        answer.clear();
        for (locus_ix, counts) in matches.iter() {
            let counts = unsafe { counts.long };
            let total_wo_unusable = total_minims - counts.unusable;
            if counts.usable >= self.params.long_read_threshold(min(self.params.stretch_minims, total_wo_unusable)) {
                if total_wo_unusable < self.params.stretch_minims || self.has_matching_stretch(minimizers, locus_ix) {
                    answer.push(locus_ix);
                }
            }
        }
    }

    /// Record one paired-end read to one or more loci.
    /// The read is recruited when both read mates satisfy the recruitment thresholods.
    fn recruit_read_pair<M: MatchesBuffer>(
        &self,
        seq1: &[u8],
        seq2: &[u8],
        answer: &mut Answer,
        minimizers: &mut Vec<Minimizer>,
        matches: &mut M,
    ) {
        // Matches have 8-byte values. Reinterpret them as four u16s: (usable1, unusable1, usable2, unusable2).
        matches.clear();
        // First mate.
        minimizers.clear();
        kmers::canon_minimizers(seq1, self.params.minimizer_k, self.params.minimizer_w, minimizers);
        let total1 = u16::try_from(minimizers.len()).expect("Paired end read has too many minimizers");
        for minimizer in minimizers.iter() {
            let Some(loci_ixs) = self.minim_to_loci.get(minimizer) else { continue };
            for &(locus_ix, is_usable) in loci_ixs.iter() {
                let counts = matches.get_or_insert(locus_ix, SHORT_ZERO_MATCHES);
                let counts = unsafe { &mut counts.short };
                if is_usable { counts.usable1 += 1 } else { counts.unusable1 += 1 };
            }
        }

        // Second mate.
        minimizers.clear();
        kmers::canon_minimizers(seq2, self.params.minimizer_k, self.params.minimizer_w, minimizers);
        let total2 = u16::try_from(minimizers.len()).expect("Paired end read has too many minimizers");
        for minimizer in minimizers.iter() {
            let Some(loci_ixs) = self.minim_to_loci.get(minimizer) else { continue };
            for &(locus_ix, is_usable) in loci_ixs.iter() {
                // No reason to insert new loci if they did not match the first read end.
                if let Some(counts) = matches.get_mut(locus_ix) {
                    let counts = unsafe { &mut counts.short };
                    if is_usable { counts.usable2 += 1 } else { counts.unusable2 += 1 };
                }
            }
        }

        answer.clear();
        for (locus_ix, counts) in matches.iter() {
            let counts = unsafe { counts.short };
            if self.params.short_read_passes(counts.usable1, counts.unusable1, total1)
                    && self.params.short_read_passes(counts.usable2, counts.unusable2, total2) {
                answer.push(locus_ix);
            }
        }
    }

    fn recruit_single_thread<T, M>(
        &self,
        mut reader: impl FastxRead<Record = T>,
        writers: &mut [impl io::Write],
    ) -> Result<(), Error>
    where T: RecruitableRecord,
          M: MatchesBuffer,
    {
        let mut record = T::default();
        let mut answer = Answer::new();
        let mut buffer1 = Default::default();
        let mut buffer2 = M::default();

        let mut stats = Stats::new();
        while reader.read_next(&mut record)? {
            record.recruit(self, &mut answer, &mut buffer1, &mut buffer2);
            for &locus_ix in answer.iter() {
                record.write_to(&mut writers[usize::from(locus_ix)]).map_err(add_path!(!))?;
            }
            stats.recruited += u64::from(!answer.is_empty());
            stats.processed += 1;
            if stats.processed % 10000 == 0 {
                stats.timed_print_log();
            }
        }
        stats.finish();
        Ok(())
    }

    fn recruit_multi_thread<T, M>(
        &self,
        reader: impl FastxRead<Record = T>,
        writers: Vec<impl io::Write>,
        threads: u16,
        chunk_size: usize,
    ) -> Result<(), Error>
    where T: RecruitableRecord,
          M: MatchesBuffer,
    {
        let n_workers = usize::from(threads - 1);
        log::info!("Starting read recruitment with 1 read/write thread and {} recruitment threads", n_workers);
        let mut main_worker = MainWorker::<T, _, _>::new::<M>(self, reader, writers, n_workers, chunk_size);
        main_worker.run()?;
        main_worker.finish()
    }

    /// How many targets are there in the set?
    pub fn n_targets(&self) -> usize {
        self.locus_minimizers.len()
    }

    /// Recruit reads to the targets, possibly in multiple threads.
    pub(crate) fn recruit(
        &self,
        reader: impl FastxRead<Record = impl RecruitableRecord>,
        mut writers: Vec<impl io::Write>,
        threads: u16,
        chunk_size: usize,
    ) -> Result<(), Error>
    {
        assert_eq!(writers.len(), self.locus_minimizers.len(), "Unexpected number of writers");
        let single_target = self.n_targets() == 1;
        if threads <= 1 {
            if single_target {
                self.recruit_single_thread::<_, SingleMatch>(reader, &mut writers)
            } else {
                self.recruit_single_thread::<_, MatchesMap>(reader, &mut writers)
            }
        } else {
            if single_target {
                self.recruit_multi_thread::<_, SingleMatch>(reader, writers, threads, chunk_size)
            } else {
                self.recruit_multi_thread::<_, MatchesMap>(reader, writers, threads, chunk_size)
            }
        }
    }
}

/// Vector of records and corresponding answers.
/// Is send between threads: main thread reads records and sends the vector to workers.
/// In the meantime, workers receive records and fills corresponding answers (recruitment targets for the record),
/// and then send shipments back to the main thread.
type Shipment<T> = Vec<(T, Answer)>;

struct Worker<T, M> {
    targets: Targets,
    buffer1: Vec<Minimizer>,
    buffer2: M,
    /// Receives records that need to be recruited.
    receiver: Receiver<Shipment<T>>,
    /// Sends already recruited reads back to the main thread.
    sender: Sender<Shipment<T>>,
}

impl<T, M: Default> Worker<T, M> {
    fn new(targets: Targets, receiver: Receiver<Shipment<T>>, sender: Sender<Shipment<T>>) -> Self {
        Self {
            buffer1: Default::default(),
            buffer2: Default::default(),
            targets, receiver, sender,
        }
    }
}

impl<T: RecruitableRecord, M: MatchesBuffer> Worker<T, M> {
    fn run(mut self) {
        // Block thread and wait for the shipment.
        while let Ok(mut shipment) = self.receiver.recv() {
            assert!(!shipment.is_empty());
            for (record, answer) in shipment.iter_mut() {
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
struct MainWorker<T, R: FastxRead<Record = T>, W> {
    /// Fasta/q reader.
    /// Becomes `None`, once the stream has ended.
    reader: Option<R>,
    /// Fasta/q writers for each of the loci.
    writers: Vec<W>,
    /// Senders from the main thread to the workers. Sends reads to be analyzed.
    senders: Vec<Sender<Shipment<T>>>,
    /// Receivers from workers to the main thread. Receives
    /// - analyzed reads with possible recruited loci,
    /// - bool: true if any of the reads were recruited.
    receivers: Vec<Receiver<Shipment<T>>>,
    /// Thread handles.
    handles: Vec<thread::JoinHandle<()>>,

    /// Send reads between threads in vectors of this size.
    chunk_size: usize,
    /// Does the worker currently recruit reads?
    is_busy: Vec<bool>,
    /// Recruitment statistics.
    stats: Stats,

    /// Chunks of reads that were read from the reader and are ready to be analyzed.
    to_send: Vec<Shipment<T>>,
    /// Chunks of reads that were analyzed and need to be writter to the writers.
    to_write: Vec<Shipment<T>>,
}

impl<T, R, W> MainWorker<T, R, W>
where T: RecruitableRecord,
      R: FastxRead<Record = T>,
      W: io::Write,
{
    fn new<M>(targets: &Targets, reader: R, writers: Vec<W>, n_workers: usize, chunk_size: usize) -> Self
    where M: MatchesBuffer,
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
            reader: Some(reader),
            is_busy: vec![false; n_workers],
            stats: Stats::new(),
            to_send: Vec::new(),
            to_write: Vec::new(),
        }
    }

    /// Starts the process: provides the first task to each worker.
    fn start(&mut self) -> Result<(), Error> {
        for (is_busy, sender) in self.is_busy.iter_mut().zip(&self.senders) {
            let shipment = read_new_shipment(&mut self.reader, self.chunk_size)?;
            if !shipment.is_empty() {
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

    fn write_read_iteration(&mut self) -> Result<(), Error> {
        if self.reader.is_none() {
            return Ok(())
        }
        while let Some(mut shipment) = self.to_write.pop() {
            write_shipment(&mut self.writers, &shipment, &mut self.stats)?;
            self.stats.timed_print_log();
            fill_shipment(&mut self.reader, &mut shipment)?;
            if !shipment.is_empty() {
                self.to_send.push(shipment);
            }
            if self.reader.is_none() {
                return Ok(());
            }
        }
        if self.to_send.is_empty() {
            let shipment = read_new_shipment(&mut self.reader, self.chunk_size)?;
            if !shipment.is_empty() {
                self.to_send.push(shipment);
            }
        }
        Ok(())
    }

    /// Main part of the multi-thread recruitment.
    /// Iterates until there are any shipments left to read from the input files.
    fn run(&mut self) -> Result<(), Error> {
        self.start()?;
        while self.reader.is_some() || !self.to_send.is_empty() {
            self.write_read_iteration()?;
            // There were no updates, and there are shipments ready to be sent.
            if !self.recv_send_iteration() && !self.to_send.is_empty() {
                const SLEEP: Duration = Duration::from_micros(100);
                thread::sleep(SLEEP);
            }
        }
        Ok(())
    }

    /// Finish the main thread: write all remaining shipments to the output files, and stop worker threads.
    fn finish(mut self) -> Result<(), Error> {
        assert!(self.reader.is_none() && self.to_send.is_empty());
        for shipment in self.to_write.into_iter() {
            write_shipment(&mut self.writers, &shipment, &mut self.stats)?;
        }
        for (&is_busy, receiver) in self.is_busy.iter().zip(&self.receivers) {
            if is_busy {
                // Block thread and wait for the task completion.
                let shipment = receiver.recv().expect("Recruitment worker has failed!");
                write_shipment(&mut self.writers, &shipment, &mut self.stats)?;
                self.stats.timed_print_log();
            }
        }
        std::mem::drop(self.senders);
        for handle in self.handles.into_iter() {
            handle.join().expect("Process failed for unknown reason");
        }
        self.stats.finish();
        Ok(())
    }
}

/// Fills `shipment` from the reader.
/// Output shipment may be empty, if the stream has ended.
fn fill_shipment<T, R>(opt_reader: &mut Option<R>, shipment: &mut Shipment<T>) -> Result<(), Error>
where T: Default,
      R: FastxRead<Record = T>,
{
    let reader = opt_reader.as_mut().expect("fill_shipment: reader must not be None");
    let mut new_len = 0;
    for (record, _) in shipment.iter_mut() {
        if reader.read_next(record)? {
            new_len += 1;
        } else {
            shipment.truncate(new_len);
            *opt_reader = None;
            break;
        }
    }
    Ok(())
}

fn read_new_shipment<T, R>(opt_reader: &mut Option<R>, chunk_size: usize) -> Result<Shipment<T>, Error>
where T: Clone + Default,
      R: FastxRead<Record = T>,
{
    let mut shipment = vec![Default::default(); chunk_size];
    fill_shipment(opt_reader, &mut shipment)?;
    Ok(shipment)
}

/// Writes recruited records to the output files.
fn write_shipment<T>(writers: &mut [impl io::Write], shipment: &Shipment<T>, stats: &mut Stats) -> Result<(), Error>
where T: fastx::WritableRecord,
{
    stats.processed += shipment.len() as u64;
    for (record, answer) in shipment.iter() {
        stats.recruited += u64::from(!answer.is_empty());
        for &locus_ix in answer.iter() {
            record.write_to(&mut writers[usize::from(locus_ix)]).map_err(add_path!(!))?;
        }
    }
    Ok(())
}
