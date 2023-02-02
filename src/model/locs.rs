use std::{
    rc::Rc,
    fmt,
    mem,
};
use htslib::bam::Record;
use intmap::IntMap;
use crate::{
    seq::{
        contigs::{ContigId, ContigNames},
        interv::Interval,
        aln::{Strand, ReadEnd, Alignment},
    },
    bg::{
        err_prof::ErrorProfile,
        insertsz::InsertDistr,
    },
    algo::{
        hash::fnv1a,
        bisect,
    },
    math::Ln,
};

// TODO: REMOVE LATER.
const DEBUG_ALNS: bool = false;

/// Single mate alignment: store alignment location, strand, read-end, and alignment ln-probability.
#[derive(Clone)]
struct MateAln {
    interval: Interval,
    strand: Strand,
    read_end: ReadEnd,
    /// Log-probability of the alignment.
    ln_prob: f64,
}

impl MateAln {
    /// Creates a new alignment extension from a htslib `Record`.
    fn from_record(record: &Record, contigs: Rc<ContigNames>, err_prof: &ErrorProfile) -> Self {
        let aln = Alignment::from_record(record, contigs);
        let ln_prob = err_prof.ln_prob(aln.cigar());
        let read_end = ReadEnd::from_record(record);
        if DEBUG_ALNS {
            log::debug!("        {}  {}  {:.3}", aln, read_end, ln_prob);
        }
        Self {
            strand: aln.strand(),
            interval: aln.take_interval(),
            read_end, ln_prob,
        }
    }

    /// Get contig id of the alignment.
    fn contig_id(&self) -> ContigId {
        self.interval.contig_id()
    }

    /// Returns sorted key: first sort by contig id, then by read end.
    fn sort_key(&self) -> u16 {
        (self.interval.contig_id().get() << 1) | (self.read_end.as_u16())
    }
}

impl fmt::Debug for MateAln {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "MateAln[{:?}, {:?}, prob={:.2}]", self.read_end, self.interval, self.ln_prob)
    }
}

impl fmt::Display for MateAln {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "MateAln[{}, {}, prob={:.2}]", self.read_end, self.interval, self.ln_prob)
    }
}

/// Preliminary, unpaired, read alignments.
/// Keys: read name hash, values: all alignments for the read pair.
pub struct PrelimAlignments {
    contigs: Rc<ContigNames>,
    alns: IntMap<Vec<MateAln>>,
}

impl PrelimAlignments {
    /// Creates `PrelimAlignments` from records.
    /// Only store read pairs where at least one alignment (for at least one read),
    /// has probability over `min_ln_prob`.
    pub fn from_records<I>(
        records: I,
        contigs: Rc<ContigNames>,
        err_prof: &ErrorProfile,
        min_ln_prob: f64,
    ) -> Self
    where I: Iterator<Item = Record>,
    {
        log::info!("[{}] Reading read alignments, alignment ln-probability threshold = {:.1}",
            contigs.tag(), min_ln_prob);
        assert!(min_ln_prob.is_finite(), "PrelimAlignments: min_ln_prob must be finite!");
        let mut alns = IntMap::new();
        let mut curr_hash = 0;
        let mut curr_alns = Vec::new();
        let mut best_prob = f64::NEG_INFINITY;

        for record in records {
            if record.is_unmapped() {
                continue;
            }

            let hash = fnv1a(record.qname());
            if hash != curr_hash {
                if best_prob >= min_ln_prob {
                    if DEBUG_ALNS {
                        log::debug!("    Saving {} alignments for hash {}  (prob {:.1})",
                            curr_alns.len(), curr_hash, best_prob);
                    }
                    // Assert must not be removed, or converted into debug_assert!
                    assert!(alns.insert(curr_hash, mem::replace(&mut curr_alns, Vec::new())).is_none(),
                        "Read pair alignments are separated by another read pair (see hash {})", curr_hash);
                } else {
                    if DEBUG_ALNS {
                        log::debug!("    Ignoring {} alignments for hash {}  (prob {:.1})",
                            curr_alns.len(), curr_hash, best_prob);
                    }
                    curr_alns.clear();
                }
                if DEBUG_ALNS {
                    log::debug!("    Read {},  hash {}:", String::from_utf8_lossy(record.qname()), hash);
                }
                curr_hash = hash;
                best_prob = f64::NEG_INFINITY;
            }
            let mate_aln = MateAln::from_record(&record, Rc::clone(&contigs), err_prof);
            best_prob = best_prob.max(mate_aln.ln_prob);
            curr_alns.push(mate_aln);
        }

        if best_prob >= min_ln_prob {
            if DEBUG_ALNS {
                log::debug!("    Saving {} alignments for hash {}  (prob {:.1})",
                    curr_alns.len(), curr_hash, best_prob);
            }
            // Assert must not be removed, or converted into debug_assert!
            assert!(alns.insert(curr_hash, curr_alns).is_none(),
                "Read pair alignments are separated by another read pair (see hash {})", curr_hash);
        } else if DEBUG_ALNS {
            log::debug!("    Ignoring {} alignments for hash {}  (prob {:.1})",
                curr_alns.len(), curr_hash, best_prob);
        }

        Self { contigs, alns }
    }

    /// Find paired-end alignment locations for each contig.
    /// Parameters:
    /// - unmapped_penalty for a single read mate (ln-space),
    /// - prob_diff: store alignment locations with ln-probability >= best_prob - prob_diff (ln-space).
    pub fn identify_locations<D>(&mut self, insert_distr: &D, unmapped_penalty: f64, prob_diff: f64)
        -> AllPairAlignments
    where D: InsertDistr
    {
        assert!(prob_diff >= 0.0, "Probability difference cannot be negative!");
        let n_reads = self.alns.len();
        let ln_ncontigs = (self.contigs.len() as f64).ln();
        log::info!("Identify paired alignment location and probabilities ({} read pairs)", n_reads);
        let mut res = AllPairAlignments::with_capacity(n_reads);
        for (&name_hash, alns) in self.alns.iter_mut() {
            // log::debug!("Read {}", name_hash);
            // Sort alignments first by contig id, then by read-end.
            alns.sort_by_key(MateAln::sort_key);
            let pair_alns = identify_pair_alignments(name_hash, alns, insert_distr,
                unmapped_penalty, prob_diff, ln_ncontigs);
            res.push(pair_alns);
        }
        res
    }
}

// Assume that there are at most 4 alignments of the to the same contig.
// (not critical, if this assumption fails).
const BISECT_RIGHT_STEP: usize = 4;

/// For a single read-pair, find all paired-read alignments to the same contig.
fn extend_pair_alignments<D>(
        new_alns: &mut Vec<PairAlignment>,
        buffer: &mut Vec<f64>,
        alns1: &[MateAln],
        alns2: &[MateAln],
        insert_distr: &D,
        unmapped_penalty: f64,
        prob_diff: f64,
    )
where D: InsertDistr,
{
    let thresh_prob = 2.0 * unmapped_penalty - prob_diff;
    let alns1_empty = alns1.is_empty();
    if !alns1_empty {
        buffer.clear();
        buffer.resize(alns2.len(), f64::NEG_INFINITY);
    }

    for aln1 in alns1.iter() {
        let mut best_prob1 = f64::NEG_INFINITY;
        for (j, aln2) in alns2.iter().enumerate() {
            let insert_size = aln1.interval.furthest_distance(&aln2.interval)
                .expect("Alignments must be on the same contig!");
            let prob = aln1.ln_prob + aln2.ln_prob + insert_distr.ln_prob(insert_size, aln1.strand == aln2.strand);
            if prob >= thresh_prob {
                best_prob1 = best_prob1.max(prob);
                buffer[j] = buffer[j].max(prob);
                new_alns.push(PairAlignment::new(
                    TwoIntervals::Both(aln1.interval.clone(), aln2.interval.clone()), prob));
            }
        }

        let prob = aln1.ln_prob + unmapped_penalty;
        if prob >= thresh_prob && prob + prob_diff >= best_prob1 {
            new_alns.push(PairAlignment::new(TwoIntervals::First(aln1.interval.clone()), prob));
        }
    }

    for (j, aln2) in alns2.iter().enumerate() {
        let prob = aln2.ln_prob + unmapped_penalty;
        if prob >= thresh_prob && (alns1_empty || prob + prob_diff >= buffer[j]) {
            new_alns.push(PairAlignment::new(TwoIntervals::Second(aln2.interval.clone()), prob));
        }
    }
}

/// For a single read pair, combine all single-mate alignments across all contigs.
fn identify_pair_alignments<D>(
    name_hash: u64,
    alns: &[MateAln],
    insert_distr: &D,
    unmapped_penalty: f64,
    prob_diff: f64,
    ln_ncontigs: f64,
) -> ReadPairAlignments
where D: InsertDistr
{
    let mut pair_alns = Vec::new();
    let mut buffer = Vec::with_capacity(16);
    let n = alns.len();
    // For the current contig id, first mates will be in i..j, and second mates in j..k.
    let mut i = 0;
    while i < n {
        let contig_id = alns[i].contig_id();
        // Sort key = contig_id * 2 + (first_mate ? 0 : 1).
        // sort_key1: current sort key for first mates, sort_key2: current sort key for second mates.
        let sort_key1 = contig_id.get() << 1;
        let sort_key2 = sort_key1 | 1;
        let j = if alns[i].read_end == ReadEnd::First {
            bisect::right_by_approx(alns, |aln| aln.sort_key().cmp(&sort_key1), i, n, BISECT_RIGHT_STEP)
        } else { i };
        let k = bisect::right_by_approx(alns, |aln| aln.sort_key().cmp(&sort_key2), j, n, BISECT_RIGHT_STEP);

        extend_pair_alignments(&mut pair_alns, &mut buffer, &alns[i..j], &alns[j..k],
            insert_distr, unmapped_penalty, prob_diff);
        i = k;
    }

    // Probability of both mates unmapped = unmapped_penalty^2.
    let mut unmapped_prob = 2.0 * unmapped_penalty;
    // Normalization factor for all pair-end alignment probabilities.
    // Unmapped probability is multiplied by the number of contigs because there should be an unmapped
    // possibility for every contig, but we do not store them explicitely.
    let norm_fct = Ln::map_sum_init(&pair_alns, PairAlignment::ln_prob, unmapped_prob + ln_ncontigs);
    unmapped_prob -= norm_fct;
    // log::debug!("    {} pair-end alignments. Unmapped prob: {:.2}", pair_alns.len(), unmapped_prob);
    for aln in pair_alns.iter_mut() {
        aln.ln_prob -= norm_fct;
        // log::debug!("        {}", aln);
    }
    ReadPairAlignments { name_hash, unmapped_prob, pair_alns }
}

#[derive(Clone)]
pub(crate) enum TwoIntervals {
    Both(Interval, Interval),
    First(Interval),
    Second(Interval),
    None,
}

/// Alignment of the read pair. At most one of two alignments may be missing!
/// If present, both alignments must map to the same contig and be relatively close to each other.
#[derive(Clone)]
pub struct PairAlignment {
    intervals: TwoIntervals,
    ln_prob: f64,
}

impl PairAlignment {
    fn new(intervals: TwoIntervals, ln_prob: f64) -> Self {
        assert!(ln_prob.is_finite(), "Pair-alignment probability must be finite!");
        Self { intervals, ln_prob }
    }

    /// Log-probability of the paired-read alignment.
    pub fn ln_prob(&self) -> f64 {
        self.ln_prob
    }

    /// Returns contig id of the alignment.
    pub fn contig_id(&self) -> ContigId {
        match &self.intervals {
            TwoIntervals::Both(x, _) | TwoIntervals::First(x) | TwoIntervals::Second(x) => x.contig_id(),
            TwoIntervals::None => panic!("Both mates are unmapped, contig id is unavailable!"),
        }
    }

    /// Creates a copy of `PairAlignment` with a new probability.
    pub fn clone_with_prob(&self, new_prob: f64) -> Self {
        Self {
            intervals: self.intervals.clone(),
            ln_prob: new_prob,
        }
    }

    pub(crate) fn intervals(&self) -> &TwoIntervals {
        &self.intervals
    }
}

impl fmt::Debug for PairAlignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self.intervals {
            TwoIntervals::Both(aln1, aln2) => write!(f, "Pair({:?}, {:?}, {:.2})", aln1, aln2, self.ln_prob),
            TwoIntervals::First(aln1) => write!(f, "Pair({:?}, *, {:.2})", aln1, self.ln_prob),
            TwoIntervals::Second(aln2) => write!(f, "Pair(*, {:?}, {:.2})", aln2, self.ln_prob),
            TwoIntervals::None => write!(f, "Pair(*, *, {:.2})", self.ln_prob),
        }
    }
}

impl fmt::Display for PairAlignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self.intervals {
            TwoIntervals::Both(aln1, aln2) => write!(f, "Pair({}, {}, {:.2})", aln1, aln2, self.ln_prob),
            TwoIntervals::First(aln1) => write!(f, "Pair({}, *, {:.2})", aln1, self.ln_prob),
            TwoIntervals::Second(aln2) => write!(f, "Pair(*, {}, {:.2})", aln2, self.ln_prob),
            TwoIntervals::None => write!(f, "Pair(*, *, {:.2})", self.ln_prob),
        }
    }
}

/// Structure, required for `ReadPairAlignments::multi_contig_alns`.
/// Stores multiple contigs (do not repeat), and their ln-coefficients.
/// Use-case: coeff = ln(contig multiplicity).
pub struct SeveralContigs {
    contigs: Vec<ContigId>,
    coeffs: Vec<f64>,
    /// Sum of coefficients.
    sum_coeff: f64,
}

impl SeveralContigs {
    /// Creates `SeveralContigs` from a slice of contig ids (may repeat).
    /// Assumes that the total number of contigs is small, and, therefore, is not the most sofisticated.
    pub fn new(contigs: &[ContigId]) -> Self {
        let mut dedup_contigs = Vec::new();
        let mut coeffs = Vec::new();
        for contig in contigs {
            if !dedup_contigs.contains(contig) {
                dedup_contigs.push(*contig);
                let multiplicity = contigs.iter().fold(0, |acc, contig2| acc + (contig == contig2) as u8);
                coeffs.push(f64::from(multiplicity).ln());
            }
        }
        Self {
            contigs: dedup_contigs,
            sum_coeff: Ln::sum(&coeffs),
            coeffs,
        }
    }

    /// Returns iterator over pairs `(ContigId, ln_coef)`.
    pub fn iter<'a>(&'a self) -> std::iter::Zip<std::slice::Iter<'a, ContigId>, std::slice::Iter<'a, f64>> {
        self.contigs.iter().zip(self.coeffs.iter())
    }

    /// Returns contig ids.
    pub fn contigs(&self) -> &[ContigId] {
        &self.contigs
    }
}

/// Read-pair alignments for a single read-pair.
pub struct ReadPairAlignments {
    /// Hash of the read name.
    name_hash: u64,

    /// Probability of being unmapped for one contig.
    unmapped_prob: f64,

    /// All pair-read alignments, sorted by contig id.
    pair_alns: Vec<PairAlignment>,
}

impl ReadPairAlignments {
    /// Returns the FNV1a hash of the read name.
    pub fn name_hash(&self) -> u64 {
        self.name_hash
    }

    /// Probability that both reads are unmapped for one specific contig (but same for all contigs).
    pub fn unmapped_prob(&self) -> f64 {
        self.unmapped_prob
    }

    /// For a given contig, returns all alignments to this contig.
    /// Note, that one or both mate may be unmapped.
    pub fn contig_alns(&self, contig_id: ContigId) -> &[PairAlignment] {
        let i = bisect::left_by(&self.pair_alns, |paln| paln.contig_id().cmp(&contig_id));
        let j = bisect::right_by_approx(&self.pair_alns, |paln| paln.contig_id().cmp(&contig_id),
            i, self.pair_alns.len(), BISECT_RIGHT_STEP);
        &self.pair_alns[i..j]
    }

    /// For a given several contigs and their coefficients,
    /// appends all pair-alignments to `out_alns`, and returns the number of added alignments.
    ///
    /// Alignments may include `None` alignments, where both mates do not map to the contigs.
    ///
    /// Output pair-alignments with ln-probability worse than `best_prob - prob_diff` are discarded.
    /// Remaining alignments have random order, and non-normalized probabilities.
    ///
    /// Any alignments that were in the vector before, stay as they are and in the same order.
    pub fn multi_contig_alns(&self,
        out_alns: &mut Vec<PairAlignment>,
        contigs: &SeveralContigs,
        prob_diff: f64,
    ) -> usize {
        let start_len = out_alns.len();
        // Probability of being unmapped to any of the contigs.
        let unmapped_prob = self.unmapped_prob + contigs.sum_coeff;
        // Current threshold, is updated during the for-loop.
        let mut thresh_prob = unmapped_prob - prob_diff;
        for (&contig_id, &coef) in contigs.iter() {
            for paln in self.contig_alns(contig_id) {
                let prob = paln.ln_prob + coef;
                if prob >= thresh_prob {
                    thresh_prob = thresh_prob.max(prob - prob_diff);
                    out_alns.push(paln.clone_with_prob(prob));
                }
            }
        }

        // As threshold was updated during the for-loop, some alignments in the beginning may need to removed.
        let mut i = start_len;
        while i < out_alns.len() {
            if out_alns[i].ln_prob < thresh_prob {
                out_alns.swap_remove(i);
            } else {
                i += 1;
            }
        }
        if unmapped_prob >= thresh_prob {
            out_alns.push(PairAlignment {
                intervals: TwoIntervals::None,
                ln_prob: unmapped_prob
            });
        }
        out_alns.len() - start_len
    }
}

/// All read-pair alignments for all read-pairs.
/// Key: read name hash.
pub struct AllPairAlignments(Vec<ReadPairAlignments>);

impl AllPairAlignments {
    /// Creates a new instance with given capacity.
    fn with_capacity(capacity: usize) -> Self {
        Self(Vec::with_capacity(capacity))
    }

    /// Inserts new read-pair alignments. Name hash is not checked for being in the vector before.
    fn push(&mut self, read_pair_alns: ReadPairAlignments) {
        self.0.push(read_pair_alns);
    }

    /// Calculates sum probability of a contig
    /// (product over all read-pair probabilities for that contig).
    pub fn sum_contig_prob(&self, contig_id: ContigId) -> f64 {
        let mut total_prob = 0.0;
        for paired_alns in self.0.iter() {
            let curr_paired_alns = paired_alns.contig_alns(contig_id);
            let unmapped_prob = paired_alns.unmapped_prob();
            total_prob += Ln::map_sum_init(curr_paired_alns, PairAlignment::ln_prob, unmapped_prob);
        }
        total_prob
    }

    /// Calculates sum probability of multiple contigs.
    /// For each read pair, calculates sum probability to map to any of the contigs,
    /// but all alignments less probable than `best_prob - prob_diff` are discarded.
    pub fn sum_multi_contig_prob(&self, contig_ids: &[ContigId], prob_diff: f64) -> f64 {
        let contigs = SeveralContigs::new(contig_ids);
        let mut total_prob = 0.0;
        // Buffer with current alignments.
        let mut curr_alns = Vec::new();

        for paired_alns in self.0.iter() {
            paired_alns.multi_contig_alns(&mut curr_alns, &contigs, prob_diff);
            total_prob += Ln::map_sum(&curr_alns, PairAlignment::ln_prob);
        }
        total_prob
    }

    /// Returns iterator over `ReadPairAlignments` for all read pairs.
    pub fn iter(&self) -> std::slice::Iter<'_, ReadPairAlignments> {
        self.0.iter()
    }
}