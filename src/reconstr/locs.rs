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
        aln::{ReadEnd, Alignment},
    },
    bg::{
        err_prof::ErrorProfile,
        insertsz::InsertDistr,
    },
    algo::{
        math::Ln,
        hash::fnv1a,
        bisect,
    },
};

/// Extension over alignment: store alignment, alignment probability and the read end.
#[derive(Clone)]
struct ExtAln {
    aln: Alignment,
    read_end: ReadEnd,
    /// Log-probability of the alignment.
    ln_prob: f64,
}

impl ExtAln {
    /// Creates a new alignment extension from a htslib `Record`.
    fn from_record<E: ErrorProfile>(record: &Record, contigs: Rc<ContigNames>, err_prof: &E) -> Self {
        let aln = Alignment::from_record(contigs, record);
        let read_end = ReadEnd::from_record(record);
        let ln_prob = err_prof.ln_prob(aln.cigar(), read_end);
        Self { aln, read_end, ln_prob }
    }

    /// Get contig id of the alignment.
    fn contig_id(&self) -> ContigId {
        self.aln.contig_id()
    }

    /// Returns sorted key: first sort by contig id, then by read end.
    fn sort_key(&self) -> u32 {
        (self.contig_id().get() << 1) | (self.read_end.as_u32())
    }
}

impl fmt::Debug for ExtAln {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ExtAln[{:?}, {:?}, prob={:.2}]", self.read_end, self.aln, self.ln_prob)
    }
}

impl fmt::Display for ExtAln {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ExtAln[{}, {}, prob={:.2}]", self.read_end, self.aln, self.ln_prob)
    }
}

/// Preliminary, unpaired, read alignments.
/// Keys: read name hash, values: all alignments for the read pair.
pub struct PrelimAlignments {
    contigs: Rc<ContigNames>,
    alns: IntMap<Vec<ExtAln>>,
}

impl PrelimAlignments {
    /// Creates `PrelimAlignments` from records.
    /// Only store read pairs where at least one alignment (for at least one read),
    /// has probability over `min_ln_prob`.
    pub fn from_records<I, E>(
        records: I,
        contigs: Rc<ContigNames>,
        err_prof: &E,
        min_ln_prob: f64,
        ) -> Self
    where I: Iterator<Item = Record>,
          E: ErrorProfile,
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
            let ext_aln = ExtAln::from_record(&record, Rc::clone(&contigs), err_prof);
            if hash != curr_hash {
                if best_prob >= min_ln_prob {
                    // log::debug!("    Saving {} alignments for hash {}  (prob {:.1})",
                    //     curr_alns.len(), curr_hash, best_prob);
                    // Assert must not be removed, or converted into debug_assert!
                    assert!(alns.insert(curr_hash, mem::replace(&mut curr_alns, Vec::new())).is_none(),
                        "Read pair alignments are separated by another read pair (see hash {})", curr_hash);
                } else {
                    // log::debug!("    Ignoring {} alignments for hash {}  (prob {:.1})",
                    //     curr_alns.len(), curr_hash, best_prob);
                    curr_alns.clear();
                }
                curr_hash = hash;
                best_prob = f64::NEG_INFINITY;
            }
            // log::debug!("    {}  {}  {:?}", String::from_utf8_lossy(record.qname()), hash, ext_aln);
            best_prob = best_prob.max(ext_aln.ln_prob);
            curr_alns.push(ext_aln);
        }

        if best_prob >= min_ln_prob {
            // log::debug!("    Saving {} alignments for hash {}  (prob {:.1})",
            //     curr_alns.len(), curr_hash, best_prob);
            // Assert must not be removed, or converted into debug_assert!
            assert!(alns.insert(curr_hash, curr_alns).is_none(),
                "Read pair alignments are separated by another read pair (see hash {})", curr_hash);
        } else {
            // log::debug!("    Ignoring {} alignments for hash {}  (prob {:.1})",
            //     curr_alns.len(), curr_hash, best_prob);
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
        let n_contigs = self.contigs.len();
        log::info!("Identify paired alignment location and probabilities ({} read pairs)", n_reads);
        let mut res = AllPairAlignments::with_capacity(n_reads);
        for (&name_hash, alns) in self.alns.iter_mut() {
            log::debug!("Read {}", name_hash);
            let pair_alns = identify_pair_alignments(alns, n_contigs as u32, insert_distr, unmapped_penalty, prob_diff);
            res.insert(name_hash, pair_alns);
        }
        res
    }
}

// Assume that there are at most 4 alignments of the to the same contig.
// (not critical, if this assumption fails).
const BISECT_RIGHT_STEP: usize = 4;

/// For a single read-pair, find all paired-read alignments to the same contig.
fn extend_pair_alignments<D>(
        contig_id: ContigId,
        new_alns: &mut Vec<PairAlignment>,
        buffer: &mut Vec<PairAlignment>,
        alns1: &[ExtAln],
        alns2: &[ExtAln],
        insert_distr: &D,
        unmapped_penalty: f64,
        prob_diff: f64,
    )
where D: InsertDistr,
{
    let unmapped2 = 2.0 * unmapped_penalty;
    // Current threshold probability. Always equal to current best probability - prob_diff.
    let mut thresh = unmapped2 - prob_diff;

    for aln1 in alns1.iter() {
        debug_assert_eq!(aln1.contig_id(), contig_id);
        for aln2 in alns2.iter() {
            let prob = aln1.ln_prob + aln2.ln_prob + insert_distr.pair_aln_prob(&aln1.aln, &aln2.aln);
            if prob >= thresh {
                thresh = thresh.max(prob - prob_diff);
                buffer.push(PairAlignment::new(TwoAlignments::Both(aln1.aln.clone(), aln2.aln.clone()), prob));
            }
        }

        let prob = aln1.ln_prob + unmapped_penalty;
        if prob >= thresh {
            thresh = thresh.max(prob - prob_diff);
            buffer.push(PairAlignment::new(TwoAlignments::First(aln1.aln.clone()), prob));
        }
    }

    for aln2 in alns2.iter() {
        debug_assert_eq!(aln2.contig_id(), contig_id);
        let prob = aln2.ln_prob + unmapped_penalty;
        if prob >= thresh {
            thresh = thresh.max(prob - prob_diff);
            buffer.push(PairAlignment::new(TwoAlignments::Second(aln2.aln.clone()), prob));
        }
    }

    let mut i = 0;
    while i < buffer.len() {
        if unsafe { buffer.get_unchecked(i) }.ln_prob >= thresh {
            new_alns.push(buffer.swap_remove(i));
        } else {
            i += 1;
        }
    }
    if unmapped2 >= thresh {
        new_alns.push(PairAlignment::new(TwoAlignments::None(contig_id), unmapped2));
    }
    buffer.clear();
}

fn identify_pair_alignments<D>(
    alns: &mut [ExtAln],
    n_contigs: u32,
    insert_distr: &D,
    unmapped_penalty: f64,
    prob_diff: f64,
    ) -> ReadPairAlignments
where D: InsertDistr
{
    // Sort alignments first by contig id, then by read-end.
    alns.sort_by_key(ExtAln::sort_key);
    let n = alns.len();
    // For the current contig id, first mates will be in i..j, and second mates in j..k.
    let mut i = 0;
    let mut pair_alns = Vec::new();
    let mut buffer = Vec::new();

    let mut contig_shift = 0;
    while i < n {
        let contig_id = alns[i].contig_id();
        // Add unmapped pair-alignments for all contigs in-between previous and this.
        for j in contig_shift..contig_id.get() {
            pair_alns.push(PairAlignment::new(TwoAlignments::None(ContigId::new(j)), 2.0 * unmapped_penalty));
        }

        // Sort key = contig_id * 2 + (first_mate ? 0 : 1).
        // sort_key1: current sort key for first mates, sort_key2: current sort key for second mates.
        let sort_key1 = contig_id.get() << 1;
        let sort_key2 = sort_key1 | 1;
        let j = if alns[i].read_end == ReadEnd::First {
            bisect::right_by_approx(alns, |aln| aln.sort_key().cmp(&sort_key1), i, n, BISECT_RIGHT_STEP)
        } else { i };
        let k = bisect::right_by_approx(alns, |aln| aln.sort_key().cmp(&sort_key2), j, n, BISECT_RIGHT_STEP);

        extend_pair_alignments(contig_id, &mut pair_alns, &mut buffer, &alns[i..j], &alns[j..k],
            insert_distr, unmapped_penalty, prob_diff);
        i = k;
        contig_shift = contig_id.get() + 1;
    }
    for j in contig_shift..n_contigs {
        pair_alns.push(PairAlignment::new(TwoAlignments::None(ContigId::new(j)), 2.0 * unmapped_penalty));
    }

    // Normalization factor for all pair-end alignment probabilities.
    let norm_fct = Ln::map_sum(&pair_alns, PairAlignment::ln_prob);
    log::debug!("    {} pair-end alignments", pair_alns.len());
    for aln in pair_alns.iter_mut() {
        aln.ln_prob -= norm_fct;
        log::debug!("        {}", aln);
    }
    ReadPairAlignments(pair_alns)
}

#[derive(Clone)]
enum TwoAlignments {
    Both(Alignment, Alignment),
    First(Alignment),
    Second(Alignment),
    None(ContigId),
}

/// Alignment of the read pair. At most one of two alignments may be missing!
/// If present, both alignments must map to the same contig and be relatively close to each other.
#[derive(Clone)]
pub struct PairAlignment {
    alns: TwoAlignments,
    ln_prob: f64,
}

impl PairAlignment {
    fn new(alns: TwoAlignments, ln_prob: f64) -> Self {
        Self { alns, ln_prob }
    }

    /// Log-probability of the paired-read alignment.
    pub fn ln_prob(&self) -> f64 {
        self.ln_prob
    }

    /// Returns contig id of the alignment.
    pub fn contig_id(&self) -> ContigId {
        match &self.alns {
            TwoAlignments::Both(aln, _) | TwoAlignments::First(aln) | TwoAlignments::Second(aln) => aln.contig_id(),
            TwoAlignments::None(contig_id) => *contig_id,
        }
    }
}

impl fmt::Debug for PairAlignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self.alns {
            TwoAlignments::Both(aln1, aln2) => write!(f, "Pair({:?}, {:?}, {:.2})", aln1, aln2, self.ln_prob),
            TwoAlignments::First(aln1) => write!(f, "Pair({:?}, *, {:.2})", aln1, self.ln_prob),
            TwoAlignments::Second(aln2) => write!(f, "Pair(*, {:?}, {:.2})", aln2, self.ln_prob),
            TwoAlignments::None(contig_id) => write!(f, "Pair({:?}, *, *, {:.2})", contig_id, self.ln_prob),
        }
    }
}

impl fmt::Display for PairAlignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self.alns {
            TwoAlignments::Both(aln1, aln2) => write!(f, "Pair({}, {}, {:.2})", aln1, aln2, self.ln_prob),
            TwoAlignments::First(aln1) => write!(f, "Pair({}, *, {:.2})", aln1, self.ln_prob),
            TwoAlignments::Second(aln2) => write!(f, "Pair(*, {}, {:.2})", aln2, self.ln_prob),
            TwoAlignments::None(contig_id) => write!(f, "Pair({}, *, *, {:.2})", contig_id, self.ln_prob),
        }
    }
}

/// Read-pair alignments for a single read-pair.
/// Vector is sorted by contig id.
pub struct ReadPairAlignments(Vec<PairAlignment>);

impl ReadPairAlignments {
    /// For a given contig, returns all alignments to this contig.
    /// Note, that one or both mate may be unmapped.
    pub fn contig_alns(&self, contig_id: ContigId) -> &[PairAlignment] {
        let i = bisect::left_by(&self.0, |paln| paln.contig_id().cmp(&contig_id));
        let j = bisect::right_by_approx(&self.0, |paln| paln.contig_id().cmp(&contig_id),
            i, self.0.len(), BISECT_RIGHT_STEP);
        assert_ne!(i, j, "No pair-alignments found for {}", contig_id);
        &self.0[i..j]
    }
}

/// All read-pair alignments for all read-pairs.
/// Key: read name hash.
pub struct AllPairAlignments(IntMap<ReadPairAlignments>);

impl AllPairAlignments {
    /// Creates a new instance with given capacity.
    fn with_capacity(capacity: usize) -> Self {
        Self(IntMap::with_capacity(capacity))
    }

    /// Inserts new read-pair alignments. Name hash must be new (not in the map).
    fn insert(&mut self, name_hash: u64, read_pair_alns: ReadPairAlignments) {
        // Assert must not be removed, or converted into debug_assert!
        assert!(self.0.insert(name_hash, read_pair_alns).is_none(), "Duplicate read name hash {}", name_hash);
    }

    /// Calculates sum probability of a contig
    /// (product over all read-pair probabilities for that contig).
    pub fn sum_contig_prob(&self, contig_id: ContigId) -> f64 {
        let mut total_prob = 0.0;
        for paired_alns in self.0.values() {
            total_prob += Ln::map_sum(paired_alns.contig_alns(contig_id), PairAlignment::ln_prob);
        }
        total_prob
    }
}