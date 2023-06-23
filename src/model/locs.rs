use std::{
    fmt,
    sync::Arc,
    io::Write,
};
use htslib::bam;
use nohash::IntSet;
use crate::{
    Error,
    seq::{
        Interval, ContigId, ContigNames,
        aln::{LightAlignment, ReadEnd},
        cigar::Cigar,
    },
    bg::{
        err_prof::{OpCounts, ErrorProfile},
        insertsz::InsertDistr,
    },
    algo::{fnv1a, bisect},
    math::Ln,
};

pub(crate) const CSV_HEADER: &'static str = "read_hash\tread_end\tinterval\tedit_dist\tpasses\tlik";

/// BAM reader, which stores the next record within.
/// Contains `read_mate_alns` method, which consecutively reads all records with the same name
/// until the next primary alignment.
struct FilteredReader<'a, R: bam::Read> {
    reader: R,
    record: bam::Record,
    has_more: bool,
    contigs: Arc<ContigNames>,
    err_prof: &'a ErrorProfile,
    /// Allow alignments with ln-probability at least `best_prob - max_diff` for the same read end.
    max_diff: f64,
}

impl<'a, R: bam::Read> FilteredReader<'a, R> {
    fn new(
        mut reader: R,
        contigs: Arc<ContigNames>,
        err_prof: &'a ErrorProfile,
        max_diff: f64
    ) -> Result<Self, Error> {
        let mut record = bam::Record::new();
        Ok(FilteredReader {
            // Reader would return None if there are no more records.
            has_more: reader.read(&mut record).transpose()?.is_some(),
            reader, record, contigs, err_prof, max_diff,
        })
    }

    /// Starting with `self.record` (already loaded), reads all alignments,
    /// corresponding to this read and current read end.
    /// Basically, the function continue to read until the next primary alignment is found,
    /// which is saved to `self.record`.
    ///
    /// If best edit distance is too high, does not add any new alignments.
    /// Otherwise, normalizes all alignment probabilities by the best probability (therefore max is always = 0).
    fn load_read_end_alns(
        &mut self,
        read_end: ReadEnd,
        alns: &mut Vec<LightAlignment>,
        dbg_writer: &mut impl Write,
    ) -> Result<u64, Error>
    {
        assert!(self.has_more, "Cannot read any more records from a BAM file.");
        let name_hash = fnv1a(self.record.qname());
        if self.record.is_unmapped() {
            writeln!(dbg_writer, "{:X}\t{}\t*\tNA", name_hash, read_end)?;
            // Read the next record, and return false if there are no more records left.
            self.has_more = self.reader.read(&mut self.record).transpose()?.is_some();
            return Ok(name_hash);
        }

        let mut any_passes = false;
        let mut max_prob = f64::NEG_INFINITY;
        let start_len = alns.len();
        loop {
            let cigar = Cigar::from_raw(&self.record);
            let read_prof = OpCounts::from_cigar(&cigar);
            let read_len = cigar.query_len();
            let edit_dist = read_prof.edit_distance(read_len);
            let aln_passes = edit_dist <= self.err_prof.allowed_edit_dist(read_len);
            any_passes |= aln_passes;
            let aln_prob = self.err_prof.ln_prob(&read_prof);
            let aln = LightAlignment::new(&self.record, &cigar, read_end, Arc::clone(&self.contigs), aln_prob);
            writeln!(dbg_writer, "{:X}\t{}\t{}\t{}\t{}\t{:.2}", name_hash, read_end, aln.interval(),
                edit_dist, if aln_passes { 'T' } else { 'F' }, Ln::to_log10(aln_prob))?;

            if aln_prob >= max_prob - self.max_diff {
                max_prob = max_prob.max(aln_prob);
                alns.push(aln);
            }

            if self.reader.read(&mut self.record).transpose()?.is_none() {
                self.has_more = false;
                break;
            }
            // Next primary record or unmapped read. Either way, it will be a new read or new read end.
            if !self.record.is_secondary() {
                break;
            }
            assert_eq!(fnv1a(self.record.qname()), name_hash,
                "Read {} first alignment is not primary", String::from_utf8_lossy(self.record.qname()));
        }

        if any_passes {
            let mut i = start_len;
            let thresh_prob = max_prob - self.max_diff;
            assert!(max_prob.is_finite(), "Maximum probability is infinite for read hash {:X}", name_hash);
            while i < alns.len() {
                let prob = alns[i].ln_prob();
                if prob >= thresh_prob {
                    alns[i].set_ln_prob(prob - max_prob);
                    i += 1;
                } else {
                    alns.swap_remove(i);
                }
            }
        } else {
            alns.truncate(start_len);
        }
        Ok(name_hash)
    }
}

/// Preliminary read alignments, where read mates are not yet connected.
pub(crate) struct PrelimAlignments {
    contigs: Arc<ContigNames>,
    /// Pairs: read name hash, all alignments for the read pair.
    all_alns: Vec<(u64, Vec<LightAlignment>)>,
}

impl PrelimAlignments {
    /// Creates `PrelimAlignments` from records.
    ///
    /// Only store read pairs where at least one read-end alignment
    /// (i)  has probability over `unmapped_penalty` AND
    /// (ii) middle of the alignment lies within an inner region of any contig (beyond the boundary region).
    pub(crate) fn load(
        reader: impl bam::Read,
        is_paired_end: bool,
        contigs: &Arc<ContigNames>,
        err_prof: &ErrorProfile,
        params: &super::Params,
        mut dbg_writer: impl Write,
    ) -> Result<Self, Error>
    {
        log::info!("[{}] Loading read alignments", contigs.tag());
        let max_diff = -(params.unmapped_penalty - params.prob_diff);
        let mut reader = FilteredReader::new(reader, Arc::clone(contigs), err_prof, max_diff)?;
        let left_boundary = params.boundary_size - params.tweak;
        let right_boundaries: Vec<_> = contigs.lengths().iter()
            .map(|len| len.saturating_sub(left_boundary)).collect();
        let mut all_alns = Vec::new();

        let mut hashes = IntSet::default();
        while reader.has_more {
            let read_name = String::from_utf8_lossy(reader.record.qname()).into_owned();
            write!(dbg_writer, "{}=", read_name)?;
            let mut alns = Vec::new();
            let name_hash = reader.load_read_end_alns(ReadEnd::First, &mut alns, &mut dbg_writer)?;
            if !hashes.insert(name_hash) {
                log::warn!("Read {} produced hash collision ({:X}). \
                    If many such messages, reads appear in an unordered fashion, will lead to errors",
                    read_name, name_hash);
            }

            if is_paired_end {
                let name_hash2 = reader.load_read_end_alns(ReadEnd::Second, &mut alns, &mut dbg_writer)?;
                if name_hash != name_hash2 {
                    return Err(Error::InvalidData(format!("Read {} with hash {:X} does not have a second read end",
                        read_name, name_hash)));
                }
            }

            // Does any of the read alignment lie within the main part of the locus?
            let any_central = alns.iter().any(|aln| left_boundary < aln.interval().end()
                && aln.interval().start() < right_boundaries[aln.contig_id().ix()]);
            if any_central {
                writeln!(dbg_writer, "{:X}\tsave\t{} alns", name_hash, alns.len())?;
                all_alns.push((name_hash, alns));
            } else {
                writeln!(dbg_writer, "{:X}\tignore\t{} alns", name_hash, alns.len())?;
            }
        }
        Ok(Self {
            contigs: Arc::clone(&contigs),
            all_alns
        })
    }

    /// Find paired-end alignment locations for each contig.
    pub fn identify_locations(&mut self, insert_distr: &InsertDistr, params: &super::Params) -> AllPairAlignments {
        let n_reads = self.all_alns.len();
        let ln_ncontigs = (self.contigs.len() as f64).ln();
        log::info!("    [{}] Identify paired alignment location and probabilities ({} {})",
            self.contigs.tag(), n_reads, if insert_distr.is_paired_end() { "read pairs" } else { "reads" });
        let mut res = AllPairAlignments::with_capacity(n_reads);
        for (name_hash, alns) in self.all_alns.iter_mut() {
            // log::debug!("Read {:X}", name_hash);
            // Sort alignments first by contig id, then by read-end.
            alns.sort_unstable_by_key(LightAlignment::sort_key);
            let pair_alns = identify_pair_alignments(*name_hash, alns, insert_distr, ln_ncontigs, params);
            res.push(pair_alns);
        }
        res
    }
}

/// Returns ln-probability, assigned to the unmapped read mate pair.
/// Is equal to two unmapped penalties (both mates unmapped) + insert size penalty.
#[inline]
fn unmapped_pair_prob(unmapped_penalty: f64, insert_penalty: f64) -> f64 {
    2.0 * unmapped_penalty + insert_penalty
}

// Assume that there are at most 4 alignments of the to the same contig.
// (not critical, if this assumption fails).
const BISECT_RIGHT_STEP: usize = 4;

/// For a single read-pair, find all paired-read alignments to the same contig.
fn extend_pair_alignments(
    new_alns: &mut Vec<PairAlignment>,
    buffer: &mut Vec<f64>,
    alns1: &[LightAlignment],
    alns2: &[LightAlignment],
    insert_distr: &InsertDistr,
    params: &super::Params,
) {
    let insert_penalty = insert_distr.worst_prob();
    let thresh_prob = unmapped_pair_prob(params.unmapped_penalty, insert_penalty) - params.prob_diff;
    let alns1_empty = alns1.is_empty();
    if !alns1_empty {
        buffer.clear();
        // Buffer stores best probabilities for each read2.
        buffer.resize(alns2.len(), f64::NEG_INFINITY);
    }

    for aln1 in alns1.iter() {
        let mut best_prob1 = f64::NEG_INFINITY;
        for (j, aln2) in alns2.iter().enumerate() {
            let prob = aln1.paired_prob(aln2, insert_distr);
            if prob >= thresh_prob {
                best_prob1 = best_prob1.max(prob);
                buffer[j] = buffer[j].max(prob);
                new_alns.push(PairAlignment::new(
                    TwoIntervals::Both(aln1.interval().clone(), aln2.interval().clone()), prob));
            }
        }

        // Add unpaired read1, if needed.
        let prob = aln1.ln_prob() + params.unmapped_penalty + insert_penalty;
        if prob >= thresh_prob && prob + params.prob_diff >= best_prob1 {
            new_alns.push(PairAlignment::new(TwoIntervals::First(aln1.interval().clone()), prob));
        }
    }

    // Add unpaired read2, if needed.
    for (j, aln2) in alns2.iter().enumerate() {
        let prob = aln2.ln_prob() + params.unmapped_penalty + insert_penalty;
        // Check on `alns1_empty`, as if they are empty, buffer was not cleaned.
        if prob >= thresh_prob && (alns1_empty || prob + params.prob_diff >= buffer[j]) {
            new_alns.push(PairAlignment::new(TwoIntervals::Second(aln2.interval().clone()), prob));
        }
    }
}

/// For a single read pair, combine all single-mate alignments across all contigs.
fn identify_pair_alignments(
    name_hash: u64,
    alns: &[LightAlignment],
    insert_distr: &InsertDistr,
    ln_ncontigs: f64,
    params: &super::Params,
) -> ReadPairAlignments {
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
        let j = if alns[i].read_end() == ReadEnd::First {
            bisect::right_by_approx(alns, |aln| aln.sort_key().cmp(&sort_key1), i, n, BISECT_RIGHT_STEP)
        } else { i };
        let k = bisect::right_by_approx(alns, |aln| aln.sort_key().cmp(&sort_key2), j, n, BISECT_RIGHT_STEP);
        extend_pair_alignments(&mut pair_alns, &mut buffer, &alns[i..j], &alns[j..k], insert_distr, params);
        i = k;
    }

    // Probability of both mates unmapped.
    let mut unmapped_prob = unmapped_pair_prob(params.unmapped_penalty, insert_distr.worst_prob());
    // Normalization factor for all pair-end alignment probabilities.
    // For normalization, unmapped probability is multiplied by the number of contigs because there is an unmapped
    // possibility for every contig, which we do not store explicitely.
    //
    // Unmapped probability is not multiplied by ln_ncontigs, it is just for normalization.
    let norm_fct = Ln::map_sum_init(&pair_alns, PairAlignment::ln_prob, unmapped_prob + ln_ncontigs);
    unmapped_prob -= norm_fct;
    // log::debug!("    {} pair-end alignments. Unmapped prob: {:.2}, norm_fct: {:.2}",
    //     pair_alns.len(), Ln::to_log10(unmapped_prob), Ln::to_log10(norm_fct));
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
    // None,
}

impl TwoIntervals {
    /// Returns the range of the first alignment.
    pub fn range1(&self) -> Option<(u32, u32)> {
        match self {
            TwoIntervals::Both(interval1, _) | TwoIntervals::First(interval1) => Some(interval1.range()),
            TwoIntervals::Second(_) => None,
        }
    }

    /// Returns the range of the second alignment.
    pub fn range2(&self) -> Option<(u32, u32)> {
        match self {
            TwoIntervals::Both(_, interval2) | TwoIntervals::Second(interval2) => Some(interval2.range()),
            TwoIntervals::First(_) => None,
        }
    }
}

/// LightAlignment of the read pair. At most one of two alignments may be missing!
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
            // TwoIntervals::None => panic!("Both mates are unmapped, contig id is unavailable!"),
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
            // TwoIntervals::None => write!(f, "Pair(*, *, {:.2})", self.ln_prob),
        }
    }
}

impl fmt::Display for PairAlignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self.intervals {
            TwoIntervals::Both(aln1, aln2) => write!(f, "Pair({}, {}, {:.2})", aln1, aln2, self.ln_prob),
            TwoIntervals::First(aln1) => write!(f, "Pair({}, *, {:.2})", aln1, self.ln_prob),
            TwoIntervals::Second(aln2) => write!(f, "Pair(*, {}, {:.2})", aln2, self.ln_prob),
            // TwoIntervals::None => write!(f, "Pair(*, *, {:.2})", self.ln_prob),
        }
    }
}

/// Read-pair alignments for a single read-pair.
pub struct ReadPairAlignments {
    /// Hash of the read name.
    name_hash: u64,

    /// Probability of that both read mates are unmapped.
    /// Needs to be normalized by ploidy.
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

    /// For a given contig, returns
    /// - index of the first appropriate alignment,
    /// - all alignments to the contig.
    pub fn contig_alns(&self, contig_id: ContigId) -> (usize, &[PairAlignment]) {
        let i = bisect::left_by(&self.pair_alns, |paln| paln.contig_id().cmp(&contig_id));
        let j = bisect::right_by_approx(&self.pair_alns, |paln| paln.contig_id().cmp(&contig_id),
            i, self.pair_alns.len(), BISECT_RIGHT_STEP);
        (i, &self.pair_alns[i..j])
    }

    /// Return `i`-th alignment of all alignments for the read pair.
    pub fn ith_aln(&self, i: usize) -> &PairAlignment {
        &self.pair_alns[i]
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
            let curr_paired_alns = paired_alns.contig_alns(contig_id).1;
            let unmapped_prob = paired_alns.unmapped_prob();
            total_prob += Ln::map_sum_init(curr_paired_alns, PairAlignment::ln_prob, unmapped_prob);
        }
        total_prob
    }

    /// Returns iterator over `ReadPairAlignments` for all read pairs.
    pub fn iter(&self) -> std::slice::Iter<'_, ReadPairAlignments> {
        self.0.iter()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }
}
