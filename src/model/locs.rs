use std::{
    sync::Arc,
    io::{self, Write},
    fmt,
};
use htslib::bam;
use nohash::IntSet;
use crate::{
    Error,
    seq::{
        Interval, ContigId, ContigNames,
        aln::{Strand, ReadEnd, Alignment},
    },
    bg::{
        err_prof::ErrorProfile,
        insertsz::InsertDistr,
    },
    algo::{fnv1a, bisect},
    math::Ln,
};

/// Write debug information about read locations.
pub(crate) trait DbgWrite {
    /// Write single-mate alignment.
    fn write_mate_aln(&mut self, qname: &[u8], hash: u64, aln: &MateAln) -> io::Result<()>;

    /// Save/ignore `count` mate alignments for read `hash`.
    fn finalize_mate_alns(&mut self, hash: u64, save: bool, count: usize) -> io::Result<()>;
}

impl DbgWrite for () {
    fn write_mate_aln(&mut self, _qname: &[u8], _hash: u64, _aln: &MateAln) -> io::Result<()> { Ok(()) }
    fn finalize_mate_alns(&mut self, _hash: u64, _save: bool, _count: usize,) -> io::Result<()> { Ok(()) }
}

impl<W: io::Write> DbgWrite for &mut W {
    fn write_mate_aln(&mut self, qname: &[u8], hash: u64, aln: &MateAln) -> io::Result<()> {
        write!(self, "{:X}", hash)?;
        if !qname.is_empty() {
            write!(self, "=")?;
            self.write_all(qname)?;
        }
        writeln!(self, "\t{}\t{}\t{:.3}", aln.read_end.as_u16() + 1, aln.interval, Ln::to_log10(aln.ln_prob))
    }

    fn finalize_mate_alns(&mut self, hash: u64, save: bool, count: usize) -> io::Result<()> {
        writeln!(self, "{:X}\t{}\t{} alns", hash, if save { "save" } else { "ignore" }, count)
    }
}

pub(crate) const CSV_HEADER: &'static str = "read_hash\tread_end\tinterval\tlik";

/// Single mate alignment: store alignment location, strand, read-end, and alignment ln-probability.
#[derive(Clone)]
pub(crate) struct MateAln {
    interval: Interval,
    strand: Strand,
    read_end: ReadEnd,
    /// Log-probability of the alignment.
    ln_prob: f64,
}

impl MateAln {
    /// Creates a new alignment extension from a htslib `Record`.
    fn from_record(
        record: &bam::Record,
        contigs: Arc<ContigNames>,
        err_prof: &ErrorProfile,
        read_end: ReadEnd,
    ) -> Self
    {
        let aln = Alignment::from_record(record, contigs);
        Self {
            strand: aln.strand(),
            ln_prob: err_prof.ln_prob(aln.cigar()),
            interval: aln.take_interval(),
            read_end,
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


struct FilteredReader<'a, R: bam::Read> {
    reader: R,
    record: bam::Record,
    has_more: bool,
    contigs: Arc<ContigNames>,
    err_prof: &'a ErrorProfile,
    thresh_prob: f64,
}

impl<'a, R: bam::Read> FilteredReader<'a, R> {
    fn new(
        mut reader: R,
        contigs: Arc<ContigNames>,
        err_prof: &'a ErrorProfile,
        thresh_prob: f64
    ) -> Result<Self, Error> {
        let mut record = bam::Record::new();
        Ok(FilteredReader {
            // Reader would return None if there are no more records.
            has_more: reader.read(&mut record).transpose()?.is_some(),
            reader, record, contigs, err_prof, thresh_prob,
        })
    }

    /// Starting with `self.record` (already loaded), reads all alignments,
    /// corresponding to this read and current read end.
    /// Basically, the function continue to read until the next primary alignment is found,
    /// which is saved to `self.record`.
    ///
    /// Function adds all alignments with probability over the threshold to `alns` and returns name hash.
    fn read_mate_alns(
        &mut self,
        read_end: ReadEnd,
        alns: &mut Vec<MateAln>,
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

        loop {
            let mate_aln = MateAln::from_record(
                &self.record, Arc::clone(&self.contigs), self.err_prof, read_end);
            writeln!(dbg_writer, "{:X}\t{}\t{}\t{:.2}",
                name_hash, read_end, mate_aln.interval, Ln::to_log10(mate_aln.ln_prob))?;
            if mate_aln.ln_prob >= self.thresh_prob {
                alns.push(mate_aln);
            }

            if self.reader.read(&mut self.record).transpose()?.is_none() {
                self.has_more = false;
                return Ok(name_hash);
            }
            // Next primary record or unmapped read. Either way, it will be a new read or new read end.
            if !self.record.is_secondary() {
                return Ok(name_hash);
            }
            assert_eq!(fnv1a(self.record.qname()), name_hash,
                "Read {} first alignment is not primary", String::from_utf8_lossy(self.record.qname()));
        }
    }
}

/// Preliminary read alignments, where read mates are not yet connected.
pub(crate) struct PrelimAlignments {
    contigs: Arc<ContigNames>,
    /// Pairs: read name hash, all alignments for the read pair.
    all_alns: Vec<(u64, Vec<MateAln>)>,
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
        let thresh_prob = params.unmapped_penalty - params.prob_diff;
        log::info!("[{}] Reading read alignments, alignment likelihood threshold = {:.1}",
            contigs.tag(), Ln::to_log10(thresh_prob));
        let mut reader = FilteredReader::new(reader, Arc::clone(contigs), err_prof, thresh_prob)?;
        let left_boundary = params.boundary_size;
        let right_boundaries: Vec<_> = contigs.lengths().iter()
            .map(|len| len.saturating_sub(left_boundary)).collect();
        let mut all_alns = Vec::new();

        let mut hashes = IntSet::default();
        while reader.has_more {
            let read_name = String::from_utf8_lossy(reader.record.qname()).into_owned();
            write!(dbg_writer, "{}=", read_name)?;
            let mut alns = Vec::new();
            let name_hash = reader.read_mate_alns(ReadEnd::First, &mut alns, &mut dbg_writer)?;
            if !hashes.insert(name_hash) {
                log::warn!("Read {} produced hash collision ({:X}). \
                    If many such messages, reads appear in an unordered fashion, will lead to errors",
                    read_name, name_hash);
            }

            if is_paired_end {
                let name_hash2 = reader.read_mate_alns(ReadEnd::Second, &mut alns, &mut dbg_writer)?;
                if name_hash != name_hash2 {
                    return Err(Error::InvalidData(format!("Read {} with hash {:X} does not have a second read end",
                        read_name, name_hash)));
                }
            }

            let any_central = alns.iter().any(|aln| {
                let middle = aln.interval.middle();
                left_boundary <= middle && middle <= right_boundaries[aln.interval.contig_id().ix()]
            });
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
            alns.sort_unstable_by_key(MateAln::sort_key);
            let pair_alns = identify_pair_alignments(*name_hash, alns, insert_distr, ln_ncontigs, params);
            res.push(pair_alns);
        }
        res
    }
}

// Assume that there are at most 4 alignments of the to the same contig.
// (not critical, if this assumption fails).
const BISECT_RIGHT_STEP: usize = 4;

/// For a single read-pair, find all paired-read alignments to the same contig.
fn extend_pair_alignments(
    new_alns: &mut Vec<PairAlignment>,
    buffer: &mut Vec<f64>,
    alns1: &[MateAln],
    alns2: &[MateAln],
    insert_distr: &InsertDistr,
    params: &super::Params,
) {
    let thresh_prob = 2.0 * params.unmapped_penalty - params.prob_diff;
    let alns1_empty = alns1.is_empty();
    if !alns1_empty {
        buffer.clear();
        // Buffer stores best probabilities for each read2.
        buffer.resize(alns2.len(), f64::NEG_INFINITY);
    }

    for aln1 in alns1.iter() {
        let mut best_prob1 = f64::NEG_INFINITY;
        for (j, aln2) in alns2.iter().enumerate() {
            let insert_size = aln1.interval.furthest_distance(&aln2.interval)
                .expect("Alignments must be on the same contig");
            let prob = aln1.ln_prob + aln2.ln_prob + insert_distr.ln_prob(insert_size, aln1.strand == aln2.strand);
            if prob >= thresh_prob {
                best_prob1 = best_prob1.max(prob);
                buffer[j] = buffer[j].max(prob);
                new_alns.push(PairAlignment::new(
                    TwoIntervals::Both(aln1.interval.clone(), aln2.interval.clone()), prob));
            }
        }

        // Add unpaired read1, if needed.
        let prob = aln1.ln_prob + params.unmapped_penalty;
        if prob >= thresh_prob && prob + params.prob_diff >= best_prob1 {
            new_alns.push(PairAlignment::new(TwoIntervals::First(aln1.interval.clone()), prob));
        }
    }

    // Add unpaired read2, if needed.
    for (j, aln2) in alns2.iter().enumerate() {
        let prob = aln2.ln_prob + params.unmapped_penalty;
        // Check on `alns1_empty`, as if they are empty, buffer was not cleaned.
        if prob >= thresh_prob && (alns1_empty || prob + params.prob_diff >= buffer[j]) {
            new_alns.push(PairAlignment::new(TwoIntervals::Second(aln2.interval.clone()), prob));
        }
    }
}

/// For a single read pair, combine all single-mate alignments across all contigs.
fn identify_pair_alignments(
    name_hash: u64,
    alns: &[MateAln],
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
        let j = if alns[i].read_end == ReadEnd::First {
            bisect::right_by_approx(alns, |aln| aln.sort_key().cmp(&sort_key1), i, n, BISECT_RIGHT_STEP)
        } else { i };
        let k = bisect::right_by_approx(alns, |aln| aln.sort_key().cmp(&sort_key2), j, n, BISECT_RIGHT_STEP);
        extend_pair_alignments(&mut pair_alns, &mut buffer, &alns[i..j], &alns[j..k], insert_distr, params);
        i = k;
    }

    // Probability of both mates unmapped = unmapped_penalty^2.
    let mut unmapped_prob = 2.0 * params.unmapped_penalty;
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
