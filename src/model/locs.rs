use std::{
    sync::Arc,
    io::Write,
    collections::hash_map::Entry,
};
use htslib::bam;
use nohash::{IntSet, IntMap};
use crate::{
    err::{Error, add_path},
    seq::{
        Interval, ContigId, ContigNames,
        aln::{Alignment, ReadEnd},
        cigar::Cigar,
    },
    bg::{
        BgDistr,
        err_prof::ErrorProfile,
        insertsz::InsertDistr,
    },
    algo::{fnv1a, bisect},
    math::Ln,
    model::windows::ContigWindows,
};

struct MateSummary {
    name_hash: u64,
    weight: f64,
    min_prob: f64,
    any_central: bool,
    mapped: bool,
}

impl MateSummary {
    fn unmapped(name_hash: u64) -> Self {
        Self {
            name_hash,
            weight: 0.0,
            min_prob: 0.0,
            any_central: false,
            mapped: false,
        }
    }
}

/// Read information: read name, two sequences and qualities (for both mates).
#[derive(Default, Clone)]
struct ReadData {
    name: String,
    sequences: [Vec<u8>; 2],
    qualities: [Vec<u8>; 2],
}

impl ReadData {
    fn clone_and_reset(&mut self) -> Self {
        let clone = self.clone();
        self.name.clear();
        for i in 0..2 {
            self.sequences[i].clear();
            self.qualities[i].clear();
        }
        clone
    }
}

/// BAM reader, which stores the next record within.
/// Contains `read_mate_alns` method, which consecutively reads all records with the same name
/// until the next primary alignment.
struct FilteredReader<'a, R: bam::Read> {
    reader: R,
    record: bam::Record,
    has_more: bool,
    contigs: Arc<ContigNames>,
    boundary: u32,
    err_prof: &'a ErrorProfile,
    contig_windows: &'a [ContigWindows],
    /// Buffer, used for discrard alignments with the same positions.
    /// Key (contig id << 32 | alignment start), value: index of the previous position.
    found_alns: IntMap<u64, usize>,
}

#[inline]
fn is_primary(record: &bam::Record) -> bool {
    (record.flags() & 2304) == 0
}

impl<'a, R: bam::Read> FilteredReader<'a, R> {
    fn new(
        mut reader: R,
        contigs: Arc<ContigNames>,
        err_prof: &'a ErrorProfile,
        contig_windows: &'a [ContigWindows],
        boundary: u32,
    ) -> Result<Self, Error> {
        let mut record = bam::Record::new();
        // Reader would return None if there are no more records.
        let has_more = reader.read(&mut record).transpose()?.is_some();
        assert!(is_primary(&record), "First record in the BAM file is secondary/supplementary");
        Ok(FilteredReader {
            found_alns: IntMap::default(),
            reader, record, contigs, err_prof, contig_windows, boundary, has_more,
        })
    }

    /// Starting with `self.record` (already loaded), reads all alignments,
    /// corresponding to this read and current read end.
    /// Basically, the function continue to read until the next primary alignment is found,
    /// which is saved to `self.record`.
    ///
    /// If read is unmapped, or best edit distance is too high, does not add any new alignments.
    fn next_alns(
        &mut self,
        read_end: ReadEnd,
        alns: &mut Vec<Alignment>,
        read_data: &mut ReadData,
        dbg_writer: &mut impl Write,
    ) -> Result<MateSummary, Error>
    {
        assert!(self.has_more, "Cannot read any more records from a BAM file");
        let name_hash = fnv1a(self.record.qname());
        let name = std::str::from_utf8(self.record.qname()).map_err(|_| Error::InvalidInput(
            format!("Read name is not UTF-8: {:?}", String::from_utf8_lossy(self.record.qname()))))?;
        if read_end == ReadEnd::First {
            read_data.name = name.to_owned();
        } else if name != read_data.name {
            return Err(Error::InvalidData(format!("Read {} does not have a second read end (old {})", name, read_data.name)));
        }
        if self.record.is_unmapped() {
            writeln!(dbg_writer, "{:X}\t{}\t*\tF\tNA\t-\tNA\tNA\t{}", name_hash, read_end, name)
                .map_err(add_path!(!))?;
            // Read the next record, and save false to `has_more` if there are no more records left.
            self.has_more = self.reader.read(&mut self.record).transpose()?.is_some();
            return Ok(MateSummary::unmapped(name_hash));
        }

        let start_len = alns.len();
        // Is any of the alignments in the central region of a contig?
        let mut any_central = false;
        // Does any alignment have good edit distance?
        let mut any_good = false;
        let mut min_prob = f64::INFINITY;
        let mut weight = f64::NAN;
        self.found_alns.clear();
        loop {
            let primary = weight.is_nan();
            let mut cigar = Cigar::from_raw(self.record.raw_cigar());
            if primary {
                assert!(!cigar.has_hard_clipping(), "Primary alignment has hard clipping");
                if self.record.seq().is_empty() {
                    return Err(Error::InvalidData(format!("Read {} does not have read sequence", read_data.name)));
                }
                read_data.sequences[read_end.ix()] = self.record.seq().as_bytes();
                read_data.qualities[read_end.ix()] = self.record.qual().to_vec();
            } else {
                cigar.hard_to_soft();
            }
            let mut aln = Alignment::new(&self.record, cigar, read_end, Arc::clone(&self.contigs), f64::NAN);
            let aln_interval = aln.interval();
            let contig_len = self.contigs.get_len(aln_interval.contig_id());
            let central = self.boundary < aln_interval.end() && aln_interval.start() < contig_len - self.boundary;
            any_central |= central;

            let read_prof = aln.count_region_operations_fast(contig_len);
            let aln_prob = self.err_prof.ln_prob(&read_prof);
            let (edit_dist, read_len) = read_prof.edit_and_read_len();
            let (good_dist, passable_dist) = self.err_prof.allowed_edit_dist(read_len);
            let is_good_dist = edit_dist <= good_dist;
            // Discard all alignments with edit distance over half of the read length.
            let is_passable_dist = edit_dist <= passable_dist;
            any_good |= is_good_dist;

            write!(dbg_writer, "{:X}\t{}\t{}\t{}\t{}/{}\t{}\t{:.2}",
                name_hash, read_end, aln_interval, if central { 'T' } else { 'F' },
                edit_dist, read_len, if is_good_dist { '+' } else if is_passable_dist { '~' } else { '-' },
                Ln::to_log10(aln_prob)).map_err(add_path!(!))?;
            if primary {
                weight = self.contig_windows[aln_interval.contig_id().ix()].get_window_weight(aln_interval.middle());
                write!(dbg_writer, "\t{:.5}\t{}", weight, read_data.name).map_err(add_path!(!))?;
            }
            writeln!(dbg_writer).map_err(add_path!(!))?;

            let pos_key = u64::from(aln.interval().contig_id().get()) << 32 | u64::from(aln.interval().start());
            aln.set_ln_prob(aln_prob);
            if is_passable_dist {
                match self.found_alns.entry(pos_key) {
                    Entry::Occupied(entry) => {
                        let aln_ix = *entry.get();
                        // Already seen this alignment.
                        if aln_prob > alns[aln_ix].ln_prob() {
                            alns[aln_ix] = aln;
                        }
                    }
                    Entry::Vacant(entry) => {
                        entry.insert(alns.len());
                        alns.push(aln);
                    }
                }
                min_prob = aln_prob.min(min_prob);
            }
            if self.reader.read(&mut self.record).transpose()?.is_none() {
                self.has_more = false;
                break;
            }
            // Next primary record or unmapped read. Either way, it will be a new read or new read end.
            if is_primary(&self.record) {
                break;
            }
            assert_eq!(self.record.qname(), read_data.name.as_bytes(),
                "Read {} first alignment is not primary", String::from_utf8_lossy(self.record.qname()));
        }

        Ok(if any_good {
            MateSummary { name_hash, weight, min_prob, any_central, mapped: true }
        } else {
            alns.truncate(start_len);
            MateSummary::unmapped(name_hash)
        })
    }

    fn has_more(&self) -> bool {
        self.has_more
    }
}

/// Enum, that encodes two events, one of them can be missing, but not both.
pub enum Pair<T> {
    Both(T, T),
    First(T),
    Second(T),
}

impl<T> Pair<T> {
    /// Returns any of the values (first if availble, second otherwise).
    pub fn any(&self) -> &T {
        match self {
            Self::Both(a, _) | Self::First(a) => a,
            Self::Second(b) => b,
        }
    }

    /// Returns first element, if available.
    pub fn first(&self) -> Option<&T> {
        match self {
            Self::Both(a, _) | Self::First(a) => Some(a),
            _ => None,
        }
    }

    /// Returns second element, if available.
    pub fn second(&self) -> Option<&T> {
        match self {
            Self::Both(_, b) | Self::Second(b) => Some(b),
            _ => None,
        }
    }
}

pub struct PairAlignment {
    intervals: Pair<Interval>,
    ln_prob: f64,
}

impl PairAlignment {
    fn new(intervals: Pair<Interval>, ln_prob: f64) -> Self {
        Self { intervals, ln_prob }
    }

    pub fn intervals(&self) -> &Pair<Interval> {
        &self.intervals
    }

    /// Combined probability of the paired alignment.
    pub fn ln_prob(&self) -> f64 {
        self.ln_prob
    }

    pub fn contig_id(&self) -> ContigId {
        self.intervals.any().contig_id()
    }
}

/// Find all paired-read alignments to a specific contig for a single paired read.
/// All alignments that are worse than `best_prob - prob_diff` are discarded.
fn identify_contig_pair_alns(
    new_alns: &mut Vec<PairAlignment>,
    alns1: &[Alignment],
    min_prob1: f64,
    alns2: &[Alignment],
    min_prob2: f64,
    weight: f64,
    max_alns: usize,
    buffer: &mut Vec<f64>,
    insert_distr: &InsertDistr,
    params: &super::Params,
) {
    let start_len = new_alns.len();
    let insert_penalty = insert_distr.insert_penalty();
    // Factor, added to `aln1.ln_prob()` to get the cost of the pair `aln1, unmapped2`.
    let term1 = insert_penalty + min_prob2 + params.unmapped_penalty;

    buffer.clear();
    buffer.resize(alns2.len(), f64::NEG_INFINITY);
    for aln1 in alns1.iter() {
        let mut max_prob1 = f64::NEG_INFINITY;
        for (aln2, max_prob2) in alns2.iter().zip(buffer.iter_mut()) {
            let prob = aln1.paired_prob(aln2, insert_distr);
            max_prob1 = max_prob1.max(prob);
            *max_prob2 = max_prob2.max(prob);
            new_alns.push(PairAlignment::new(Pair::Both(aln1.interval().clone(), aln2.interval().clone()),
                weight * prob));
        }

        // Only add alignment `aln1,unmapped2` if it is better than any existing paired alignment to the same contig.
        let alone_prob1 = aln1.ln_prob() + term1;
        if alone_prob1 >= max_prob1 {
            new_alns.push(PairAlignment::new(Pair::First(aln1.interval().clone()), weight * alone_prob1));
        }
    }

    let term2 = insert_penalty + min_prob1 + params.unmapped_penalty;
    for (aln2, &max_prob2) in alns2.iter().zip(buffer.iter()) {
        let alone_prob2 = aln2.ln_prob() + term2;
        // Only add alignment `unmapped1,aln2` if it is better than existing `aln1,aln2` pairs.
        if alone_prob2 >= max_prob2 {
            new_alns.push(PairAlignment::new(Pair::Second(aln2.interval().clone()), weight * alone_prob2));
        }
    }

    let slice = &mut new_alns[start_len..];
    // Decreasing sort by ln-probability.
    slice.sort_unstable_by(|a, b| b.ln_prob.total_cmp(&a.ln_prob));
    let thresh_prob = slice[0].ln_prob - params.prob_diff;
    let keep_alns = slice[..slice.len().min(max_alns)].partition_point(|aln| aln.ln_prob >= thresh_prob);
    new_alns.truncate(start_len + keep_alns);
}

/// All appropriate alignments for one read pair / single read.
pub struct GrouppedAlignments {
    /// Read name.
    read_data: ReadData,
    /// Hash of the read name.
    name_hash: u64,
    /// Probability that both mates are unmapped.
    unmapped_prob: f64,
    /// Vector of alignments & their corresponding ln-probabilities, sorted by contig.
    alns: Vec<PairAlignment>,
}

impl GrouppedAlignments {
    /// Returns read name.
    pub fn read_name(&self) -> &str {
        &self.read_data.name
    }

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
        let i = bisect::left_by(&self.alns, |paln| paln.contig_id().cmp(&contig_id));
        let j = bisect::right_boundary(&self.alns, |paln| contig_id == paln.contig_id(), i, self.alns.len());
        (i, &self.alns[i..j])
    }

    /// Returns best alignment probability for each contig (must be sorted).
    pub fn best_for_each_contig<'a>(&'a self, contig_ids: &'a [ContigId]) -> impl Iterator<Item = f64> + 'a {
        let n = self.alns.len();
        let mut j = 0;
        contig_ids.iter().map(move |&id| {
            let i = j; // Do this way so that we can immediately return value two lines later.
            j = bisect::right_boundary(&self.alns, |aln| aln.contig_id() == id, i, n);
            if i == j { self.unmapped_prob } else { self.alns[i].ln_prob() }
        })
    }

    /// Return `i`-th alignment of all alignments for the read pair.
    pub fn ith_aln(&self, i: usize) -> &PairAlignment {
        &self.alns[i]
    }

    pub fn len(&self) -> usize {
        self.alns.len()
    }

    pub fn is_empty(&self) -> bool {
        self.alns.is_empty()
    }
}

#[inline]
fn store_by_prob(alns: &mut [Alignment]) {
    alns.sort_unstable_by(|a, b| b.ln_prob().total_cmp(&a.ln_prob()));
}

/// Store at most 3 alignments for cases when the read/read pair will not be used later.
const MAX_UNUSED_ALNS: usize = 3;
/// Store at most 10 alignments for cases when the read/read pair will not be used later.
/// NOTE: If we put usize::MAX here, we need to replace `i + max_alns` into `i.saturating_add(max_alns)`.
const MAX_USED_ALNS: usize = 10;

/// For a paired-end read, combine and pair all mate alignments across all contigs.
/// Input alignments are sorted first by contig, and then by read end.
/// Returns None if the read pair is ignored.
/// Returns groupped alignments and true if the alignments should be used for read assignment.
fn identify_paired_end_alignments(
    read_data: ReadData,
    alns: &mut [Alignment],
    summary1: &MateSummary,
    summary2: &MateSummary,
    buffer: &mut Vec<f64>,
    insert_distr: &InsertDistr,
    ln_ncontigs: f64,
    params: &super::Params,
) -> (GrouppedAlignments, bool)
{
    // If both mates available: mean weight, otherwise: half weight of the available alignment.
    let weight = 0.5 * (summary1.weight + summary2.weight);
    let use_pair = weight >= params.min_weight && (params.use_unpaired || (summary1.mapped && summary2.mapped));
    let max_alns = if use_pair { MAX_USED_ALNS } else { MAX_UNUSED_ALNS };

    let mut groupped_alns = Vec::new();
    let n = alns.len();
    alns.sort_unstable_by_key(Alignment::sort_key);
    let mut i = 0;
    while i < n {
        // For the current contig id, first mates will be in i..j, and second mates in j..k.
        let contig_id = alns[i].contig_id();
        // Sort key = contig_id * 2 + (first_mate ? 0 : 1).
        // sort_key1: current sort key for first mates, sort_key2: current sort key for second mates.
        let sort_key1 = contig_id.get() << 1;
        let sort_key2 = sort_key1 | 1;
        let j = if alns[i].read_end() == ReadEnd::First {
            bisect::right_boundary(alns, |aln| sort_key1 == aln.sort_key(), i, n)
        } else { i };
        let k = bisect::right_boundary(alns, |aln| sort_key2 == aln.sort_key(), j, n);

        store_by_prob(&mut alns[i..j]);
        store_by_prob(&mut alns[j..k]);
        identify_contig_pair_alns(&mut groupped_alns,
            &alns[i..j.min(i + max_alns)], summary1.min_prob,
            &alns[j..k.min(j + max_alns)], summary2.min_prob,
            weight, max_alns, buffer, insert_distr, params);
        i = k;
    }

    let both_unmapped = weight
        * (summary1.min_prob + summary2.min_prob + 2.0 * params.unmapped_penalty + insert_distr.insert_penalty());
    // Only for normalization, unmapped probability is multiplied by the number of contigs
    // because there is an unmapped possibility for every contig, which we do not store explicitely.
    let norm_fct = Ln::map_sum_init(&groupped_alns, PairAlignment::ln_prob, both_unmapped + ln_ncontigs);
    for aln in groupped_alns.iter_mut() {
        aln.ln_prob -= norm_fct;
    }
    (GrouppedAlignments {
        read_data,
        name_hash: summary1.name_hash,
        unmapped_prob: both_unmapped - norm_fct,
        alns: groupped_alns,
    }, use_pair)
}

/// For a single-end read, sort alignment across contigs, discard improbable alignments, and normalize probabilities.
/// Input alignments are sorted first by contig.
/// Returns groupped alignments and true if the alignments should be used for read assignment.
fn identify_single_end_alignments(
    read_data: ReadData,
    alns: &mut [Alignment],
    summary: &MateSummary,
    ln_ncontigs: f64,
    params: &super::Params,
) -> (GrouppedAlignments, bool)
{
    let weight = summary.weight;
    let use_read = weight >= params.min_weight;
    let max_alns = if use_read { MAX_USED_ALNS } else { MAX_UNUSED_ALNS };

    let mut groupped_alns = Vec::new();
    alns.sort_unstable_by_key(Alignment::contig_id);
    let n = alns.len();
    let mut i = 0;
    while i < n {
        let contig_id = alns[i].contig_id();
        let j = bisect::right_boundary(alns, |aln| contig_id == aln.contig_id(), i, n);
        store_by_prob(&mut alns[i..j]);
        let max_prob = alns[i].ln_prob();
        let thresh_prob = weight * max_prob - params.prob_diff;
        for aln in alns[i..j.min(i + max_alns)].iter() {
            let wprob = weight * aln.ln_prob();
            if wprob >= thresh_prob {
                groupped_alns.push(PairAlignment::new(Pair::First(aln.interval().clone()), wprob));
            } else {
                break;
            }
        }
        i = j;
    }

    let unmapped_prob = weight * (summary.min_prob + params.unmapped_penalty);
    // Only for normalization, unmapped probability is multiplied by the number of contigs
    // because there is an unmapped possibility for every contig, which we do not store explicitely.
    let norm_fct = Ln::map_sum_init(&groupped_alns, PairAlignment::ln_prob, unmapped_prob + ln_ncontigs);
    for aln in groupped_alns.iter_mut() {
        aln.ln_prob -= norm_fct;
    }
    (GrouppedAlignments {
        read_data,
        name_hash: summary.name_hash,
        unmapped_prob: unmapped_prob - norm_fct,
        alns: groupped_alns,
    }, use_read)
}

/// Paired-end/single-end read alignments, sorted by contig.
pub struct AllAlignments {
    /// Each element: one read pair.
    /// This vector stores reads with non-0 weight.
    reads: Vec<GrouppedAlignments>,
    /// This vector stores reads with 0 (or close to 0) weight.
    w0_reads: Vec<GrouppedAlignments>,
}

impl AllAlignments {
    /// Loads all paired-end or single-end alignments.
    ///
    /// Alignments are only kept if
    /// (i)  there is an alignment with small edit distance (under `err_prof.allowed_edit_dist`)  AND
    /// (ii) there is an alignment that overlaps an inner region of any contig (beyond the boundary region).
    pub fn load(
        reader: impl bam::Read,
        contigs: &Arc<ContigNames>,
        bg_distr: &BgDistr,
        contig_windows: &[ContigWindows],
        params: &super::Params,
        mut dbg_writer: impl Write,
    ) -> Result<Self, Error>
    {
        log::info!("    Loading read alignments");
        let boundary = params.boundary_size.checked_sub(params.tweak.unwrap()).unwrap();
        assert!(contigs.lengths().iter().all(|&len| len > 2 * boundary),
            "[{}] Some contigs are too short (must be over {})", contigs.tag(), 2 * boundary);

        writeln!(dbg_writer, "read_hash\tread_end\tinterval\tcentral\tedit_dist\tedit_status\
            \tlik\tweight\tread_name").map_err(add_path!(!))?;
        let mut reader = FilteredReader::new(reader, Arc::clone(contigs), bg_distr.error_profile(),
            contig_windows, boundary)?;

        let insert_distr = bg_distr.insert_distr();
        let is_paired_end = insert_distr.is_paired_end();
        let ln_ncontigs = (contigs.len() as f64).ln();
        let mut hashes = IntSet::default();
        let mut tmp_alns = Vec::new();
        let mut buffer = Vec::with_capacity(16);
        let mut collisions = 0;

        let mut total_reads = 0;
        let mut reads = Vec::new();
        let mut w0_reads = Vec::new();
        let mut read_data = ReadData::default();
        while reader.has_more() {
            total_reads += 1;
            tmp_alns.clear();
            let summary1 = reader.next_alns(ReadEnd::First, &mut tmp_alns, &mut read_data, &mut dbg_writer)?;
            if !hashes.insert(summary1.name_hash) {
                log::warn!("Read {} produced hash collision ({:X})", read_data.name, summary1.name_hash);
                collisions += 1;
            }
            let mut central = summary1.any_central;
            let summary2 = if is_paired_end {
                let summary2 = reader.next_alns(ReadEnd::Second, &mut tmp_alns, &mut read_data, &mut dbg_writer)?;
                central |= summary2.any_central;
                Some(summary2)
            } else {
                None
            };
            if !central {
                continue;
            }

            let (groupped_alns, use_read) = if is_paired_end {
                identify_paired_end_alignments(read_data.clone_and_reset(), &mut tmp_alns,
                    &summary1, summary2.as_ref().unwrap(), &mut buffer, insert_distr, ln_ncontigs, params)
            } else {
                identify_single_end_alignments(read_data.clone_and_reset(), &mut tmp_alns, &summary1,
                    ln_ncontigs, params)
            };
            if use_read {
                reads.push(groupped_alns);
            } else {
                w0_reads.push(groupped_alns);
            }
        }
        log::info!("    Total {} read{},  {} good,  {} low weight{},  {} out of bounds",
            total_reads, if is_paired_end { " pairs" } else { "s" }, reads.len(),
            w0_reads.len(), if is_paired_end && !params.use_unpaired { "/unpaired" } else { "" },
            total_reads - reads.len() - w0_reads.len());
        if collisions > 2 && collisions * 100 > total_reads {
            return Err(Error::RuntimeError(format!("Too many hash collisions ({})", collisions)));
        }
        Ok(Self { reads, w0_reads })
    }

    /// Reads that have non-0 weight.
    pub fn reads(&self) -> &[GrouppedAlignments] {
        &self.reads
    }

    /// Reads that have 0 or close to 0 weight.
    pub fn w0_reads(&self) -> &[GrouppedAlignments] {
        &self.w0_reads
    }

    /// Returns matrix `contig_id -> read pair -> best aln prob`.
    pub fn best_aln_matrix(&self, contig_ids: &[ContigId]) -> Vec<Vec<f64>> {
        let mut best_aln_probs = vec![Vec::<f64>::with_capacity(self.reads().len()); contig_ids.len()];
        for read_alns in self.reads.iter() {
            for (contig_probs, best_prob) in best_aln_probs.iter_mut()
                    .zip(read_alns.best_for_each_contig(contig_ids)) {
                contig_probs.push(best_prob);
            }
        }
        best_aln_probs
    }
}
