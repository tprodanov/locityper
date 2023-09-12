use std::{
    cmp::min,
    sync::Arc,
    io::{self, Write},
    collections::hash_map::Entry,
};
use htslib::bam;
use nohash::{IntSet, IntMap};
use once_cell::sync::OnceCell;
use crate::{
    err::{Error, add_path},
    seq::{
        self,
        Interval, ContigId, ContigNames,
        aln::{Alignment, ReadEnd, Strand},
        cigar::Cigar,
    },
    bg::{
        BgDistr,
        err_prof::ErrorProfile,
        insertsz::InsertDistr,
    },
    algo::{bisect, TwoU32, get_hash},
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

#[derive(Clone)]
pub struct MateData {
    strand: Strand,
    sequence: Vec<u8>,
    opp_sequence: OnceCell<Vec<u8>>,
    qualities: Vec<u8>,
    opp_qualities: OnceCell<Vec<u8>>,
}

impl MateData {
    fn new(record: &bam::Record) -> Self {
        Self {
            strand: Strand::from_record(record),
            sequence: record.seq().as_bytes(),
            opp_sequence: OnceCell::new(),
            // Need to do this as otherwise htslib will panic.
            qualities: if record.qual().is_empty() { vec![255; record.seq().len()] } else { record.qual().to_vec() },
            opp_qualities: OnceCell::new(),
        }
    }

    /// Returns sequence and quality for the corresponding strand.
    pub fn get_seq_and_qual(&self, strand: Strand) -> (&[u8], &[u8]) {
        if strand == self.strand {
            (&self.sequence, &self.qualities)
        } else {
            (
                self.opp_sequence.get_or_init(|| seq::reverse_complement(&self.sequence)),
                self.opp_qualities.get_or_init(|| self.qualities.iter().rev().copied().collect()),
            )
        }
    }
}

/// Read information: read name, two sequences and qualities (for both mates).
#[derive(Default, Clone)]
pub struct ReadData {
    name: String,
    mates: [Option<MateData>; 2],
}

impl ReadData {
    pub fn name(&self) -> &String {
        &self.name
    }

    pub fn mate_data(&self, read_end: ReadEnd) -> &MateData {
        self.mates[read_end.ix()].as_ref().expect("Mate data undefined")
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
    /// Key (contig id, alignment start), value: index of the previous position.
    found_alns: IntMap<TwoU32, usize>,
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
        let name_hash = get_hash(self.record.qname());
        let name = std::str::from_utf8(self.record.qname()).map_err(|_| Error::InvalidInput(
            format!("Read name is not UTF-8: {:?}", String::from_utf8_lossy(self.record.qname()))))?;
        if read_end == ReadEnd::First {
            assert!(read_data.name.is_empty());
            read_data.name.push_str(name);
        } else if name != read_data.name {
            return Err(Error::InvalidData(format!("Read {} does not have a second read end", name)));
        }
        if self.record.seq().is_empty() {
            return Err(Error::InvalidData(format!("Read {} does not have read sequence", read_data.name)));
        }
        let old_option = read_data.mates[read_end.ix()].replace(MateData::new(&self.record));
        assert!(old_option.is_none(), "Mate data defined twice");

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

            aln.set_ln_prob(aln_prob);
            if is_passable_dist {
                match self.found_alns.entry(TwoU32(aln.contig_id().get().into(), aln.interval().start())) {
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

#[derive(Clone)]
pub struct PairAlignment {
    ln_prob: f64,
    contig_id: ContigId,
    /// Index (across all alignment of this read) and middle of the first alignment (if mapped).
    aln1: Option<(u32, u32)>,
    /// Index and middle of the second alignment (if mapped).
    aln2: Option<(u32, u32)>,
}

impl PairAlignment {
    fn new_first(ix: usize, interv: &Interval, ln_prob: f64) -> Self {
        Self {
            ln_prob,
            contig_id: interv.contig_id(),
            aln1: Some((ix as u32, interv.middle())),
            aln2: None,
        }
    }

    fn new_second(ix: usize, interv: &Interval, ln_prob: f64) -> Self {
        Self {
            ln_prob,
            contig_id: interv.contig_id(),
            aln1: None,
            aln2: Some((ix as u32, interv.middle())),
        }
    }

    fn new_both(ix1: usize, interv1: &Interval, ix2: usize, interv2: &Interval, ln_prob: f64) -> Self {
        assert_eq!(interv1.contig_id(), interv2.contig_id());
        Self {
            ln_prob,
            contig_id: interv1.contig_id(),
            aln1: Some((ix1 as u32, interv1.middle())),
            aln2: Some((ix2 as u32, interv2.middle())),
        }
    }

    /// Combined probability of the paired alignment.
    pub fn ln_prob(&self) -> f64 {
        self.ln_prob
    }

    pub fn contig_id(&self) -> ContigId {
        self.contig_id
    }

    /// Returns the middle of the first alignment.
    pub fn middle1(&self) -> Option<u32> {
        self.aln1.map(|(_, middle)| middle)
    }

    /// Returns the middle of the second alignment.
    pub fn middle2(&self) -> Option<u32> {
        self.aln2.map(|(_, middle)| middle)
    }

    /// Returns the index of the first alignment (across all alignments for this read pair).
    pub fn ix1(&self) -> Option<u32> {
        self.aln1.map(|(ix, _)| ix)
    }

    /// Returns the index of the second alignment (across all alignments for this read pair).
    pub fn ix2(&self) -> Option<u32> {
        self.aln2.map(|(ix, _)| ix)
    }
}

/// Groups first end alignments `alignments[i..j]` and second end alignments `alignments[j..]` into `aln_pairs`.
fn identify_contig_pair_alns(
    alignments: &[Alignment],
    i: usize,
    j: usize,
    aln_pairs: &mut Vec<PairAlignment>,
    min_prob1: f64,
    min_prob2: f64,
    max_alns: usize,
    buffer: &mut Vec<f64>,
    insert_distr: &InsertDistr,
    params: &super::Params,
) {
    let start_len = aln_pairs.len();
    let insert_penalty = insert_distr.insert_penalty();
    // Factor, added to `aln1.ln_prob()` to get the cost of the pair `aln1, unmapped2`.
    let term1 = insert_penalty + min_prob2 + params.unmapped_penalty;

    let k = alignments.len();
    assert!(i <= j && j <= k);
    buffer.clear();
    buffer.resize(k - j, f64::NEG_INFINITY);
    for ix1 in i..j {
        let aln1 = unsafe { alignments.get_unchecked(ix1) };
        let mut max_prob1 = f64::NEG_INFINITY;
        for (ix2, max_prob2) in (j..k).zip(buffer.iter_mut()) {
            let aln2 = unsafe { alignments.get_unchecked(ix2) };
            let prob = aln1.paired_prob(aln2, insert_distr);
            max_prob1 = max_prob1.max(prob);
            *max_prob2 = max_prob2.max(prob);
            aln_pairs.push(PairAlignment::new_both(ix1, aln1.interval(), ix2, aln2.interval(), prob));
        }

        // Only add alignment `aln1,unmapped2` if it is better than any existing paired alignment to the same contig.
        let alone_prob1 = aln1.ln_prob() + term1;
        if alone_prob1 >= max_prob1 {
            aln_pairs.push(PairAlignment::new_first(ix1, aln1.interval(), alone_prob1));
        }
    }

    let term2 = insert_penalty + min_prob1 + params.unmapped_penalty;
    for (ix2, &max_prob2) in (j..k).zip(buffer.iter()) {
        let aln2 = unsafe { alignments.get_unchecked(ix2) };
        let alone_prob2 = aln2.ln_prob() + term2;
        // Only add alignment `unmapped1,aln2` if it is better than existing `aln1,aln2` pairs.
        if alone_prob2 >= max_prob2 {
            aln_pairs.push(PairAlignment::new_second(ix2, aln2.interval(), alone_prob2));
        }
    }

    let slice = &mut aln_pairs[start_len..];
    // Decreasing sort by ln-probability.
    slice.sort_unstable_by(|a, b| b.ln_prob.total_cmp(&a.ln_prob));
    let thresh_prob = slice[0].ln_prob - params.prob_diff;
    let keep_alns = slice[..slice.len().min(max_alns)].partition_point(|aln| aln.ln_prob >= thresh_prob);
    aln_pairs.truncate(start_len + keep_alns);
}

/// All appropriate alignments for one read pair / single read.
pub struct GrouppedAlignments {
    /// Read name.
    read_data: ReadData,
    /// Hash of the read name.
    name_hash: u64,
    /// Probability that both mates are unmapped.
    unmapped_prob: f64,
    /// All read alignments.
    alignments: Vec<Alignment>,
    /// Read alignments, groupped into pairs.
    aln_pairs: Vec<PairAlignment>,
    /// Weight of the read pair.
    weight: f64,
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

    /// For a given contig, returns alignment pairs, corresponding to this contig.
    pub fn contig_aln_pairs(&self, contig_id: ContigId) -> &[PairAlignment] {
        let i = bisect::left_by(&self.aln_pairs, |paln| paln.contig_id().cmp(&contig_id));
        let j = bisect::right_boundary(&self.aln_pairs, |paln| contig_id == paln.contig_id(), i, self.aln_pairs.len());
        &self.aln_pairs[i..j]
    }

    /// Returns best alignment probability for each contig (must be sorted).
    pub fn best_for_each_contig<'a>(&'a self, contig_ids: &'a [ContigId]) -> impl Iterator<Item = f64> + 'a {
        let n = self.aln_pairs.len();
        let mut j = 0;
        contig_ids.iter().map(move |&id| {
            let i = j; // Do this way so that we can immediately return value two lines later.
            j = bisect::right_boundary(&self.aln_pairs, |aln| aln.contig_id() == id, i, n);
            if i == j { self.unmapped_prob } else { self.aln_pairs[i].ln_prob() }
        })
    }

    /// Return `i`-th alignment of all alignments for the read pair.
    pub fn ith_aln(&self, i: u32) -> &Alignment {
        &self.alignments[i as usize]
    }

    /// Returns the total number of alignment pairs.
    pub fn alignment_pairs(&self) -> usize {
        self.aln_pairs.len()
    }

    /// Returns read information, such as sequences, qualities, and read name.
    pub fn read_data(&self) -> &ReadData {
        &self.read_data
    }

    /// Returns the weight of the read pair.
    pub fn weight(&self) -> f64 {
        self.weight
    }

    fn write_read_pair_info(&self, f: &mut impl Write, contigs: &ContigNames) -> io::Result<()> {
        for pair in self.aln_pairs.iter() {
            write!(f, "{:X}\t{}\t", self.name_hash, contigs.get_name(pair.contig_id))?;
            if let Some((i, _)) = pair.aln1 {
                write!(f, "{}\t", self.alignments[i as usize].interval().start() + 1)?;
            } else {
                write!(f, "*\t")?;
            }
            if let Some((j, _)) = pair.aln2 {
                write!(f, "{}\t", self.alignments[j as usize].interval().start() + 1)?;
            } else {
                write!(f, "*\t")?;
            }
            writeln!(f, "{:.4}", Ln::to_log10(pair.ln_prob))?;
        }
        Ok(())
    }
}

/// Store at most 3 alignments for cases when the read/read pair will not be used later.
const MAX_UNUSED_ALNS: usize = 2;
/// Store at most 10 alignments for cases when the read/read pair will not be used later.
/// NOTE: If we put usize::MAX here, we need to replace `i + max_alns` into `i.saturating_add(max_alns)`.
const MAX_USED_ALNS: usize = 10;

/// For a paired-end read, combine and pair all mate alignments across all contigs.
/// Input alignments are sorted first by contig, and then by read end.
/// Returns None if the read pair is ignored.
/// Returns groupped alignments and true if the alignments should be used for read assignment.
fn identify_paired_end_alignments(
    read_data: ReadData,
    tmp_alns: &mut Vec<Alignment>,
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

    // First: by contig (decr.), second: by read end (decr.), third: by aln probability (incr.).
    tmp_alns.sort_unstable_by(|a, b|
        (b.contig_id(), b.read_end(), a.ln_prob()).partial_cmp(&(a.contig_id(), a.read_end(), b.ln_prob())).unwrap());
    let mut alignments = Vec::new();
    let mut aln_pairs = Vec::new();

    let mut curr_contig = ContigId::new(0);
    // First-end alignments for this contig are stored in i..j,
    let mut i = 0;
    // Second-end alignments are stored in j..k.
    // j = MAX if there were no second-end alignments.
    let mut j = usize::MAX;
    while let Some(aln) = tmp_alns.pop() {
        let k = alignments.len();
        if curr_contig != aln.contig_id() {
            if i < k {
                // There are some alignments that need to be saved.
                identify_contig_pair_alns(&alignments, i, min(j, k), &mut aln_pairs,
                    summary1.min_prob, summary2.min_prob, max_alns, buffer, insert_distr, params);
            }
            curr_contig = aln.contig_id();
            i = k;
            j = usize::MAX;
        }
        if aln.read_end() == ReadEnd::First {
            if k - i < max_alns {
                alignments.push(aln);
            }
        } else {
            j = min(j, k);
            if k - j < max_alns {
                alignments.push(aln);
            }
        }
    }
    let k = alignments.len();
    if i < k {
        // There are some alignments that need to be saved.
        identify_contig_pair_alns(&alignments, i, min(j, k), &mut aln_pairs,
            summary1.min_prob, summary2.min_prob, max_alns, buffer, insert_distr, params);
    }

    let both_unmapped = summary1.min_prob + summary2.min_prob
        + 2.0 * params.unmapped_penalty + insert_distr.insert_penalty();
    // Only for normalization, unmapped probability is multiplied by the number of contigs
    // because there is an unmapped possibility for every contig, which we do not store explicitely.
    let norm_fct = Ln::map_sum_init(&aln_pairs, PairAlignment::ln_prob, both_unmapped + ln_ncontigs);
    for aln in aln_pairs.iter_mut() {
        aln.ln_prob = weight * (aln.ln_prob - norm_fct);
    }
    (GrouppedAlignments {
        read_data, weight,
        name_hash: summary1.name_hash,
        unmapped_prob: weight * (both_unmapped - norm_fct),
        alignments, aln_pairs,
    }, use_pair)
}

/// For a single-end read, sort alignment across contigs, discard improbable alignments, and normalize probabilities.
/// Input alignments are sorted first by contig.
/// Returns groupped alignments and true if the alignments should be used for read assignment.
fn identify_single_end_alignments(
    read_data: ReadData,
    tmp_alns: &mut Vec<Alignment>,
    summary: &MateSummary,
    ln_ncontigs: f64,
    params: &super::Params,
) -> (GrouppedAlignments, bool)
{
    let weight = summary.weight;
    let use_read = weight >= params.min_weight;
    let max_alns = if use_read { MAX_USED_ALNS } else { MAX_UNUSED_ALNS };

    let mut alignments = Vec::new();
    let mut aln_pairs = Vec::new();
    // First: by contig (decreasing), then decreasing by ln-probability (increasing).
    tmp_alns.sort_unstable_by(|a, b| (b.contig_id(), a.ln_prob()).partial_cmp(&(a.contig_id(), b.ln_prob())).unwrap());

    let mut curr_contig = None;
    let mut thresh_prob = f64::NAN;
    let mut curr_saved = 0;
    // Go from end to start.
    while let Some(aln) = tmp_alns.pop() {
        if curr_contig != Some(aln.contig_id()) {
            curr_contig = Some(aln.contig_id());
            thresh_prob = aln.ln_prob() - params.prob_diff;
            curr_saved = 0;
        }
        if aln.ln_prob() >= thresh_prob && curr_saved < max_alns {
            aln_pairs.push(PairAlignment::new_first(alignments.len(), aln.interval(), aln.ln_prob()));
            alignments.push(aln);
            curr_saved += 1;
        }
    }

    let unmapped_prob = summary.min_prob + params.unmapped_penalty;
    // Only for normalization, unmapped probability is multiplied by the number of contigs
    // because there is an unmapped possibility for every contig, which we do not store explicitely.
    let norm_fct = Ln::map_sum_init(&aln_pairs, PairAlignment::ln_prob, unmapped_prob + ln_ncontigs);
    for aln in aln_pairs.iter_mut() {
        aln.ln_prob = weight * (aln.ln_prob - norm_fct);
    }
    (GrouppedAlignments {
        read_data, weight,
        name_hash: summary.name_hash,
        unmapped_prob: weight * (unmapped_prob - norm_fct),
        alignments, aln_pairs,
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
        while reader.has_more() {
            let mut read_data = ReadData::default();
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
                identify_paired_end_alignments(read_data, &mut tmp_alns, &summary1, summary2.as_ref().unwrap(),
                    &mut buffer, insert_distr, ln_ncontigs, params)
            } else {
                identify_single_end_alignments(read_data, &mut tmp_alns, &summary1, ln_ncontigs, params)
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

    /// Writes debug information about all read pairs.
    pub fn write_read_pair_info(&self, mut f: impl Write, contigs: &ContigNames, include_w0: bool) -> io::Result<()> {
        writeln!(f, "read_hash\tcontig\tpos1\tpos2\tlik")?;
        for read in self.reads.iter() {
            read.write_read_pair_info(&mut f, contigs)?;
        }
        if include_w0 {
            for read in self.w0_reads.iter() {
                read.write_read_pair_info(&mut f, contigs)?;
            }
        }
        Ok(())
    }
}
