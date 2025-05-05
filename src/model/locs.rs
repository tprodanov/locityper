use std::{
    cmp::min,
    fmt::{self, Write as FmtWrite},
    sync::{Arc, OnceLock},
    io::{self, Write},
};
use htslib::bam;
use crate::{
    err::{Error, error, add_path},
    seq::{
        self,
        Interval, ContigId, ContigNames, ContigSet,
        kmers::{self, CANONICAL},
        aln::{Alignment, ReadEnd, Strand},
        cigar::Cigar,
    },
    bg::{
        BgDistr,
        err_prof::{ErrorProfile, EditDistCache},
        insertsz::InsertDistr,
    },
    algo::{bisect, TwoU32, get_hash, HashSet, IntSet, IntMap, hash_map::Entry},
    math::Ln,
    model::windows::ContigInfos,
    ext::{
        sys::GzFile,
        vec::VecExt,
    },
};

// ------------------------- Read-end and pair-end data, such as sequences and qualities -------------------------

#[derive(Default, Clone, Copy, PartialEq, Eq)]
pub struct NameHash(u64);

impl NameHash {
    pub fn new(name: &[u8]) -> Self {
        Self(get_hash(name))
    }
}

impl fmt::Display for NameHash {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        use base64ct::Encoding;
        // Each base64 letter encodes 6 bits, therefore 64-bit number will take 11 letters.
        let mut buf = [0_u8; 11];
        f.write_str(base64ct::Base64UrlUnpadded::encode(&self.0.to_le_bytes(), &mut buf).unwrap())
    }
}

impl std::hash::Hash for NameHash {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        hasher.write_u64(self.0)
    }
}

impl nohash::IsEnabled for NameHash {}

#[derive(Clone)]
pub struct MateData {
    strand: Strand,
    sequence: Vec<u8>,
    opp_sequence: OnceLock<Vec<u8>>,
    qualities: Vec<u8>,
    opp_qualities: OnceLock<Vec<u8>>,
    unique_kmers: u16,
}

impl MateData {
    fn new(record: &bam::Record) -> Self {
        Self {
            strand: Strand::from_record(record),
            sequence: record.seq().as_bytes(),
            opp_sequence: OnceLock::new(),
            // Need to do this as otherwise htslib will panic.
            qualities: if record.qual().is_empty() { vec![255; record.seq().len()] } else { record.qual().to_vec() },
            opp_qualities: OnceLock::new(),
            unique_kmers: 0,
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

    /// Returns the number of unique k-mers in the read.
    pub fn unique_kmers(&self) -> u16 {
        self.unique_kmers
    }
}

/// Read information: read name, two sequences and qualities (for both mates).
#[derive(Default, Clone)]
pub struct ReadData {
    name: String,
    name_hash: NameHash,
    mates: [Option<MateData>; 2],
}

impl ReadData {
    pub fn name(&self) -> &String {
        &self.name
    }

    pub fn name_hash(&self) -> NameHash {
        self.name_hash
    }

    pub fn mate_data(&self, read_end: ReadEnd) -> &MateData {
        self.mates[read_end.ix()].as_ref().expect("Mate data undefined")
    }
}

// ------------------------- BAM file reader, that skips bad alignments -------------------------

/// BAM reader.
/// Contains `next_alns` method, which consecutively reads all records with the same name
/// until the next primary alignment.
struct FilteredReader<'a, R: bam::Read> {
    reader: R,
    /// Next record is stored but not yet returned.
    record: bam::Record,
    has_more: bool,
    contigs: Arc<ContigNames>,
    err_prof: &'a ErrorProfile,
    are_short_reads: bool,
    edit_dist_cache: &'a EditDistCache,
    /// Buffer, used for discrard alignments with the same positions.
    /// Key (contig id, alignment start), value: index of the previous position.
    found_alns: IntMap<TwoU32, usize>,
    buffer_set: IntSet<u16>,
}

#[inline]
fn is_primary(record: &bam::Record) -> bool {
    (record.flags() & 2304) == 0
}

impl<'a, R: bam::Read> FilteredReader<'a, R> {
    fn new(
        mut reader: R,
        contigs: Arc<ContigNames>,
        bg_distr: &'a BgDistr,
        edit_dist_cache: &'a EditDistCache,
    ) -> crate::Result<Self> {
        let mut record = bam::Record::new();
        // Reader would return None if there are no more records.
        let has_more = reader.read(&mut record).transpose()?.is_some();
        assert!(is_primary(&record), "First record in the BAM file is secondary/supplementary");
        Ok(FilteredReader {
            found_alns: IntMap::default(),
            buffer_set: IntSet::default(),
            are_short_reads: bg_distr.seq_info().technology().are_short_reads(),
            err_prof: bg_distr.error_profile(),
            reader, record, contigs, edit_dist_cache, has_more,
        })
    }

    /// Starting with `self.record` (already loaded), reads all alignments,
    /// corresponding to this read and current read end.
    /// Basically, the function continue to read until the next primary alignment is found,
    /// which is saved to `self.record`.
    ///
    /// If read is unmapped, or best edit distance is too high, does not add any new alignments.
    ///
    /// Multiplies input weight by the weight of the best alignment.
    /// Returns true if any of the alignments were good, otherwise the input vector did not change.
    fn next_alns(
        &mut self,
        read_end: ReadEnd,
        alns: &mut Vec<Alignment>,
        weight: &mut f64,
        read_data: &mut ReadData,
        dbg_writer: &mut impl Write,
    ) -> crate::Result<bool>
    {
        assert!(self.has_more, "Cannot read any more records from a BAM file");
        let name = std::str::from_utf8(self.record.qname())
            .map_err(|_| Error::Utf8("read name", self.record.qname().to_vec()))?;
        let name_hash = NameHash::new(self.record.qname());

        // Check if everything is correct.
        if read_end == ReadEnd::First {
            assert!(read_data.name.is_empty());
            read_data.name.push_str(name);
            read_data.name_hash = name_hash;
        } else if name_hash != read_data.name_hash {
            return Err(error!(InvalidData, "Read {} does not have a second read end", name));
        }
        if self.record.seq().is_empty() {
            if self.record.qname().is_empty() {
                return Err(error!(InvalidData,
                    "Alignment file contains absolutely empty read (no read name or sequence)"));
            }
            return Err(error!(InvalidData, "Read {} does not have read sequence", read_data.name));
        }
        let old_option = read_data.mates[read_end.ix()].replace(MateData::new(&self.record));
        assert!(old_option.is_none(), "Mate data defined twice");

        if self.record.is_unmapped() {
            writeln!(dbg_writer, "{}\t{}\t*\tNA\t-\t{}", name_hash, read_end, name).map_err(add_path!(!))?;
            // Read the next record, and save false to `has_more` if there are no more records left.
            self.has_more = self.reader.read(&mut self.record).transpose()?.is_some();
            return Ok(false);
        }

        let start_len = alns.len();
        // Does any alignment have good edit distance?
        let mut best_edit = u32::MAX;
        let mut max_prob = f64::NEG_INFINITY;
        self.found_alns.clear();
        let mut primary = true;

        loop {
            let mut cigar = Cigar::from_raw(self.record.raw_cigar());
            if primary {
                assert!(!cigar.has_hard_clipping(), "Primary alignment has hard clipping");
            } else {
                cigar.hard_to_soft();
            }
            let mut aln = Alignment::new(&self.record, cigar, read_end, Arc::clone(&self.contigs), f64::NAN);
            let contig_len = self.contigs.get_len(aln.contig_id());
            let read_prof = aln.count_region_operations_fast(contig_len);
            let dist = read_prof.edit_distance();
            best_edit = best_edit.min(dist.edit());
            let aln_prob = self.err_prof.ln_prob(&read_prof);
            max_prob = aln_prob.max(max_prob);
            aln.set_distance(dist);
            aln.set_ln_prob(aln_prob);

            write!(dbg_writer, "{}\t{}\t{}\t{}\t{:.2}", name_hash, read_end, aln.interval(),
                dist, Ln::to_log10(aln_prob)).map_err(add_path!(!))?;
            if primary {
                write!(dbg_writer, "\t{}", read_data.name).map_err(add_path!(!))?;
                primary = false;
            }
            writeln!(dbg_writer).map_err(add_path!(!))?;

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

        let read_len = alns[start_len].distance().unwrap().read_len();
        let (good_dist, passable_dist) = self.edit_dist_cache.get(read_len);
        // TODO: Use parameter.
        let threshold_dist = if self.are_short_reads { read_len / 2 } else { good_dist };
        if best_edit > threshold_dist {
            alns.truncate(start_len);
            return Ok(false);
        }

        alns[start_len..].iter_mut().for_each(|aln| aln.set_ln_prob(aln.ln_prob() - max_prob));
        // TODO: Use parameter.
        let prob_thresh = Ln::from_log10(-20.0);
        if self.are_short_reads {
            VecExt::unstable_retain(alns, start_len, |aln| aln.ln_prob() >= prob_thresh);
        } else {
            VecExt::unstable_retain(alns, start_len, |aln| aln.distance().unwrap().edit() >= passable_dist);
        }

        *weight *= if best_edit <= good_dist { 1.0 } else { (f64::from(good_dist) / f64::from(best_edit)).sqrt() };
        Ok(true)
    }

    fn has_more(&self) -> bool {
        self.has_more
    }
}

// ------------------------- Structures for storing all alignments for one read pair -------------------------

/// All appropriate alignments for one read pair / single read.
pub struct GrouppedAlignments {
    /// Read name.
    read_data: ReadData,
    /// Read weight.
    weight: f64,
    /// Highest likelihood across all possible locations.
    max_lik: f64,
    /// Probability that both mates are unmapped.
    unmapped_prob: f64,
    /// All read alignments.
    alignments: Vec<Alignment>,
    /// Read alignments, groupped into pairs.
    aln_pairs: Vec<PairAlignment>,
}

impl GrouppedAlignments {
    /// Returns read name.
    #[inline(always)]
    pub fn read_name(&self) -> &str {
        &self.read_data.name
    }

    /// Read pair weight.
    #[inline(always)]
    pub fn weight(&self) -> f64 {
        self.weight
    }

    /// Maximum likelihood across all contigs.
    #[inline(always)]
    pub fn max_lik(&self) -> f64 {
        self.max_lik
    }

    /// Probability that both reads are unmapped for one specific contig (but same for all contigs).
    #[inline(always)]
    pub fn unmapped_prob(&self) -> f64 {
        self.unmapped_prob
    }

    /// Returns the highest probability at given contig.
    pub fn best_at_contig(&self, contig_id: ContigId) -> f64 {
        let i = bisect::left_by(&self.aln_pairs, |paln| paln.contig_id().cmp(&contig_id));
        if let Some(aln) = self.aln_pairs.get(i) {
            if aln.contig_id() == contig_id {
                return aln.ln_prob();
            }
        }
        self.unmapped_prob
    }

    /// For a given contig, returns alignment pairs, corresponding to this contig.
    pub fn contig_alns(&self, contig_id: ContigId) -> &[PairAlignment] {
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
    #[inline(always)]
    pub fn ith_aln(&self, i: u32) -> &Alignment {
        &self.alignments[i as usize]
    }

    /// Returns the total number of alignment pairs.
    #[inline(always)]
    pub fn alignment_pairs(&self) -> usize {
        self.aln_pairs.len()
    }

    /// Returns read information, such as sequences, qualities, and read name.
    pub fn read_data(&self) -> &ReadData {
        &self.read_data
    }

    fn write_read_pair_info(&self, f: &mut impl Write, contigs: &ContigNames) -> io::Result<()> {
        for pair in self.aln_pairs.iter() {
            write!(f, "{}\t{}\t", self.read_data.name_hash, contigs.get_name(pair.contig_id))?;
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
        writeln!(f, "{}\t*\t*\t*\t{:.4}", self.read_data.name_hash, Ln::to_log10(self.unmapped_prob))?;
        Ok(())
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

// ---------------------- Loading alignments, groupping them by contigs and organizing into pairs ----------------------

/// Store at most 3 alignments for cases when the read/read pair will not be used later.
const MAX_UNUSED_ALNS: usize = 2;
/// Store at most 10 alignments for cases when the read/read pair will not be used later.
/// NOTE: If we put usize::MAX here, we need to replace `i + max_alns` into `i.saturating_add(max_alns)`.
const MAX_USED_ALNS: usize = 10;

/// Groups first end alignments `alignments[i..j]` and second end alignments `alignments[j..]` into `aln_pairs`.
fn identify_contig_pair_alns(
    alignments: &[Alignment],
    i: usize,
    j: usize,
    aln_pairs: &mut Vec<PairAlignment>,
    max_alns: usize,
    buffer: &mut Vec<f64>,
    insert_distr: &InsertDistr,
    params: &super::Params,
) {
    let start_len = aln_pairs.len();
    let insert_penalty = insert_distr.insert_penalty();
    let one_unmapped = insert_penalty + params.unmapped_penalty;

    let k = alignments.len();
    assert!(i <= j && j <= k);
    buffer.clear();
    buffer.resize(k - j, f64::NEG_INFINITY);
    for ix1 in i..j {
        let aln1 = unsafe { alignments.get_unchecked(ix1) };
        let mut max_prob1 = f64::NEG_INFINITY;
        for (ix2, max_prob2) in (j..k).zip(buffer.iter_mut()) {
            let aln2 = unsafe { alignments.get_unchecked(ix2) };
            if aln1.strand() != aln2.strand() {
                let prob = aln1.paired_prob(aln2, insert_distr);
                if prob.is_finite() {
                    max_prob1 = max_prob1.max(prob);
                    *max_prob2 = max_prob2.max(prob);
                    aln_pairs.push(PairAlignment::new_both(ix1, aln1.interval(), ix2, aln2.interval(), prob));
                }
            }
        }

        // Only add alignment `aln1,unmapped2` if it is better than any existing paired alignment to the same contig.
        let alone_prob1 = aln1.ln_prob() + one_unmapped;
        if alone_prob1 >= max_prob1 {
            aln_pairs.push(PairAlignment::new_first(ix1, aln1.interval(), alone_prob1));
        }
    }

    for (ix2, &max_prob2) in (j..k).zip(buffer.iter()) {
        let aln2 = unsafe { alignments.get_unchecked(ix2) };
        let alone_prob2 = aln2.ln_prob() + one_unmapped;
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

/// For a paired-end read, combine and pair all mate alignments across all contigs.
/// Input alignments are sorted first by contig, and then by read end.
/// Returns None if the read pair is ignored.
/// Returns groupped alignments and true if the alignments should be used for read assignment.
fn identify_paired_end_alignments(
    read_data: ReadData,
    tmp_alns: &mut Vec<Alignment>,
    mut weight: f64,
    max_alns: usize,
    buffer: &mut Vec<f64>,
    insert_distr: &InsertDistr,
    contig_infos: &ContigInfos,
    params: &super::Params,
) -> GrouppedAlignments
{
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
                identify_contig_pair_alns(&alignments, i, min(j, k), &mut aln_pairs, max_alns,
                    buffer, insert_distr, params);
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
        identify_contig_pair_alns(&alignments, i, min(j, k), &mut aln_pairs, max_alns,
            buffer, insert_distr, params);
    }

    weight *= contig_infos.explicit_read_weight(&aln_pairs);
    let mut max_lik = f64::NEG_INFINITY;
    for aln in aln_pairs.iter_mut() {
        aln.ln_prob *= weight;
        max_lik = max_lik.max(aln.ln_prob);
    }
    GrouppedAlignments {
        read_data, max_lik, alignments, aln_pairs, weight,
        unmapped_prob: weight * (2.0 * params.unmapped_penalty + insert_distr.insert_penalty()),
    }
}

/// For a single-end read, sort alignment across contigs, discard improbable alignments, and normalize probabilities.
/// Input alignments are sorted first by contig.
/// Returns groupped alignments and true if the alignments should be used for read assignment.
fn identify_single_end_alignments(
    read_data: ReadData,
    tmp_alns: &mut Vec<Alignment>,
    mut weight: f64,
    max_alns: usize,
    contig_infos: &ContigInfos,
    params: &super::Params,
) -> GrouppedAlignments
{
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

    weight *= contig_infos.explicit_read_weight(&aln_pairs);
    let mut max_lik = f64::NEG_INFINITY;
    for aln in aln_pairs.iter_mut() {
        aln.ln_prob *= weight;
        max_lik = max_lik.max(aln.ln_prob);
    }
    GrouppedAlignments {
        read_data, max_lik, alignments, aln_pairs, weight,
        unmapped_prob: weight * params.unmapped_penalty,
    }
}

// ------------------- Checks reads for the number of unique non-overlapping k-mers -------------------

struct UniqueKmers {
    k: u8,
    /// k - 2
    k_2: usize,
    unique_kmers: HashSet<u128>,
    kmers_buf: Vec<u128>,
    /// Calculate read weight as (x - Th + 1) / (Ts - Th + 1).
    /// Weight multiplier = 1 / (Ts - Th + 1),
    weight_mult: f64,
    /// Intercept = (1 - Th) * weight_mult.
    weight_interc: f64,
}

impl UniqueKmers {
    /// Stores all k-mers, unique to the current locus.
    fn new(
        contig_set: &ContigSet,
        hard_threshold: u16,
        soft_threshold: u16,
    ) -> Self {
        let kmer_counts = contig_set.kmer_counts();
        let k = u8::try_from(kmer_counts.k()).unwrap();
        assert!(k > 1);
        let mut kmers_buf = Vec::new();
        let mut unique_kmers = HashSet::default();
        let mut off_target_kmers = HashSet::default();
        for (seq, counts) in contig_set.seqs().iter().zip(kmer_counts.iter()) {
            kmers_buf.clear();
            kmers::kmers::<u128, _, CANONICAL>(seq, k, &mut kmers_buf);
            assert_eq!(kmers_buf.len(), counts.len());
            for (&kmer, &count) in kmers_buf.iter().zip(counts) {
                if count == 0 {
                    unique_kmers.insert(kmer);
                } else {
                    off_target_kmers.insert(kmer);
                }
            }
        }
        log::info!("    {} k-mers unique to this locus, {} k-mers appear off target",
            unique_kmers.len(), off_target_kmers.len());

        assert!(hard_threshold <= soft_threshold);
        let weight_mult = 1.0 / f64::from(soft_threshold + 1 - hard_threshold);
        let weight_interc = (1.0 - f64::from(hard_threshold)) * weight_mult;
        Self {
            k, unique_kmers, kmers_buf, weight_mult, weight_interc,
            k_2: usize::from(k - 2),
        }
    }

    /// Counts the number of unique k-mers in both read mates and
    // returns weight the read/read pair (can be zero).
    fn read_weight(
        &mut self,
        read_data: &mut ReadData,
        dbg_writer: &mut impl Write,
    ) -> io::Result<f64>
    {
        let mut paired_count = 0_u16;
        write!(dbg_writer, "{}\t", read_data.name_hash)?;
        for mate_data in read_data.mates.iter_mut() {
            if let Some(data) = mate_data.as_mut() {
                self.kmers_buf.clear();
                kmers::kmers::<u128, _, CANONICAL>(&data.sequence, self.k, &mut self.kmers_buf);

                let mut kmers_iter = self.kmers_buf.iter();
                let mut count = 0_u16;
                while let Some(kmer) = kmers_iter.next() {
                    if self.unique_kmers.contains(kmer) {
                        count = count.saturating_add(1);
                        // Skip several k-mers because we are interested in non-overlapping k-mers.
                        // NOTE: Should instead use `advance_by(k - 1)`, but it is not yet stable.
                        kmers_iter.nth(self.k_2);
                    }
                }
                write!(dbg_writer, "{}\t", count)?;
                data.unique_kmers = count;
                paired_count += count;
            } else {
                write!(dbg_writer, "*\t")?;
            }
        }
        let weight = (self.weight_interc + f64::from(paired_count) * self.weight_mult).clamp(0.0, 1.0);
        writeln!(dbg_writer, "{:.2}", weight)?;
        Ok(weight)
    }
}

// ------------------------- Reading alignments for all read pairs -------------------------

/// Checks if any of the alignments is within the "central" region: not in the boundary.
fn in_bounds(alns: &[Alignment], boundary: u32, contigs: &ContigNames) -> bool {
    alns.iter().any(|aln| {
        let contig_len = contigs.get_len(aln.interval().contig_id());
        let aln_middle = aln.interval().middle();
        boundary <= aln_middle && aln_middle < contig_len - boundary
    })
}

#[derive(Default)]
struct ReadCounts {
    total: u32,
    out_of_bounds: u32,
    both_unmapped: u32,
    pair_unmapped: u32,
    few_kmers: u32,
}

impl ReadCounts {
    fn to_string(&self, use_reads: usize, is_paired_end: bool) -> String {
        let mut s = format!("Use {} read{}s. Discard ", use_reads, if is_paired_end { " pair" } else { "" });
        if is_paired_end {
            write!(s, "{} + {} one|two mates", self.pair_unmapped, self.both_unmapped).unwrap();
        } else {
            write!(s, "{}", self.both_unmapped).unwrap();
        }
        write!(s, " poorly mapped, {} out of bounds and {} with few unique k-mers",
            self.out_of_bounds, self.few_kmers).unwrap();
        s
    }
}

/// Paired-end/single-end read alignments, sorted by contig.
pub struct AllAlignments {
    /// Each element: one read pair.
    reads: Vec<GrouppedAlignments>,
    /// This vector stores unused reads (but which passed necessary thresholds).
    unused_reads: Vec<GrouppedAlignments>,
}

impl AllAlignments {
    /// Loads all paired-end or single-end alignments.
    ///
    /// Alignments are only kept if
    /// (i)  there is an alignment with small edit distance (under `err_prof.allowed_edit_dist`)  AND
    /// (ii) there is an alignment that overlaps an inner region of any contig (beyond the boundary region).
    pub fn load(
        reader: impl bam::Read,
        contig_set: &ContigSet,
        bg_distr: &BgDistr,
        edit_dist_cache: &EditDistCache,
        contig_infos: &ContigInfos,
        params: &super::Params,
        mut dbg_writer1: impl Write,
        mut dbg_writer2: impl Write,
        mut dbg_writer3: Option<GzFile>,
    ) -> crate::Result<Self>
    {
        let contigs = contig_set.contigs();
        let boundary = params.boundary_size.checked_sub(params.tweak.unwrap()).unwrap();
        assert!(contigs.lengths().iter().all(|&len| len > 2 * boundary),
            "[{}] Some contigs are too short (must be over twice boundary size = {})", contigs.tag(), 2 * boundary);
        // TODO: Provide hard threshold from CLI.
        let mut unique_kmers = UniqueKmers::new(contig_set, params.kmer_hard_thresh, params.kmer_soft_thresh);

        log::info!("    Loading read alignments");
        writeln!(dbg_writer1, "read_hash\tread_end\tinterval\tedit_dist\tedit_status\tlik\tu5\tread_name")
            .map_err(add_path!(!))?;
        writeln!(dbg_writer2, "read_hash\tuniq_kmers1\tuniq_kmers2\tweight").map_err(add_path!(!))?;
        if let Some(w) = &mut dbg_writer3 {
            writeln!(w, "read_hash\tread_end\tinterval\told_lik\tnew_lik\toverall_weight\tweights")
                .map_err(add_path!(!))?;
        }
        let mut reader = FilteredReader::new(reader, Arc::clone(contigs), bg_distr, edit_dist_cache)?;

        let is_paired_end = bg_distr.insert_distr().is_paired_end();
        let mut hashes = IntSet::default();
        let mut tmp_alns = Vec::new();
        let mut buffer = Vec::with_capacity(16);
        let mut collisions = 0;

        let mut counts = ReadCounts::default();
        let mut reads = Vec::new();
        let mut unused_reads = Vec::new();
        while reader.has_more() {
            let mut read_data = ReadData::default();
            counts.total += 1;
            tmp_alns.clear();
            let mut weight = 1.0;
            let has_first = reader.next_alns(
                ReadEnd::First, &mut tmp_alns, &mut weight, &mut read_data, &mut dbg_writer1)?;
            if !hashes.insert(read_data.name_hash) {
                log::debug!("Read {} produced hash collision ({})", read_data.name, read_data.name_hash);
                collisions += 1;
            }
            let has_second = is_paired_end && reader.next_alns(
                ReadEnd::Second, &mut tmp_alns, &mut weight, &mut read_data, &mut dbg_writer1)?;

            if !has_first && !has_second {
                counts.both_unmapped += 1;
                continue;
            } else if is_paired_end && !(has_first && has_second) {
                counts.pair_unmapped += 1;
                continue;
            } else if !in_bounds(&tmp_alns, boundary, contigs) {
                counts.out_of_bounds += 1;
                continue;
            }

            weight *= unique_kmers.read_weight(&mut read_data, &mut dbg_writer2).map_err(add_path!(!))?;
            let max_alns = if weight >= params.min_weight { MAX_USED_ALNS } else { MAX_UNUSED_ALNS };
            let groupped_alns = if is_paired_end {
                identify_paired_end_alignments(read_data, &mut tmp_alns, weight,
                    max_alns, &mut buffer, bg_distr.insert_distr(), contig_infos, params)
            } else {
                identify_single_end_alignments(read_data, &mut tmp_alns, weight,
                    max_alns, contig_infos, params)
            };
            // TODO: Rethink this, reads may be needed for read depth!
            if groupped_alns.weight() >= params.min_weight {
                reads.push(groupped_alns);
            } else {
                unused_reads.push(groupped_alns);
                counts.few_kmers += 1;
            }
        }
        log::debug!("    {}", counts.to_string(reads.len(), is_paired_end));
        if collisions > 2 && collisions * 100 > counts.total {
            return Err(error!(RuntimeError, "Too many read name collisions ({}). \
                Possibly, paired-end reads are processed as single-end reads.", collisions))
        }
        Ok(Self { reads, unused_reads })
    }

    /// Reads that are used in the model.
    pub fn reads(&self) -> &[GrouppedAlignments] {
        &self.reads
    }

    /// Reads that are not used.
    pub fn unused_reads(&self) -> &[GrouppedAlignments] {
        &self.unused_reads
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
    pub fn write_read_pair_info<const INCLUDE_UNUSED: bool>(
        &self,
        mut f: impl Write,
        contigs: &ContigNames,
    ) -> io::Result<()>
    {
        writeln!(f, "read_hash\tcontig\tpos1\tpos2\tlik")?;
        for read in self.reads.iter() {
            read.write_read_pair_info(&mut f, contigs)?;
        }
        if INCLUDE_UNUSED {
            for read in self.unused_reads.iter() {
                read.write_read_pair_info(&mut f, contigs)?;
            }
        }
        Ok(())
    }
}
