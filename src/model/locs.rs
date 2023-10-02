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
    algo::{bisect, TwoU32, get_hash, HashSet},
    math::Ln,
};

// ------------------------- Read-end and pair-end data, such as sequences and qualities -------------------------

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
    name_hash: u64,
    mates: [Option<MateData>; 2],
}

impl ReadData {
    pub fn name(&self) -> &String {
        &self.name
    }

    pub fn name_hash(&self) -> u64 {
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
    edit_dist_cache: &'a EditDistCache,
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
        edit_dist_cache: &'a EditDistCache,
    ) -> Result<Self, Error> {
        let mut record = bam::Record::new();
        // Reader would return None if there are no more records.
        let has_more = reader.read(&mut record).transpose()?.is_some();
        assert!(is_primary(&record), "First record in the BAM file is secondary/supplementary");
        Ok(FilteredReader {
            found_alns: IntMap::default(),
            reader, record, contigs, err_prof, edit_dist_cache, has_more,
        })
    }

    /// Starting with `self.record` (already loaded), reads all alignments,
    /// corresponding to this read and current read end.
    /// Basically, the function continue to read until the next primary alignment is found,
    /// which is saved to `self.record`.
    ///
    /// If read is unmapped, or best edit distance is too high, does not add any new alignments.
    ///
    /// Returns smallest saved alignment probability (None if unmapped).
    fn next_alns(
        &mut self,
        read_end: ReadEnd,
        alns: &mut Vec<Alignment>,
        read_data: &mut ReadData,
        dbg_writer: &mut impl Write,
    ) -> Result<Option<f64>, Error>
    {
        assert!(self.has_more, "Cannot read any more records from a BAM file");
        let name_hash = get_hash(self.record.qname());
        let name = std::str::from_utf8(self.record.qname()).map_err(|_| Error::InvalidInput(
            format!("Read name is not UTF-8: {:?}", String::from_utf8_lossy(self.record.qname()))))?;

        // Check if everything is correct.
        if read_end == ReadEnd::First {
            assert!(read_data.name.is_empty());
            read_data.name.push_str(name);
            read_data.name_hash = name_hash;
        } else if name_hash != read_data.name_hash {
            return Err(Error::InvalidData(format!("Read {} does not have a second read end", name)));
        }
        if self.record.seq().is_empty() {
            return Err(Error::InvalidData(format!("Read {} does not have read sequence", read_data.name)));
        }
        let old_option = read_data.mates[read_end.ix()].replace(MateData::new(&self.record));
        assert!(old_option.is_none(), "Mate data defined twice");

        if self.record.is_unmapped() {
            writeln!(dbg_writer, "{:X}\t{}\t*\tNA\t-\t{}", name_hash, read_end, name).map_err(add_path!(!))?;
            // Read the next record, and save false to `has_more` if there are no more records left.
            self.has_more = self.reader.read(&mut self.record).transpose()?.is_some();
            return Ok(None);
        }

        let start_len = alns.len();
        // Does any alignment have good edit distance?
        let mut any_good = false;
        let mut min_prob = f64::INFINITY;
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
            let read_prof = aln.count_region_operations_fast(self.contigs.get_len(aln.interval().contig_id()));
            let aln_prob = self.err_prof.ln_prob(&read_prof);
            let (edit_dist, read_len) = read_prof.edit_and_read_len();
            let (good_dist, passable_dist) = self.edit_dist_cache.get(read_len);
            let is_good_dist = edit_dist <= good_dist;
            // Discard all alignments with edit distance over half of the read length.
            let is_passable_dist = edit_dist <= passable_dist;
            any_good |= is_good_dist;

            write!(dbg_writer, "{:X}\t{}\t{}\t{}/{}\t{}\t{:.2}", name_hash, read_end, aln.interval(),
                edit_dist, read_len, if is_good_dist { '+' } else if is_passable_dist { '~' } else { '-' },
                Ln::to_log10(aln_prob)).map_err(add_path!(!))?;
            if primary {
                write!(dbg_writer, "\t{}", read_data.name).map_err(add_path!(!))?;
                primary = false;
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

        if any_good {
            Ok(Some(min_prob))
        } else {
            alns.truncate(start_len);
            Ok(None)
        }
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
            write!(f, "{:X}\t{}\t", self.read_data.name_hash, contigs.get_name(pair.contig_id))?;
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

/// For a paired-end read, combine and pair all mate alignments across all contigs.
/// Input alignments are sorted first by contig, and then by read end.
/// Returns None if the read pair is ignored.
/// Returns groupped alignments and true if the alignments should be used for read assignment.
fn identify_paired_end_alignments(
    read_data: ReadData,
    tmp_alns: &mut Vec<Alignment>,
    opt_min_prob1: Option<f64>,
    opt_min_prob2: Option<f64>,
    mut weight: f64,
    buffer: &mut Vec<f64>,
    insert_distr: &InsertDistr,
    ln_ncontigs: f64,
    params: &super::Params,
) -> (GrouppedAlignments, bool)
{
    let mut use_pair = true;
    if opt_min_prob2.is_none() || opt_min_prob2.is_none() {
        weight *= 0.5; // Halve weight if one of the mates is unmapped.
        use_pair = params.use_unpaired;
    }
    use_pair &= weight >= params.min_weight;
    let max_alns = if use_pair { MAX_USED_ALNS } else { MAX_UNUSED_ALNS };
    // Will be useful for in some calculations below.
    let min_prob1 = opt_min_prob1.unwrap_or(0.0);
    let min_prob2 = opt_min_prob2.unwrap_or(0.0);

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
                identify_contig_pair_alns(&alignments, i, min(j, k), &mut aln_pairs, min_prob1, min_prob2, max_alns,
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
        identify_contig_pair_alns(&alignments, i, min(j, k), &mut aln_pairs, min_prob1, min_prob2, max_alns,
            buffer, insert_distr, params);
    }

    let both_unmapped = min_prob1 + min_prob2
        + 2.0 * params.unmapped_penalty + insert_distr.insert_penalty();
    // Only for normalization, unmapped probability is multiplied by the number of contigs
    // because there is an unmapped possibility for every contig, which we do not store explicitely.
    let norm_fct = Ln::map_sum_init(&aln_pairs, PairAlignment::ln_prob, both_unmapped + ln_ncontigs);
    for aln in aln_pairs.iter_mut() {
        aln.ln_prob = weight * (aln.ln_prob - norm_fct);
    }
    (GrouppedAlignments {
        read_data, weight,
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
    min_prob: f64,
    weight: f64,
    ln_ncontigs: f64,
    params: &super::Params,
) -> (GrouppedAlignments, bool)
{
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

    let unmapped_prob = min_prob + params.unmapped_penalty;
    // Only for normalization, unmapped probability is multiplied by the number of contigs
    // because there is an unmapped possibility for every contig, which we do not store explicitely.
    let norm_fct = Ln::map_sum_init(&aln_pairs, PairAlignment::ln_prob, unmapped_prob + ln_ncontigs);
    for aln in aln_pairs.iter_mut() {
        aln.ln_prob = weight * (aln.ln_prob - norm_fct);
    }
    (GrouppedAlignments {
        read_data, weight,
        unmapped_prob: weight * (unmapped_prob - norm_fct),
        alignments, aln_pairs,
    }, use_read)
}

// ------------------- Assign weights to reads based on the number of unique non-overlapping k-mers -------------------

struct UniqueKmers {
    k: u8,
    unique_kmers: HashSet<u128>,
    kmers_buf: Vec<u128>,
}

impl UniqueKmers {
    /// Stores all k-mers, unique to the current locus.
    fn new(contig_set: &ContigSet) -> Self {
        let kmer_counts = contig_set.kmer_counts();
        let k = u8::try_from(kmer_counts.k()).unwrap();
        assert!(k > 1);
        let mut kmers_buf = Vec::new();
        let mut unique_kmers = HashSet::default();
        for (seq, counts) in contig_set.seqs().iter().zip(kmer_counts.iter()) {
            kmers_buf.clear();
            kmers::kmers::<u128, CANONICAL>(seq, k, &mut kmers_buf);
            assert_eq!(kmers_buf.len(), counts.len());
            // Add all k-mers, for which off-target count is 0.
            unique_kmers.extend(kmers_buf.iter().zip(counts).filter(|(_, &count)| count == 0).map(|(&kmer, _)| kmer));
        }
        Self { k, unique_kmers, kmers_buf }
    }

    /// Count the number of unique k-mers in both read mates and returns read pair weight.
    fn read_weight(&mut self, read_data: &ReadData, dbg_writer: &mut impl Write) -> io::Result<f64> {
        let mut paired_count = 0;
        write!(dbg_writer, "{:X}\t", read_data.name_hash)?;
        for mate_data in &read_data.mates {
            if let Some(data) = mate_data {
                self.kmers_buf.clear();
                kmers::kmers::<u128, CANONICAL>(&data.sequence, self.k, &mut self.kmers_buf);

                let mut kmers_iter = self.kmers_buf.iter();
                let mut count = 0;
                while let Some(kmer) = kmers_iter.next() {
                    if self.unique_kmers.contains(kmer) {
                        count += 1;
                        // Skip several k-mers because we are interested in non-overlapping k-mers.
                        // NOTE: Should instead use `advance_by(k - 1)`, but it is not yet stable.
                        kmers_iter.nth(usize::from(self.k) - 2);
                    }
                }
                write!(dbg_writer, "{}\t", count)?;
                paired_count += count;
            } else {
                write!(dbg_writer, "*\t")?;
            }
        }
        let weight = 1.0; // TODO: More sophisticated.
        writeln!(dbg_writer, "{:.5}", weight)?;
        Ok(weight)
    }
}

// ------------------------- All alignments for all read pairs -------------------------

/// Checks if any of the alignments is within the "central" region: not in the boundary.
fn in_bounds(alns: &[Alignment], boundary: u32, contigs: &ContigNames) -> bool {
    alns.iter().any(|aln| {
        let contig_len = contigs.get_len(aln.interval().contig_id());
        let aln_middle = aln.interval().middle();
        boundary <= aln_middle && aln_middle < contig_len - boundary
    })
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
        contig_set: &ContigSet,
        bg_distr: &BgDistr,
        edit_dist_cache: &EditDistCache,
        params: &super::Params,
        mut dbg_writer1: impl Write,
        mut dbg_writer2: impl Write,
    ) -> Result<Self, Error>
    {
        log::info!("    Loading read alignments");
        let contigs = contig_set.contigs();
        let boundary = params.boundary_size.checked_sub(params.tweak.unwrap()).unwrap();
        assert!(contigs.lengths().iter().all(|&len| len > 2 * boundary),
            "[{}] Some contigs are too short (must be over {})", contigs.tag(), 2 * boundary);
        let mut unique_kmers = UniqueKmers::new(contig_set);

        writeln!(dbg_writer1, "read_hash\tread_end\tinterval\tedit_dist\tedit_status\tlik\tread_name")
            .map_err(add_path!(!))?;
        writeln!(dbg_writer2, "read_hash\tuniq_kmers1\tuniq_kmers2\tweight").map_err(add_path!(!))?;

        let mut reader = FilteredReader::new(reader, Arc::clone(contigs), bg_distr.error_profile(), edit_dist_cache)?;

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
        let mut out_of_bounds_reads = 0;
        let mut unmapped_reads = 0;
        while reader.has_more() {
            let mut read_data = ReadData::default();
            total_reads += 1;
            tmp_alns.clear();
            let opt_min_prob1 = reader.next_alns(ReadEnd::First, &mut tmp_alns, &mut read_data, &mut dbg_writer1)?;
            if !hashes.insert(read_data.name_hash) {
                log::warn!("Read {} produced hash collision ({:X})", read_data.name, read_data.name_hash);
                collisions += 1;
            }
            let opt_min_prob2 = if is_paired_end {
                reader.next_alns(ReadEnd::Second, &mut tmp_alns, &mut read_data, &mut dbg_writer1)?
            } else { None };

            if opt_min_prob1.is_none() && opt_min_prob2.is_none() {
                unmapped_reads += 1;
                continue;
            } else if !in_bounds(&tmp_alns, boundary, contigs) {
                out_of_bounds_reads += 1;
                continue;
            }

            let weight = unique_kmers.read_weight(&read_data, &mut dbg_writer2).map_err(add_path!(!))?;
            let (groupped_alns, use_read) = if is_paired_end {
                identify_paired_end_alignments(read_data, &mut tmp_alns, opt_min_prob1, opt_min_prob2, weight,
                    &mut buffer, insert_distr, ln_ncontigs, params)
            } else {
                identify_single_end_alignments(read_data, &mut tmp_alns, opt_min_prob1.unwrap(), weight,
                    ln_ncontigs, params)
            };
            if use_read {
                reads.push(groupped_alns);
            } else {
                w0_reads.push(groupped_alns);
            }
        }
        log::info!("    Total {} read{},  {} good,  {} low weight{},  {} out of bounds,  {} unmapped",
            total_reads, if is_paired_end { " pairs" } else { "s" }, reads.len(),
            w0_reads.len(), if is_paired_end && !params.use_unpaired { "/unpaired" } else { "" },
            out_of_bounds_reads, unmapped_reads);
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
