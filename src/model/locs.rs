use std::{
    sync::Arc,
    io::Write,
};
use htslib::bam;
use nohash::IntSet;
use crate::{
    Error,
    seq::{
        Interval, ContigId, ContigNames,
        aln::{LightAlignment, Alignment, ReadEnd},
        cigar::Cigar,
    },
    bg::{
        BgDistr,
        err_prof::ErrorProfile,
        insertsz::InsertDistr,
    },
    algo::{fnv1a, bisect},
    math::Ln,
};

pub(crate) const CSV_HEADER: &'static str = "read_hash\tread_end\tinterval\tcentral\tedit_dist\tedit_status\tlik";

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
}

impl<'a, R: bam::Read> FilteredReader<'a, R> {
    fn new(
        mut reader: R,
        contigs: Arc<ContigNames>,
        err_prof: &'a ErrorProfile,
        boundary: u32,
    ) -> Result<Self, Error> {
        let mut record = bam::Record::new();
        Ok(FilteredReader {
            // Reader would return None if there are no more records.
            has_more: reader.read(&mut record).transpose()?.is_some(),
            reader, record, contigs, err_prof, boundary,
        })
    }

    /// Starting with `self.record` (already loaded), reads all alignments,
    /// corresponding to this read and current read end.
    /// Basically, the function continue to read until the next primary alignment is found,
    /// which is saved to `self.record`.
    ///
    /// If best edit distance is too high, does not add any new alignments.
    /// Otherwise, normalizes all alignment probabilities by the best probability (therefore max is always = 0).
    ///
    /// Return triple
    /// - name_hash,
    /// - min_prob: smallest ln-probability of all non-discarded probabilities [0 if there are no alignments],
    /// - any_central: true if at least one alignment is in the central region of the corresponding contig.
    fn next_alns(
        &mut self,
        read_end: ReadEnd,
        alns: &mut Vec<LightAlignment>,
        dbg_writer: &mut impl Write,
    ) -> Result<(u64, f64, bool), Error>
    {
        assert!(self.has_more, "Cannot read any more records from a BAM file");
        let name_hash = fnv1a(self.record.qname());
        if self.record.is_unmapped() {
            writeln!(dbg_writer, "{:X}\t{}\t*\tF\tNA\t-\tNA", name_hash, read_end)?;
            // Read the next record, and save false to `has_more` if there are no more records left.
            self.has_more = self.reader.read(&mut self.record).transpose()?.is_some();
            return Ok((name_hash, 0.0, false));
        }

        let start_len = alns.len();
        // Is any of the alignments in the central region of a contig?
        let mut any_central = false;
        // Does any alignment have good edit distance?
        let mut any_good = false;
        let mut min_prob = f64::INFINITY;
        loop {
            let mut aln = Alignment::new_without_name(
                &self.record, Cigar::from_raw(&self.record), read_end, Arc::clone(&self.contigs), f64::NAN);
            let aln_interval = aln.interval();
            let contig_len = self.contigs.get_len(aln_interval.contig_id());
            let central = self.boundary < aln_interval.end() && aln_interval.start() < contig_len - self.boundary;
            any_central |= central;

            let read_prof = aln.count_region_operations_fast(contig_len);
            let aln_prob = self.err_prof.ln_prob(&read_prof);

            let (edit_dist, read_len) = read_prof.edit_and_read_len();
            let allowed_edit_dist = self.err_prof.allowed_edit_dist(read_len);
            let good_dist = edit_dist <= allowed_edit_dist;
            // Discard all alignments with edit distance over half of the read length.
            let medium_dist = 2 * edit_dist < read_len;
            any_good |= good_dist;

            writeln!(dbg_writer, "{:X}\t{}\t{}\t{}\t{}/{}\t{}\t{:.2}",
                name_hash, read_end, aln_interval, if central { 'T' } else { 'F' },
                edit_dist, read_len, if good_dist { '+' } else if medium_dist { '~' } else { '-' },
                Ln::to_log10(aln_prob))?;
            aln.set_ln_prob(aln_prob);

            if medium_dist {
                alns.push(aln.take_light_aln());
                min_prob = aln_prob.min(min_prob);
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

        Ok(if any_good {
            (name_hash, min_prob, any_central)
        } else {
            alns.truncate(start_len);
            (name_hash, 0.0, false)
        })
    }

    /// Returns the name of the next record.
    fn next_name(&self) -> Result<&str, Error> {
        std::str::from_utf8(self.record.qname()).map_err(|_|
            Error::InvalidInput(format!("Read name is not UTF-8: {:?}",
                String::from_utf8_lossy(self.record.qname()))))
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
    alns1: &[LightAlignment],
    min_prob1: f64,
    alns2: &[LightAlignment],
    min_prob2: f64,
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
    let mut max_prob = f64::NEG_INFINITY;
    for aln1 in alns1.iter() {
        let mut max_prob1 = f64::NEG_INFINITY;
        for (aln2, max_prob2) in alns2.iter().zip(buffer.iter_mut()) {
            let prob = aln1.paired_prob(aln2, insert_distr);
            max_prob1 = max_prob1.max(prob);
            *max_prob2 = max_prob2.max(prob);
            new_alns.push(PairAlignment::new(Pair::Both(aln1.interval().clone(), aln2.interval().clone()), prob));
        }

        // Only add alignment `aln1,unmapped2` if it is better than any existing paired alignment to the same contig.
        let alone_prob1 = aln1.ln_prob() + term1;
        if alone_prob1 >= max_prob1 {
            max_prob1 = alone_prob1;
            new_alns.push(PairAlignment::new(Pair::First(aln1.interval().clone()), alone_prob1));
        }
        max_prob = max_prob.max(max_prob1);
    }

    let term2 = insert_penalty + min_prob1 + params.unmapped_penalty;
    for (aln2, &max_prob2) in alns2.iter().zip(buffer.iter()) {
        let alone_prob2 = aln2.ln_prob() + term2;
        // Only add alignment `unmapped1,aln2` if it is better than existing `aln1,aln2` pairs.
        if alone_prob2 >= max_prob2 {
            max_prob = max_prob.max(alone_prob2);
            new_alns.push(PairAlignment::new(Pair::Second(aln2.interval().clone()), alone_prob2));
        }
    }

    let thresh_prob = max_prob - params.prob_diff;
    let mut i = start_len;
    while i < new_alns.len() {
        if new_alns[i].ln_prob >= thresh_prob {
            i += 1;
        } else {
            new_alns.swap_remove(i);
        }
    }
}

/// All appropriate alignments for one read pair / single read.
pub struct GrouppedAlignments {
    /// Read name.
    read_name: String,
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
        &self.read_name
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
        let j = bisect::right_linear(&self.alns, |paln| paln.contig_id().cmp(&contig_id), i, self.alns.len());
        (i, &self.alns[i..j])
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

/// For a paired-end read, combine and pair all single-mate alignments across all contigs.
///
/// Input alignments are sorted first by contig, and then by read end.
fn identify_paired_end_alignments(
    read_name: String,
    name_hash: u64,
    alns: &[LightAlignment],
    min_prob1: f64,
    min_prob2: f64,
    buffer: &mut Vec<f64>,
    insert_distr: &InsertDistr,
    ln_ncontigs: f64,
    params: &super::Params,
) -> GrouppedAlignments
{
    let mut groupped_alns = Vec::new();
    let n = alns.len();
    let mut i = 0;
    while i < n {
        // For the current contig id, first mates will be in i..j, and second mates in j..k.
        let contig_id = alns[i].contig_id();
        // Sort key = contig_id * 2 + (first_mate ? 0 : 1).
        // sort_key1: current sort key for first mates, sort_key2: current sort key for second mates.
        let sort_key1 = contig_id.get() << 1;
        let sort_key2 = sort_key1 | 1;
        let j = if alns[i].read_end() == ReadEnd::First {
            bisect::right_linear(alns, |aln| aln.sort_key().cmp(&sort_key1), i, n)
        } else { i };
        let k = bisect::right_linear(alns, |aln| aln.sort_key().cmp(&sort_key2), j, n);
        identify_contig_pair_alns(&mut groupped_alns, &alns[i..j], min_prob1, &alns[j..k], min_prob2,
            buffer, insert_distr, params);
        i = k;
    }

    let both_unmapped = min_prob1 + min_prob2 + 2.0 * params.unmapped_penalty + insert_distr.insert_penalty();
    // Only for normalization, unmapped probability is multiplied by the number of contigs
    // because there is an unmapped possibility for every contig, which we do not store explicitely.
    let norm_fct = Ln::map_sum_init(&groupped_alns, PairAlignment::ln_prob, both_unmapped + ln_ncontigs);
    for aln in groupped_alns.iter_mut() {
        aln.ln_prob -= norm_fct;
    }
    GrouppedAlignments {
        read_name, name_hash,
        unmapped_prob: both_unmapped - norm_fct,
        alns: groupped_alns,
    }
}

/// For a single-end read, sort alignment across contigs, discard improbable alignments, and normalize probabilities.
/// Input alignments are sorted first by contig.
fn identify_single_end_alignments(
    read_name: String,
    name_hash: u64,
    alns: &[LightAlignment],
    min_prob: f64,
    ln_ncontigs: f64,
    params: &super::Params,
) -> GrouppedAlignments
{
    let mut groupped_alns = Vec::new();
    let n = alns.len();
    let mut i = 0;
    while i < n {
        // For the current contig id, first mates will be in i..j, and second mates in j..k.
        let contig_id = alns[i].contig_id();
        let j = bisect::right_linear(alns, |aln| aln.contig_id().cmp(&contig_id), i, n);
        let thresh_prob = alns[i..j].iter()
            .fold(f64::NEG_INFINITY, |m, aln| m.max(aln.ln_prob())) - params.prob_diff;
        groupped_alns.extend(alns[i..j].iter()
            .filter(|aln| aln.ln_prob() >= thresh_prob)
            .map(|aln| PairAlignment::new(Pair::First(aln.interval().clone()), aln.ln_prob())));
        i = j;
    }

    let unmapped_prob = min_prob + params.unmapped_penalty;
    // Only for normalization, unmapped probability is multiplied by the number of contigs
    // because there is an unmapped possibility for every contig, which we do not store explicitely.
    let norm_fct = Ln::map_sum_init(&groupped_alns, PairAlignment::ln_prob, unmapped_prob + ln_ncontigs);
    for aln in groupped_alns.iter_mut() {
        aln.ln_prob -= norm_fct;
    }
    GrouppedAlignments {
        read_name, name_hash,
        unmapped_prob: unmapped_prob - norm_fct,
        alns: groupped_alns,
    }
}

pub struct AllAlignments(Vec<GrouppedAlignments>);

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
        params: &super::Params,
        mut dbg_writer: impl Write,
    ) -> Result<Self, Error>
    {
        log::info!("[{}] Loading read alignments", contigs.tag());
        let boundary = params.boundary_size.checked_sub(params.tweak).unwrap();
        assert!(contigs.lengths().iter().all(|&len| len > 2 * boundary),
            "[{}] Some contigs are too short (must be over {})", contigs.tag(), 2 * boundary);
        let mut reader = FilteredReader::new(reader, Arc::clone(contigs), bg_distr.error_profile(), boundary)?;

        let insert_distr = bg_distr.insert_distr();
        let is_paired_end = insert_distr.is_paired_end();
        let ln_ncontigs = (contigs.len() as f64).ln();
        let mut all_alns = Vec::new();
        let mut hashes = IntSet::default();
        let mut tmp_alns = Vec::new();
        let mut buffer = Vec::with_capacity(16);
        while reader.has_more {
            tmp_alns.clear();
            let read_name = reader.next_name()?.to_owned();
            write!(dbg_writer, "{}=", read_name)?;
            let (hash, min_prob1, mut central) = reader.next_alns(ReadEnd::First, &mut tmp_alns, &mut dbg_writer)?;
            if !hashes.insert(hash) {
                log::warn!("Read {} produced hash collision ({:X}). \
                    If many such messages, reads appear in an unordered fashion, will lead to errors",
                    read_name, hash);
            }
            let min_prob2 = if is_paired_end {
                let (hash2, min_prob2, central2) = reader.next_alns(ReadEnd::Second, &mut tmp_alns, &mut dbg_writer)?;
                if hash != hash2 {
                    return Err(Error::InvalidData(format!("Read {} with hash {:X} does not have a second read end",
                        read_name, hash)));
                }
                central |= central2;
                min_prob2
            } else { 0.0 };
            if !central {
                continue;
            }

            let groupped_alns = if is_paired_end {
                tmp_alns.sort_unstable_by_key(LightAlignment::sort_key);
                identify_paired_end_alignments(read_name, hash,
                    &tmp_alns, min_prob1, min_prob2, &mut buffer, insert_distr, ln_ncontigs, params)
            } else {
                tmp_alns.sort_unstable_by_key(LightAlignment::contig_id);
                identify_single_end_alignments(read_name, hash, &tmp_alns, min_prob1, ln_ncontigs, params)
            };
            all_alns.push(groupped_alns);
        }
        Ok(Self(all_alns))
    }

    /// Returns the total number of saved reads/read pairs.
    pub fn n_reads(&self) -> usize {
        self.0.len()
    }

    /// Returns iterator over `GrouppedAlignments` for all read pairs.
    pub fn iter(&self) -> std::slice::Iter<'_, GrouppedAlignments> {
        self.0.iter()
    }
}
