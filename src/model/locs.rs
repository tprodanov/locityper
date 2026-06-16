use std::{
    cmp::{min, max},
    fmt,
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
        counts::KmerCounts,
        transfer::HapAlns,
        wfa::Aligner,
    },
    bg::{
        BgDistr,
        err_prof::{ErrorProfile, EditDistCache},
        insertsz::InsertDistr,
    },
    algo::{bisect, get_hash, HashSet, IntSet, IntMap, hash_map::Entry},
    math::Ln,
    model::{
        windows::ContigInfos,
    },
    ext::{
        sys::GzFile,
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

    /// Returns sequence for the corresponding strand.
    pub fn get_seq(&self, strand: Strand) -> &[u8] {
        if strand == self.strand {
            &self.sequence
        } else {
            self.opp_sequence.get_or_init(|| seq::reverse_complement(&self.sequence))
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

    /// Read length.
    pub fn len(&self) -> u32 {
        self.sequence.len() as u32
    }
}

/// Read information: read name, two sequences and qualities (for both mates).
#[derive(Clone)]
pub struct ReadData {
    name: String,
    name_hash: NameHash,
    mates: [Option<MateData>; 2],
    /// Read pair weight, taken from edit distance.
    weight: f64,
}

impl Default for ReadData {
    fn default() -> Self {
        Self {
            name: String::new(),
            name_hash: NameHash::default(),
            mates: Default::default(),
            weight: 1.0,
        }
    }
}

impl ReadData {
    pub fn set_name(&mut self, record: &bam::Record, read_end: ReadEnd) -> crate::Result<()> {
        let qname = record.qname();
        if read_end == ReadEnd::Second {
            return if qname != self.name.as_bytes() {
                Err(error!(InvalidData, "Read {} does not have a second read end", self.name))
            } else {
                Ok(())
            }
        }
        assert!(self.name.is_empty());
        self.name = String::from_utf8(qname.to_vec())
            .map_err(|_| Error::Utf8("read name", record.qname().to_vec()))?;
        self.name_hash = NameHash::new(qname);
        Ok(())
    }

    #[inline(always)]
    pub fn name(&self) -> &String {
        &self.name
    }

    #[inline(always)]
    pub fn name_hash(&self) -> NameHash {
        self.name_hash
    }

    #[inline(always)]
    pub fn mate_data(&self, read_end: ReadEnd) -> &MateData {
        self.mates[read_end.ix()].as_ref().expect("Mate data undefined")
    }
}

/// Store start positions divided by 128.
const STEP_PWR: u32 = 7;
/// During query, check the position div 128 and one of the neighboring positions (two queries in total).
/// This way we will for certain find the position by checking queries within 64 bp.
const DIST_PWR: u32 = STEP_PWR - 1;
/// To check if the value is closer to zero or to 128, check (val & UP_DOWN_MASK).
const UP_DOWN_MASK: u32 = 1 << DIST_PWR;

#[inline(always)]
fn encode(read_end: ReadEnd, contig: ContigId, pos: u32) -> u64 {
    (read_end.as_int::<u64>() << 48) | (u64::from(contig.get()) << 32) | u64::from(pos >> STEP_PWR)
}

#[derive(Clone, Copy)]
pub(crate) struct PosCollectionValue {
    index: u32,
    pos: u32,
}

impl PosCollectionValue {
    #[inline(always)]
    pub(crate) fn new(index: u32, pos: u32) -> Self {
        Self { index, pos }
    }

    #[inline(always)]
    pub(crate) fn index(self) -> u32 {
        self.index
    }

    // #[inline(always)]
    // pub(crate) fn pos(self) -> u32 {
    //     self.pos
    // }
}

/// Use this value for indices of alignments that are not saved due to poor edit distance.
const NOT_SAVED: u32 = u32::MAX;

/// Store previously observed alignment starts.
#[derive(Default)]
pub(crate) struct PosCollection {
    /// Key = encode(contig, pos); value = PosCollectionValue { index, pos }.
    map: IntMap<u64, PosCollectionValue>,
}

impl PosCollection {
    // /// Adds position to the collection. Returns false if the 128-bin is already in use.
    // /// Returns `Some(index)` and does not update collection if a similar position (within 128 bp) is already stored.
    // #[inline(always)]
    // pub(crate) fn add(&mut self, read_end: ReadEnd, contig: ContigId, pos: u32, index: u32) -> bool {
    //     match self.entry(read_end, contig, pos) {
    //         Entry::Occupied(_) => false,
    //         Entry::Vacant(entry) => {
    //             entry.insert(PosCollectionValue::new(index, pos));
    //             true
    //         }
    //     }
    // }

    /// Returns entry corresponding to a given contig and the 128 bp bin around the position.
    /// When modifying the entry, be careful that start position is still inside the bin
    /// (won't panic if not, but can screw up the distance check in `get`).
    #[inline(always)]
    pub(crate) fn entry<'a>(
        &'a mut self,
        read_end: ReadEnd,
        contig: ContigId,
        pos: u32,
    ) -> Entry<'a, u64, PosCollectionValue> {
        self.map.entry(encode(read_end, contig, pos))
    }

    /// Returns previously observed index if the position is similar (within 64 bp or inside the same 128 bp bin).
    #[inline(always)]
    pub(crate) fn get(&self, read_end: ReadEnd, contig: ContigId, pos: u32) -> Option<u32> {
        let key = encode(read_end, contig, pos);
        if let Some(&val) = self.map.get(&key) {
            return Some(val.index());
        }
        // +1 if start is closer to 128, -1 if start is closer to 0.
        let key_increment = 1 - (i64::from((pos & UP_DOWN_MASK) == 0) << 1);
        // Don't care about wrapping since we won't get any matches in the map.
        // Could also use `key.checked_add_signed()?`, but it may be slower.
        let key2 = key.wrapping_add_signed(key_increment);
        match self.map.get(&key2) {
            Some(&val) if (val.pos.abs_diff(pos) >> DIST_PWR) != 0 => Some(val.index()),
            _ => None,
        }
    }
}

/// Structure for storing initial alignments, before any subsequent pairing.
pub(crate) struct PrelimAlignments {
    alns: Vec<Alignment>,
    pos_collection: PosCollection,
    passable_dist: [u32; 2],
    best_edit: [u32; 2],
    best_lik: [f64; 2],
}

impl PrelimAlignments {
    fn new() -> Self {
        Self {
            alns: Vec::new(),
            pos_collection: Default::default(),
            passable_dist: [u32::MAX; 2],
            best_edit: [u32::MAX; 2],
            best_lik: [f64::NEG_INFINITY; 2],
        }
    }

    #[inline]
    fn set_passable_dist(&mut self, read_end: ReadEnd, dist: u32) {
        self.passable_dist[read_end.ix()] = dist;
    }

    /// Adds new alignment. Returns false if edit distance is under `passable_dist`.
    pub(crate) fn push(
        &mut self,
        mut aln: Alignment,
        contigs: &ContigNames,
        err_prof: &ErrorProfile,
    ) -> bool {
        let read_end = aln.read_end();
        let read_prof = aln.count_region_operations_fast(contigs.get_len(aln.contig_id()));
        let dist = read_prof.edit_distance();
        let aln_prob = err_prof.ln_prob(&read_prof);
        aln.set_distance(dist);
        aln.set_ln_prob(aln_prob);
        self.best_edit[read_end.ix()] = self.best_edit[read_end.ix()].min(dist.edit());
        self.best_lik[read_end.ix()] = self.best_lik[read_end.ix()].max(aln_prob);

        let new_aln_ix = self.alns.len() as u32;
        let save = dist.edit() <= self.passable_dist[read_end.ix()];
        // Return false now before we update pos collection.
        if new_aln_ix == 0 && !save {
            return false;
        }

        let aln_start = aln.interval().start();
        match (self.pos_collection.entry(read_end, aln.contig_id(), aln_start), save) {
            (Entry::Occupied(mut entry), true) => {
                let entry = entry.get_mut();
                let aln_ix = entry.index();
                if aln_ix == NOT_SAVED {
                    entry.index = new_aln_ix;
                    entry.pos = aln_start;
                    self.alns.push(aln);
                } else if aln_prob > self.alns[aln_ix as usize].ln_prob() {
                    self.alns[aln_ix as usize] = aln;
                    entry.pos = aln_start;
                }
            }
            (Entry::Occupied(_entry), false) => {}
            (Entry::Vacant(entry), true) => {
                entry.insert(PosCollectionValue::new(new_aln_ix, aln_start));
                self.alns.push(aln);
            }
            (Entry::Vacant(entry), false) => {
                entry.insert(PosCollectionValue::new(NOT_SAVED, aln_start));
            }
        }
        save
    }

    fn debug_write(&self, read_data: &ReadData, dbg_writer: &mut impl Write) -> crate::Result<()> {
        for (i, aln) in self.alns.iter().enumerate() {
            write!(dbg_writer, "{}\t{}\t{}\t{}\t{:.2}", read_data.name_hash, aln.read_end(), aln.interval(),
                aln.distance().unwrap(), Ln::to_log10(aln.ln_prob())).map_err(add_path!(!))?;
            if i == 0 {
                write!(dbg_writer, "\t{}", read_data.name).map_err(add_path!(!))?;
            }
            writeln!(dbg_writer).map_err(add_path!(!))?;
        }
        Ok(())
    }

    fn finalize(&mut self) {
        self.alns.iter_mut().for_each(|aln| aln.set_ln_prob(aln.ln_prob() - self.best_lik[aln.read_end().ix()]));
    }

    /// Returns the number of stored alignments.
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.alns.len()
    }

    #[inline(always)]
    pub fn pos_collection(&self) -> &PosCollection {
        &self.pos_collection
    }

    #[inline(always)]
    pub fn passable_dist(&self, read_end: ReadEnd) -> u32 {
        self.passable_dist[read_end.ix()]
    }
}

impl std::ops::Index<usize> for PrelimAlignments {
    type Output = Alignment;

    fn index(&self, i: usize) -> &Alignment {
        self.alns.index(i)
    }
}

// ------------------------- BAM file reader, that skips bad alignments -------------------------

/// For each BAM contig, store its index across `ContigNames`.
/// BAM contigs must be a subset (possibly equal) of the contig names since if there are new contigs,
/// we cannot transfer alignments from them.
fn construct_tid_to_contig_map(header: &bam::HeaderView, contigs: &ContigNames) -> crate::Result<Vec<ContigId>> {
    let mut conversion = Vec::with_capacity(header.target_count() as usize);
    for i in 0..header.target_count() {
        let name_bytes = header.tid2name(i);
        let name = std::str::from_utf8(name_bytes).map_err(|_| Error::Utf8("read name", name_bytes.to_vec()))?;
        let id = contigs.try_get_id(name)
            .ok_or_else(|| error!(InvalidData, "Intermediate BAM file contains unexpected contigs (e.g. {})", name))?;
        conversion.push(id);
    }
    Ok(conversion)
}

#[inline(always)]
fn is_primary(record: &bam::Record) -> bool {
    (record.flags() & 2304) == 0
}

/// BAM reader that allows to peek the last BAM record.
struct LaggedReader<R> {
    reader: R,
    record: bam::Record,
    has_more: bool,
}

impl<R: bam::Read> LaggedReader<R> {
    fn new(mut reader: R) -> crate::Result<Self> {
        let mut record = bam::Record::new();
        // Reader would return None if there are no more records.
        let has_more = reader.read(&mut record).transpose()?.is_some();
        assert!(is_primary(&record), "First record in the BAM file is secondary/supplementary");
        Ok(Self { reader, record, has_more })
    }

    /// Reads next read into `self.record`, returns true if there are more reads left.
    #[inline]
    fn proceed(&mut self) -> crate::Result<bool> {
        self.has_more = self.reader.read(&mut self.record).transpose()?.is_some();
        Ok(self.has_more)
    }

    /// Returns the next record unless it has primary alignment.
    #[inline]
    fn next_unless_primary<'a>(&'a mut self) -> crate::Result<Option<&'a bam::Record>> {
        if !self.proceed()? || is_primary(&self.record) {
            Ok(None)
        } else {
            Ok(Some(&self.record))
        }
    }

    #[inline]
    fn skip_until_primary(&mut self) -> crate::Result<()> {
        while let Some(_) = self.next_unless_primary()? {}
        Ok(())
    }

    #[inline]
    fn current(&self) -> Option<&bam::Record> {
        if self.has_more {
            Some(&self.record)
        } else {
            None
        }
    }

    #[inline]
    fn has_more(&self) -> bool {
        self.has_more
    }
}

/// BAM reader.
/// Contains `next_alns` method, which consecutively reads all records with the same name
/// until the next primary alignment.
struct Data<'a> {
    tid2contig: Vec<ContigId>,
    contig_set: &'a ContigSet,
    contig_infos: &'a ContigInfos,
    err_prof: &'a ErrorProfile,
    are_short_reads: bool,
    edit_dist_cache: &'a EditDistCache,
    params: &'a super::Params,
}

impl<'a> Data<'a> {
    fn new(
        contig_set: &'a ContigSet,
        contig_infos: &'a ContigInfos,
        tid2contig: Vec<ContigId>,
        bg_distr: &'a BgDistr,
        edit_dist_cache: &'a EditDistCache,
        params: &'a super::Params,
    ) -> crate::Result<Self> {
        Ok(Self {
            are_short_reads: bg_distr.seq_info().technology().are_short_reads(),
            err_prof: bg_distr.error_profile(),
            contig_set, contig_infos, tid2contig, edit_dist_cache, params,
        })
    }
}

/// Starting with `self.record` (already loaded), reads all alignments,
/// corresponding to this read and current read end.
/// Basically, the function continue to read until the next primary alignment is found,
/// which is saved to `self.record`.
///
/// Returns false if the read is unmapped, or best edit distance is too high.
fn read_next_alns<'a>(
    data: &Data<'a>,
    reader: &mut LaggedReader<impl bam::Read>,
    read_end: ReadEnd,
    read_data: &mut ReadData,
    prelim_alignments: &mut PrelimAlignments,
) -> crate::Result<bool> {
    let record = reader.current().expect("Cannot read any more records from a BAM file");
    read_data.set_name(&record, read_end)?;
    if record.seq().is_empty() {
        if record.qname().is_empty() {
            return Err(
                error!(InvalidData, "Alignment file contains absolutely empty read (no read name or sequence)"));
        }
        return Err(error!(InvalidData, "Read {} does not have read sequence", read_data.name));
    }
    let old_option = read_data.mates[read_end.ix()].replace(MateData::new(&record));
    assert!(old_option.is_none(), "Mate data defined twice");
    if record.is_unmapped() {
        reader.proceed()?;
        return Ok(false);
    }

    let cigar = Cigar::from_raw(record.raw_cigar());
    assert!(!cigar.has_hard_clipping(), "Primary alignment has hard clipping");
    let aln = Alignment::from_record_w_contig_id(record, data.tid2contig[record.tid() as usize],
        cigar, read_end, Arc::clone(data.contig_set.contigs()));
    let neighb_complexity = if data.are_short_reads { data.contig_infos.neighb_complexity(&aln) } else { 1.0 };
    let read_len = record.seq().len() as u32;
    let (good_dist, mut passable_dist) = data.edit_dist_cache.get(read_len);
    let mut threshold_dist = good_dist;
    if neighb_complexity <= data.params.poor_compl {
        threshold_dist = max(good_dist, (data.params.poor_compl_edit * read_len as f64) as u32);
        passable_dist += threshold_dist - good_dist;
    }

    prelim_alignments.set_passable_dist(read_end, passable_dist);
    if !prelim_alignments.push(aln, data.contig_set.contigs(), data.err_prof) {
        // Primary alignment is not good enough.
        reader.skip_until_primary()?;
        return Ok(false);
    }

    while let Some(record) = reader.next_unless_primary()? {
        assert_eq!(record.qname(), read_data.name.as_bytes(),
            "Read {} first alignment is not primary", String::from_utf8_lossy(record.qname()));
        let mut cigar = Cigar::from_raw(record.raw_cigar());
        cigar.hard_to_soft();
        if cigar.is_empty() {
            log::warn!("    Read {} is mapped (flag {}) and has empty CIGAR; skipping this alignment",
                read_data.name, record.flags());
            continue;
        }
        let aln = Alignment::from_record_w_contig_id(record, data.tid2contig[record.tid() as usize],
            cigar, read_end, Arc::clone(data.contig_set.contigs()));
        prelim_alignments.push(aln, data.contig_set.contigs(), data.err_prof);
    }

    let best_edit = prelim_alignments.best_edit[read_end.ix()];
    if best_edit > threshold_dist {
        return Ok(false);
    }
    read_data.weight *= if best_edit <= good_dist { 1.0 } else { (f64::from(good_dist) / f64::from(best_edit)).sqrt() };
    Ok(true)
}

// ------------------------- Structures for storing all alignments for one read pair -------------------------

/// All appropriate alignments for one read pair / single read.
pub struct GrouppedAlignments {
    /// Read name.
    read_data: ReadData,
    /// Read weight.
    weight: f64,
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
    unm_ins_penalty: f64, // Unmapped + insert size penalties
    prob_diff: f64,
) {
    let start_len = aln_pairs.len();
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
        let alone_prob1 = aln1.ln_prob() + unm_ins_penalty;
        if alone_prob1 >= max_prob1 {
            aln_pairs.push(PairAlignment::new_first(ix1, aln1.interval(), alone_prob1));
        }
    }

    for (ix2, &max_prob2) in (j..k).zip(buffer.iter()) {
        let aln2 = unsafe { alignments.get_unchecked(ix2) };
        let alone_prob2 = aln2.ln_prob() + unm_ins_penalty;
        // Only add alignment `unmapped1,aln2` if it is better than existing `aln1,aln2` pairs.
        if alone_prob2 >= max_prob2 {
            aln_pairs.push(PairAlignment::new_second(ix2, aln2.interval(), alone_prob2));
        }
    }

    let slice = &mut aln_pairs[start_len..];
    // Decreasing sort by ln-probability.
    slice.sort_unstable_by(|a, b| b.ln_prob.total_cmp(&a.ln_prob));
    let thresh_prob = slice[0].ln_prob - prob_diff;
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
    max_alns: usize,
    buffer: &mut Vec<f64>,
    insert_distr: &InsertDistr,
    contig_infos: &ContigInfos,
    params: &super::Params,
) -> GrouppedAlignments
{
    let insert_penalty = insert_distr.insert_penalty();
    let unm_ins_penalty = params.unmapped_penalty + insert_penalty;

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
                    buffer, insert_distr, unm_ins_penalty, params.prob_diff);
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
            buffer, insert_distr, unm_ins_penalty, params.prob_diff);
    }

    let weight = read_data.weight * contig_infos.explicit_read_weight(&aln_pairs);
    for aln in aln_pairs.iter_mut() {
        aln.ln_prob *= weight;
    }
    GrouppedAlignments {
        read_data, alignments, aln_pairs, weight,
        unmapped_prob: weight * (2.0 * params.unmapped_penalty + insert_penalty),
    }
}

/// For a single-end read, sort alignment across contigs, discard improbable alignments, and normalize probabilities.
/// Input alignments are sorted first by contig.
/// Returns groupped alignments and true if the alignments should be used for read assignment.
fn identify_single_end_alignments(
    read_data: ReadData,
    tmp_alns: &mut Vec<Alignment>,
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

    let weight = read_data.weight * contig_infos.explicit_read_weight(&aln_pairs);
    for aln in aln_pairs.iter_mut() {
        aln.ln_prob *= weight;
    }
    GrouppedAlignments {
        unmapped_prob: weight * params.unmapped_penalty,
        read_data, alignments, aln_pairs, weight,
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
        kmer_counts: &KmerCounts,
        hard_threshold: u16,
        soft_threshold: u16,
    ) -> Self {
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
    poorly_mapped: u32,
    out_of_bounds: u32,
    few_kmers: u32,
}

impl ReadCounts {
    fn to_string(&self, use_reads: usize, is_paired_end: bool) -> String {
        format!("Use {} read{}s. Discard {} poorly mapped, {} out of bounds and {} with few unique k-mers",
            use_reads, if is_paired_end { " pair" } else { "" },
            self.poorly_mapped, self.out_of_bounds, self.few_kmers)
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
        kmer_counts: &KmerCounts,
        bg_distr: &BgDistr,
        edit_dist_cache: &EditDistCache,
        contig_infos: &ContigInfos,
        params: &super::Params,
        mut dbg_writer1: impl Write,
        mut dbg_writer2: impl Write,
        mut dbg_writer3: Option<GzFile>,
        opt_hap_alns: Option<&HapAlns>,
    ) -> crate::Result<Self>
    {
        let contigs = contig_set.contigs();
        let boundary = params.boundary_size.strict_sub(params.tweak.unwrap());
        assert!(contigs.lengths().iter().all(|&len| len > 2 * boundary),
            "[{}] Some contigs are too short (must be over twice boundary size = {})", contigs.tag(), 2 * boundary);
        let mut unique_kmers = UniqueKmers::new(contig_set, kmer_counts,
            params.kmer_hard_thresh, params.kmer_soft_thresh);

        log::info!("    Loading read alignments");
        writeln!(dbg_writer1, "read_hash\tread_end\tinterval\tedit_dist\tlik\tread_name")
            .map_err(add_path!(!))?;
        writeln!(dbg_writer2, "read_hash\tuniq_kmers1\tuniq_kmers2\tweight").map_err(add_path!(!))?;
        if let Some(w) = &mut dbg_writer3 {
            writeln!(w, "read_hash\tread_end\tinterval\told_lik\tnew_lik\toverall_weight\tweights")
                .map_err(add_path!(!))?;
        }
        // [TODO] Explain params.
        let aligner = Aligner::new(Default::default(), 6, Some(10), true);
        let tid2contig = construct_tid_to_contig_map(reader.header(), contig_set.contigs())?;
        let mut reader = LaggedReader::new(reader)?;
        let data = Data::new(contig_set, contig_infos, tid2contig, bg_distr, edit_dist_cache, params)?;
        let is_paired_end = bg_distr.insert_distr().is_paired_end();
        let mut hashes = IntSet::default();
        let mut buffer = Vec::with_capacity(16);
        let mut collisions = 0;

        let mut counts = ReadCounts::default();
        let mut reads = Vec::new();
        let mut unused_reads = Vec::new();

        let mut total_alns = 0;
        let mut n_recovered = 0;
        while reader.has_more() {
            let mut read_data = ReadData::default();
            counts.total += 1;
            let mut prelim_alignments = PrelimAlignments::new();
            let mut well_mapped = read_next_alns(&data, &mut reader, ReadEnd::First,
                &mut read_data, &mut prelim_alignments)?;
            if !hashes.insert(read_data.name_hash) {
                log::debug!("Read {} produced hash collision ({})", read_data.name, read_data.name_hash);
                collisions += 1;
            }
            if is_paired_end {
                if well_mapped {
                    well_mapped = read_next_alns(&data, &mut reader, ReadEnd::Second,
                        &mut read_data, &mut prelim_alignments)?;
                } else {
                    reader.skip_until_primary()?;
                }
            }
            if !well_mapped {
                counts.poorly_mapped += 1;
                continue;
            } else if !in_bounds(&prelim_alignments.alns, boundary, contigs) {
                counts.out_of_bounds += 1;
                continue;
            }

            prelim_alignments.debug_write(&read_data, &mut dbg_writer1)?;
            read_data.weight *= unique_kmers.read_weight(&mut read_data, &mut dbg_writer2).map_err(add_path!(!))?;
            let max_alns = if read_data.weight >= params.min_weight { MAX_USED_ALNS } else { MAX_UNUSED_ALNS };

            if let Some(hap_alns) = opt_hap_alns {
                n_recovered += hap_alns.transfer_alignments(
                    &mut prelim_alignments, &read_data, contig_set, &aligner, &data.err_prof);
            }
            total_alns += prelim_alignments.len();
            prelim_alignments.finalize();
            let groupped_alns = if is_paired_end {
                identify_paired_end_alignments(read_data, &mut prelim_alignments.alns,
                    max_alns, &mut buffer, bg_distr.insert_distr(), contig_infos, params)
            } else {
                identify_single_end_alignments(read_data, &mut prelim_alignments.alns,
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
        if n_recovered > 0 {
            log::debug!("    Loaded {} alignments, recovered additional {}", total_alns - n_recovered, total_alns)
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
