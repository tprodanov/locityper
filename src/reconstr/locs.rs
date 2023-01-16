use std::{
    rc::Rc,
    fmt,
};
use htslib::bam::Record;
use intmap::{IntMap, Entry};
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
        self.aln.ref_interval().contig_id()
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

/// Preliminary read locations for multiple reads.
pub struct PrelimLocations<'a, E: ErrorProfile> {
    contigs: Rc<ContigNames>,
    err_prof: &'a E,
    min_ln_prob: f64,
    /// All alignments with probability over `min_ln_prob`.
    /// Keys: read name hashes, values: all alignments for the read.
    alns: IntMap<Vec<ExtAln>>,
}

impl<'a, E: ErrorProfile> PrelimLocations<'a, E> {
    /// Creates a new, empty, `PrelimLocations`.
    pub fn new(contigs: Rc<ContigNames>, err_prof: &'a E, min_ln_prob: f64) -> Self {
        Self {
            contigs, err_prof, min_ln_prob,
            alns: IntMap::new(),
        }
    }

    /// Create a new `ExtAln` from the record, and push it in the vector,
    /// if alignment probability is over `min_ln_prob`.
    pub fn push(&mut self, record: &Record) {
        let ext_aln = ExtAln::from_record(record, Rc::clone(&self.contigs), self.err_prof);
        if ext_aln.ln_prob >= self.min_ln_prob {
            let hash = fnv1a(record.qname());

            log::debug!("Read {}  {}  {}", String::from_utf8_lossy(&record.qname()), hash, ext_aln);
            match self.alns.entry(hash) {
                Entry::Occupied(mut entry) => entry.get_mut().push(ext_aln),
                Entry::Vacant(entry) => entry.insert(Vec::new()).push(ext_aln),
            };
        }
    }

    /// Push multiple records in the `PrelimLocations`.
    pub fn extend(&mut self, records: impl Iterator<Item = Record>) {
        for record in records {
            if !record.is_unmapped() {
                self.push(&record);
            }
        }
    }

    pub fn identify_locations<D>(&mut self, unmapped_penalty: f64, insert_distr: &D) -> AllPairAlignments
    where D: InsertDistr
    {
        log::info!("Identify paired alignment location and probabilities ({} read pairs)", self.alns.len());
        let mut res = AllPairAlignments {
            alns: IntMap::with_capacity(self.alns.len()),
        };
        for (&key, values) in self.alns.iter_mut() {
            log::debug!("Read {}", key);
            res.alns.insert(key, identify_pair_alignments(values, unmapped_penalty, insert_distr));
        }
        res
    }
}

// Assume that there are at most 4 alignments of the to the same contig.
// (not critical, if this assumption fails).
const BISECT_RIGHT_STEP: usize = 4;

/// For a single read-pair, find all paired-read alignments to the same contig.
fn extend_pair_alignments<D>(new_alns: &mut Vec<PairAlignment>,
        alns1: &[ExtAln], alns2: &[ExtAln], unmapped_penalty: f64, insert_distr: &D,
    )
where D: InsertDistr,
{
    let thresh_prob = unmapped_penalty * 2.0;
    if alns1.len() > 0 && alns2.len() > 0 {
        for aln1 in alns1.iter() {
            for aln2 in alns2.iter() {
                let new_prob = aln1.ln_prob + aln2.ln_prob + insert_distr.pair_aln_prob(&aln1.aln, &aln2.aln);
                if new_prob > thresh_prob {
                    new_alns.push(PairAlignment::new(Some(aln1.aln.clone()), Some(aln2.aln.clone()), new_prob));
                }
            }
        }
    } else if alns1.len() > 0 {
        for aln1 in alns1.iter() {
            let new_prob = aln1.ln_prob + unmapped_penalty;
            if new_prob > thresh_prob {
                new_alns.push(PairAlignment::new(Some(aln1.aln.clone()), None, new_prob));
            }
        }
    } else if alns2.len() > 0 {
        for aln2 in alns2.iter() {
            let new_prob = aln2.ln_prob + unmapped_penalty;
            if new_prob > thresh_prob {
                new_alns.push(PairAlignment::new(None, Some(aln2.aln.clone()), new_prob));
            }
        }
    } else {
        panic!("Both read mates have no alignments to a certain contig.");
    }
}

fn identify_pair_alignments<D: InsertDistr>(alns: &mut [ExtAln], unmapped_penalty: f64, insert_distr: &D)
        -> ReadPairAlignments {
    // Sort alignments first by contig id, then by read-end.
    alns.sort_by_key(ExtAln::sort_key);
    let n = alns.len();
    // For the current contig id, first mates will be in i..j, and second mates in j..k.
    let mut i = 0;
    let mut pair_alns = Vec::new();
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
        extend_pair_alignments(&mut pair_alns, &alns[i..j], &alns[j..k], unmapped_penalty, insert_distr);
        i = k;
    }
    // Normalization factor for all pair-end alignment probabilities.
    let norm_fct = Ln::add(2.0 * unmapped_penalty, Ln::map_sum(&pair_alns, PairAlignment::ln_prob));
    log::debug!("    {} pair-end alignments, unmapped prob. = {:.2}",
        pair_alns.len(), 2.0 * unmapped_penalty - norm_fct);
    for aln in pair_alns.iter_mut() {
        aln.ln_prob -= norm_fct;
        log::debug!("        {:?}", aln);
    }
    ReadPairAlignments {
        unmapped_prob: 2.0 * unmapped_penalty - norm_fct,
        pair_alns,
    }
}

/// Alignment of the read pair. At most one of two alignments may be missing!
/// If present, both alignments must map to the same contig and be relatively close to each other.
#[derive(Clone)]
pub struct PairAlignment {
    first: Option<Alignment>,
    second: Option<Alignment>,
    ln_prob: f64,
}

impl PairAlignment {
    pub fn new(first: Option<Alignment>, second: Option<Alignment>, ln_prob: f64) -> Self {
        assert!(first.is_some() || second.is_some(), "PairAlignment: at least alignment must be present!");
        Self { first, second, ln_prob }
    }

    /// Returns first-mate alignment, if present.
    pub fn first(&self) -> &Option<Alignment> {
        &self.first
    }

    /// Returns true if the first mate is mapped.
    pub fn first_mapped(&self) -> bool {
        self.first.is_some()
    }

    /// Returns second-mate alignment, if present.
    pub fn second(&self) -> &Option<Alignment> {
        &self.second
    }

    /// Returns true if the second mate is mapped.
    pub fn second_mapped(&self) -> bool {
        self.second.is_some()
    }

    /// Log-probability of the paired-read alignment.
    pub fn ln_prob(&self) -> f64 {
        self.ln_prob
    }

    pub fn contig_id(&self) -> ContigId {
        self.first.as_ref().or(self.second.as_ref())
            .expect("One of the alignments must be present!")
            .ref_interval().contig_id()
    }
}

impl fmt::Debug for PairAlignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Pair({:?}, {:?}, {:.2})", self.first, self.second, self.ln_prob)
    }
}

impl fmt::Display for PairAlignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Pair(")?;
        if let Some(aln) = &self.first {
            write!(f, "{}, ", aln)?;
        } else {
            write!(f, "*, ")?;
        }
        if let Some(aln) = &self.second {
            write!(f, "{}, ", aln)?;
        } else {
            write!(f, "*, ")?;
        }
        write!(f, "prob={:.2})", self.ln_prob)
    }
}

/// Read-pair alignments for a single read-pair.
pub struct ReadPairAlignments {
    /// Probability that both read mates are unaligned.
    unmapped_prob: f64,

    /// All pair-alignments for the read. One of the mates may be unaligned, but not both.
    /// Vector is sorted by contig id.
    pair_alns: Vec<PairAlignment>,
}

impl ReadPairAlignments {
    /// For a given contig, returns:
    /// - all alignments to this contig (may be empty),
    /// - probability that the read pair is not mapped to this contig (may be -INF).
    pub fn contig_alns(&self, contig_id: ContigId) -> (&[PairAlignment], f64) {
        let i = bisect::left_by(&self.pair_alns, |paln| paln.contig_id().cmp(&contig_id));
        let j = bisect::right_by_approx(&self.pair_alns, |paln| paln.contig_id().cmp(&contig_id),
            i, self.pair_alns.len(), BISECT_RIGHT_STEP);
        (&self.pair_alns[i..j], if i == j { self.unmapped_prob } else { f64::NEG_INFINITY })
    }
}

/// All read-pair alignments for all read-pairs.
pub struct AllPairAlignments {
    /// Key: read name hash.
    alns: IntMap<ReadPairAlignments>,
}

impl AllPairAlignments {
}