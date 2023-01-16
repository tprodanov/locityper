use std::{
    rc::Rc,
    fmt,
    cmp::min,
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

    /// Returns log-probability of the alignment.
    fn ln_prob(&self) -> f64 {
        self.ln_prob
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
fn extend_pair_alignments<D>(
        new_alns: &mut Vec<PairAlignment>,
        alns1: &[ExtAln], alns2: &[ExtAln],
        norm1: f64, norm2: f64, unmapped_penalty: f64, thresh_prob: f64,
        insert_distr: &D,
    )
where D: InsertDistr,
{
    if alns1.len() > 0 && alns2.len() > 0 {
        for aln1 in alns1.iter() {
            for aln2 in alns2.iter() {
                let new_prob = aln1.ln_prob - norm1 + aln2.ln_prob - norm2
                    + insert_distr.pair_aln_prob(&aln1.aln, &aln2.aln);
                log::debug!("        prob1 {:.1}   prob2 {:.1}   insert {:.1}  -> {:.1}",
                    aln1.ln_prob, aln2.ln_prob, insert_distr.pair_aln_prob(&aln1.aln, &aln2.aln), new_prob);
                if new_prob > thresh_prob {
                    new_alns.push(PairAlignment::new(Some(aln1.aln.clone()), Some(aln2.aln.clone()), new_prob));
                    log::debug!("        Push {}", new_alns.last().unwrap());
                }
            }
        }
    } else if alns1.len() > 0 {
        for aln1 in alns1.iter() {
            let new_prob = aln1.ln_prob - norm1 + unmapped_penalty - norm2;
            if new_prob > thresh_prob {
                new_alns.push(PairAlignment::new(Some(aln1.aln.clone()), None, new_prob));
                log::debug!("        Push {}", new_alns.last().unwrap());
            }
        }
    } else if alns2.len() > 0 {
        for aln2 in alns2.iter() {
            let new_prob = unmapped_penalty - norm1 + aln2.ln_prob - norm2;
            if new_prob > thresh_prob {
                new_alns.push(PairAlignment::new(None, Some(aln2.aln.clone()), new_prob));
                log::debug!("        Push {}", new_alns.last().unwrap());
            }
        }
    } else {
        panic!("Both read mates have no alignments to a certain contig.");
    }
}

fn identify_pair_alignments<D: InsertDistr>(alns: &mut [ExtAln], unmapped_penalty: f64, insert_distr: &D)
        -> ReadPairAlignments {
    // Sort alignments first by read-end, then by contig id.
    alns.sort_by_key(|ext_aln| (ext_aln.read_end, ext_aln.contig_id()));
    // Find the first alignment corresponding to the second read mate.
    let m = bisect::left_by(alns, |aln| aln.read_end.cmp(&ReadEnd::Second));
    let n = alns.len();
    // All first mate alignments are in 0..m, second mate alignments are in m..n.
    debug_assert!(m == 0 || alns[m - 1].read_end == ReadEnd::First);
    debug_assert!(m == n || alns[m].read_end == ReadEnd::Second);
    log::debug!("    {} first alns:", m);
    for i in 0..m {
        log::debug!("        {:3} {:?}", i, alns[i]);
    }
    log::debug!("    {} second alns:", n - m);
    for i in m..n {
        log::debug!("        {:3} {:?}", i, alns[i]);
    }

    // Normalizing factors for first-mate and second-mate alignment probabilities.
    let norm1 = if m > 0 {
        Ln::add(Ln::map_sum(&alns[..m], ExtAln::ln_prob), unmapped_penalty)
    } else { unmapped_penalty };
    let norm2 = if m < n {
        Ln::add(Ln::map_sum(&alns[m..], ExtAln::ln_prob), unmapped_penalty)
    } else { unmapped_penalty };
    // TODO: Should we penalize read pairs where one mate does not map anywhere at all?
    
    let unmapped_prob = 2.0 * unmapped_penalty - norm1 - norm2;
    log::debug!("    Norm factors {:.2}   {:.2}.   Unmapped prob: {:.1}", norm1, norm2, unmapped_prob);
    let mut pair_alns = Vec::new();
    // For the current contigs, alns[i1..j1]: first mate alignments, alns[i2..j2]: sceond mate alignments.
    let mut i1 = 0;
    let mut i2 = m;
    loop {
        let contig_id = match (i1 < m, i2 < n) {
            (false, false) => break,
            (true,  false) => alns[i1].contig_id(),
            (false, true ) => alns[i2].contig_id(),
            (true,  true ) => min(alns[i1].contig_id(), alns[i2].contig_id()),
        };
        log::debug!("    Searching for contig {},  i1 = {}, i2 = {}", contig_id, i1, i2);
        let j1 = bisect::right_by_approx(alns, |aln| aln.contig_id().cmp(&contig_id), i1, m, BISECT_RIGHT_STEP);
        debug_assert!(i1 <= j1 && j1 <= m);
        let j2 = bisect::right_by_approx(alns, |aln| aln.contig_id().cmp(&contig_id), i2, n, BISECT_RIGHT_STEP);
        debug_assert!(i2 <= j2 && j2 <= n);
        log::debug!("    Contig #{}:   {}..{}   {}..{}", contig_id, i1, j1, i2, j2);

        extend_pair_alignments(&mut pair_alns, &alns[i1..j1], &alns[i2..j2],
            norm1, norm2, unmapped_penalty, unmapped_prob, insert_distr);
        i1 = j1;
        i2 = j2;
    }
    ReadPairAlignments { unmapped_prob, pair_alns }
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