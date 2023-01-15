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
        // interv::Interval,
        // cigar::Cigar,
        aln::{READ_ENDS, ReadEnd, Alignment},
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
#[derive(Clone, Debug)]
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
            log::debug!("Read {}  {}", String::from_utf8_lossy(&record.qname()), hash);
            log::debug!("    {}", ext_aln);
            match self.alns.entry(hash) {
                Entry::Occupied(mut entry) => entry.get_mut().push(ext_aln),
                Entry::Vacant(entry) => entry.insert(Vec::new()).push(ext_aln),
            };
        }
    }

    /// Push multiple records in the `PrelimLocations`.
    pub fn extend<'b>(&mut self, records: impl Iterator<Item = &'b Record>) {
        for record in records {
            self.push(record);
        }
    }

    // /// Sorts alignments, first by read hash, then by contig id, then by read mate.
    // pub fn sort(&mut self) {
    //     self.alns.sort_by_key(|ext_aln| (ext_aln.name_hash, ext_aln.ix()))
    // }

    // pub fn identify_locations(&mut self) -
}

fn extend_pair_alignments<D: InsertDistr>(new_alns: &mut Vec<PairAlignment>, alns1: &[ExtAln], alns2: &[ExtAln],
        norm1: f64, norm2: f64, unmapped_penalty: f64, insert_distr: &D) {
    let thresh = new_alns[0].ln_prob;
    if alns1.len() > 0 && alns2.len() > 0 {
        for aln1 in alns1.iter() {
            for aln2 in alns2.iter() {
                let new_prob = aln1.ln_prob - norm1 + aln2.ln_prob - norm2
                    + insert_distr.pair_aln_prob(&aln1.aln, &aln2.aln);
                if new_prob > thresh {
                    new_alns.push(PairAlignment::new(Some(aln1.aln.clone()), Some(aln2.aln.clone()), new_prob));
                }
            }
        }
    } else if alns1.len() > 0 {
        for aln1 in alns1.iter() {
            let new_prob = aln1.ln_prob - norm1 + unmapped_penalty - norm2;
            if new_prob > thresh {
                new_alns.push(PairAlignment::new(Some(aln1.aln.clone()), None, new_prob));
            }
        }
    } else if alns2.len() > 0 {
        for aln2 in alns2.iter() {
            let new_prob = unmapped_penalty - norm1 + aln2.ln_prob - norm2;
            if new_prob > thresh {
                new_alns.push(PairAlignment::new(None, Some(aln2.aln.clone()), new_prob));
            }
        }
    } else {
        panic!("Read alignments are unavailable to a contig for both read mates.");
    }
}

/// Find all pair-alignments for a single read pair.
/// `unmapped_penalty` is in the ln-space.
fn identify_pair_alignments<D: InsertDistr>(alns: &mut [ExtAln], unmapped_penalty: f64, insert_distr: &D)
        -> Vec<PairAlignment> {
    // Sort alignments first by read-end, then by contig id.
    alns.sort_by_key(|ext_aln| (ext_aln.read_end, ext_aln.contig_id()));
    // Find the first alignment corresponding to the second read mate.
    let m = bisect::left_by(alns, |aln| aln.read_end.cmp(&ReadEnd::Second));
    let n = alns.len();
    // All first mate alignments are in 0..m, second mate alignments are in m..n.
    debug_assert!(m == 0 || alns[m - 1].read_end == ReadEnd::First);
    debug_assert!(m == n || alns[m].read_end == ReadEnd::Second);

    // Normalizing factors for first-mate and second-mate alignment probabilities.
    let norm1 = if m > 0 {
        Ln::add(Ln::map_sum(&alns[..m], ExtAln::ln_prob), unmapped_penalty)
    } else { unmapped_penalty };
    let norm2 = if m < n {
        Ln::add(Ln::map_sum(&alns[m..], ExtAln::ln_prob), unmapped_penalty)
    } else { unmapped_penalty };
    // TODO: Should we penalize read pairs where one mate does not map anywhere at all?

    // Probability of pair-read alignment where both reads are unmapped.
    let mut pair_alns = vec![PairAlignment::new_unmapped(2.0 * unmapped_penalty - norm1 - norm2)];
    // For the current contigs, alns[i1..j1]: first mate alignments, alns[i2..j2]: sceond mate alignments.
    let mut i1 = 0;
    let mut i2 = m;
    loop {
        // Assume that there are at most 4 alignments of the same read mate to the same contig.
        // (not critical, if this assumption fails).
        const STEP: usize = 4;
        let contig_id = match (i1 < m, i2 < n) {
            (false, false) => break,
            (true,  false) => alns[i1].contig_id(),
            (false, true ) => alns[i2].contig_id(),
            (true,  true ) => min(alns[i1].contig_id(), alns[i2].contig_id()),
        };
        let j1 = bisect::right_by_approx(alns, |aln| aln.contig_id().cmp(&contig_id), i1, i1 + STEP, m);
        debug_assert!(i1 <= j1 && j1 <= m);
        let j2 = bisect::right_by_approx(alns, |aln| aln.contig_id().cmp(&contig_id), i2, i2 + STEP, n);
        debug_assert!(i2 <= j2 && j2 <= m);
        
        extend_pair_alignments(&mut pair_alns, &alns[i1..j1], &alns[i2..j2],
            norm1, norm2, unmapped_penalty, insert_distr);
    }
    pair_alns
}

/// Alignment of the read pair.
/// One or both of the alignments may be missing.
/// If present, both alignments must map to the same contig and be relatively close to each other.
#[derive(Clone, Debug)]
pub struct PairAlignment {
    first: Option<Alignment>,
    second: Option<Alignment>,
    ln_prob: f64,
}

impl PairAlignment {
    pub fn new(first: Option<Alignment>, second: Option<Alignment>, ln_prob: f64) -> Self {
        Self { first, second, ln_prob }
    }

    /// Creates a pair alignment, where both mates are unmapped.
    pub fn new_unmapped(ln_prob: f64) -> Self {
        Self {
            first: None,
            second: None,
            ln_prob,
        }
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
}

pub struct ReadPairLocations {
    /// Keys: read hashes, values: all pair-alignments for the read.
    /// Vectors are sorted by contig id.
    /// First vector entry: both read mates are unmapped.
    alns: IntMap<Vec<PairAlignment>>,
}

impl ReadPairLocations {
    pub fn new(capacity: usize) -> Self {
        Self {
            alns: IntMap::with_capacity(capacity),
        }
    }


}