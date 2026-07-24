use std::{
    thread, io,
    sync::{Arc, OnceLock},
    cmp::Ordering,
    fmt::Write as FmtWrite,
    ops::Deref,
};
use arc_swap::ArcSwapOption;
use ruint::aliases::U256;
use smallvec::SmallVec;
use crate::{
    seq::{
        wfa, minim_div,
        contigs::{ContigId, ContigSet},
        kmers::{self, Kmer},
        wfa::{Aligner, Penalties},
        cigar::{Cigar, Operation},
    },
    err::{error, validate_param, add_path},
    math::RoundDiv,
    ext::TriangleMatrix,
};

// ========== Parameters ==========

/// Alignment/divergence calculation parameters.
#[derive(Clone)]
pub struct Params {
    pub skip_div: bool,
    pub div_k: u8,
    pub div_w: u8,
    /// Do not align sequences with higher minimizer divergence.
    pub thresh_div: f64,
    /// Same as thresh_div, but for specified `--against` targets.
    pub against_div: f64,
    pub penalties: Penalties,
    pub backbone_ks: Vec<u8>,
    pub accuracy: u8,
    pub max_gap: u32,

    /// For all-v-all alignments, transfer alignments
    /// if one of the two constructed alignments has divergence <= this value.
    pub transitive_div: f64,
    /// When transferring alignments, use regions of this or larger size where one or both of the CIGARs have "=".
    pub transitive_anchor: u32,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            skip_div: false,
            div_k: 15,
            div_w: 15,
            thresh_div: 0.5,
            against_div: 1.0,
            penalties: Default::default(),
            backbone_ks: vec![25, 51, 101],
            accuracy: 9,
            max_gap: 500,
            transitive_div: 0.02,
            transitive_anchor: 101,
        }
    }
}

impl Params {
    pub fn validate(&mut self) -> crate::Result<()> {
        validate_param!(0 < self.div_k && self.div_k <= u64::MAX_KMER_SIZE,
            "k-mer size ({}) must be between 1 and {}", self.div_k, u64::MAX_KMER_SIZE);
        validate_param!(0 < self.div_w && self.div_w <= kmers::MAX_MINIMIZER_W,
            "Minimizer window ({}) must be between 1 and {}", self.div_w, kmers::MAX_MINIMIZER_W);
        validate_param!(0.0 <= self.thresh_div && self.thresh_div <= 1.0,
            "Maximum divergence ({}) must be within [0, 1]", self.thresh_div);
        validate_param!(1 <= self.accuracy && self.accuracy <= wfa::MAX_ACCURACY,
            "Alignment accuracy level ({}) must be between 0 and {}.", self.accuracy, wfa::MAX_ACCURACY);

        if self.thresh_div == 0.0 {
            // Never run alignment.
            self.thresh_div = -1.0;
            self.backbone_ks.clear();
        } else {
            validate_param!(!self.backbone_ks.is_empty(), "Expect at least one backbone k-mer");
            validate_param!(self.backbone_ks.iter().all(|&k| 5 <= k && k <= U256::MAX_KMER_SIZE),
                "Backbone k-mer sizes must be between 5 and {} ({})", U256::MAX_KMER_SIZE,
                self.backbone_str());
        }
        self.penalties.validate()?;
        Ok(())
    }

    /// Combine backbone k via comma.
    pub fn backbone_str(&self) -> String {
        self.backbone_ks.iter().map(u8::to_string).collect::<Vec<_>>().join(",")
    }
}

// ========== Backbone k-mers ==========

const CAPACITY: usize = 4;
type SeqKmers = Vec<(U256, SmallVec<[u32; CAPACITY]>)>;

fn precompute_kmers(seq: &[u8], k: u8, buf: &mut Vec<(u32, U256)>) -> SeqKmers {
    buf.clear();
    kmers::kmers::<U256, _, { kmers::NON_CANONICAL }>(seq, k, buf);
    buf.sort_by(|a, b| a.1.cmp(&b.1));

    // k-mers with combined positions.
    let mut compressed_kmers: SeqKmers = Vec::with_capacity(buf.len());
    let mut curr_kmer = buf[0].1;
    let mut curr_positions = SmallVec::new();
    for &(pos, kmer) in buf.iter() {
        if kmer != curr_kmer {
            compressed_kmers.push((curr_kmer, std::mem::take(&mut curr_positions)));
            curr_kmer = kmer;
        }
        curr_positions.push(pos);
    }
    compressed_kmers.push((curr_kmer, curr_positions));
    compressed_kmers
}

fn fill_kmers_singlethread(
    contig_set: &ContigSet,
    in_use: &[bool],
    start: usize,
    end: usize,
    params: &Params,
) -> (Vec<Vec<u64>>, Vec<SeqKmers>) {
    let n = end - start;
    let m = params.backbone_ks.len();
    let mut minimizers = Vec::with_capacity(if params.skip_div { 0 } else { n });
    // Backbone ks will be empty if no alignments need to be calculated.
    let mut kmers = Vec::with_capacity(n * m);

    let mut kmer_buf = Vec::new();
    for (use_seq, seq) in itertools::izip!(&in_use[start..end], &contig_set.seqs()[start..end]) {
        if !use_seq {
            if !params.skip_div {
                minimizers.push(Vec::new());
            }
            kmers.extend((0..m).map(|_| Vec::new()));
            continue;
        }
        if !params.skip_div {
            // Expected num of minimizers = 2L / (w + 1), here we take 5/2 * ... more to be safe.
            let mut buf = Vec::with_capacity((5 * seq.len()).fast_round_div(2 * usize::from(params.div_w) + 2));
            kmers::minimizers::<u64, _, { kmers::NON_CANONICAL }>(seq, params.div_k, params.div_w, &mut buf);
            buf.sort_unstable();
            minimizers.push(buf);
        }
        for &k in &params.backbone_ks {
            kmers.push(precompute_kmers(seq, k, &mut kmer_buf));
        }
    }
    (minimizers, kmers)
}

fn fill_kmers(
    contig_set: &Arc<ContigSet>,
    in_use: &[bool],
    params: &Params,
    threads: usize,
) -> (Vec<Vec<u64>>, Vec<SeqKmers>) {
    let total = in_use.iter().copied().map(u32::from).sum::<u32>();
    let n_contigs = contig_set.len();
    let mut end0 = n_contigs;
    let mut thread_ranges = Vec::with_capacity(threads - 1);
    if threads != 1 && total >= 16 {
        let mut start = 0;
        for worker_ix in 0..threads {
            if start == n_contigs {
                break;
            }
            let end = start + (n_contigs - start).fast_ceil_div(threads - worker_ix);
            assert!(start < end);
            if worker_ix == 0 {
                end0 = end;
            } else {
                thread_ranges.push((start, end));
            }
            start = end;
        }
        assert!(start == n_contigs);
    }
    let mut handles = Vec::with_capacity(thread_ranges.len());
    for (start, end) in thread_ranges {
        let contig_set = Arc::clone(contig_set);
        let params = params.clone();
        let in_use = in_use.to_vec();
        handles.push(thread::spawn(move || fill_kmers_singlethread(&contig_set, &in_use, start, end, &params)));
    }
    let (mut minimizers, mut kmers) = fill_kmers_singlethread(&contig_set, &in_use, 0, end0, &params);

    for handle in handles {
        let (curr_minimizers, curr_kmers) = handle.join().expect("Thread unexpectedly terminated");
        minimizers.extend(curr_minimizers.into_iter());
        kmers.extend(curr_kmers.into_iter());
    }
    (minimizers, kmers)
}

fn get_kmer_matches(kmers1: &SeqKmers, kmers2: &SeqKmers, buf: &mut Vec<(u32, u32)>) {
    buf.clear();
    let mut iter1 = kmers1.iter();
    let mut iter2 = kmers2.iter();
    let mut opt_x = iter1.next();
    let mut opt_y = iter2.next();
    while let (Some(x), Some(y)) = (opt_x, opt_y) {
        match x.0.cmp(&y.0) {
            Ordering::Equal => {
                for &pos1 in &x.1 {
                    for &pos2 in &y.1 {
                        buf.push((pos1, pos2));
                    }
                }
                opt_x = iter1.next();
                opt_y = iter2.next();
            }
            Ordering::Less => opt_x = iter1.next(),
            Ordering::Greater => opt_y = iter2.next(),
        }
    }
    buf.sort_unstable();
}

// ========== Direct alignment from k-mer backbones ==========

#[derive(Clone, Copy)]
struct RefEntry<'a> {
    id: ContigId,
    name: &'a str,
    seq: &'a [u8],
}

impl<'a> RefEntry<'a> {
    #[inline]
    fn from_set(contig_set: &'a ContigSet, id: ContigId) -> Self {
        Self {
            id,
            name: contig_set.contigs().get_name(id),
            seq: contig_set.get_seq(id),
        }
    }
}

fn align_from_backbone(
    aligner: &Aligner,
    entry1: RefEntry<'_>,
    entry2: RefEntry<'_>,
    kmer_matches: &[(u32, u32)],
    backbone_k: u32,
    max_gap: u32,
) -> crate::Result<(Cigar, i32)>
{
    let sparse_aln = bio::alignment::sparse::lcskpp(&kmer_matches, backbone_k as usize);
    let mut cigar = Cigar::new();
    let mut score = 0;
    let mut i1 = 0;
    let mut j1 = 0;
    let mut curr_match = 0;

    for ix in sparse_aln.path.into_iter() {
        let (i2, j2) = kmer_matches[ix];
        if i1 > i2 {
            curr_match += 1;
            i1 += 1;
            j1 += 1;
            continue;
        }

        if curr_match > 0 {
            cigar.push_unchecked(Operation::Equal, curr_match);
            curr_match = 0;
        }
        score += aligner.smart_align(entry1.seq, i1, i2, entry2.seq, j1, j2, max_gap, &mut cigar);
        curr_match += backbone_k;
        i1 = i2 + backbone_k;
        j1 = j2 + backbone_k;
    }

    if curr_match > 0 {
        cigar.push_unchecked(Operation::Equal, curr_match);
    }
    let n1 = entry1.seq.len() as u32;
    let n2 = entry2.seq.len() as u32;
    score += aligner.smart_align(entry1.seq, i1, n1, entry2.seq, j1, n2, max_gap, &mut cigar);
    assert_eq!(cigar.ref_len(), entry1.seq.len() as u32,
        "Alignment {} - {} produced incorrect CIGAR {}", entry1.name, entry2.name, cigar);
    assert_eq!(cigar.query_len(), entry2.seq.len() as u32,
        "Alignment {} - {} produced incorrect CIGAR {}", entry1.name, entry2.name, cigar);
    Ok((cigar, score))
}

fn align_multik(
    aligner: &Aligner,
    entry1: RefEntry<'_>,
    entry2: RefEntry<'_>,
    kmers: &[SeqKmers],
    params: &Params,
    buf: &mut Vec<(u32, u32)>,
) -> crate::Result<(Cigar, i32)>
{
    let mut best_cigar = None;
    let mut best_score = i32::MIN;
    let n_backbones = params.backbone_ks.len();
    for (&k, kmers1, kmers2) in itertools::izip!(
            &params.backbone_ks, &kmers[n_backbones * entry1.id.ix()..], &kmers[n_backbones * entry2.id.ix()..])
    {
        get_kmer_matches(kmers1, kmers2, buf);
        let (cigar, score) = align_from_backbone(aligner, entry1, entry2, buf, u32::from(k), params.max_gap)?;
        if score > best_score {
            best_score = score;
            best_cigar = Some(cigar);
        }
    }
    best_cigar.ok_or_else(|| error!(RuntimeError, "No alignment found between {} and {}", entry1.name, entry2.name))
        .map(|cigar| (cigar, best_score))
}

// ========== Alignment strategy: backbone/transitive alignment ==========

#[derive(Default, Clone)]
struct Counts {
    aligned: u32,
    accelerated: u32,
    skipped: u32,
}

impl Counts {
    fn add(&mut self, oth: &Counts) {
        self.aligned += oth.aligned;
        self.accelerated += oth.accelerated;
        self.skipped += oth.skipped;
    }

    fn summarize(&self) {
        let mut s = format!("Constructed {} alignments", self.aligned);
        if self.accelerated > 0 {
            write!(s, ", accelerated {:.1}% of them",
                100.0 * f64::from(self.accelerated) / f64::from(self.aligned)).unwrap();
        }
        if self.skipped > 0 {
            write!(s, ", skipped {} pairs due to high minimizer divergence", self.skipped).unwrap();
        }
        log::debug!("{}", s);
    }
}

trait Strategy {
    /// Aligns two sequences (entry1 = reference, entry2 = query).
    /// Returns CIGAR, its score, and flag indicating whether an alignment was accelerated (e.g. transitive alignment).
    fn align<'a>(
        &mut self,
        aligner: &Aligner,
        entry1: RefEntry<'a>,
        entry2: RefEntry<'a>,
        params: &Params,
    ) -> crate::Result<(Cigar, i32, bool)>;

    fn save_cigar(
        &mut self,
        _ref_id: ContigId,
        _query_id: ContigId,
        _cigar: &Cigar,
        _div: f64,
        _params: &Params,
    ) {}
}

trait StrategyFactory {
    type Strategy: Strategy + Send + 'static;

    fn new_strategy(&self) -> Self::Strategy;
}

// ========== Backbone-based alignment strategy ==========

struct BackboneStrategy {
    kmers: Arc<Vec<SeqKmers>>,
    buf: Vec<(u32, u32)>,
}

impl BackboneStrategy {
    fn new(kmers: Arc<Vec<SeqKmers>>) -> Self {
        Self {
            kmers,
            buf: Vec::new(),
        }
    }
}

impl Strategy for BackboneStrategy {
    #[inline(always)]
    fn align<'a>(
        &mut self,
        aligner: &Aligner,
        entry1: RefEntry<'a>,
        entry2: RefEntry<'a>,
        params: &Params,
    ) -> crate::Result<(Cigar, i32, bool)> {
        align_multik(aligner, entry1, entry2, &self.kmers, params, &mut self.buf)
            .map(|(cigar, score)| (cigar, score, false))
    }
}

struct BackboneStrategyFactory {
    kmers: Arc<Vec<SeqKmers>>,
}

impl BackboneStrategyFactory {
    #[inline]
    fn new(kmers: Arc<Vec<SeqKmers>>) -> Self {
        Self { kmers }
    }
}

impl StrategyFactory for BackboneStrategyFactory {
    type Strategy = BackboneStrategy;

    #[inline]
    fn new_strategy(&self) -> Self::Strategy {
        BackboneStrategy::new(Arc::clone(&self.kmers))
    }
}

// ========== Transitive alignment strategy (single thread) ==========

/// CIGAR for contigs a < b (indices not stored).
/// `ref_smaller` is true when `a` is reference and `b` is query.
#[derive(Clone)]
struct DirectedCigar {
    cigar: Cigar,
    ref_smaller: bool,
}

impl DirectedCigar {
    #[inline(always)]
    fn new(cigar: Cigar, ref_id: ContigId, query_id: ContigId) -> Self {
        Self {
            cigar,
            ref_smaller: ref_id < query_id,
        }
    }

    /// Assuming both i and j are true contig indices for this CIGAR, returns true if `j` is reference and `i` is query.
    #[inline(always)]
    fn second_is_ref(&self, i: ContigId, j: ContigId) -> bool {
        (j < i) == self.ref_smaller
    }
}

struct TransitiveStrategy {
    backbone_strategy: BackboneStrategy,
    /// For each contig, store its closest other contig, corresponding CIGAR and diversity.
    /// Value is only stored if diversity is not greater than params.transitive_div
    /// Closest other contig will always have smaller index than i.
    closest: Vec<Option<(ContigId, DirectedCigar, f64)>>,
    cigars: TriangleMatrix<Option<DirectedCigar>>,
}

impl TransitiveStrategy {
    fn new(kmers: Arc<Vec<SeqKmers>>, n_contigs: usize) -> Self {
        Self {
            backbone_strategy: BackboneStrategy::new(kmers),
            closest: vec![None; n_contigs],
            cigars: TriangleMatrix::new(n_contigs, None),
        }
    }
}

impl Strategy for TransitiveStrategy {
    fn align<'a>(
        &mut self,
        aligner: &Aligner,
        entry_k: RefEntry<'a>,
        entry_i: RefEntry<'a>,
        params: &Params,
    ) -> crate::Result<(Cigar, i32, bool)> {
        let k = entry_k.id;
        let i = entry_i.id;
        let transitive_edge = if let Some((j, cigar_jk, _)) = &self.closest[k.ix()]
            && let Some(cigar_ij) = &self.cigars.get_symmetric(i.ix(), j.ix())
        {
            Some((cigar_ij, *j, cigar_jk))
        } else if let Some((j, cigar_ij, _)) = &self.closest[i.ix()]
            && let Some(cigar_jk) = &self.cigars.get_symmetric(k.ix(), j.ix())
        {
            Some((cigar_ij, *j, cigar_jk))
        } else {
            None
        };

        if let Some((cigar_ij, j, cigar_jk)) = transitive_edge {
            let cigar_ik = Cigar::find_transitive_alignment(
                &cigar_ij.cigar, cigar_ij.second_is_ref(i, j), &cigar_jk.cigar, cigar_jk.second_is_ref(j, k),
                entry_i.seq, entry_k.seq, aligner, params.max_gap, params.transitive_anchor);
            let score = aligner.penalties().calculate_score(&cigar_ik);
            Ok((cigar_ik, score, true))
        } else {
            self.backbone_strategy.align(aligner, entry_k, entry_i, params)
        }
    }

    fn save_cigar(&mut self, ref_id: ContigId, query_id: ContigId, cigar: &Cigar, div: f64, params: &Params) {
        let dir_cigar = DirectedCigar::new(cigar.to_owned(), ref_id, query_id);
        if div <= params.transitive_div {
            let prev_closest = &mut self.closest[query_id.ix()];
            if let Some((_, _, prev_div)) = prev_closest && *prev_div <= div { /* Do nothing */ } else {
                *prev_closest = Some((ref_id, dir_cigar.clone(), div));
            }
        }
        *self.cigars.get_symmetric_mut(ref_id.ix(), query_id.ix()) = Some(dir_cigar);
    }
}

// ========== Transitive alignment strategy (parallel) ==========

struct ParallelTransitiveStrategy {
    backbone_strategy: BackboneStrategy,
    /// For each contig, store its closest other contig, corresponding CIGAR and diversity.
    /// Value is only stored if diversity is not greater than params.transitive_div
    /// Closest other contig will always have smaller index than i.
    closest: Arc<Vec<ArcSwapOption<(ContigId, DirectedCigar, f64)>>>,
    cigars: Arc<TriangleMatrix<OnceLock<DirectedCigar>>>,
}

impl ParallelTransitiveStrategy {
    fn new(
        kmers: Arc<Vec<SeqKmers>>,
        closest: Arc<Vec<ArcSwapOption<(ContigId, DirectedCigar, f64)>>>,
        cigars: Arc<TriangleMatrix<OnceLock<DirectedCigar>>>,
    ) -> Self {
        Self {
            backbone_strategy: BackboneStrategy::new(kmers),
            closest, cigars,
        }
    }
}

impl Strategy for ParallelTransitiveStrategy {
    fn align<'a>(
        &mut self,
        aligner: &Aligner,
        entry_k: RefEntry<'a>,
        entry_i: RefEntry<'a>,
        params: &Params,
    ) -> crate::Result<(Cigar, i32, bool)> {
        let k = entry_k.id;
        let i = entry_i.id;
        let closest_k = self.closest[k.ix()].load();
        let closest_i = self.closest[i.ix()].load();

        let transitive_edge = if let Some((j, cigar_jk, _)) = closest_k.as_ref().map(Deref::deref)
            && let Some(cigar_ij) = &self.cigars.get_symmetric(i.ix(), j.ix()).get()
        {
            Some((*cigar_ij, *j, cigar_jk))
        } else if let Some((j, cigar_ij, _)) = closest_i.as_ref().map(Deref::deref)
            && let Some(cigar_jk) = &self.cigars.get_symmetric(k.ix(), j.ix()).get()
        {
            Some((cigar_ij, *j, *cigar_jk))
        } else {
            None
        };

        if let Some((cigar_ij, j, cigar_jk)) = transitive_edge {
            let cigar_ik = Cigar::find_transitive_alignment(
                &cigar_ij.cigar, cigar_ij.second_is_ref(i, j), &cigar_jk.cigar, cigar_jk.second_is_ref(j, k),
                entry_i.seq, entry_k.seq, aligner, params.max_gap, params.transitive_anchor);
            let score = aligner.penalties().calculate_score(&cigar_ik);
            Ok((cigar_ik, score, true))
        } else {
            self.backbone_strategy.align(aligner, entry_k, entry_i, params)
        }
    }

    fn save_cigar(&mut self, ref_id: ContigId, query_id: ContigId, cigar: &Cigar, div: f64, params: &Params) {
        let dir_cigar = DirectedCigar::new(cigar.to_owned(), ref_id, query_id);
        if div <= params.transitive_div {
            let query_closest = &self.closest[query_id.ix()];
            if query_closest.load().as_ref().map(|arc_tuple| div < arc_tuple.deref().2).unwrap_or(true) {
                // It is possible that there are two concurrent threads that check divergence and
                // write after that. So, we could end up with suboptimal closest alignment.
                // However, that should not happen too often and is not the end of the world.
                query_closest.store(Some(Arc::new((ref_id, dir_cigar.clone(), div))));
            }
        }
        let _old_val = self.cigars.get_symmetric(ref_id.ix(), query_id.ix()).set(dir_cigar);
    }
}

struct TransitiveStrategyFactory {
    kmers: Arc<Vec<SeqKmers>>,
    /// For each contig, store its closest other contig, corresponding CIGAR and diversity.
    /// Value is only stored if diversity is not greater than params.transitive_div
    /// Closest other contig will always have smaller index than i.
    closest: Arc<Vec<ArcSwapOption<(ContigId, DirectedCigar, f64)>>>,
    cigars: Arc<TriangleMatrix<OnceLock<DirectedCigar>>>,
}

impl TransitiveStrategyFactory {
    fn new(kmers: Arc<Vec<SeqKmers>>, n_contigs: usize) -> Self {
        Self {
            kmers,
            closest: Arc::new((0..n_contigs).map(|_| ArcSwapOption::const_empty()).collect()),
            cigars: Arc::new(TriangleMatrix::new(n_contigs, OnceLock::new())),
        }
    }
}

impl StrategyFactory for TransitiveStrategyFactory {
    type Strategy = ParallelTransitiveStrategy;

    fn new_strategy(&self) -> Self::Strategy {
        ParallelTransitiveStrategy::new(Arc::clone(&self.kmers), Arc::clone(&self.closest), Arc::clone(&self.cigars))
    }
}

// ========== Single-thread/parallel processing logic ==========

/// If necessary, calculates minimizer-based sequence divergence;
/// if necessary, constructs actual alignment using `aligner`;
/// writes resulting PAF entry to `out` and return CIGAR.
///
/// In the alignment, `entry1` is reference and `entry2` is query.
///
/// Close `construct_alignment` takes two entries and returns CIGAR and its score.
fn process_pair<'a>(
    strategy: &mut impl Strategy,
    aligner: &Aligner,
    entry1: RefEntry<'a>,
    entry2: RefEntry<'a>,
    against_contig: &[bool],
    minimizers: &[Vec<u64>],
    params: &Params,
    out: &mut impl io::Write,
    counts: &mut Counts,
) -> crate::Result<()>
{
    write!(out, "{}\t{len}\t0\t{len}\t+\t", entry2.name, len = entry2.seq.len()).map_err(add_path!(!))?;
    write!(out, "{}\t{len}\t0\t{len}\t", entry1.name, len = entry1.seq.len()).map_err(add_path!(!))?;

    let opt_div = if params.skip_div { None } else {
        Some(minim_div::jaccard_distance(&minimizers[entry1.id.ix()], &minimizers[entry2.id.ix()]))
    };
    let mut opt_cigar = None::<Cigar>;
    let thresh_div = if against_contig[entry1.id.ix()] || against_contig[entry2.id.ix()]
        { params.against_div } else { params.thresh_div };
    if opt_div.map(|(_, dv)| dv <= thresh_div).unwrap_or(true) {
        let (cigar, score, accelerated) = strategy.align(aligner, entry1, entry2, params)?;
        counts.aligned += 1;
        counts.accelerated += u32::from(accelerated);
        let mut nmatches = 0;
        let mut nerrs = 0;
        for item in cigar.iter() {
            match item.operation() {
                Operation::Equal => nmatches += item.len(),
                _ => nerrs += item.len(),
            }
        }
        let aln_len = nmatches + nerrs;
        write!(out, "{}\t{}\t255", nmatches, aln_len).map_err(add_path!(!))?;
        let dv = f64::from(nerrs) / f64::from(aln_len);
        let qv = if dv.is_finite() { -10.0 * dv.log10() } else { f64::INFINITY };
        write!(out, "\tNM:i:{}\tAS:i:{}\tdv:f:{:.9}\tqv:f:{:.6}", nerrs, score, dv, qv).map_err(add_path!(!))?;
        strategy.save_cigar(entry1.id, entry2.id, &cigar, dv, params);
        opt_cigar = Some(cigar);
    } else {
        counts.skipped += 1;
        write!(out, "0\t0\t255").map_err(add_path!(!))?;
    }
    if let Some((uniq_minims, minim_dv)) = opt_div {
        write!(out, "\tum:i:{}\tmd:f:{:.9}", uniq_minims, minim_dv).map_err(add_path!(!))?;
    }
    if let Some(cigar) = opt_cigar {
        write!(out, "\tcg:Z:{}", cigar).map_err(add_path!(!))?;
    }
    writeln!(out).map_err(add_path!(!))?;
    Ok(())
}

fn align_pairs_singlethread(
    mut strategy: impl Strategy,
    contig_set: &ContigSet,
    pairs: &[(ContigId, ContigId)],
    against_contig: &[bool],
    minimizers: &[Vec<u64>],
    params: &Params,
    out: &mut impl io::Write,
    verbose: bool,
) -> crate::Result<Counts>
{
    let aligner = Aligner::new(params.penalties.clone(), params.accuracy, None, false);
    let mut counts = Counts::default();
    let mult = 100.0 / pairs.len() as f64;
    // Power of 2 minus 1.
    const LOG_FREQ: usize = 511;
    for (ix, &(i, j)) in pairs.iter().enumerate() {
        process_pair(&mut strategy, &aligner, RefEntry::from_set(contig_set, i), RefEntry::from_set(contig_set, j),
            against_contig, minimizers, params, out, &mut counts)?;
        if verbose && (ix & LOG_FREQ) == LOG_FREQ {
            log::debug!("    Aligned ≈{:5.1}% pairs", mult * ix as f64);
        }
    }
    Ok(counts)
}

fn align_pairs_parallel(
    factory: impl StrategyFactory,
    contig_set: Arc<ContigSet>,
    pairs: Vec<(ContigId, ContigId)>,
    against_contig: Vec<bool>,
    minimizers: Vec<Vec<u64>>,
    params: &Params,
    threads: usize,
    outputs: Vec<impl io::Write + Send + 'static>,
) -> crate::Result<Counts> {
    let pairs = Arc::new(pairs);
    let minimizers = Arc::new(minimizers);
    let mut handles = Vec::with_capacity(threads);
    let n_pairs = pairs.len();
    let mut start = 0;
    debug_assert_eq!(outputs.len(), threads);
    for (worker_ix, mut out) in outputs.into_iter().enumerate() {
        if start == n_pairs {
            break;
        }
        let end = start + (n_pairs - start).fast_ceil_div(threads - worker_ix);
        // Closure with cloned data.
        {
            let contig_set = Arc::clone(&contig_set);
            let pairs = Arc::clone(&pairs);
            let against_contig = against_contig.clone();
            let minimizers = Arc::clone(&minimizers);
            let params = params.clone();
            let verbose = worker_ix == 0;
            let strategy = factory.new_strategy();
            handles.push(thread::spawn(move || {
                align_pairs_singlethread(strategy, &contig_set, &pairs[start..end],
                    &against_contig, &minimizers, &params, &mut out, verbose)
            }));
        }
        start = end;
    }
    assert_eq!(start, n_pairs);
    let mut counts = Counts::default();
    for handle in handles {
        let curr_counts = handle.join().expect("Worker process failed")?;
        counts.add(&curr_counts);
    }
    Ok(counts)
}

// ========== Main public function ==========

/// Aligns sequences to each other.
pub fn align_sequences(
    contig_set: Arc<ContigSet>,
    pairs: Vec<(ContigId, ContigId)>,
    against_contig: Vec<bool>,
    params: &Params,
    threads: u16,
    mut outputs: Vec<impl io::Write + Send + 'static>,
) -> crate::Result<()>
{
    let threads = usize::from(threads);
    let n_contigs = contig_set.len();
    // Find sequences that actually appear.
    let mut in_use = vec![false; n_contigs];
    let mut n_rem = n_contigs;
    for &(i, j) in &pairs {
        // Will do -1 if old value was false.
        n_rem -= usize::from(std::mem::replace(&mut in_use[i.ix()], true));
        n_rem -= usize::from(std::mem::replace(&mut in_use[j.ix()], true));
        if n_rem == 0 {
            break;
        }
    }
    log::debug!("    Preprocessing k-mers and minimizers");
    let (minimizers, kmers) = fill_kmers(&contig_set, &in_use, params, threads);
    let kmers = Arc::new(kmers);

    log::debug!("    Started alignment");
    // Assume that pairs do not repeat.
    let use_transitivity = params.transitive_div > 0.0
        && pairs.len() * 10 >= TriangleMatrix::calc_linear_len(contig_set.len());
    let counts = if threads == 1 {
        if use_transitivity {
            align_pairs_singlethread(
                TransitiveStrategy::new(kmers, n_contigs),
                &contig_set, &pairs, &against_contig, &minimizers, params, &mut outputs[0], true)?
        } else {
            align_pairs_singlethread(
                BackboneStrategy::new(kmers),
                &contig_set, &pairs, &against_contig, &minimizers, params, &mut outputs[0], true)?
        }
    } else {
        if use_transitivity {
            align_pairs_parallel(
                TransitiveStrategyFactory::new(kmers, n_contigs),
                contig_set, pairs, against_contig, minimizers, params, threads, outputs)?
        } else {
            align_pairs_parallel(
                BackboneStrategyFactory::new(kmers),
                contig_set, pairs, against_contig, minimizers, params, threads, outputs)?
        }
    };
    counts.summarize();
    Ok(())
}
