use std::{
    thread, io,
    sync::Arc,
    cmp::Ordering,
    fmt::Write as FmtWrite,
};
use ruint::aliases::U256;
use smallvec::SmallVec;
use crate::{
    seq::{
        wfa, div,
        contigs::{ContigId, ContigSet},
        kmers::{self, Kmer},
        wfa::{Aligner, Penalties},
        cigar::{Cigar, Operation},
    },
    err::{error, validate_param, add_path},
    math::RoundDiv,
    ext::TriangleMatrix,
};

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
    let (minimizers, kmers) = fill_kmers_singlethread(&contig_set, &in_use, params);
    let kmers = Arc::new(kmers);

    // Assume that pairs do not repeat.
    let use_transitivity = params.transitive_div > 0.0
        && pairs.len() * 10 >= TriangleMatrix::calc_linear_len(contig_set.len());
    if threads == 1 {
        if use_transitivity {
            transitive_align_pairs_singlethread(
                &contig_set, &pairs, &against_contig, &minimizers, kmers, params, &mut outputs[0])
        } else {
            align_pairs_singlethread(
                &contig_set, &pairs, &against_contig, &minimizers, kmers, params, &mut outputs[0], true)
        }
    } else {
        align_pairs_parallel(contig_set, pairs, against_contig, minimizers, kmers, params, threads, outputs)
    }
}

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
    params: &Params,
) -> (Vec<Vec<u64>>, Vec<SeqKmers>)
{
    let mut minimizers = Vec::with_capacity(if params.skip_div { 0 } else { contig_set.len() });
    // Backbone ks will be empty if no alignments need to be calculated.
    let mut kmers = Vec::with_capacity(contig_set.len() * params.backbone_ks.len());

    let mut kmer_buf = Vec::new();
    for (seq, &in_use) in contig_set.seqs().iter().zip(in_use) {
        if !in_use {
            minimizers.push(Vec::new());
            for _ in 0..params.backbone_ks.len() {
                kmers.push(Vec::new());
            }
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

#[derive(Clone, Copy)]
struct RefEntry<'a> {
    id: ContigId,
    name: &'a str,
    seq: &'a [u8],
}

impl<'a> RefEntry<'a> {
    // #[inline]
    // fn new(id: ContigId, name: &'a str, seq: &'a [u8]) -> Self {
    //     Self { id, name, seq}
    // }

    #[inline]
    fn from_set(contig_set: &'a ContigSet, id: ContigId) -> Self {
        Self {
            id,
            name: contig_set.contigs().get_name(id),
            seq: contig_set.get_seq(id),
        }
    }
}

fn align(
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
        let (cigar, score) = align(aligner, entry1, entry2, buf, u32::from(k), params.max_gap)?;
        if score > best_score {
            best_score = score;
            best_cigar = Some(cigar);
        }
    }
    best_cigar.ok_or_else(|| error!(RuntimeError, "No alignment found between {} and {}", entry1.name, entry2.name))
        .map(|cigar| (cigar, best_score))
}

trait AlignerWrapper {
    fn align<'a>(
        &mut self,
        entry1: RefEntry<'a>,
        entry2: RefEntry<'a>,
        params: &Params,
    ) -> crate::Result<(Cigar, i32)>;
}

struct DirectAligner {
    aligner: Aligner,
    kmers: Arc<Vec<SeqKmers>>,
    buf: Vec<(u32, u32)>,
}

impl DirectAligner {
    fn new(aligner: Aligner, kmers: Arc<Vec<SeqKmers>>) -> Self {
        Self {
            aligner, kmers,
            buf: Vec::new(),
        }
    }
}

impl AlignerWrapper for DirectAligner {
    #[inline(always)]
    fn align<'a>(
        &mut self,
        entry1: RefEntry<'a>,
        entry2: RefEntry<'a>,
        params: &Params,
    ) -> crate::Result<(Cigar, i32)> {
        align_multik(&self.aligner, entry1, entry2, &self.kmers, params, &mut self.buf)
    }
}

/// CIGAR for contigs a < b (indices not stored).
/// `ref_smaller` is true when `a` is reference and `b` is query.
#[derive(Clone)]
struct DirectedCigar {
    cigar: Cigar,
    ref_smaller: bool,
}

impl DirectedCigar {
    #[inline(always)]
    fn new(cigar: Cigar, query_id: ContigId, ref_id: ContigId) -> Self {
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

struct TransitiveAligner {
    direct_aligner: DirectAligner,
    /// For each contig, store its closest other contig, corresponding CIGAR and diversity.
    /// Value is only stored if diversity is not greater than params.transitive_div
    /// Closest other contig will always have smaller index than i.
    closest: Vec<Option<(ContigId, DirectedCigar, f64)>>,
    cigars: TriangleMatrix<Option<DirectedCigar>>,
    /// Number of transitive and the total number of constructed alignments.
    n_transitive: u32,
    total: u32,
}

impl TransitiveAligner {
    fn new(direct_aligner: DirectAligner, n_contigs: usize) -> Self {
        Self {
            direct_aligner,
            closest: vec![None; n_contigs],
            cigars: TriangleMatrix::new(n_contigs, None),
            n_transitive: 0,
            total: 0,
        }
    }
}

impl AlignerWrapper for TransitiveAligner {
    fn align<'a>(
        &mut self,
        entry_k: RefEntry<'a>,
        entry_i: RefEntry<'a>,
        params: &Params,
    ) -> crate::Result<(Cigar, i32)> {
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

        self.total += 1;
        if let Some((cigar_ij, j, cigar_jk)) = transitive_edge {
            let cigar_ik = Cigar::find_transitive_alignment(
                &cigar_ij.cigar, cigar_ij.second_is_ref(i, j), &cigar_jk.cigar, cigar_jk.second_is_ref(j, k),
                entry_i.seq, entry_k.seq, &self.direct_aligner.aligner, params.max_gap, params.transitive_anchor);
            let score = self.direct_aligner.aligner.penalties().calculate_score(&cigar_ik);
            self.n_transitive += 1;
            Ok((cigar_ik, score))
        } else {
            self.direct_aligner.align(entry_k, entry_i, params)
        }
    }
}

/// If necessary, calculates minimizer-based sequence divergence;
/// if necessary, constructs actual alignment using `aligner`;
/// writes resulting PAF entry to `out` and return CIGAR.
///
/// In the alignment, `entry1` is reference and `entry2` is query.
///
/// Close `construct_alignment` takes two entries and returns CIGAR and its score.
///
/// If CIGAR was constructed, returns pair (CIGAR, sequence divergence).
fn process_pair<'a>(
    aligner: &mut impl AlignerWrapper,
    entry1: RefEntry<'a>,
    entry2: RefEntry<'a>,
    against_contig: &[bool],
    minimizers: &[Vec<u64>],
    params: &Params,
    buf_cigar: &mut String,
    out: &mut impl io::Write,
) -> crate::Result<Option<(Cigar, f64)>>
{
    write!(out, "{}\t{len}\t0\t{len}\t+\t", entry2.name, len = entry2.seq.len()).map_err(add_path!(!))?;
    write!(out, "{}\t{len}\t0\t{len}\t", entry1.name, len = entry1.seq.len()).map_err(add_path!(!))?;

    let opt_div = if params.skip_div { None } else {
        Some(div::jaccard_distance(&minimizers[entry1.id.ix()], &minimizers[entry2.id.ix()]))
    };
    buf_cigar.clear();
    let mut ret_val = None;
    let thresh_div = if against_contig[entry1.id.ix()] || against_contig[entry2.id.ix()]
        { params.against_div } else { params.thresh_div };
    if opt_div.map(|(_, dv)| dv <= thresh_div).unwrap_or(true) {
        let (cigar, score) = aligner.align(entry1, entry2, params)?;
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
        write!(buf_cigar, "{}", cigar).unwrap();
        ret_val = Some((cigar, dv));
    } else {
        // Skip alignment.
        write!(out, "0\t0\t255").map_err(add_path!(!))?;
    }
    if let Some((uniq_minims, minim_dv)) = opt_div {
        write!(out, "\tum:i:{}\tmd:f:{:.9}", uniq_minims, minim_dv).map_err(add_path!(!))?;
    }
    if !buf_cigar.is_empty() {
        write!(out, "\tcg:Z:{}", buf_cigar).map_err(add_path!(!))?;
    }
    writeln!(out).map_err(add_path!(!))?;
    Ok(ret_val)
}

fn align_pairs_singlethread(
    contig_set: &ContigSet,
    pairs: &[(ContigId, ContigId)],
    against_contig: &[bool],
    minimizers: &[Vec<u64>],
    kmers: Arc<Vec<SeqKmers>>,
    params: &Params,
    out: &mut impl io::Write,
    verbose: bool,
) -> crate::Result<()>
{
    let mut buf = Default::default();
    let aligner = Aligner::new(params.penalties.clone(), params.accuracy, None, false);
    let mut aligner_wrapper = DirectAligner::new(aligner, kmers);
    let mult = 100.0 / pairs.len() as f64;
    // Power of 2 minus 1.
    const LOG_FREQ: usize = 255;
    for (ix, &(i, j)) in pairs.iter().enumerate() {
        process_pair(&mut aligner_wrapper, RefEntry::from_set(contig_set, i), RefEntry::from_set(contig_set, j),
            against_contig, minimizers, params, &mut buf, out)?;
        if verbose && (ix & LOG_FREQ) == LOG_FREQ {
            log::debug!("    Aligned ≈{:5.1}% pairs", mult * ix as f64);
        }
    }
    Ok(())
}

fn align_pairs_parallel(
    contig_set: Arc<ContigSet>,
    pairs: Vec<(ContigId, ContigId)>,
    against_contig: Vec<bool>,
    minimizers: Vec<Vec<u64>>,
    kmers: Arc<Vec<SeqKmers>>,
    params: &Params,
    threads: u16,
    outputs: Vec<impl io::Write + Send + 'static>,
) -> crate::Result<()> {
    let pairs = Arc::new(pairs);
    let minimizers = Arc::new(minimizers);

    let threads = usize::from(threads);
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
            let kmers = Arc::clone(&kmers);
            let params = params.clone();
            let verbose = worker_ix == 0;
            handles.push(thread::spawn(move || {
                align_pairs_singlethread(&contig_set, &pairs[start..end],
                    &against_contig, &minimizers, kmers, &params, &mut out, verbose)
            }));
        }
        start = end;
    }
    assert_eq!(start, n_pairs);
    handles.into_iter().map(|handle| handle.join().expect("Worker process failed")).collect()
}

/// Align all sequence pairs to each other, speeding up things using transitive
fn transitive_align_pairs_singlethread(
    contig_set: &ContigSet,
    pairs: &[(ContigId, ContigId)],
    against_contig: &[bool],
    minimizers: &[Vec<u64>],
    kmers: Arc<Vec<SeqKmers>>,
    params: &Params,
    out: &mut impl io::Write,
) -> crate::Result<()> {
    let aligner = Aligner::new(params.penalties.clone(), params.accuracy, None, false);
    let mut transitive_aligner = TransitiveAligner::new(DirectAligner::new(aligner, kmers), contig_set.len());
    let mut buf = Default::default();

    const LOG_FREQ: usize = 255;
    let mult = 100.0 / pairs.len() as f64;
    for (ix, &(i, j)) in pairs.iter().enumerate() {
        let res = process_pair(&mut transitive_aligner,
            RefEntry::from_set(contig_set, i), RefEntry::from_set(contig_set, j),
            against_contig, minimizers, params, &mut buf, out)?;
        if let Some((cigar, div)) = res {
            let dir_cigar = DirectedCigar::new(cigar, j, i);
            if div <= params.transitive_div {
                let prev_closest = &mut transitive_aligner.closest[j.ix()];
                if let Some((_, _, prev_div)) = prev_closest && *prev_div < div {} else {
                    *prev_closest = Some((i, dir_cigar.clone(), div));
                }
            }
            transitive_aligner.cigars[(i.ix(), j.ix())] = Some(dir_cigar);
        }
        if (ix & LOG_FREQ) == LOG_FREQ {
            log::debug!("    Aligned ≈{:5.1}% pairs ({:.1}% transitive)", mult * ix as f64,
                100.0 * f64::from(transitive_aligner.n_transitive) / f64::from(transitive_aligner.total));
        }
    }
    log::debug!("    Sped up alignment for {}/{} pairs ({:.1}%)",
        transitive_aligner.n_transitive, transitive_aligner.total,
        100.0 * f64::from(transitive_aligner.n_transitive) / f64::from(transitive_aligner.total));
    Ok(())
}
