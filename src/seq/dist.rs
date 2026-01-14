use std::{
    thread, io,
    sync::Arc,
    cmp::Ordering,
    fmt::Write as FmtWrite,
};
use ruint::aliases::U256;
use smallvec::SmallVec;
use rand::{Rng, seq::SliceRandom};
use crate::{
    seq::{
        wfa, div,
        NamedSeq,
        kmers::{self, Kmer},
        wfa::{Aligner, Penalties},
        cigar::{Cigar, Operation},
    },
    err::{error, validate_param, add_path},
    math::RoundDiv,
};

/// Alignment/divergence calculation parameters.
#[derive(Clone)]
pub struct Params {
    pub skip_div: bool,
    pub div_k: u8,
    pub div_w: u8,
    pub thresh_div: f64,
    pub penalties: Penalties,
    pub backbone_ks: Vec<u8>,
    pub accuracy: u8,
    pub max_gap: u32,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            skip_div: false,
            div_k: 15,
            div_w: 15,
            thresh_div: 0.5,
            penalties: Default::default(),
            backbone_ks: vec![25, 51, 101],
            accuracy: 9,
            max_gap: 500,
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
    entries: Vec<NamedSeq>,
    mut pairs: Vec<(u32, u32)>,
    params: &Params,
    threads: u16,
    mut outputs: Vec<impl io::Write + Send + 'static>,
    rng: &mut impl Rng,
) -> crate::Result<()>
{
    let n_entries = entries.len();
    // Find sequences that actually appear.
    let mut entry_in_use = vec![false; n_entries];
    let mut n_rem = n_entries;
    for &(i, j) in &pairs {
        // Will do -1 if old value was false.
        n_rem -= usize::from(std::mem::replace(&mut entry_in_use[i as usize], true));
        n_rem -= usize::from(std::mem::replace(&mut entry_in_use[j as usize], true));
        if n_rem == 0 {
            break;
        }
    }
    let (minimizers, kmers) = fill_kmers_singlethread(&entries, &entry_in_use, params);

    if threads == 1 {
        align_all_singlethread(&entries, &pairs, &minimizers, &kmers, params, &mut outputs[0])
    } else {
        if params.thresh_div < 1.0 {
            // Shuffle pairs for better job distribution since some alignments may be skipped.
            pairs.shuffle(rng);
        }
        align_all_parallel(entries, pairs, minimizers, kmers, params, usize::from(threads), outputs)
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
    entries: &[NamedSeq],
    entry_in_use: &[bool],
    params: &Params,
) -> (Vec<Vec<u64>>, Vec<SeqKmers>)
{
    let mut minimizers = Vec::with_capacity(if params.skip_div { 0 } else { entries.len() });
    // Backbone ks will be empty if no alignments need to be calculated.
    let mut kmers = Vec::with_capacity(entries.len() * params.backbone_ks.len());

    let mut kmer_buf = Vec::new();
    for (entry, &in_use) in entries.iter().zip(entry_in_use) {
        if !in_use {
            minimizers.push(Vec::new());
            for _ in 0..params.backbone_ks.len() {
                kmers.push(Vec::new());
            }
            continue;
        }

        let seq = entry.seq();
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

#[inline(always)]
fn align_gap(
    seq1: &[u8],
    seq2: &[u8],
    i1: u32,
    i2: u32,
    j1: u32,
    j2: u32,
    aligner: &Aligner,
    max_gap: u32,
    cigar: &mut Cigar,
) -> crate::Result<i32>
{
    let penalties = aligner.penalties();
    debug_assert!(i1 <= i2 && j1 <= j2);
    let jump1 = i2 - i1;
    let jump2 = j2 - j1;
    if jump1 > 0 && jump2 > 0 {
        let subseq1 = &seq1[i1 as usize..i2 as usize];
        let subseq2 = &seq2[j1 as usize..j2 as usize];
        if jump1 > max_gap || jump2 > max_gap {
            Ok(penalties.align_simple(subseq1, subseq2, cigar))
        } else {
            aligner.align(subseq1, subseq2, cigar)
        }
    } else if jump1 > 0 {
        cigar.push_unchecked(Operation::Del, jump1);
        Ok(-penalties.gap_open - jump1 as i32 * penalties.gap_extend)
    } else if jump2 > 0 {
        cigar.push_unchecked(Operation::Ins, jump2);
        Ok(-penalties.gap_open - jump2 as i32 * penalties.gap_extend)
    } else {
        Ok(0)
    }
}

fn align(
    aligner: &Aligner,
    entry1: &NamedSeq,
    entry2: &NamedSeq,
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

    let seq1 = entry1.seq();
    let seq2 = entry2.seq();
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
        score += align_gap(seq1, seq2, i1, i2, j1, j2, aligner, max_gap, &mut cigar)?;
        curr_match += backbone_k;
        i1 = i2 + backbone_k;
        j1 = j2 + backbone_k;
    }

    if curr_match > 0 {
        cigar.push_unchecked(Operation::Equal, curr_match);
    }
    let n1 = seq1.len() as u32;
    let n2 = seq2.len() as u32;
    score += align_gap(seq1, seq2, i1, n1, j1, n2, aligner, max_gap, &mut cigar)?;
    assert_eq!(cigar.ref_len(), seq1.len() as u32,
        "Alignment {} - {} produced incorrect CIGAR {}", entry1.name(), entry2.name(), cigar);
    assert_eq!(cigar.query_len(), seq2.len() as u32,
        "Alignment {} - {} produced incorrect CIGAR {}", entry1.name(), entry2.name(), cigar);
    Ok((cigar, score))
}

fn align_multik(
    aligner: &Aligner,
    i: usize,
    entry1: &NamedSeq,
    j: usize,
    entry2: &NamedSeq,
    kmers: &[SeqKmers],
    params: &Params,
    buf: &mut Vec<(u32, u32)>,
) -> crate::Result<(Cigar, i32)>
{
    let mut best_cigar = None;
    let mut best_score = i32::MIN;
    let n_backbones = params.backbone_ks.len();
    for (&k, kmers1, kmers2) in itertools::izip!(
            &params.backbone_ks, &kmers[n_backbones * i..], &kmers[n_backbones * j..])
    {
        get_kmer_matches(kmers1, kmers2, buf);
        let (cigar, score) = align(aligner, entry1, entry2, buf, u32::from(k), params.max_gap)?;
        if score > best_score {
            best_score = score;
            best_cigar = Some(cigar);
        }
    }
    best_cigar.ok_or_else(|| error!(RuntimeError, "No alignment found between {} and {}", entry1.name(), entry2.name()))
        .map(|cigar| (cigar, best_score))
}

fn process_pair(
    entries: &[NamedSeq],
    i: usize,
    j: usize,
    minimizers: &[Vec<u64>],
    kmers: &[SeqKmers],
    params: &Params,
    aligner: &Aligner,
    buf1: &mut Vec<(u32, u32)>,
    buf2: &mut String,
    out: &mut impl io::Write,
) -> crate::Result<()>
{
    let entry1 = &entries[i as usize];
    let entry2 = &entries[j as usize];
    write!(out, "{}\t{len}\t0\t{len}\t+\t", entry2.name(), len = entry2.seq().len()).map_err(add_path!(!))?;
    write!(out, "{}\t{len}\t0\t{len}\t", entry1.name(), len = entry1.seq().len()).map_err(add_path!(!))?;

    let opt_div = if params.skip_div { None } else { Some(div::jaccard_distance(&minimizers[i], &minimizers[j])) };
    buf2.clear();
    if opt_div.map(|(_, dv)| dv <= params.thresh_div).unwrap_or(false) {
        let (cigar, score) = align_multik(aligner, i, entry1, j, entry2, kmers, params, buf1)?;
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
        write!(buf2, "{}", cigar).unwrap();
    } else {
        // Skip alignment.
        write!(out, "0\t0\t255").map_err(add_path!(!))?;
    }
    if let Some((uniq_minims, minim_dv)) = opt_div {
        write!(out, "\tum:i:{}\tmd:f:{:.9}", uniq_minims, minim_dv).map_err(add_path!(!))?;
    }
    if !buf2.is_empty() {
        write!(out, "\tcg:Z:{}", buf2).map_err(add_path!(!))?;
    }
    writeln!(out).map_err(add_path!(!))?;
    Ok(())
}

fn align_all_singlethread(
    entries: &[NamedSeq],
    pairs: &[(u32, u32)],
    minimizers: &[Vec<u64>],
    kmers: &[SeqKmers],
    params: &Params,
    out: &mut impl io::Write,
) -> crate::Result<()>
{
    let mut buf1 = Default::default();
    let mut buf2 = Default::default();
    let aligner = Aligner::new(params.penalties.clone(), params.accuracy);
    for &(i, j) in pairs {
        process_pair(entries, i as usize, j as usize, minimizers, kmers, params, &aligner,
            &mut buf1, &mut buf2, out)?;
    }
    Ok(())
}

fn align_all_parallel(
    entries: Vec<NamedSeq>,
    pairs: Vec<(u32, u32)>,
    minimizers: Vec<Vec<u64>>,
    kmers: Vec<SeqKmers>,
    params: &Params,
    threads: usize,
    outputs: Vec<impl io::Write + Send + 'static>,
) -> crate::Result<()>
{
    let entries = Arc::new(entries);
    let pairs = Arc::new(pairs);
    let minimizers = Arc::new(minimizers);
    let kmers = Arc::new(kmers);

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
            let pairs = Arc::clone(&pairs);
            let entries = Arc::clone(&entries);
            let minimizers = Arc::clone(&minimizers);
            let kmers = Arc::clone(&kmers);
            let params = params.clone();
            handles.push(thread::spawn(move || {
                let aligner = Aligner::new(params.penalties.clone(), params.accuracy);
                let mut buf1 = Default::default();
                let mut buf2 = Default::default();
                pairs[start..end].iter()
                    .map(|&(i, j)| process_pair(&entries, i as usize, j as usize, &minimizers, &kmers, &params,
                        &aligner, &mut buf1, &mut buf2, &mut out))
                    .collect::<crate::Result<()>>()
            }));
        }
        start = end;
    }
    assert_eq!(start, n_pairs);
    handles.into_iter().map(|handle| handle.join().expect("Worker process failed")).collect()
}
