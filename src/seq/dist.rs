use std::{
    thread,
    path::Path,
    sync::Arc,
    rc::Rc,
};
use smallvec::SmallVec;
use htslib::bam::{
    self,
    Read as BamRead,
    record::Aux,
};
use crate::{
    ext::{self, TriangleMatrix},
    seq::{
        NamedSeq,
        contigs::{ContigNames, ContigId},
        wfa::{Aligner, Penalties},
        cigar::{Cigar, CigarItem, Operation},
    },
    err::{Error, error},
    algo::HashMap,
    math::RoundDiv,
};

const CAPACITY: usize = 4;
type KmerCache<'a> = HashMap<&'a [u8], SmallVec<[u32; CAPACITY]>>;

/// Creates a HashMap containing all the k-mers in the sequence.
/// A good rolling hash function should speed up the code.
fn cache_kmers(seq: &[u8], k: u32) -> KmerCache {
    let mut map = KmerCache::default();
    for (i, kmer) in seq.windows(k as usize).enumerate() {
        map.entry(kmer)
            .or_insert_with(SmallVec::new)
            .push(i as u32);
    }
    map
}

/// Finds all matches between k-mers in two caches.
fn find_kmer_matches(cache1: &KmerCache, cache2: &KmerCache) -> Vec<(u32, u32)> {
    let mut matches = Vec::new();
    for (kmer1, positions1) in cache1.iter() {
        if let Some(positions2) = cache2.get(kmer1) {
            for &pos1 in positions1 {
                for &pos2 in positions2 {
                    matches.push((pos1, pos2));
                }
            }
        }
    }
    matches.sort_unstable();
    matches
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

/// Aligns two sequences between k-mer matches.
/// First sequence: reference, sequence: query.
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
    entry1: &NamedSeq,
    entry2: &NamedSeq,
    kmer_matches: &[Vec<(u32, u32)>],
    backbone_ks: &[u32],
    max_gap: u32,
) -> crate::Result<(Cigar, i32)>
{
    let mut best_cigar = None;
    // Instead of i32::MIN, use ~ i32::MIN / 2, so that we notice errors if the score is too low.
    let mut best_score = -0x3fffffff_i32;
    for (&k, matches) in backbone_ks.iter().zip(kmer_matches) {
        let (cigar, score) = align(aligner, entry1, entry2, matches, k, max_gap)?;
        if score > best_score {
            best_score = score;
            best_cigar = Some(cigar);
        }
    }
    best_cigar.ok_or_else(|| error!(RuntimeError, "No alignment found between {} and {}", entry1.name(), entry2.name()))
        .map(|cigar| (cigar, best_score))
}

/// Aligns sequences to each other.
pub fn align_sequences(
    entries: &[NamedSeq],
    pairs: Vec<(u32, u32)>,
    penalties: &Penalties,
    backbone_ks: &[u32],
    accuracy: u8,
    max_gap: u32,
    threads: u16,
) -> crate::Result<Vec<(Cigar, i32)>>
{
    let n_entries = entries.len();
    let n_pairs = pairs.len();
    let mut all_kmer_matches = vec![Vec::with_capacity(backbone_ks.len()); n_pairs];
    for &backbone_k in backbone_ks {
        let mut caches = vec![None; n_entries];
        for (&(i, j), kmer_matches) in pairs.iter().zip(all_kmer_matches.iter_mut()) {
            let i = i as usize;
            let j = j as usize;
            caches[i].get_or_insert_with(|| cache_kmers(entries[i].seq(), backbone_k));
            caches[j].get_or_insert_with(|| cache_kmers(entries[j].seq(), backbone_k));
            kmer_matches.push(find_kmer_matches(caches[i].as_ref().unwrap(), caches[j].as_ref().unwrap()));
        }
    }

    if threads == 1 {
        let mut alns = Vec::with_capacity(n_pairs);
        let aligner = Aligner::new(penalties.clone(), accuracy);
        for (&(i, j), matches) in pairs.iter().zip(&all_kmer_matches) {
            let i = i as usize;
            let j = j as usize;
            alns.push(align_multik(&aligner, &entries[i], &entries[j], matches, backbone_ks, max_gap)?);
        }
        Ok(alns)
    } else {
        multithread_align(entries, pairs, all_kmer_matches, penalties, backbone_ks, accuracy, max_gap, threads)
    }
}

fn multithread_align(
    entries: &[NamedSeq],
    pairs: Vec<(u32, u32)>,
    // [pair_index][backbone_index] -> indices of matches
    all_kmer_matches: Vec<Vec<Vec<(u32, u32)>>>,
    penalties: &Penalties,
    backbone_ks: &[u32],
    accuracy: u8,
    max_gap: u32,
    threads: u16,
) -> crate::Result<Vec<(Cigar, i32)>>
{
    let pairs = Arc::new(pairs);
    let entries = Arc::new(entries.to_vec());
    let all_kmer_matches = Arc::new(all_kmer_matches);

    let threads = usize::from(threads);
    let mut handles = Vec::with_capacity(threads);
    let n_pairs = pairs.len();
    let mut start = 0;
    for worker_ix in 0..threads {
        if start == n_pairs {
            break;
        }
        let rem_workers = threads - worker_ix;
        let end = start + (n_pairs - start).fast_ceil_div(rem_workers);
        // Closure with cloned data.
        {
            let pairs = Arc::clone(&pairs);
            let entries = Arc::clone(&entries);
            let all_kmer_matches = Arc::clone(&all_kmer_matches);
            let penalties = penalties.clone();
            let curr_backbone_ks = backbone_ks.to_vec();
            handles.push(thread::spawn(move || {
                assert!(start < end);
                let aligner = Aligner::new(penalties, accuracy);
                pairs[start..end].iter().zip(all_kmer_matches[start..end].iter())
                    .map(|(&(i, j), matches)|
                        align_multik(&aligner, &entries[i as usize], &entries[j as usize],
                            matches, &curr_backbone_ks, max_gap))
                    .collect::<crate::Result<Vec<_>>>()
            }));
        }
        start = end;
    }
    assert_eq!(start, n_pairs);

    // Collect all results so that even if there is a problem in one of the sequences, all handles are closed.
    let results: Vec<_> = handles.into_iter().map(|handle| handle.join().expect("Worker process failed")).collect();
    let mut alns = Vec::with_capacity(n_pairs);
    for res in results.into_iter() {
        alns.extend(res?.into_iter());
    }
    Ok(alns)
}

fn create_bam_header(entries: &[NamedSeq], backbone_k: u32, accuracy: u8) -> bam::header::Header {
    let mut header = bam::header::Header::new();
    let mut prg = bam::header::HeaderRecord::new(b"PG");
    prg.push_tag(b"ID", crate::command::PROGRAM);
    prg.push_tag(b"PN", crate::command::PROGRAM);
    prg.push_tag(b"VN", crate::command::VERSION);
    prg.push_tag(b"CL", &std::env::args().collect::<Vec<_>>().join(" "));
    header.push_record(&prg);
    header.push_comment(format!("backbone-k={};accuracy-lvl={}", backbone_k, accuracy).as_bytes());

    for entry in entries {
        let mut record = bam::header::HeaderRecord::new(b"SQ");
        record.push_tag(b"SN", entry.name());
        record.push_tag(b"LN", entry.len());
        header.push_record(&record);
    }
    header
}

fn fill_cigar_buffer(buffer: &mut bam::record::CigarString, iter: impl Iterator<Item = CigarItem>) {
    buffer.0.clear();
    buffer.0.extend(iter.map(CigarItem::to_htslib));
}

fn create_record(
    header: &Rc<bam::HeaderView>,
    query: &NamedSeq,
    refid: usize,
    cigar_view: &bam::record::CigarString,
    edit_dist: u32,
    score: i32,
    divergence: f64,
) -> bam::Record
{
    let mut record = bam::Record::new();
    record.set_header(Rc::clone(header));
    record.set_tid(refid as i32);
    record.set_pos(0);
    record.set_mtid(-1);
    record.set_mpos(-1);
    record.set(query.name().as_bytes(), Some(cigar_view), &[], &[]);
    record.push_aux(b"NM", bam::record::Aux::U32(edit_dist)).expect("Cannot set NM tag");
    record.push_aux(b"AS", bam::record::Aux::I32(score)).expect("Cannot set AS tag");
    record.push_aux(b"dv", bam::record::Aux::Double(divergence)).expect("Cannot set `dv` tag");
    record
}

/// Creates BAM file, writes all pairwise records, and returns linear vector of divergences.
pub fn write_all(
    bam_path: &Path,
    entries: &[NamedSeq],
    alns: TriangleMatrix<(Cigar, i32)>,
    backbone_k: u32,
    accuracy: u8,
) -> crate::Result<TriangleMatrix<f64>> {
    let header = create_bam_header(entries, backbone_k, accuracy);
    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam)?;
    writer.set_compression_level(bam::CompressionLevel::Maximum)?;
    let mut cigar_buffer = bam::record::CigarString(Vec::new());

    let n = entries.len();
    let mut pairwise_records = vec![None; n * n];
    let mut divergences = TriangleMatrix::new(n, 0.0);
    let header_view = Rc::new(bam::HeaderView::from_header(&header));
    for (i, refer) in entries.iter().enumerate() {
        let rlen = refer.len();
        fill_cigar_buffer(&mut cigar_buffer, std::iter::once(CigarItem::new(Operation::Equal, rlen)));
        pairwise_records[i * n + i] = Some(create_record(&header_view, refer, i, &cigar_buffer, 0, 0, 0.0));

        for (j, query) in (i + 1..).zip(&entries[i + 1..]) {
            let (cigar, score) = &alns[(i, j)];
            if query.len() != cigar.query_len() || rlen != cigar.ref_len() {
                return Err(error!(RuntimeError,
                    "Generated invalid alignment between {} ({} bp) and {} ({} bp), CIGAR qlen {}, rlen {}",
                    query.name(), query.len(), refer.name(), rlen, cigar.query_len(), cigar.ref_len()));
            }

            let (nmatches, total_size) = cigar.frac_matches();
            let edit_dist = total_size - nmatches;
            let divergence = f64::from(edit_dist) / f64::from(total_size);
            divergences[(i, j)] = divergence;
            fill_cigar_buffer(&mut cigar_buffer, cigar.iter().copied());
            pairwise_records[i * n + j] = Some(
                create_record(&header_view, query, i, &cigar_buffer, edit_dist, *score, divergence));

            fill_cigar_buffer(&mut cigar_buffer, cigar.iter().map(CigarItem::invert));
            pairwise_records[j * n + i] = Some(
                create_record(&header_view, refer, j, &cigar_buffer, edit_dist, *score, divergence));
        }
    }

    for record in pairwise_records.into_iter() {
        writer.write(&record.expect("Alignment record is not set"))?;
    }
    Ok(divergences)
}

/// Loads edit distances between all contigs based on a BAM file.
/// All contigs must be present in the BAM file.
pub fn load_edit_distances(path: impl AsRef<Path>, contigs: &ContigNames) -> crate::Result<TriangleMatrix<u32>> {
    let path = path.as_ref();
    let mut reader = bam::Reader::from_path(&path)?;
    let mut record = bam::Record::new();
    let mut matrix = TriangleMatrix::new(contigs.len(), u32::MAX);
    while reader.read(&mut record).transpose()?.is_some() {
        let qname = std::str::from_utf8(record.qname())
            .map_err(|_| Error::Utf8("contig name", record.qname().to_vec()))?;
        let i = match contigs.try_get_id(&qname) {
            Some(val) => val.ix(),
            None => continue,
        };
        let rname = reader.header().tid2name(record.tid() as u32);
        let rname = std::str::from_utf8(rname).map_err(|_| Error::Utf8("contig name", rname.to_vec()))?;
        let j = match contigs.try_get_id(&rname) {
            Some(val) => val.ix(),
            None => continue,
        };
        if i >= j {
            continue;
        }
        let edit_dist: u32 = match record.aux(b"NM")
            .map_err(|_| error!(InvalidData, "BAM file {} does not contain NM field", ext::fmt::path(&path)))?
        {
            Aux::I8(val) => val.try_into().unwrap(),
            Aux::U8(val) => val.into(),
            Aux::I16(val) => val.try_into().unwrap(),
            Aux::U16(val) => val.into(),
            Aux::I32(val) => val.try_into().unwrap(),
            Aux::U32(val) => val,
            _ => return Err(error!(InvalidData, "Invalid value for NM field in {}", ext::fmt::path(&path))),
        };
        matrix[(i, j)] = edit_dist;
    }
    if let Some(k) = matrix.iter().position(|&val| val == u32::MAX) {
        let (i, j) = matrix.from_linear_index(k);
        Err(error!(InvalidData, "BAM file {} does not contain alignment between contigs {} and {}",
            ext::fmt::path(&path), contigs.get_name(ContigId::new(i)), contigs.get_name(ContigId::new(j))))
    } else {
        Ok(matrix)
    }
}
