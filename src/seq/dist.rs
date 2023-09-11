use std::{
    thread,
    path::Path,
    sync::Arc,
    rc::Rc,
};
use fx::FxHashMap;
use smallvec::SmallVec;
use htslib::bam;
use crate::{
    seq::{
        NamedSeq,
        wfa::{Aligner, Penalties},
        cigar::{Cigar, CigarItem, Operation},
    },
    Error,
};

const CAPACITY: usize = 4;
type KmerCache<'a> = FxHashMap<&'a [u8], SmallVec<[u32; CAPACITY]>>;

/// Creates a HashMap containing all the k-mers in the sequence.
/// A good rolling hash function should speed up the code.
fn cache_kmers(seq: &[u8], k: usize) -> KmerCache {
    let mut map = KmerCache::default();
    for (i, kmer) in seq.windows(k).enumerate() {
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

/// Aligns two sequences between k-mer matches.
/// First sequence: reference, sequence: query.
fn align(
    aligner: &Aligner,
    seq1: &[u8],
    seq2: &[u8],
    kmer_matches: &[(u32, u32)],
    k: usize,
) -> Result<(Cigar, i32), Error>
{
    let penalties = aligner.penalties();
    let sparse_aln = bio::alignment::sparse::lcskpp(&kmer_matches, k);
    let mut cigar = Cigar::new();
    let mut score = 0;
    let mut i1 = 0;
    let mut j1 = 0;
    let k = k as u32;
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
        let jump1 = i2 - i1;
        let jump2 = j2.checked_sub(j1).unwrap();
        if jump1 > 0 && jump2 > 0 {
            let subseq1 = &seq1[i1 as usize..i2 as usize];
            let subseq2 = &seq2[j1 as usize..j2 as usize];
            score += aligner.align(subseq1, subseq2, &mut cigar)?;
        } else if jump1 > 0 {
            score -= penalties.gap_open + jump1 as i32 * penalties.gap_extend;
            cigar.push_unchecked(Operation::Del, jump1);
        } else if jump2 > 0 {
            score -= penalties.gap_open + jump2 as i32 * penalties.gap_extend;
            cigar.push_unchecked(Operation::Ins, jump2);
        }

        curr_match += k;
        i1 = i2 + k;
        j1 = j2 + k;
    }

    if curr_match > 0 {
        cigar.push_unchecked(Operation::Equal, curr_match);
    }
    let n = seq1.len() as u32;
    let m = seq2.len() as u32;
    let jump1 = n.checked_sub(i1).unwrap();
    let jump2 = m.checked_sub(j1).unwrap();
    if jump1 > 0 && jump2 > 0 {
        score += aligner.align(&seq1[i1 as usize..n as usize], &seq2[j1 as usize..m as usize], &mut cigar)?;
    } else if jump1 > 0 {
        score -= penalties.gap_open + jump1 as i32 * penalties.gap_extend;
        cigar.push_unchecked(Operation::Del, jump1);
    } else if jump2 > 0 {
        score -= penalties.gap_open + jump2 as i32 * penalties.gap_extend;
        cigar.push_unchecked(Operation::Ins, jump2);
    }
    Ok((cigar, score))
}

/// Calculates all pairwise divergences between all sequences, writes alignments to PAF file,
/// and returns condensed distance matrix (see `kodama` crate).
pub fn pairwise_divergences(
    bam_path: &Path,
    entries: &[NamedSeq],
    penalties: &Penalties,
    k: usize,
    threads: u16,
) -> Result<Vec<f64>, Error> {
    let caches: Vec<_> = entries.iter().map(|entry| cache_kmers(entry.seq(), k)).collect();
    let alns = if threads == 1 {
        log::debug!("        Aligning sequences in 1 thread");
        let n = entries.len();
        let mut alns = Vec::with_capacity(n * (n - 1) / 2);
        let aligner = Aligner::new(penalties.clone());
        for (i, (entry1, cache1)) in entries.iter().zip(&caches).enumerate() {
            for (entry2, cache2) in entries[i + 1..].iter().zip(&caches[i + 1..]) {
                let kmer_matches = find_kmer_matches(cache1, cache2);
                alns.push(align(&aligner, entry1.seq(), entry2.seq(), &kmer_matches, k)?);
            }
        }
        alns
    } else {
        divergences_multithread(entries, &caches, penalties, k, threads)?
    };
    write_all(bam_path, entries, alns, k, penalties.accuracy)
}

fn divergences_multithread(
    entries: &[NamedSeq],
    caches: &[KmerCache],
    penalties: &Penalties,
    k: usize,
    threads: u16,
) -> Result<Vec<(Cigar, i32)>, Error> {
    let threads = usize::from(threads);
    log::debug!("        Aligning sequences in {} threads", threads);
    let n = entries.len() as u32;
    let pairs: Arc<Vec<(u32, u32)>> = Arc::new(
        (0..n - 1).flat_map(|i| (i + 1..n).map(move |j| (i, j))).collect());
    let n_pairs = pairs.len();
    let mut handles = Vec::with_capacity(threads);

    let entries = Arc::new(entries.to_vec());
    // Need to precompute all kmer matches as `caches` cannot be passed between threads.
    let all_kmer_matches = Arc::new(pairs.iter()
        .map(|&(i, j)| find_kmer_matches(&caches[i as usize], &caches[j as usize])).collect::<Vec<_>>());
    let mut start = 0;
    for i in 0..threads {
        if start == n_pairs {
            break;
        }
        let rem_workers = threads - i;
        // Ceiling division.
        let end = start + ((n_pairs - start) + rem_workers - 1) / rem_workers;
        // Closure with cloned data.
        {
            let pairs = Arc::clone(&pairs);
            let entries = Arc::clone(&entries);
            let all_kmer_matches = Arc::clone(&all_kmer_matches);
            let penalties = penalties.clone();
            handles.push(thread::spawn(move || {
                assert!(start < end);
                let aligner = Aligner::new(penalties);
                pairs[start..end].iter().zip(all_kmer_matches[start..end].iter()).map(|(&(i, j), kmer_matches)|
                        align(&aligner, &entries[i as usize].seq(), entries[j as usize].seq(), kmer_matches, k))
                    .collect::<Result<Vec<_>, Error>>()
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

fn create_bam_header(entries: &[NamedSeq], k: usize, accuracy: u8) -> bam::header::Header {
    let mut header = bam::header::Header::new();
    let mut prg = bam::header::HeaderRecord::new(b"PG");
    prg.push_tag(b"ID", crate::command::PROGRAM);
    prg.push_tag(b"PN", crate::command::PROGRAM);
    prg.push_tag(b"VN", crate::command::VERSION);
    prg.push_tag(b"CL", &std::env::args().collect::<Vec<_>>().join(" "));
    header.push_record(&prg);
    header.push_comment(format!("backbone-k={};accuracy-lvl={}", k, accuracy).as_bytes());

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
fn write_all(
    bam_path: &Path,
    entries: &[NamedSeq],
    alns: Vec<(Cigar, i32)>,
    k: usize,
    accuracy: u8,
) -> Result<Vec<f64>, Error> {
    let header = create_bam_header(entries, k, accuracy);
    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam)?;
    writer.set_compression_level(bam::CompressionLevel::Maximum)?;
    let mut cigar_buffer = bam::record::CigarString(Vec::new());

    let n = entries.len();
    let mut pairwise_records = vec![None; n * n];
    let mut divergences = Vec::new();
    let header_view = Rc::new(bam::HeaderView::from_header(&header));
    let mut alns_iter = alns.into_iter();
    for (i, refer) in entries.iter().enumerate() {
        let rlen = refer.len();
        fill_cigar_buffer(&mut cigar_buffer, std::iter::once(CigarItem::new(Operation::Equal, rlen)));
        pairwise_records[i * n + i] = Some(create_record(&header_view, refer, i, &cigar_buffer, 0, 0, 0.0));

        for (j, query) in (i + 1..).zip(&entries[i + 1..]) {
            let (cigar, score) = alns_iter.next().expect("Too few pairwse alignments");
            if query.len() != cigar.query_len() || rlen != cigar.ref_len() {
                return Err(Error::RuntimeError(format!(
                    "Generated invalid alignment between {} ({} bp) and {} ({} bp), CIGAR qlen {}, rlen {}",
                    query.name(), query.len(), refer.name(), rlen, cigar.query_len(), cigar.ref_len())));
            }

            let (nmatches, total_size) = cigar.frac_matches();
            let edit_dist = total_size - nmatches;
            let divergence = f64::from(edit_dist) / f64::from(total_size);
            divergences.push(divergence);
            fill_cigar_buffer(&mut cigar_buffer, cigar.iter().copied());
            pairwise_records[i * n + j] = Some(
                create_record(&header_view, query, i, &cigar_buffer, edit_dist, score, divergence));

            fill_cigar_buffer(&mut cigar_buffer, cigar.iter().map(CigarItem::invert));
            pairwise_records[j * n + i] = Some(
                create_record(&header_view, refer, j, &cigar_buffer, edit_dist, score, divergence));
        }
    }
    assert!(alns_iter.next().is_none(), "Too many pairwise alignments");

    for record in pairwise_records.into_iter() {
        writer.write(&record.expect("Alignment record is not set"))?;
    }
    Ok(divergences)
}
