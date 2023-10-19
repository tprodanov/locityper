use std::{
    thread,
    cmp::{min, Ordering},
    path::Path,
    sync::Arc,
    rc::Rc,
    ops::{Index, IndexMut},
};
use smallvec::SmallVec;
use htslib::bam::{
    self,
    Read as BamRead,
    record::Aux,
};
use crate::{
    ext,
    seq::{
        NamedSeq,
        contigs::{ContigNames, ContigId, Genotype},
        wfa::{Aligner, Penalties},
        cigar::{Cigar, CigarItem, Operation},
    },
    Error,
    algo::HashMap,
    math::RoundDiv,
};

/// Upper triangle matrix, excluding diagonal (i < j).
pub struct TriangleMatrix<T> {
    side: usize,
    data: Vec<T>,
}

impl TriangleMatrix<()> {
    /// Returns iterator over all pairs `(i, j)` such that `0 <= i < j < side`.
    pub fn indices(side: usize) -> impl Iterator<Item = (usize, usize)> {
        (0..side - 1).flat_map(move |i| (i + 1..side).map(move |j| (i, j)))
    }
}

impl<T> TriangleMatrix<T> {
    /// Returns 0-based index (i,j) in the regular matrix n×n based on linear index k.
    pub fn from_linear_index(&self, k: usize) -> (usize, usize) {
        assert!(k < self.data.len());
        let under_root = (8 * self.data.len()).checked_sub(8 * k + 7).unwrap();
        let i = self.side.checked_sub(2 + (0.5 * (under_root as f64).sqrt() - 0.5).floor() as usize).unwrap();
        let j = (k + i * (i + 3) / 2 + 1).checked_sub(self.side * i).unwrap();
        (i, j)
    }

    #[inline]
    fn expected_len(side: usize) -> usize {
        side.saturating_sub(1) * side / 2
    }

    #[inline]
    fn to_linear_index(&self, i: usize, j: usize) -> usize {
        assert!(i < j && j < self.side, "Incorrect indices ({}, {}) to triangle matrix", i, j);
        (2 * self.side - 3 - i) * i / 2 + j - 1
    }

    /// Creates triangle matrix from linear storage (must have correct order: sorted first by row, then by column).
    pub fn from_linear(side: usize, data: Vec<T>) -> Self {
        assert_eq!(data.len(), Self::expected_len(side), "Incorrect triangle matrix size");
        Self { side, data }
    }

    /// Creates triangle matrix by running `f(i, j)` for corresponding indices.
    pub fn create(side: usize, f: impl FnMut((usize, usize)) -> T) -> Self {
        Self {
            side,
            data: TriangleMatrix::indices(side).map(f).collect(),
        }
    }

    /// Total number of elements in the triangle matrix.
    pub fn linear_len(&self) -> usize {
        self.data.len()
    }

    pub fn linear_data(&self) -> &[T] {
        &self.data
    }

    pub fn take_linear(self) -> Vec<T> {
        self.data
    }

    /// Size of the matrix side.
    pub fn side(&self) -> usize {
        self.side
    }

    /// Linear iterator over elements (first by row, second by column).
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.data.iter()
    }

    /// Returns `self[(i, j)]` if `i < j`, `self[(j, i)] if `j < i`, and `None` when `i == j`.
    pub fn get_symmetric(&self, i: usize, j: usize) -> Option<&T> {
        match i.cmp(&j) {
            Ordering::Less => Some(self.index((i, j))),
            Ordering::Equal => None,
            Ordering::Greater => Some(self.index((j, i))),
        }
    }
}

impl<T: Clone> TriangleMatrix<T> {
    pub fn new(side: usize, val: T) -> Self {
        Self {
            side,
            data: vec![val; Self::expected_len(side)],
        }
    }
}

impl<T> Index<(usize, usize)> for TriangleMatrix<T> {
    type Output = T;

    #[inline]
    fn index(&self, (i, j): (usize, usize)) -> &T {
        self.data.index(self.to_linear_index(i, j))
    }
}

impl<T> IndexMut<(usize, usize)> for TriangleMatrix<T> {
    #[inline]
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut T {
        self.data.index_mut(self.to_linear_index(i, j))
    }
}

const CAPACITY: usize = 4;
type KmerCache<'a> = HashMap<&'a [u8], SmallVec<[u32; CAPACITY]>>;

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
    accuracy: u8,
    threads: u16,
) -> Result<TriangleMatrix<f64>, Error> {
    let caches: Vec<_> = entries.iter().map(|entry| cache_kmers(entry.seq(), k)).collect();
    let alns = if threads == 1 {
        log::debug!("        Aligning sequences in 1 thread");
        let n = entries.len();
        let mut alns = Vec::with_capacity(n * (n - 1) / 2);
        let aligner = Aligner::new(penalties.clone(), accuracy);
        for (i, (entry1, cache1)) in entries.iter().zip(&caches).enumerate() {
            for (entry2, cache2) in entries[i + 1..].iter().zip(&caches[i + 1..]) {
                let kmer_matches = find_kmer_matches(cache1, cache2);
                alns.push(align(&aligner, entry1.seq(), entry2.seq(), &kmer_matches, k)?);
            }
        }
        TriangleMatrix::from_linear(n, alns)
    } else {
        divergences_multithread(entries, &caches, penalties, k, accuracy, threads)?
    };
    write_all(bam_path, entries, alns, k, accuracy)
}

fn divergences_multithread(
    entries: &[NamedSeq],
    caches: &[KmerCache],
    penalties: &Penalties,
    k: usize,
    accuracy: u8,
    threads: u16,
) -> Result<TriangleMatrix<(Cigar, i32)>, Error> {
    let threads = usize::from(threads);
    log::debug!("        Aligning sequences in {} threads", threads);
    let n = entries.len();
    let pairs: Arc<Vec<(u32, u32)>> = Arc::new(TriangleMatrix::indices(n).map(|(i, j)| (i as u32, j as u32)).collect());
    let n_pairs = pairs.len();
    let mut handles = Vec::with_capacity(threads);

    let entries = Arc::new(entries.to_vec());
    // Need to precompute all kmer matches as `caches` cannot be passed between threads.
    let all_kmer_matches = Arc::new(pairs.iter()
        .map(|&(i, j)| find_kmer_matches(&caches[i as usize], &caches[j as usize])).collect::<Vec<_>>());
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
            handles.push(thread::spawn(move || {
                assert!(start < end);
                let aligner = Aligner::new(penalties, accuracy);
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
    Ok(TriangleMatrix::from_linear(n, alns))
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
    alns: TriangleMatrix<(Cigar, i32)>,
    k: usize,
    accuracy: u8,
) -> Result<TriangleMatrix<f64>, Error> {
    let header = create_bam_header(entries, k, accuracy);
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
                return Err(Error::RuntimeError(format!(
                    "Generated invalid alignment between {} ({} bp) and {} ({} bp), CIGAR qlen {}, rlen {}",
                    query.name(), query.len(), refer.name(), rlen, cigar.query_len(), cigar.ref_len())));
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
pub fn load_edit_distances(path: impl AsRef<Path>, contigs: &ContigNames) -> Result<TriangleMatrix<u32>, Error> {
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
        let edit_dist: u32 = match record.aux(b"NM").map_err(|_| Error::InvalidData(format!(
                "BAM file {} does not contain NM field", ext::fmt::path(&path))))? {
            Aux::I8(val) => val.try_into().unwrap(),
            Aux::U8(val) => val.into(),
            Aux::I16(val) => val.try_into().unwrap(),
            Aux::U16(val) => val.into(),
            Aux::I32(val) => val.try_into().unwrap(),
            Aux::U32(val) => val,
            _ => return Err(Error::InvalidData(format!("Invalid value for NM field in {}", ext::fmt::path(&path)))),
        };
        matrix[(i, j)] = edit_dist;
    }
    if let Some(k) = matrix.iter().position(|&val| val == u32::MAX) {
        let (i, j) = matrix.from_linear_index(k);
        Err(Error::InvalidData(format!("BAM file {} does not contain alignment between contigs {} and {}",
            ext::fmt::path(&path), contigs.get_name(ContigId::new(i)), contigs.get_name(ContigId::new(j)))))
    } else {
        Ok(matrix)
    }
}

/// Calculates distance between two genotypes as minimum sum distance between permutations of contigs within genotypes.
pub fn genotype_distance(gt1: &Genotype, gt2: &Genotype, distances: &TriangleMatrix<u32>) -> u32 {
    let mut min_dist = u32::MAX;
    ext::vec::gen_permutations(gt1.ids(), |perm_ids1| {
        let dist: u32 = perm_ids1.iter().zip(gt2.ids())
            .map(|(i, j)| distances.get_symmetric(i.ix(), j.ix()).copied().unwrap_or(0))
            .sum();
        min_dist = min(min_dist, dist);
    });
    min_dist
}
