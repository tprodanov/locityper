use std::{
    fmt,
    io::{self, BufRead},
    fs::File,
    sync::Arc,
    path::Path,
};
use smallvec::SmallVec;
use fnv::FnvHashMap;
use bio::io::fasta::{self, FastaRead};
use crate::{
    err::{Error, add_path},
    ext::{self, vec::VecOrNone},
    seq::kmers::KmerCounts,
};

/// Contig identificator - newtype over u16.
/// Can be converted to `usize` using `id.ix()` method.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub struct ContigId(u16);

impl ContigId {
    /// Creates a new ContigId.
    pub fn new<T>(val: T) -> ContigId
    where T: TryInto<u16>,
          T::Error: fmt::Debug,
    {
        ContigId(val.try_into().expect("Contig ID too large"))
    }

    /// Get `u16` value of the contig id.
    #[inline]
    pub fn get(self) -> u16 {
        self.0
    }

    /// Converts `ContigId` into `usize`.
    #[inline]
    pub fn ix(self) -> usize {
        self.0 as usize
    }
}

impl fmt::Display for ContigId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Id({})", self.0)
    }
}

/// Structure that stores contig names.
/// Contig name and lengths can be accessed using `ContigId`.
/// Additionally, contig id can be extracting knowing contig name using `id()` method.
#[derive(Clone)]
pub struct ContigNames {
    tag: String,
    names: Vec<String>,
    lengths: Vec<u32>,
    name_to_id: FnvHashMap<String, ContigId>,
}

impl ContigNames {
    /// Create contig names from an iterator over pairs (name, length).
    /// Names must not repeat.
    /// First argument: overall name of the contig set.
    pub fn new(tag: impl Into<String>, it: impl Iterator<Item = (String, u32)>) -> Self {
        let mut names = Vec::new();
        let mut lengths = Vec::new();
        let mut name_to_id = FnvHashMap::default();

        for (name, length) in it {
            let contig_id = ContigId::new(names.len());
            if let Some(prev_id) = name_to_id.insert(name.clone(), contig_id) {
                panic!("Contig {} appears twice in the contig list ({} and {})!", name, prev_id, contig_id);
            }
            names.push(name);
            lengths.push(length);
        }

        const MAX_CONTIGS: usize = (std::u16::MAX as usize - 1) / 2;
        assert!(names.len() < MAX_CONTIGS,
            "Cannot support {} contigs, maximum number is {}", names.len(), MAX_CONTIGS);

        names.shrink_to_fit();
        lengths.shrink_to_fit();
        name_to_id.shrink_to_fit();
        Self {
            tag: tag.into(),
            names, lengths, name_to_id,
        }
    }

    /// Creates contig names from FASTA index.
    /// First argument: overall name of the contig name set.
    pub fn from_index(tag: impl Into<String>, index: &fasta::Index) -> Self {
        Self::new(tag, index.sequences().into_iter().map(|seq| (seq.name, u32::try_from(seq.len).unwrap())))
    }

    /// Reads all entries from fasta and saves them into memory.
    ///
    /// Descriptions are stored in either `&mut Vector<Option<String>>` or `()`, if they are not needed.
    ///
    /// Returns tuple
    /// - ContigNames,
    /// - Vector of contig sequences.
    pub fn load_fasta(
        tag: impl Into<String>,
        stream: impl BufRead,
        mut descriptions: impl VecOrNone<Option<String>>,
    ) -> io::Result<(Arc<Self>, Vec<Vec<u8>>)>
    {
        let mut reader = fasta::Reader::from_bufread(stream);
        let mut names_lengths = Vec::new();
        let mut seqs = Vec::new();
        let mut record = fasta::Record::new();
        loop {
            reader.read(&mut record)?;
            if record.is_empty() {
                let contigs = Self::new(tag, names_lengths.into_iter());
                return Ok((Arc::new(contigs), seqs));
            }

            let mut ref_seq = record.seq().to_vec();
            super::standardize(&mut ref_seq);
            names_lengths.push((record.id().to_string(), u32::try_from(ref_seq.len()).unwrap()));
            seqs.push(ref_seq.to_vec());
            descriptions.push(record.desc().map(str::to_string));
        }
    }

    /// Loads indexed fasta and stored contig names and lengths.
    pub fn load_indexed_fasta(
        tag: impl Into<String>,
        filename: &(impl AsRef<Path> + fmt::Debug),
    ) -> Result<(Arc<Self>, fasta::IndexedReader<File>), Error>
    {
        let fasta = fasta::IndexedReader::from_file(&filename)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))
            .map_err(add_path!(filename))?;
        let contigs = ContigNames::from_index(tag, &fasta.index);
        Ok((Arc::new(contigs), fasta))
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// Get the number of contigs.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Returns the tag of the contig set.
    pub fn tag(&self) -> &str {
        &self.tag
    }

    /// Get contig name from an id.
    pub fn get_name(&self, id: ContigId) -> &str {
        &self.names[id.ix()]
    }

    /// Get contig length from an id.
    pub fn get_len(&self, id: ContigId) -> u32 {
        self.lengths[id.ix()]
    }

    /// Get contig id from its name.
    pub fn get_id(&self, name: &str) -> ContigId {
        self.name_to_id[name]
    }

    /// Returns contig id, if it is available.
    pub fn try_get_id(&self, name: &str) -> Result<ContigId, Error> {
        self.name_to_id.get(name).copied()
            .ok_or_else(|| Error::ParsingError(format!("Unknown contig {:?}", name)))
    }

    /// Returns iterator over all contig IDs.
    pub fn ids(&self) -> impl Iterator<Item = ContigId> + std::iter::ExactSizeIterator {
        (0..u16::try_from(self.len()).unwrap()).map(ContigId::new)
    }

    /// Returns iterator over all contig names.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Returns iterator over all contig lengths.
    pub fn lengths(&self) -> &[u32] {
        &self.lengths
    }

    /// Returns true if the interval is in bounds for its contig (interval end does not exceed contig length).
    pub fn in_bounds(&self, interval: &super::Interval) -> bool {
        interval.end() <= self.get_len(interval.contig_id())
    }

    /// Returns sum length of all contigs.
    pub fn genome_size(&self) -> u64 {
        self.lengths.iter().copied().map(u64::from).sum()
    }
}

impl fmt::Debug for ContigNames {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ContigNames({}, {} entries)", self.tag, self.names.len())
    }
}

/// Contigs, their complete sequences, and k-mer counts.
pub struct ContigSet {
    contigs: Arc<ContigNames>,
    seqs: Vec<Vec<u8>>,
    kmer_counts: KmerCounts,
}

impl ContigSet {
    /// Loads contigs, their sequences and k-mer counts from two files.
    /// `descriptions` can either be `&mut Vec<String>`, or `()`.
    pub fn load(
        tag: impl Into<String>,
        fasta_filename: &Path,
        kmers_filename: &Path,
        descriptions: impl VecOrNone<Option<String>>,
    ) -> Result<Self, Error>
    {
        let (contigs, seqs) = ContigNames::load_fasta(tag, ext::sys::open(fasta_filename)?, descriptions)
            .map_err(add_path!(fasta_filename))?;
        let kmer_counts = KmerCounts::load(ext::sys::open(kmers_filename)?, contigs.lengths())?;
        Ok(Self { contigs, seqs, kmer_counts })
    }

    /// Returns the number of contigs in the set.
    pub fn len(&self) -> usize {
        self.seqs.len()
    }

    /// Returns true if there are no contigs in the set.
    pub fn is_empty(&self) -> bool {
        self.seqs.is_empty()
    }

    /// Returns inner Contig names.
    pub fn contigs(&self) -> &Arc<ContigNames> {
        &self.contigs
    }

    /// Returns the tag of the contig set.
    pub fn tag(&self) -> &str {
        self.contigs.tag()
    }

    /// Returns all sequences.
    pub fn seqs(&self) -> &[Vec<u8>] {
        &self.seqs
    }

    /// Returns the sequence of the corresponding contig.
    pub fn get_seq(&self, contig_id: ContigId) -> &[u8] {
        &self.seqs[contig_id.ix()]
    }

    /// Returns k-mer counts.
    pub fn kmer_counts(&self) -> &KmerCounts {
        &self.kmer_counts
    }
}

type GtStorage = SmallVec<[ContigId; 4]>;

/// Genotype: a tuple of contigs.
#[derive(Clone)]
pub struct Genotype {
    ids: GtStorage,
    name: String,
}

impl Genotype {
    pub fn new(ids_iter: impl Iterator<Item = ContigId>, contigs: &ContigNames) -> Self {
        let mut ids = GtStorage::new();
        let mut name = String::new();
        for id in ids_iter {
            if !name.is_empty() {
                name.push(',');
            }
            name.push_str(contigs.get_name(id));
            ids.push(id);
        }
        assert!(!ids.is_empty(), "Empty genotypes are not allowed");
        Self { ids, name }
    }

    /// Returns contigs within the genotype (they can repeat).
    pub fn ids(&self) -> &[ContigId] {
        &self.ids
    }

    /// Returns genotype ploidy (number of contigs).
    pub fn ploidy(&self) -> usize {
        self.ids.len()
    }

    /// Returns genotype name (all contig names through comma).
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Parses genotype (contig names through comma).
    pub fn parse(s: &str, contigs: &ContigNames) -> Result<Self, Error> {
        let ids = s.split(',').map(|name| contigs.try_get_id(name))
            .collect::<Result<GtStorage, Error>>()?;
        assert!(!ids.is_empty(), "Empty genotypes are not allowed");
        Ok(Self {
            ids,
            name: s.to_owned(),
        })
    }

    /// Generates all genotype of given ploidy.
    pub fn generate_all(avail_ids: &[ContigId], contigs: &ContigNames, ploidy: usize) -> Vec<Self> {
        let mut res = Vec::with_capacity(ext::vec::count_combinations_with_repl(avail_ids.len(), ploidy));
        ext::vec::gen_combinations_with_repl(avail_ids, ploidy,
            |ids| res.push(Self::new(ids.iter().copied(), contigs)));
        res
    }
}

impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(&self.name)
    }
}
