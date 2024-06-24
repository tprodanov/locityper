use std::{
    fmt, io,
    fs::File,
    sync::Arc,
    path::Path,
};
use smallvec::SmallVec;
use htslib::bam;
use bio::io::fasta;
use crate::{
    err::{Error, error, add_path},
    ext,
    seq::{
        fastx,
        NamedSeq,
        kmers::KmerCounts,
    },
    algo::HashMap,
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
    name_to_id: HashMap<String, ContigId>,
}

impl ContigNames {
    /// Create new, empty, contig names.
    pub fn empty() -> Self {
        ContigNames {
            tag: "UNINIT".to_string(),
            names: Vec::new(),
            lengths: Vec::new(),
            name_to_id: HashMap::default(),
        }
    }

    /// Create contig names from an iterator over pairs (name, length).
    /// Names must not repeat.
    /// First argument: overall name of the contig set.
    pub fn new<T: TryInto<u32> + std::fmt::Display + Copy>(
        tag: impl Into<String>,
        mut it: impl Iterator<Item = crate::Result<(String, T)>>
    ) -> crate::Result<Self>
    {
        let mut names = Vec::new();
        let mut lengths = Vec::new();
        let mut name_to_id = HashMap::default();

        while let Some((name, length)) = it.next().transpose()? {
            let length: u32 = length.try_into()
                .map_err(|_| error!(InvalidData, "Contig {} is too long ({} bp)", name, length))?;
            let contig_id = ContigId::new(names.len());
            if let Some(prev_id) = name_to_id.insert(name.clone(), contig_id) {
                panic!("Contig {} appears twice in the contig list ({} and {})!", name, prev_id, contig_id);
            }
            names.push(name);
            lengths.push(length);
        }

        const MAX_CONTIGS: usize = (std::u16::MAX as usize - 1) / 2;
        if names.is_empty() {
            return Err(error!(InvalidData, "Contig set is empty"));
        } else if names.len() >= MAX_CONTIGS {
            return Err(error!(InvalidData, "Too many contigs ({}), can support at most {}",
                names.len(), MAX_CONTIGS));
        }

        names.shrink_to_fit();
        lengths.shrink_to_fit();
        name_to_id.shrink_to_fit();
        Ok(Self {
            tag: tag.into(),
            names, lengths, name_to_id,
        })
    }

    /// Creates contig names from FASTA index.
    /// First argument: overall name of the contig name set.
    pub fn from_index(tag: impl Into<String>, index: &fasta::Index) -> crate::Result<Self> {
        Self::new(tag, index.sequences().into_iter().map(|seq| Ok((seq.name, seq.len))))
    }

    /// Creates contig names from BAM header.
    pub fn from_bam_header(tag: impl Into<String>, header: &bam::HeaderView) -> crate::Result<Self> {
        Self::new(tag, (0..header.target_count()).map(|i| {
            let byte_name = header.tid2name(i);
            let name = String::from_utf8(byte_name.to_vec())
                .map_err(|_| Error::Utf8("contig name", byte_name.to_vec()))?;
            Ok((name, header.target_len(i).unwrap()))
        }))
    }

    /// Loads indexed fasta and stored contig names and lengths.
    pub fn load_indexed_fasta(
        tag: impl Into<String>,
        filename: &Path,
    ) -> crate::Result<(Self, fasta::IndexedReader<File>)>
    {
        let fasta = fasta::IndexedReader::from_file(&filename)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))
            .map_err(add_path!(filename))?;
        let contigs = ContigNames::from_index(tag, &fasta.index)?;
        Ok((contigs, fasta))
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
    pub fn try_get_id(&self, name: &str) -> Option<ContigId> {
        self.name_to_id.get(name).copied()
    }

    /// Returns true if the contig exists.
    pub fn exists(&self, name: &str) -> bool {
        self.name_to_id.contains_key(name)
    }

    /// Returns iterator over all contig IDs.
    pub fn ids(&self) -> impl Iterator<Item = ContigId> + std::iter::ExactSizeIterator {
        (0..u16::try_from(self.len()).unwrap()).map(ContigId::new)
    }

    /// Returns iterator over all contig names.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn take_names(self) -> Vec<String> {
        self.names
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

/// Version of the reference genome, supported by this program.
#[derive(Clone, Copy, Debug)]
#[repr(u8)]
pub enum GenomeVersion {
    Chm13,
    GRCh38,
    GRCh37,
}

impl GenomeVersion {
    pub fn to_str(self) -> &'static str {
        match self {
            Self::Chm13 => "CHM13",
            Self::GRCh38 => "GRCh38",
            Self::GRCh37 => "GRCh37",
        }
    }

    /// Based on the length of the first chromosome, tries to identify reference version.
    pub fn guess(contigs: &ContigNames) -> Option<Self> {
        let chr1_len = contigs.try_get_id("chr1")
            .or_else(|| contigs.try_get_id("1"))
            .map(|id| contigs.get_len(id))?;
        match chr1_len {
            248_387_328 => Some(Self::Chm13),
            248_956_422 => Some(Self::GRCh38),
            249_250_621 => Some(Self::GRCh37),
            _ => None,
        }
    }
}

impl fmt::Display for GenomeVersion {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(self.to_str())
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
    ) -> crate::Result<Self>
    {
        let mut fasta_reader = fastx::Reader::from_path(fasta_filename)?;
        let mut named_seqs = fasta_reader.read_all()?;
        let contigs = ContigNames::new(tag,
            named_seqs.iter_mut().map(|entry| Ok((std::mem::replace(entry.name_mut(), String::new()), entry.len()))))?;
        let seqs: Vec<_> = named_seqs.into_iter().map(NamedSeq::take_seq).collect();
        // k-mer counts file contains both off-target and regular k-mer counts.
        // For now, we only need off-target counts, so we only read k-mer counts once.
        let mut kmers_file = ext::sys::open(kmers_filename)?;
        let kmer_counts = KmerCounts::load(&mut kmers_file).map_err(add_path!(kmers_filename))?;
        kmer_counts.validate(&contigs)?;
        Ok(Self {
            contigs: Arc::new(contigs),
            seqs, kmer_counts,
        })
    }

    pub fn new(contigs: Arc<ContigNames>, seqs: Vec<Vec<u8>>, kmer_counts: KmerCounts) -> Self {
        Self { contigs, seqs, kmer_counts }
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
    pub fn new(ids: &[ContigId], contigs: &ContigNames) -> Self {
        let ids: GtStorage = ids.iter().copied().collect();
        let mut name = String::new();
        for &id in ids.iter() {
            if !name.is_empty() {
                name.push(',');
            }
            name.push_str(contigs.get_name(id));
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
    pub fn parse(s: &str, contigs: &ContigNames) -> crate::Result<Self> {
        let ids = s.split(',').map(|name| contigs.try_get_id(name)
            .ok_or_else(|| error!(ParsingError, "Unknown contig {:?}", name)))
            .collect::<crate::Result<GtStorage>>()?;
        assert!(!ids.is_empty(), "Empty genotypes are not allowed");
        Ok(Self {
            ids,
            name: s.to_owned(),
        })
    }
}

impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(&self.name)
    }
}
