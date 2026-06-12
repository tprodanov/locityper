use std::{
    fmt,
    io::{self, BufRead},
    fs::File,
    sync::Arc,
    path::Path,
};
use smallvec::SmallVec;
use htslib::bam;
use bio::io::fasta;
use itertools::izip;
use crate::{
    err::{Error, error, add_path},
    ext,
    seq::{
        fastx,
        counts::KmerCounts,
    },
    algo::{HashMap, HashSet, IntMap},
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
    #[inline(always)]
    pub fn get(self) -> u16 {
        self.0
    }

    /// Converts `ContigId` into `usize`.
    #[inline(always)]
    pub fn ix(self) -> usize {
        usize::from(self.0)
    }
}

impl fmt::Display for ContigId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Id({})", self.0)
    }
}

impl std::hash::Hash for ContigId {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        hasher.write_u16(self.0)
    }
}

impl nohash::IsEnabled for ContigId {}

impl Into<usize> for ContigId {
    #[inline(always)]
    fn into(self) -> usize {
        self.ix()
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
    /// Creates new contig names.
    pub fn new(tag: impl Into<String>) -> Self {
        Self {
            tag: tag.into(),
            names: Vec::new(),
            lengths: Vec::new(),
            name_to_id: HashMap::default(),
        }
    }

    pub fn add<T>(&mut self, name: String, length: T) -> crate::Result<ContigId>
    where T: TryInto<u32> + std::fmt::Display + Copy,
    {
        let contig_id = TryInto::<u16>::try_into(self.names.len())
            .map_err(|_| error!(InvalidData, "[{}] Too many contigs (at most {} supported)", self.tag, u16::MAX))
            .map(ContigId::new)?;
        let length: u32 = length.try_into()
            .map_err(|_| error!(InvalidData, "[{}] Contig {} is too long ({} bp do not fit into 32 bits)",
                self.tag, name, length))?;
        if let Some(old_id) = self.name_to_id.insert(name.clone(), contig_id) {
            *self.name_to_id.get_mut(&name).unwrap() = old_id;
            panic!("[{}] Contig {} appears twice in the contig list", self.tag, name);
        }
        self.names.push(name);
        self.lengths.push(length);
        Ok(contig_id)
    }

    /// Creates contig names from FASTA index.
    /// First argument: overall name of the contig name set.
    pub fn from_index(tag: impl Into<String>, index: &fasta::Index) -> crate::Result<Self> {
        let mut contigs = Self::new(tag);
        for seq in index.sequences().into_iter() {
            contigs.add(seq.name, seq.len)?;
        }
        Ok(contigs)
    }

    /// Creates contig names from BAM header.
    pub fn from_bam_header(tag: impl Into<String>, header: &bam::HeaderView) -> crate::Result<Self> {
        let mut contigs = Self::new(tag);
        for i in 0..header.target_count() {
            let byte_name = header.tid2name(i);
            let name = String::from_utf8(byte_name.to_vec())
                .map_err(|_| Error::Utf8("contig name", byte_name.to_vec()))?;
            contigs.add(name, header.target_len(i).expect("Contig length undefined in BAM"))?;
        }
        Ok(contigs)
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
    pub fn contains(&self, name: &str) -> bool {
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

    /// Adds a different name that points to an existing contig. Only relevant for `get_id`/`try_get_id`.
    pub fn add_synonym(&mut self, name: String, id: ContigId) {
        let old_val = self.name_to_id.insert(name, id);
        assert!(old_val.is_none(), "Cannot add synonym (duplicated name)");
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
}

impl ContigSet {
    pub fn new(contigs: Arc<ContigNames>, seqs: Vec<Vec<u8>>) -> Self {
        Self { contigs, seqs }
    }

    /// Loads contigs, their sequences and k-mer counts from two files.
    /// `descriptions` can either be `&mut Vec<String>`, or `()`.
    pub fn load(
        tag: impl Into<String>,
        fasta_filename: &Path,
    ) -> crate::Result<Self>
    {
        let mut fasta_reader = fastx::Reader::from_path(fasta_filename)?;
        let mut contigs = ContigNames::new(tag);
        let mut seqs = Vec::new();
        fasta_reader.read_all(|name, seq| {
            contigs.add(name, seq.len())?;
            seqs.push(seq.to_owned());
            Ok(())
        })?;

        Ok(Self::new(Arc::new(contigs), seqs))
    }

    pub fn load_with_kmer_counts(
        tag: impl Into<String>,
        fasta_filename: &Path,
        kmers_filename: &Path,
    ) -> crate::Result<(Self, KmerCounts)> {
        let set = Self::load(tag, fasta_filename)?;
        // k-mer counts file contains both off-target and regular k-mer counts.
        // For now, we only need off-target counts, so we only read k-mer counts once.
        let mut kmers_file = ext::sys::open(kmers_filename)?;
        let kmer_counts = KmerCounts::load(&mut kmers_file).map_err(add_path!(kmers_filename))?;
        kmer_counts.validate(&set.contigs)?;
        Ok((set, kmer_counts))
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

    /// Returns sequence by contig name, if found.
    pub fn get_seq_by_name(&self, name: &str) -> Option<&[u8]> {
        self.contigs.try_get_id(name).map(|id| self.get_seq(id))
    }

    /// Extract a subset of this contig set counts removing leave-out haplotype names.
    /// Returns subset indices and a new contig set.
    /// Subset indices will follow be sorted, so in the same order as this contigs set.
    pub fn extract_subset(
        &self,
        leave_out: &HashSet<String>,
        disc_filename: &Path,
    ) -> crate::Result<(Vec<usize>, Self)> {
        let disc_haps = if disc_filename.exists() {
            load_discarded_haplotypes(ext::sys::open(&disc_filename)?, &self.contigs)?
        } else {
            Default::default()
        };
        if !discarded_all_identical(&disc_haps) {
            log::warn!("[{}] Haplotypes were previously pruned (~ in discarded haplotypes). \
                Leave-out genotyping may lose relevant previously discarded haplotypes",
                self.contigs.tag());
        }

        let mut ixs = Vec::new();
        let mut contigs = ContigNames::new(self.contigs.tag());
        let mut discarded = Vec::new();
        let mut replaced = Vec::new();
        for (i, (name, length)) in izip!(self.contigs.names(), self.contigs.lengths()).enumerate() {
            let mut save_name = name.to_string();
            let mut found = false;
            let mut replaced_now = false;
            if sample_or_haplotype_in_set(&name, &leave_out) {
                if let Some(discarded) = disc_haps.get(&ContigId::new(i)) {
                    for (oth_name, is_identical) in discarded {
                        if *is_identical && !sample_or_haplotype_in_set(oth_name, &leave_out) {
                            replaced_now = true;
                            replaced.push(format!("{}->{}", name, oth_name));
                            save_name = oth_name.to_string();
                            found = true;
                            break;
                        }
                    }
                }
                if !found {
                    discarded.push(name.to_string());
                    continue;
                }
            }
            ixs.push(i);
            let contig_id = contigs.add(save_name, *length)?;
            if replaced_now {
                contigs.add_synonym(name.to_string(), contig_id);
            }
        }

        log::debug!("    [{}] Leave-out: discarded {} haplotype(s) [{}], replaced {} haplotype(s) with identical [{}]",
            self.contigs.tag(),
            discarded.len(), ext::vec::join_up_to(&discarded, 5),
            replaced.len(), ext::vec::join_up_to(&replaced, 5));
        let seqs = ixs.iter().map(|&i| self.seqs[i].clone()).collect();
        Ok((ixs, ContigSet::new(Arc::new(contigs), seqs)))
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

pub type DiscardedHaplotypes = IntMap<ContigId, Vec<(String, bool)>>;

/// Loads `discarded_haplotypes.txt` file, with lines `haplotype (= or ~) haplotype2, haplotype3, ...`.
/// Returns map (retained contig id) -> Vec of (discarded name, is_identical).
pub fn load_discarded_haplotypes(
    f: impl BufRead,
    contigs: &ContigNames,
) -> crate::Result<DiscardedHaplotypes>
{
    // Map { retained -> (discarded, identical) }, where retained is not found in `contigs`.
    let mut corresp_unknown: HashMap<String, Vec<(String, bool)>> = HashMap::default();
    // Same, but retained haplotype if found, and replaced with its id.
    let mut corresp = IntMap::default();

    for line in f.lines() {
        let line = line.map_err(add_path!(!))?;
        let split: Vec<_> = line.split_whitespace().collect();
        if split.len() < 3 {
            return Err(error!(InvalidInput, "Each line in discarded haplotypes must have at least 3 columns"));
        }
        let identical = split[1] == "=";

        let mut rhs = Vec::with_capacity(split.len() - 2);
        for contig in &split[2..] {
            let contig = contig.strip_suffix(',').unwrap_or(contig);
            if contigs.contains(contig) {
                log::warn!("{}Haplotype {} is marked as discarded, but present in the haplotypes fasta",
                    if contigs.tag().is_empty() { String::new() } else { format!("[{}] ", contigs.tag()) }, contig);
                continue;
            }
            rhs.push((contig.to_string(), identical));
            if let Some(v) = corresp_unknown.remove(contig) {
                rhs.extend(v.into_iter().map(|(contig2, identical2)| (contig2, identical && identical2)));
            }
        }

        match contigs.try_get_id(split[0]) {
            Some(id) => corresp.insert(id, rhs),
            None => corresp_unknown.insert(split[0].to_string(), rhs),
        };
    }

    if let Some(contig) = corresp_unknown.keys().next() {
        log::warn!("{}Haplotype {} is on the left side of the discarded haplotypes, but missing from the fasta",
            if contigs.tag().is_empty() { String::new() } else { format!("[{}] ", contigs.tag()) }, contig);
    }
    Ok(corresp)
}

/// Returns true if all discarded haplotypes are identical to one of the retained haplotypes.
/// Will be false if `locityper prune` was used.
pub fn discarded_all_identical(disc_haplotypes: &DiscardedHaplotypes) -> bool {
    disc_haplotypes.values().flat_map(|v| v.iter()).all(|(_, identical)| *identical)
}

/// Returns true if `name` is in `set`, or, if `name` ends with `.N` / `_N` / `-N` (N is a single digit)
/// and corresponding prefix is in the set.
pub fn sample_or_haplotype_in_set(name: &str, set: &HashSet<String>) -> bool {
    let bytes = name.as_bytes();
    let n = bytes.len();
    set.contains(name) ||
        (n > 2 && set.contains(&std::str::from_utf8(&bytes[..n - 2]).unwrap() as &str)
        && b'0' <= bytes[n - 1] && bytes[n - 1] <= b'9'
        && (bytes[n - 2] == b'.' || bytes[n - 2] == b'_' || bytes[n - 2] == b'-'))
}
