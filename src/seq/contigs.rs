use std::{
    collections::HashMap,
    fmt, io,
};
use bio::io::fasta::{FastaRead, Index, Record};

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

    /// Get `u32` value of the contig id.
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
    /// Create contig names from an iterator over pairs (name, length).
    /// Names must not repeat.
    /// First argument: overall name of the contig set.
    pub fn new(tag: String, it: impl Iterator<Item = (String, u32)>) -> Self {
        let mut names = Vec::new();
        let mut lengths = Vec::new();
        let mut name_to_id = HashMap::new();

        for (name, length) in it {
            let contig_id = ContigId::new(names.len());
            if let Some(prev_id) = name_to_id.insert(name.clone(), contig_id) {
                panic!("Contig {} appears twice in the contig list ({} and {})!", name, prev_id, contig_id);
            }
            names.push(name);
            lengths.push(length);
        }

        assert!(names.len() < (std::u16::MAX as usize - 1) / 2,
            "Cannot supports {} contigs, maximum number is {}", names.len(), (std::u16::MAX as usize - 1) / 2);

        names.shrink_to_fit();
        lengths.shrink_to_fit();
        name_to_id.shrink_to_fit();
        Self { tag, names, lengths, name_to_id }
    }

    /// Creates contig names from FASTA index.
    /// First argument: overall name of the contig name set.
    pub fn from_index(tag: String, index: &Index) -> Self {
        Self::new(tag, index.sequences().into_iter().map(|seq| (seq.name, u32::try_from(seq.len).unwrap())))
    }

    /// Reads all entries from fasta and saves them into memory.
    /// Returns pair (ContigNames, Vec<nt sequence>).
    pub fn from_fasta<R: FastaRead>(tag: String, reader: &mut R) -> io::Result<(Self, Vec<Vec<u8>>)> {
        let mut names_lengths = Vec::new();
        let mut seqs = Vec::new();
        let mut record = Record::new();
        loop {
            reader.read(&mut record)?;
            if record.is_empty() {
                return Ok((Self::new(tag, names_lengths.into_iter()), seqs));
            }

            let mut ref_seq = record.seq().to_vec();
            super::standardize(&mut ref_seq);
            names_lengths.push((record.id().to_string(), u32::try_from(ref_seq.len()).unwrap()));
            seqs.push(ref_seq.to_vec());
        }
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// Get the number of contigs.
    pub fn len(&self) -> usize {
        self.names.len()
    }

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

    /// Returns iterator over all contig IDs.
    pub fn ids(&self) -> impl Iterator<Item = ContigId> + std::iter::ExactSizeIterator {
        (0..u16::try_from(self.len()).unwrap()).map(ContigId::new)
    }

    /// Returns iterator over all contig names.
    pub fn names(&self) -> std::slice::Iter<'_, String> {
        self.names.iter()
    }

    /// Returns iterator over all contig lengths.
    pub fn lengths(&self) -> std::iter::Cloned<std::slice::Iter<'_, u32>> {
        self.lengths.iter().cloned()
    }

    /// Returns true if the interval is in bounds for its contig (interval end does not exceed contig length).
    pub fn in_bounds(&self, interval: &super::Interval) -> bool {
        interval.end() <= self.get_len(interval.contig_id())
    }
}

impl fmt::Debug for ContigNames {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ContigNames({}, {} entries)", self.tag, self.names.len())
    }
}

impl fmt::Display for ContigNames {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}
