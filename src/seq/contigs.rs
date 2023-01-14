use std::collections::HashMap;
use std::fmt;

/// Contig identificator - newtype over u32.
/// Can be converted to `usize` using `id.ix()` method.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub struct ContigId(u32);

impl ContigId {
    /// Creates a new ContigId.
    #[inline]
    pub fn new(val: u32) -> ContigId {
        ContigId(val)
    }

    /// Get `u32` value of the contig id.
    #[inline]
    pub fn get(self) -> u32 {
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
    set_name: String,
    names: Vec<String>,
    lengths: Vec<u32>,
    name_to_id: HashMap<String, ContigId>,
}

impl ContigNames {
    /// Create contig names from an iterator over pairs (name, length).
    /// Names must not repeat.
    /// First argument: overall name of the contig name set.
    pub fn new(set_name: String, it: impl Iterator<Item = (String, u32)>) -> Self {
        let mut names = Vec::new();
        let mut lengths = Vec::new();
        let mut name_to_id = HashMap::new();

        for (name, length) in it {
            let contig_id = ContigId(names.len() as u32);
            if let Some(prev_id) = name_to_id.insert(name.clone(), contig_id) {
                panic!("Contig {} appears twice in the contig list ({} and {})!", name, prev_id, contig_id);
            }
            names.push(name);
            lengths.push(length);
        }

        names.shrink_to_fit();
        lengths.shrink_to_fit();
        name_to_id.shrink_to_fit();
        Self { set_name, names, lengths, name_to_id }
    }

    /// Creates contig names from FASTA index.
    /// First argument: overall name of the contig name set.
    pub fn from_index(set_name: String, index: &bio::io::fasta::Index) -> Self {
        Self::new(set_name, index.sequences().into_iter().map(|seq| (seq.name, seq.len as u32)))
    }

    /// Get the number of contigs.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn set_name(&self) -> &str {
        &self.set_name
    }

    /// Get contig name from an id.
    pub fn name(&self, id: ContigId) -> &str {
        &self.names[id.ix()]
    }

    /// Get contig length from an id.
    pub fn length(&self, id: ContigId) -> u32 {
        self.lengths[id.ix()]
    }

    /// Get contig id from its name.
    pub fn id(&self, name: &str) -> ContigId {
        self.name_to_id[name]
    }

    /// Returns contig id, if it is available.
    pub fn get_id(&self, name: &str) -> Option<ContigId> {
        self.name_to_id.get(name).copied()
    }
}

impl fmt::Debug for ContigNames {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ContigNames({}, {} entries)", self.set_name, self.names.len())
    }
}

impl fmt::Display for ContigNames {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}
