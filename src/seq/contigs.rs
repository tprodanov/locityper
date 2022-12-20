use std::collections::HashMap;
use std::fmt;

/// Contig identificator - newtype over u32.
/// Can be converted to `usize` using `id.ix()` method.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub struct ContigId(u32);

impl ContigId {
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
    names: Vec<String>,
    lengths: Vec<u32>,
    name_to_id: HashMap<String, ContigId>,
}

impl ContigNames {
    /// Create contig names from an iterator over pairs (name, length).
    /// Names must not repeat.
    pub fn new(it: impl Iterator<Item = (String, u32)>) -> Self {
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
        Self { names, lengths, name_to_id }
    }

    /// Get the number of contigs.
    pub fn len(&self) -> usize {
        self.names.len()
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
}

impl fmt::Debug for ContigNames {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ContigNames({:p})", &self)
    }
}
