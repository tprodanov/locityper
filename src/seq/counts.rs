use std::{
    io::{self, BufRead},
    fs,
    path::PathBuf,
    cmp::{min, max, Ord},
    process::{Stdio, Command},
};
use varint_rs::{VarintWriter, VarintReader};
use crate::{
    seq::{
        self,
        kmers::{kmers, Kmer, CANONICAL},
        ContigNames, ContigId
    },
    algo::HashMap,
    err::{Error, error, add_path},
    bg::ser::json_get,
    ext::sys::PipeGuard,
};

/// Store k-mer counts as u16.
pub type KmerCount = u16;

/// Parses integer, fails if it is not integer, and returns at most `max_value` (must fit in `KmerCount`).
#[inline]
fn parse_count(s: &str, max_value: u64) -> Result<KmerCount, std::num::ParseIntError> {
    s.parse::<u64>().map(|x| KmerCount::try_from(min(x, max_value)).unwrap())
}

/// Stores k-mer counts for each input k-mer across a set of sequences.
#[derive(Clone)]
pub struct KmerCounts {
    k: u32,
    /// Maximum value that can appear.
    /// Will be equal to `KmerCount::MAX` if Jellyfish was run with appropriate counter-len,
    /// and will be smaller otherwise.
    max_value: KmerCount,
    /// k-mer counts for every contig.
    /// If value == KmerCount::MAX, k-mer count may be too high (cannot be stored in N bytes).
    counts: Vec<Vec<KmerCount>>,
}

impl KmerCounts {
    /// Returns k-mer size.
    pub fn k(&self) -> u32 {
        self.k
    }

    /// Returns the maximum value, that can be stored.
    pub fn max_value(&self) -> KmerCount {
        self.max_value
    }

    /// Creates new k-mer counts, where all values are zero.
    pub fn new_zeroed(k: u32, lengths: impl Iterator<Item = usize>) -> Self {
        let k_usize = k as usize;
        let mut counts = Vec::new();
        for len in lengths {
            counts.push(vec![0; (len + 1).saturating_sub(k_usize)]);
        }
        Self {
            k, counts,
            max_value: KmerCount::MAX,
        }
    }

    /// Returns iterator over k-mer counts for every contig.
    pub fn iter(&self) -> std::slice::Iter<'_, Vec<KmerCount>> {
        self.counts.iter()
    }

    /// Returns the number of entries (one per each sequences) in the collection.
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// Returns all k-mer counts for the corresponding contig.
    /// Length: contig len - k + 1.
    pub fn get(&self, ix: usize) -> &[KmerCount] {
        &self.counts[ix]
    }

    /// Returns all k-mer counts for the first contig (must be the only contig).
    /// Length: contig len - k + 1.
    pub fn get_first(&self) -> &[KmerCount] {
        assert_eq!(self.counts.len(), 1,
            "Cannot call KmerCounts.get_first when the number of contigs is different than one");
        &self.counts[0]
    }

    /// Takes k-mer counts for the last contig out of this collection.
    pub fn pop(&mut self) -> Vec<KmerCount> {
        self.counts.pop().expect("No k-mer counts present in the collection")
    }

    /// Returns all k-mer counts for the first contig (must be the only contig).
    /// Length: contig len - k + 1.
    pub fn take_first(mut self) -> Vec<KmerCount> {
        assert_eq!(self.counts.len(), 1,
            "Cannot call KmerCounts.take_first when the number of contigs is different than one");
        self.counts.pop().unwrap()
    }

    /// Saves collection into binary format:
    /// First two bytes: k-mer size and counter length in bytes.
    /// Next: N = number of contigs as 4-byte varint.
    /// Next: N times number of k-mers as 4-byte varint and
    ///     all k-mer counts consecutively in as 8-byte varints (for future compatibility).
    pub fn save<W>(&self, f: &mut W) -> io::Result<()>
    where W: VarintWriter<Error = io::Error>,
    {
        f.write(u8::try_from(self.k).unwrap())?;
        assert!(self.max_value.leading_zeros() + self.max_value.count_ones() == KmerCount::BITS,
            "Invalid maximum value for k-mer counts: {}", self.max_value);
        f.write(u8::try_from(self.max_value.count_ones() / 8).unwrap())?;
        f.write_u32_varint(self.counts.len() as u32)?;

        for counts in self.counts.iter() {
            f.write_u32_varint(counts.len() as u32)?;
            for &count in counts {
                f.write_u64_varint(u64::from(count))?;
            }
        }
        Ok(())
    }

    /// Loads k-mer counts from a binary file.
    pub fn load<R>(f: &mut R) -> io::Result<Self>
    where R: VarintReader<Error = io::Error>,
    {
        let k = u32::from(f.read()?);
        let byte_len = u32::from(f.read()?);
        assert!(byte_len <= 8);
        let max_value = min(u64::from(KmerCount::MAX), if byte_len == 8 { u64::MAX } else { (1 << byte_len * 8) - 1 });

        let n_contigs = f.read_u32_varint()?;
        let mut counts = Vec::with_capacity(n_contigs as usize);

        for _ in 0..n_contigs {
            let n_kmers = f.read_u32_varint()?;
            let curr_counts = (0..n_kmers)
                .map(|_| f.read_u64_varint().map(|x| min(x, max_value) as KmerCount))
                .collect::<io::Result<Vec<_>>>()?;
            counts.push(curr_counts);
        }

        Ok(Self {
            k, counts,
            max_value: max_value as KmerCount,
        })
    }

    /// Checks the number of and sizes of contigs.
    pub fn validate(&self, contigs: &ContigNames) -> crate::Result<()> {
        if self.counts.len() != contigs.len() {
            return Err(error!(InvalidData, "k-mer counts contain {} contigs, while there should be {} contigs",
                self.counts.len(), contigs.len()));
        }
        for (i, (counts, len)) in self.counts.iter().zip(contigs.lengths()).enumerate() {
            let expected = (len + 1).saturating_sub(self.k);
            if expected != counts.len() as u32 {
                return Err(error!(InvalidData, "k-mer counts contain {} k-mers for contig {} (expected {})",
                    counts.len(), contigs.get_name(ContigId::new(i)), expected));
            }
        }
        Ok(())
    }

    /// Counts off-target k-mers across all contigs.
    /// This is done by subtracting target k-mers from the whole reference k-mer counts,
    /// and then for each contig we keep k-mer counts as they are for k-mers that don't appear in the target,
    /// while we replace other k-mer counts with off-target counts.
    ///
    /// If one of the counts is equal to `max_count`, keep it as it is.
    pub fn off_target_counts(
        &self,
        seqs: &[Vec<u8>],
        target_seq: &[u8],
        target_counts: &[KmerCount],
        check_negatives: bool,
    ) -> Self
    {
        let k = u8::try_from(self.k).expect("k-mer size is too large");
        assert!(k <= u128::MAX_KMER_SIZE, "Cannot subtract k-mer counts: k is too high ({})", self.k);

        let mut buffer = Vec::new();
        kmers::<u128, CANONICAL>(target_seq, k, &mut buffer);
        assert_eq!(buffer.len(), target_counts.len(), "Unexpected number of subtract k-mers");

        // Value: off-target count.
        let mut off_target_map = HashMap::default();
        // Insert max value for the undefined k-mer (with Ns), just in case.
        off_target_map.insert(u128::UNDEF, self.max_value);
        let mut have_negatives = false;
        for (&kmer, &count) in buffer.iter().zip(target_counts) {
            let val = off_target_map.entry(kmer).or_insert(count);
            if *val != self.max_value {
                have_negatives |= *val == 0;
                // Decrease off-target count by one.
                *val = val.saturating_sub(1);
            }
        }
        if check_negatives && have_negatives {
            log::error!("Reference sequence does not match completely with Jellyfish k-mer counts.\n    \
                Perhaps Jellyfish counts were calculated for another reference?");
        }

        let mut all_new_counts = Vec::with_capacity(seqs.len());
        for (seq, old_counts) in seqs.iter().zip(&self.counts) {
            buffer.clear();
            kmers::<u128, CANONICAL>(seq, k, &mut buffer);
            let new_counts: Vec<_> = buffer.iter().zip(old_counts)
                .map(|(kmer, &old_count)| off_target_map.get(kmer).copied().unwrap_or(old_count))
                .collect();
            assert_eq!(old_counts.len(), new_counts.len(), "Unexpected number of k-mers for one of the sequences");
            all_new_counts.push(new_counts);
        }
        assert_eq!(self.counts.len(), all_new_counts.len(),
            "Number of sequences does not match the number of k-mer counts");
        Self {
            k: self.k,
            max_value: self.max_value,
            counts: all_new_counts,
        }
    }

    /// Create k-mer counts for a subregion.
    /// Only works for a single sequence.
    pub fn subregion(&self, seq_start: u32, seq_end: u32) -> KmerCounts {
        assert_eq!(self.counts.len(), 1, "k-mer subregion counts can only be calculated for a single sequence");
        let full_counts = &self.counts[0];
        let n = full_counts.len() as u32;
        let end = max(seq_start, (seq_end + 1).saturating_sub(self.k));
        assert!(end <= n);
        let sub_counts = full_counts[seq_start as usize..end as usize].to_vec();
        Self {
            k: self.k,
            max_value: self.max_value,
            counts: vec![sub_counts],
        }
    }
}

/// Structure, that runs `jellyfish` and queries k-mers per sequence or from the reference file.
pub struct JfKmerGetter {
    /// Jellyfish executable.
    jf_exe: PathBuf,
    /// Jellyfish database.
    jf_counts: PathBuf,
    /// k-mer size. Extracted from the filename of the jellyfish database.
    k: u32,
    /// Maximum available k-mer count (limited by `KmerCount::MAX` as well).
    max_value: u64,
}

impl JfKmerGetter {
    /// Creates a new jellyfish k-mer getter.
    /// Arguments: `jellyfish` executable and database. Database filename must be in the form `*/<kmer-size>.*`.
    pub fn new(jf_exe: PathBuf, jf_counts: PathBuf) -> crate::Result<Self> {
        let header = parse_jellyfish_header(&jf_counts)?;
        json_get!(&header => canonical (as_bool), key_len (as_u32), counter_len (as_u32));
        if !canonical {
            log::warn!("Jellyfish counts file contains non-canonical k-mer counts. \
                Consider using canonical counts (jellyfish --canonical)");
        }
        if key_len % 2 != 0 {
            return Err(error!(InvalidData, "Cannot parse jellyfish counts file: key length is odd"));
        }
        let k = key_len / 2;
        if k % 2 != 1 {
            return Err(error!(InvalidData,
                "Jellyfish counts are calculated for an even k size ({}), odd k is required", k));
        } else if k > u32::from(u128::MAX_KMER_SIZE) {
            return Err(error!(InvalidData,
                "Jellyfish counts are calculated for too big k ({}), at most {} can be used", k, u128::MAX_KMER_SIZE));
        }
        log::info!("Detected jellyfish k-mer size: {}", k);
        if counter_len > 8 {
            return Err(error!(InvalidData,
                "Jellyfish was run with {}-byte counters, at most 8-byte counters are allowed", counter_len));
        }

        // Calculate maximum k-mer count based on (i) counter length (in bytes) in the counts file,
        // and (ii) maximum value that can be stored in `KmerCount`.
        let max_value = if counter_len * 8 >= KmerCount::BITS {
            KmerCount::MAX as u64
        } else {
            (1 << counter_len * 8) - 1
        };

        Ok(Self { jf_exe, jf_counts, k, max_value })
    }

    /// Returns k-mer size.
    pub fn k(&self) -> u32 {
        self.k
    }

    /// Runs jellyfish to return k-mer counts across all sequences.
    /// Sequences are consumed, as they are passed to another thread.
    pub fn fetch(&self, seqs: impl IntoIterator<Item = Vec<u8>> + Send + 'static) -> crate::Result<KmerCounts> {
        let mut child = Command::new(&self.jf_exe)
            .args(&["query", "-s", "/dev/stdin"]).arg(&self.jf_counts)
            .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped())
            .spawn().map_err(add_path!(self.jf_exe))?;
        let mut child_stdin = io::BufWriter::new(child.stdin.take().unwrap());
        let pipe_guard = PipeGuard::new(self.jf_exe.clone(), child);
        let k_usize = self.k as usize;
        let handle = std::thread::spawn(move || -> crate::Result<Vec<usize>> {
            let mut n_kmers = Vec::new();
            for seq in seqs.into_iter() {
                if seq::has_n(&seq) {
                    return Err(error!(RuntimeError, "Cannot count k-mers for sequence with Ns"));
                }
                n_kmers.push((seq.len() + 1).saturating_sub(k_usize));
                seq::write_fasta(&mut child_stdin, b"", &seq).map_err(add_path!(!))?;
            }
            Ok(n_kmers)
        });

        let output = pipe_guard.wait()?;
        let n_kmers = handle.join().expect("Process failed for unknown reason")?;

        let mut jf_lines = match std::str::from_utf8(&output.stdout) {
            Ok(val) => val.split('\n'),
            Err(_) => return Err(Error::Utf8("Jellyfish output",
                output.stdout[..output.stdout.len().min(500)].to_vec())),
        };
        let mut counts = Vec::with_capacity(n_kmers.len());
        for n in n_kmers.into_iter() {
            let mut curr_counts = Vec::with_capacity(n as usize);
            for _ in 0..n {
                let line = jf_lines.next()
                    .ok_or_else(|| error!(ParsingError, "Not enough k-mer counts!"))?;
                let count = line.split_once(' ')
                    .map(|tup| tup.1)
                    .and_then(|s| parse_count(s, self.max_value).ok())
                    .ok_or_else(||
                        error!(ParsingError, "Failed to parse Jellyfish output line '{}'", line))?;
                // Assume that each k-mer appears at least 1.
                curr_counts.push(max(count, 1));
            }
            counts.push(curr_counts);
        }

        match jf_lines.next() {
            Some("") | None => Ok(KmerCounts {
                counts,
                k: self.k,
                max_value: KmerCount::try_from(self.max_value).unwrap(),
            }),
            _ => Err(error!(InvalidData, "Too many k-mer counts!")),
        }
    }
}

/// Parses Jellyfish counts file header and returns json object.
fn parse_jellyfish_header(filename: &std::path::Path) -> crate::Result<json::JsonValue> {
    let mut reader = io::BufReader::new(fs::File::open(filename).map_err(add_path!(filename))?);
    let mut buffer = Vec::new();
    // Skip until JSON starts.
    if reader.read_until(b'{', &mut buffer).map_err(add_path!(filename))? == 0 {
        return Err(error!(InvalidData, "Cannot parse jellyfish counts file: empty file"));
    }
    buffer.clear();
    buffer.push(b'{');

    let mut prev_len = buffer.len();
    let mut depth = 1;
    while depth > 0 {
        if reader.read_until(b'}', &mut buffer).map_err(add_path!(filename))? == 0 {
            return Err(error!(InvalidData, "Cannot parse jellyfish counts file: Header stops unexpectedly"));
        }
        depth += buffer[prev_len..].iter().filter(|&&ch| ch == b'{').count() - 1;
        prev_len = buffer.len();
    }
    let json_str = std::str::from_utf8(&buffer).map_err(|_| Error::Utf8("Jellyfish header", buffer.clone()))?;
    json::parse(json_str).map_err(|_|
        error!(InvalidData, "Cannot parse jellyfish counts file: header contains invalid JSON"))
}
