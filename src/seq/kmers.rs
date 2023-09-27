use std::{
    fmt::Display,
    io::{self, Write, BufRead},
    fs,
    path::PathBuf,
    cmp::{min, max, Ord},
    process::{Stdio, Command},
    ops::{Shl, Shr, BitOr, BitAnd},
};
use crate::{
    err::{Error, add_path},
    bg::ser::json_get,
    seq::{self, ContigId},
    ext::sys::PipeGuard,
};

/// Store k-mers in integers of different length (15-mers in u32, 31-mers in u64, 63-mers in u128).
pub trait Kmer: From<u8>
    + Copy + Ord + Eq + Display
    + Shl<u8, Output = Self>
    + Shr<u8, Output = Self>
    + BitOr<Self, Output = Self>
    + BitAnd<Self, Output = Self>
{
    /// Zero value.
    const ZERO: Self;

    /// Maximum k-mer size `(BITS / 2 - 1)`.
    /// It is actually possible to store `BITS / 2` k-mers in case of minimizers and canonical k-mers, but not in case
    /// of non-canonical k-mers.
    /// To be on the safe side, we forbid this k-mer size.
    const MAX_KMER_SIZE: u8;

    /// Undefined k-mer/minimizer.
    /// In case of minimizers, it is possible that hash function would produce `UNDEF` value inadvertently.
    /// Nevertheless, such minimizers would have the highest hash function anyway, and should not go to the output.
    const UNDEF: Self;

    fn create_mask(k: u8) -> Self;
}

macro_rules! impl_kmer {
    ($prim:ty) => {
        impl Kmer for $prim {
            const ZERO: Self = 0;

            const MAX_KMER_SIZE: u8 = (Self::BITS / 2 - 1) as u8;

            const UNDEF: Self = Self::MAX;

            #[inline]
            fn create_mask(k: u8) -> Self {
                // if k == 16 { -1_i32 as u32 } else { (1_u32 << 2 * k) - 1 }
                // Already know that k is not `Self::BITS / 2`.
                (1 << 2 * k) - 1
            }
        }
    }
}

impl_kmer!(u8);
impl_kmer!(u16);
impl_kmer!(u32);
impl_kmer!(u64);
impl_kmer!(u128);

pub trait Minimizer: Kmer {
    /// Quickly calculate hash function.
    /// For now, use Murmur3 function, this can change later.
    fn fast_hash(self) -> Self;
}

impl Minimizer for u32 {
    /// Murmur3 for uint32.
    #[inline]
    fn fast_hash(mut self) -> u32 {
        self ^= self >> 16;
        self = self.wrapping_mul(0x85ebca6b);
        self ^= self >> 13;
        // self = self.wrapping_mul(0xc2b2ae35);
        // self ^= self >> 16;
        // 0 will always produce 0. Here, we take bit-inversion, so that 0 will produce u32::MAX = UNDEF.
        // This corresponds to AA..AA k-mer, which is a very bad minimizer.
        !self
    }
}

impl Minimizer for u64 {
    /// fasthash (https://github.com/rurban/smhasher/blob/master/fasthash.cpp) mix function.
    #[inline]
    fn fast_hash(mut self) -> u64 {
        self ^= self >> 23;
        self = self.wrapping_mul(0x2127599bf4325c37);
        self ^= self >> 47;
        !self
    }
}

impl Minimizer for u128 {
    #[inline]
    fn fast_hash(self) -> Self {
        let [h1, h2] = unsafe { std::mem::transmute::<u128, [u64; 2]>(self) };
        let h1 = h1.fast_hash();
        let h2 = h2.fast_hash();
        unsafe { std::mem::transmute::<[u64; 2], u128>([h1, h2]) }
    }
}

pub const CANONICAL: bool = true;
pub const NON_CANONICAL: bool = false;

/// Returns all k-mers in the sequence.
/// Returns `K::UNDEF` for k-mers containing N.
///
/// If `CANONICAL` is true, returns canonical k-mers (min between forward and reverse-complement k-mer).
pub fn kmers<K: Kmer, const CANON: bool>(seq: &[u8], k: u8, buffer: &mut Vec<K>) {
    debug_assert!(k <= K::MAX_KMER_SIZE, "k-mer size ({}) can be at most {}", k, K::MAX_KMER_SIZE);
    let mask = K::create_mask(k);
    let rv_shift = 2 * k - 2;
    let mut fw_kmer = K::ZERO;
    let mut rv_kmer = K::ZERO;

    let k = u32::from(k);
    let mut reset = k - 1;

    for (i, &nt) in seq.iter().enumerate() {
        let i = i as u32;
        let fw_enc: u8 = match nt {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                reset = i + k;
                if i + 1 >= k {
                    buffer.push(K::UNDEF);
                }
                continue;
            },
        };
        fw_kmer = (fw_kmer << 2) | K::from(fw_enc);
        if CANON {
            rv_kmer = (K::from(3 - fw_enc) << rv_shift) | (rv_kmer >> 2);
        }

        if i >= reset {
            if CANON {
                buffer.push(min(fw_kmer & mask, rv_kmer));
            } else {
                buffer.push(fw_kmer & mask);
            }
        } else if i + 1 >= k {
            buffer.push(K::UNDEF);
        }
    }
}

/// Biggest allowed window size (number of consecutive k-mer). Must be the power of two.
pub const MAX_MINIMIZER_W: u8 = 64;
const MOD_MAXW: u32 = MAX_MINIMIZER_W as u32 - 1;

/// Function that goes over indices `start..end`, and returns new `start`.
/// Additionally, the function pushes the new minimizer to the buffer, if it is not `UNDEF`.
#[inline]
fn select_minimizer<K: Kmer>(
    buffer: &mut Vec<K>,
    hashes: &mut [K; MAX_MINIMIZER_W as usize],
    start: u32, end: u32,
) -> u32 {
    let mut minimizer = K::UNDEF;
    let mut new_start = end;
    for j in start..end {
        let h = hashes[(j & MOD_MAXW) as usize];
        if h < minimizer {
            minimizer = h;
            new_start = j + 1;
        }
    }
    if minimizer != K::UNDEF {
        buffer.push(minimizer);
    }
    new_start
}

/// Finds sequence minimizers.
/// Minimizer is a k-mer with the smallest hash value across `w` consecutive k-mers.
pub fn minimizers<K: Minimizer>(seq: &[u8], k: u8, w: u8, buffer: &mut Vec<K>) {
    debug_assert!(k <= K::MAX_KMER_SIZE, "k-mer size ({}) can be at most {}", k, K::MAX_KMER_SIZE);
    debug_assert!(w <= MAX_MINIMIZER_W, "Minimizer window ({}) can be at most {}", w, MAX_MINIMIZER_W);
    let mask = K::create_mask(k);
    let rv_shift = 2 * k - 2;
    // Hashes in a window, stored in a cycling array.
    let mut hashes = [K::UNDEF; MAX_MINIMIZER_W as usize];

    let mut fw_kmer = K::ZERO;
    let mut rv_kmer = K::ZERO;

    let k = u32::from(k);
    let w = u32::from(w);
    // At what index will the first k-mer be available.
    let mut reset = k - 1;
    // Start of the window with consecutive k-mers.
    let mut start = reset;

    for (i, &nt) in seq.iter().enumerate() {
        let i = i as u32;
        let (fw_enc, rv_enc): (u8, u8) = match nt {
            b'A' => (0, 3),
            b'C' => (1, 2),
            b'G' => (2, 1),
            b'T' => (3, 0),
            _ => {
                reset = i + k;
                if reset >= start + w {
                    if i > start {
                        select_minimizer(buffer, &mut hashes, start, i);
                    }
                    start = reset;
                }
                hashes[(i & MOD_MAXW) as usize] = K::UNDEF;
                continue;
            },
        };
        fw_kmer = (fw_kmer << 2) | K::from(fw_enc);
        rv_kmer = (K::from(rv_enc) << rv_shift) | (rv_kmer >> 2);
        if i < reset {
            hashes[(i & MOD_MAXW) as usize] = K::UNDEF;
            continue;
        }

        hashes[(i & MOD_MAXW) as usize] = min(fw_kmer & mask, rv_kmer).fast_hash();
        if i == start + w - 1 {
            start = select_minimizer(buffer, &mut hashes, start, i + 1);
        }
    }
    let l = seq.len() as u32;
    if l >= start {
        debug_assert!(l <= start + w - 1);
        select_minimizer(buffer, &mut hashes, start, l);
    }
}

/// Store k-mer counts as u16.
pub type KmerCount = u8;

#[inline]
fn parse_count(s: &str) -> Result<KmerCount, std::num::ParseIntError> {
    s.parse::<u64>().map(|x| KmerCount::try_from(x).unwrap_or(KmerCount::MAX))
}

/// Stores k-mer counts for each input k-mer across a set of sequences.
pub struct KmerCounts {
    k: u32,
    counts: Vec<Vec<KmerCount>>,
}

impl KmerCounts {
    /// Load k-mer counts from a file (see `save` for format).
    /// Panics if the number of k-mer counts does not match the contig lengths exactly.
    pub fn load<R: BufRead>(f: R, contig_lengths: &[u32]) -> Result<Self, Error> {
        assert!(!contig_lengths.is_empty(), "Cannot load k-mer counts for empty contigs set!");
        let mut lines = f.lines();
        let first = lines.next().ok_or_else(|| Error::InvalidData("Empty file with k-mer counts!".to_owned()))?
            .map_err(add_path!(!))?;
        let k = match first.split_once('=') {
            Some(("k", v)) => v.parse().map_err(|_| Error::ParsingError(format!("Cannot parse line {}", first)))?,
            _ => return Err(Error::InvalidData(
                format!("Incorrect k-mer counts format, first line ({}) must be in format \"k=integer\"", first))),
        };

        let mut counts = Vec::with_capacity(contig_lengths.len());
        for &contig_len in contig_lengths.iter() {
            let curr_counts = lines.next()
                .ok_or_else(|| Error::InvalidData(
                    "File with k-mer counts does not contain enough contigs".to_owned()))?.map_err(add_path!(!))?
                .split(' ')
                .map(parse_count)
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| Error::ParsingError(format!("Cannot parse k-mer counts: {}", e)))?;

            if curr_counts.len() as u32 != contig_len - k + 1 {
                return Err(Error::InvalidData("Incorrect number of k-mers counts".to_owned()));
            }
            counts.push(curr_counts);
        }
        match lines.next() {
            Some(Err(e)) => Err(Error::Io(e, vec![]))?,
            Some(Ok(s)) if s.is_empty() => Ok(Self { k, counts }),
            None => Ok(Self { k, counts }),
            _ => Err(Error::InvalidData("Too many k-mer counts!".to_owned())),
        }
    }

    /// Writes k-mer counts into file in the following format (for example), each line - new contig:
    /// ```
    /// k=25
    /// 0 0 0 10 24 35 23 9 0 0 0
    /// 0 0 0 0 0 0
    /// ```
    pub fn save<W: Write>(&self, f: &mut W) -> io::Result<()> {
        writeln!(f, "k={}", self.k)?;
        for counts in self.counts.iter() {
            if !counts.is_empty() {
                write!(f, "{}", counts[0])?;
                for count in &counts[1..] {
                    write!(f, " {}", count)?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }

    /// Returns k-mer size.
    pub fn k(&self) -> u32 {
        self.k
    }

    /// Returns all k-mer counts for the corresponding contig.
    /// Length: contig len - k + 1.
    pub fn get(&self, contig_id: ContigId) -> &[KmerCount] {
        &self.counts[contig_id.ix()]
    }

    /// Returns all k-mer counts for the first contig (must be the only contig).
    /// Length: contig len - k + 1.
    pub fn get_first(&self) -> &[KmerCount] {
        assert_eq!(self.counts.len(), 1,
            "Cannot call KmerCounts.get_first when the number of contigs is different than one");
        &self.counts[0]
    }

    /// Returns all k-mer counts for the first contig (must be the only contig).
    /// Length: contig len - k + 1.
    pub fn take_first(mut self) -> Vec<KmerCount> {
        assert_eq!(self.counts.len(), 1,
            "Cannot call KmerCounts.take_first when the number of contigs is different than one");
        self.counts.pop().unwrap()
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
}

impl JfKmerGetter {
    /// Creates a new jellyfish k-mer getter.
    /// Arguments: `jellyfish` executable and database. Database filename must be in the form `*/<kmer-size>.*`.
    pub fn new(jf_exe: PathBuf, jf_counts: PathBuf) -> Result<Self, Error> {
        let k = parse_jellyfish_header(&jf_counts)?;
        if k % 2 != 1 {
            return Err(Error::InvalidData(
                format!("Jellyfish counts are calculated for an even k size ({}), odd k is required", k)));
        }
        log::debug!("    Jellyfish k-mer size: {}", k);
        Ok(Self { jf_exe, jf_counts, k })
    }

    /// Returns k-mer size.
    pub fn k(&self) -> u32 {
        self.k
    }

    /// Runs jellyfish to return k-mer counts across all sequences.
    /// Sequences are consumed, as they are passed to another thread.
    pub fn fetch(&self, seqs: Vec<Vec<u8>>) -> Result<KmerCounts, Error> {
        let n_kmers: Vec<_> = seqs.iter()
            .map(|seq| (seq.len() + 1).saturating_sub(self.k as usize)).collect();
        let mut child = Command::new(&self.jf_exe)
            .args(&["query", "-s", "/dev/stdin"]).arg(&self.jf_counts)
            .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped())
            .spawn().map_err(add_path!(self.jf_exe))?;
        let mut child_stdin = io::BufWriter::new(child.stdin.take().unwrap());
        let pipe_guard = PipeGuard::new(self.jf_exe.clone(), child);
        let handle = std::thread::spawn(move || -> Result<(), Error> {
            for seq in seqs.iter() {
                seq::write_fasta(&mut child_stdin, b"", seq).map_err(add_path!(!))?;
            }
            Ok(())
        });

        let output = pipe_guard.wait()?;
        // handle.join() returns Result<Result<(), crate::Error>, Any>.
        // expect unwraps the outer Err, then ? returns the inner Err, if any.
        handle.join().expect("Process failed for unknown reason")?;

        let mut jf_lines = std::str::from_utf8(&output.stdout)?.split('\n');
        let mut counts = Vec::with_capacity(n_kmers.len());
        for n in n_kmers.into_iter() {
            let mut curr_counts = Vec::with_capacity(n as usize);
            for _ in 0..n {
                let line = jf_lines.next()
                    .ok_or_else(|| Error::ParsingError("Not enough k-mer counts!".to_owned()))?;
                let count = line.split_once(' ')
                    .map(|tup| tup.1)
                    .and_then(|s| parse_count(s).ok())
                    .ok_or_else(||
                        Error::ParsingError(format!("Failed to parse Jellyfish output line '{}'", line)))?;
                // Assume that each k-mer appears at least 1.
                curr_counts.push(max(count, 1));
            }
            counts.push(curr_counts);
        }

        match jf_lines.next() {
            Some("") | None => Ok(KmerCounts {
                counts,
                k: self.k,
            }),
            _ => Err(Error::InvalidData("Too many k-mer counts!".to_owned())),
        }
    }

    /// Runs jellyfish and returns k-mer frequencies for all k-mers in the sequence.
    /// Sequence is consumed, as it needs to be passed to another thread.
    pub fn fetch_one(&self, seq: Vec<u8>) -> Result<KmerCounts, Error> {
        self.fetch(vec![seq])
    }
}

/// Based on the jellyfish counts file, extracts k-mer size, and checks if k-mers are canonical.
fn parse_jellyfish_header(filename: &std::path::Path) -> Result<u32, Error> {
    let mut reader = io::BufReader::new(fs::File::open(filename).map_err(add_path!(filename))?);
    let mut buffer = Vec::new();
    // Skip until JSON starts.
    if reader.read_until(b'{', &mut buffer).map_err(add_path!(filename))? == 0 {
        return Err(Error::InvalidData("Cannot parse jellyfish counts file: empty file".to_owned()));
    }
    buffer.clear();
    buffer.push(b'{');

    let mut prev_len = buffer.len();
    let mut depth = 1;
    while depth > 0 {
        if reader.read_until(b'}', &mut buffer).map_err(add_path!(filename))? == 0 {
            return Err(Error::InvalidData("Cannot parse jellyfish counts file: Header stops unexpectedly".to_owned()));
        }
        depth += buffer[prev_len..].iter().filter(|&&ch| ch == b'{').count() - 1;
        prev_len = buffer.len();
    }
    let json_str = std::str::from_utf8(&buffer).map_err(|_|
        Error::InvalidData("Cannot parse jellyfish counts file: Header is not UTF-8".to_owned()))?;
    let obj = json::parse(json_str).map_err(|_|
        Error::InvalidData("Cannot parse jellyfish counts file: header contains invalid JSON".to_owned()))?;
    json_get!(&obj => canonical (as_bool), key_len (as_u32));
    if !canonical {
        log::warn!("Jellyfish counts file contains non-canonical k-mer counts. \
            Consider using canonical counts (jellyfish --canonical)");
    }
    if key_len % 2 != 0 {
        return Err(Error::InvalidData("Cannot parse jellyfish counts file: key length is odd".to_owned()));
    }
    Ok(key_len / 2)
}
