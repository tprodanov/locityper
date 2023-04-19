use std::{
    io::{self, Write, BufRead},
    path::PathBuf,
    cmp::{min, max},
    process::{Stdio, Command},
};
use crate::{
    Error,
    seq::{self, ContigId},
};

/// Store k-mer counts as u16.
pub type KmerCount = u16;

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
        let first = lines.next().ok_or_else(|| Error::InvalidData("Empty file with k-mer counts!".to_owned()))??;
        let k = match first.split_once('=') {
            Some(("k", v)) => v.parse().map_err(|_| Error::ParsingError(format!("Cannot parse line {}", first)))?,
            _ => return Err(Error::InvalidData(
                format!("Incorrect k-mer counts format, first line ({}) must be in format \"k=integer\"", first))),
        };

        let mut counts = Vec::with_capacity(contig_lengths.len());
        for &contig_len in contig_lengths.iter() {
            let curr_counts = lines.next()
                .ok_or_else(|| Error::InvalidData(
                    "File with k-mer counts does not contain enough contigs".to_owned()))??
                .split(' ')
                .map(str::parse)
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| Error::ParsingError(format!("Cannot parse k-mer counts: {}", e)))?;

            if curr_counts.len() as u32 != contig_len - k + 1 {
                return Err(Error::InvalidData("Incorrect number of k-mers counts".to_owned()));
            }
            counts.push(curr_counts);
        }
        match lines.next() {
            Some(Err(e)) => Err(e)?,
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

/// k-mer that contains Ns.
pub const N_KMER: u64 = u64::MAX;

/// Maximal allowed k-mer size.
pub const MAX_KMER: u8 = 31;

/// Iterates over canonical k-mers in the sequence.
/// Returns N_KMER for k-mers containing N.
/// K-mer size `k` must be at most 31.
///
/// Canonical k-mer is calculated as minimal between the k-mer and its rev.complement.
pub fn canonical_kmers(seq: &[u8], k: u8, buffer: &mut Vec<u64>) {
    assert!(0 < k && k <= MAX_KMER, "k-mer size must be within [1, {}].", MAX_KMER);
    let k_usize = usize::from(k);
    assert!(seq.len() >= k_usize, "Sequence is too short!");
    let mask: u64 = if k == 32 { -1_i64 as u64 } else { (1_u64 << 2 * k) - 1 };
    let rv_shift = 2 * k - 2;
    let mut fw_kmer: u64 = 0;
    let mut rv_kmer: u64 = 0;
    let mut reset = k_usize - 1;

    for (i, &nt) in seq.iter().enumerate() {
        let fw_enc: u64 = match nt {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                reset = i + k_usize;
                if i + 1 >= k_usize {
                    buffer.push(N_KMER);
                }
                continue;
            },
        };
        fw_kmer = (fw_kmer << 2) | fw_enc;
        rv_kmer = ((3 - fw_enc) << rv_shift) | (rv_kmer >> 2);

        if i >= reset {
            buffer.push(min(fw_kmer & mask, rv_kmer));
        } else if i + 1 >= k_usize {
            buffer.push(N_KMER);
        }
    }
}

// fn murmur3_hash(mut key: u64) -> u64 {
//     key ^= key >> 33;
//     key = key.wrapping_mul(0xff51afd7ed558ccd);
//     key ^= key >> 33;
//     key = key.wrapping_mul(0xc4ceb9fe1a85ec53);
//     key ^= key >> 33;
//     key
// }

// /// Fast u32 hash, taken from [here](https://github.com/skeeto/hash-prospector).
// fn u32_hash(mut key: u32) -> u32 {
//     key ^= key >> 16;
//     key = key.wrapping_mul(0x7feb352d);
//     key ^= key >> 15;
//     key = key.wrapping_mul(0x846ca68b);
//     key ^= key >> 16;
//     key
// }

/// Fast u32 hash, adapted from [here](https://github.com/skeeto/hash-prospector).
/// Additionally, u32 `FNV_OFFSET` is used to scramble input value.
fn kmer_hash(mut key: u32) -> u32 {
    const FNV_OFFSET: u32 = 0x811c9dc5;
    key ^= FNV_OFFSET;
    key ^= key >> 16;
    key = key.wrapping_mul(0x7feb352d);
    key ^= key >> 15;
    // key = key.wrapping_mul(0x846ca68b);
    // key ^= key >> 16;
    key
}

/// Finds sequence minimizers.
/// Minimizer is a k-mer with the smallest hash value across `n` consecutive k-mers.
/// Usually, letter `w` is used: in that case, minimizer is selected per each rolling `w`-window.
/// Here, `w = n + k - 1`.
///
/// `k` must be at most 16. `n` must be at most 64.
pub fn minimizers(seq: &[u8], k: u8, n: u8, buffer: &mut Vec<u32>) {
    let k = u32::from(k);
    let n = u32::from(n);
    const MAXN: usize = 64;
    const MOD_MAXN: u32 = MAXN as u32 - 1; // Must be the power of two.
    const UNDEF: u32 = u32::MAX;
    debug_assert!(0 < k && k <= 16, "k-size must be within [1, 16]");
    debug_assert!(1 < n && n <= MAXN as u32, "n must be within [2, {}]", MAXN);

    let mask: u32 = if k == 16 { -1_i32 as u32 } else { (1_u32 << 2 * k) - 1 };
    let rv_shift = 2 * k - 2;
    let mut fw_kmer: u32 = 0;
    let mut rv_kmer: u32 = 0;

    // Hashes in a window, stored in a cycling array.
    let mut hashes = [UNDEF; MAXN];
    /// At what index will the first k-mer be available.
    let mut reset = k - 1;
    // Start of the window with consecutive k-mers.
    let mut start = reset;

    /// Function that goes over indices `start..end`, and returns new `start`.
    /// Additionally, the function pushes the new minimizer to the buffer, if it is not `UNDEF`.
    #[inline]
    fn select_minimizer(buffer: &mut Vec<u32>, hashes: &mut [u32; MAXN], start: u32, end: u32) -> u32 {
        let mut minimizer = UNDEF;
        let mut new_start = end;
        for j in start..end {
            let h = hashes[(j & MOD_MAXN) as usize];
            if h < minimizer {
                minimizer = h;
                new_start = j + 1;
            }
        }
        if minimizer != UNDEF {
            buffer.push(minimizer);
        }
        new_start
    }

    for (i, &nt) in seq.iter().enumerate() {
        let i = i as u32;
        let (fw_enc, rv_enc): (u32, u32) = match nt {
            b'A' => (0, 3),
            b'C' => (1, 2),
            b'G' => (2, 1),
            b'T' => (3, 0),
            _ => {
                reset = i + k;
                if i > start && reset >= start + n {
                    select_minimizer(buffer, &mut hashes, start, i);
                    start = reset;
                }
                hashes[(i & MOD_MAXN) as usize] = UNDEF;
                continue;
            },
        };
        fw_kmer = (fw_kmer << 2) | fw_enc;
        rv_kmer = (rv_enc << rv_shift) | (rv_kmer >> 2);
        if i < reset {
            hashes[(i & MOD_MAXN) as usize] = UNDEF;
            continue;
        }

        hashes[(i & MOD_MAXN) as usize] = kmer_hash(min(fw_kmer & mask, rv_kmer));
        if i == start + n - 1 {
            start = select_minimizer(buffer, &mut hashes, start, i + 1);
        }
    }
    let l = seq.len() as u32;
    if l >= start {
        debug_assert!(l <= start + n - 1);
        select_minimizer(buffer, &mut hashes, start, l);
    }
}

/// Structure, that runs `jellyfish` and queries k-mers per sequence or from the reference file.
pub struct JfKmerGetter {
    /// Jellyfish executable.
    jf_exe: PathBuf,
    /// Jellyfish database.
    jf_db: PathBuf,
    /// k-mer size. Extracted from the filename of the jellyfish database.
    k: u32,
}

impl JfKmerGetter {
    /// Creates a new jellyfish k-mer getter.
    /// Arguments: `jellyfish` executable and database. Database filename must be in the form `*/<kmer-size>.*`.
    pub fn new(jf_exe: PathBuf, jf_db: PathBuf) -> Result<Self, Error> {
        let k = jf_db.file_stem()
            .and_then(std::ffi::OsStr::to_str)
            .and_then(|s| s.parse().ok())
            .ok_or_else(||
                Error::ParsingError(format!("Cannot parse jellyfish database filename {:?}", jf_db)))?;
        assert!(k % 2 == 1, "k-mer size ({}) must be odd!", k);
        Ok(Self { jf_exe, jf_db, k })
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
            .args(&["query", "-s", "/dev/stdin"]).arg(&self.jf_db)
            .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped())
            .spawn()?;
        let mut child_stdin = io::BufWriter::new(child.stdin.take().unwrap());
        let handle = std::thread::spawn(move || -> Result<(), Error> {
            for seq in seqs.iter() {
                seq::write_fasta(&mut child_stdin, b"", None, seq)?;
            }
            Ok(())
        });

        let output = child.wait_with_output()?;
        if !output.status.success() {
            return Err(Error::SubprocessFail(output));
        }
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
                    .and_then(|s| s.parse().ok())
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