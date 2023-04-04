use std::{
    io::{self, Write, BufRead},
    path::PathBuf,
    cmp::{min, max},
    process::{Stdio, Command},
};
use crate::{
    Error,
    seq::{ContigId, ContigNames},
};

/// Store k-mer counts as u16.
pub type KmerCount = u16;

/// Stores k-mer counts for each input k-mer across a set of sequences.
pub struct KmerCounts {
    k: u32,
    counts: Vec<Vec<KmerCount>>,
}

impl KmerCounts {
    /// Creates k-mer counts for each sequence in `seqs`.
    pub fn create<'a>(seqs: impl Iterator<Item = &'a [u8]>, kmer_getter: &JfKmerGetter) -> Result<Self, Error> {
        Ok(Self {
            counts: kmer_getter.fetch_multiple(seqs)?,
            k: kmer_getter.k(),
        })
    }

    /// Load k-mer counts from a file (see `save` for format).
    /// Panics if the number of k-mer counts does not match the contig lengths exactly.
    pub fn load<R: BufRead>(f: &mut R, contigs: &ContigNames) -> Result<Self, Error> {
        assert!(!contigs.is_empty(), "Cannot load k-mer counts for empty contigs set!");
        let mut lines = f.lines();
        let first = lines.next().ok_or_else(|| Error::InvalidData("Empty file with k-mer counts!".to_string()))??;
        let k = match first.split_once('=') {
            Some(("k", v)) => v.parse().map_err(|_| Error::ParsingError(format!("Cannot parse line {}", first)))?,
            _ => return Err(Error::InvalidData(
                format!("Incorrect k-mer counts format, first line ({}) must be in format \"k=integer\"", first))),
        };

        let mut counts = Vec::with_capacity(contigs.len());
        for contig_len in contigs.lengths() {
            let curr_counts = lines.next()
                .ok_or_else(|| Error::InvalidData(
                    "File with k-mer counts does not contain enough contigs".to_string()))??
                .split(' ')
                .map(str::parse)
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| Error::ParsingError(format!("Cannot parse k-mer counts: {}", e)))?;

            if curr_counts.len() as u32 != contig_len - k + 1 {
                return Err(Error::InvalidData("Incorrect number of k-mers counts".to_string()));
            }
            counts.push(curr_counts);
        }
        match lines.next() {
            Some(Err(e)) => Err(e)?,
            Some(Ok(s)) if s.is_empty() => Ok(Self { k, counts }),
            None => Ok(Self { k, counts }),
            _ => Err(Error::InvalidData("Too many k-mer counts!".to_string())),
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
}

/// k-mer that contains Ns.
pub const N_KMER: u64 = u64::MAX;

/// Iterates over canonical k-mers in the sequence.
/// Returns N_KMER for k-mers containing N.
/// K-mer size `k` must be at most 32.
///
/// Canonical k-mer is calculated as minimal between the k-mer and its rev.complement.
pub fn canonical_kmers(seq: &[u8], k: u8) -> Vec<u64> {
    assert!(0 < k && k <= 32, "k-mer size must be within [1, 32].");
    let k_usize = usize::from(k);
    assert!(seq.len() >= k_usize, "Sequence is too short!");
    let mask: u64 = if k == 32 { -1_i64 as u64 } else { (1_u64 << 2 * k) - 1 };
    let rv_shift = 2 * k - 2;
    let mut fw_kmer: u64 = 0;
    let mut rv_kmer: u64 = 0;
    let mut reset = k_usize - 1;
    let mut kmers = Vec::with_capacity(seq.len() - k_usize + 1);

    for (i, &nt) in seq.iter().enumerate() {
        let fw_enc: u64 = match nt {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                reset = i + k_usize;
                if i + 1 >= k_usize {
                    kmers.push(N_KMER);
                }
                continue;
            },
        };
        let rv_enc = 3 - fw_enc;
        fw_kmer = (fw_kmer << 2) | fw_enc;
        rv_kmer = (rv_enc << rv_shift) | (rv_kmer >> 2);

        if i >= reset {
            kmers.push(min(fw_kmer & mask, rv_kmer));
        } else if i + 1 >= k_usize {
            kmers.push(N_KMER);
        }
    }
    debug_assert_eq!(kmers.len(), seq.len() - k_usize + 1);
    kmers
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
                Error::ParsingError(format!("Cannot parse jellyfish database filename '{}'", jf_db.display())))?;
        assert!(k % 2 == 1, "k-mer size ({}) must be odd!", k);
        Ok(Self { jf_exe, jf_db, k })
    }

    /// Returns k-mer size.
    pub fn k(&self) -> u32 {
        self.k
    }

    /// Private function, that provides all sequences to Jellyfish, and returns a single vector of k-mer counts.
    /// k-mers are then split into different sequences in public functions.
    /// Second returned element: number of k-mers in each sequence.
    fn fetch_inner<'a>(&self, seqs: impl Iterator<Item = &'a [u8]>) -> Result<(Vec<KmerCount>, Vec<usize>), Error> {
        let mut child = Command::new(&self.jf_exe)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .args(&["query", "-s", "/dev/stdin"])
            .arg(&self.jf_db)
            .spawn()?;

        let mut n_kmers = Vec::new();
        // unwrap as stdin was specified above.
        let mut child_stdin = child.stdin.take().unwrap();
        // let mut child_stdin = child.stdin.as_ref().unwrap();
        log::info!("Writing to child");
        for seq in seqs {
            log::debug!("    Sequence of length {}", seq.len());
            child_stdin.write_all(b">\n")?;
            child_stdin.write_all(&seq)?;
            child_stdin.write_all(b"\n")?;
            n_kmers.push((seq.len() + 1).saturating_sub(self.k as usize));
        }
        std::mem::drop(child_stdin);

        log::info!("Finished writing, {:?} kmers", n_kmers);
        let output = child.wait_with_output()?;
        log::info!("Finished waiting");
        if !output.status.success() {
            return Err(Error::SubcommandFail(output));
        }

        let jf_out = std::str::from_utf8(&output.stdout)?;
        let exp_size = n_kmers.iter().cloned().sum();
        let mut counts: Vec<KmerCount> = Vec::with_capacity(exp_size);
        for line in jf_out.split('\n') {
            if line.is_empty() {
                break;
            }
            let count = line.split_once(' ')
                .map(|tup| tup.1)
                .and_then(|s| s.parse().ok())
                .ok_or_else(||
                    Error::ParsingError(format!("Failed to parse Jellyfish output line '{}'", line)))?;
            // Assume that each k-mer appears at least 1.
            counts.push(max(count, 1));
        }

        if counts.len() != exp_size {
            Err(Error::RuntimeError(format!("Failed to run `jellyfish query`. \
                Expected {} k-mer counts, found {}", exp_size, counts.len())))
        } else {
            Ok((counts, n_kmers))
        }
    }

    /// Returns all k-mers in the sequence.
    pub fn fetch_one(&self, seq: &[u8]) -> Result<Vec<KmerCount>, Error> {
        self.fetch_inner(std::iter::once(seq)).map(|(counts, _)| counts)
    }

    /// Returns all k-mers in the multipe sequences.
    pub fn fetch_multiple<'a>(&self, seqs: impl Iterator<Item = &'a [u8]>) -> Result<Vec<Vec<KmerCount>>, Error> {
        let mut res = if let Some(capacity) = seqs.size_hint().1 {
            Vec::with_capacity(capacity)
        } else {
            Vec::new()
        };
        let (counts, n_kmers) = self.fetch_inner(seqs)?;
        let mut slice = &counts[..];
        for n in n_kmers.into_iter() {
            res.push(slice[..n].to_vec());
            slice = &slice[n..];
        }
        Ok(res)
    }
}