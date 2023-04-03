use std::{
    io::{self, ErrorKind, Write},
    path::PathBuf,
    cmp::{min, max},
    process::{Stdio, Command},
};
use crate::{
    Error,
    seq::{Interval, ContigId, ContigNames},
};

/// Store k-mer counts as u16.
pub type KmerCount = u16;

/// Stores k-mer counts for each input k-mer.
pub struct KmerCounts {
    k: u32,
    counts: Vec<Vec<KmerCount>>,
}

impl KmerCounts {
    /// Load k-mer counts from a string.
    /// First line is "k=<number>". All consecutive lines contain just a single number.
    /// Must contain exact number of k-mer as `contigs`.
    pub fn load(full_contents: &str, contigs: &ContigNames) -> io::Result<Self> {
        assert!(!contigs.is_empty(), "Cannot load k-mer counts for empty contigs set!");
        let mut split = full_contents.split('\n');
        let first = split.next().ok_or_else(|| io::Error::new(ErrorKind::InvalidData,
            "Empty file with k-mer counts!"))?;
        let k = match first.split_once('=') {
            Some(("k", v)) => v.parse().map_err(|e: std::num::ParseIntError|
                io::Error::new(ErrorKind::InvalidData, e))?,
            _ => return Err(io::Error::new(ErrorKind::InvalidData,
                format!("Incorrect k-mer counts format, first line ({}) must be in format \"k=integer\"", first))),
        };

        let mut counts = Vec::with_capacity(contigs.len());
        for contig_len in contigs.lengths() {
            assert!(contig_len >= k, "Contig too short!");
            let n_kmers = contig_len - k + 1;
            let mut curr_counts = Vec::with_capacity(n_kmers as usize);
            for _ in 0..n_kmers {
                let val: KmerCount = split.next()
                    .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Not enough k-mer counts!"))?
                    .parse().map_err(|e: std::num::ParseIntError| io::Error::new(ErrorKind::InvalidData, e))?;
                // We assume that each k-mer appears at least once.
                curr_counts.push(max(val, 1));
            }
            counts.push(curr_counts);
        }
        match split.next() {
            Some("") | None => Ok(Self { k, counts }),
            _ => Err(io::Error::new(ErrorKind::InvalidData, "Too many k-mer counts!")),
        }
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
        Ok(Self { jf_exe, jf_db, k })
    }

    /// Returns k-mer size.
    pub fn k(&self) -> u32 {
        self.k
    }

    /// Returns all k-mers in the sequence.
    pub fn fetch(&self, seq: &[u8]) -> Result<Vec<KmerCount>, Error> {
        let mut child = Command::new(&self.jf_exe)
            .stdin(Stdio::piped())
            .args(&["query", "-s", "-"])
            .arg(&self.jf_db)
            .spawn()?;

        {
            // unwrap as stdin was specified above.
            let mut child_stdin = child.stdin.as_ref().unwrap();
            child_stdin.write_all(b">\n")?;
            child_stdin.write_all(&seq)?;
        }
        let output = child.wait_with_output()?;
        if !output.status.success() {
            return Err(Error::SubcommandFail(output));
        }

        let jf_out = std::str::from_utf8(&output.stdout)?;
        let exp_size = (seq.len() + 1).saturating_sub(self.k as usize);
        let mut counts: Vec<KmerCount> = Vec::with_capacity(exp_size);
        for line in jf_out.split('\n') {
            let count = line.split_once(' ')
                .map(|tup| tup.1)
                .and_then(|s| s.parse().ok())
                .ok_or_else(||
                    Error::ParsingError(format!("Failed to parse Jellyfish output line '{}'", line)))?;
            counts.push(count);
        }

        if counts.len() != exp_size {
            Err(Error::RuntimeError(format!("Failed to run jellyfish query on sequence of length {}. \
                Expected {} k-mer counts, found {}", seq.len(), exp_size, counts.len())))
        } else {
            Ok(counts)
        }
    }
}