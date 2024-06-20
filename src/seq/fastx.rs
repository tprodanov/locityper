use std::{
    fmt,
    ffi::OsStr,
    borrow::Cow,
    io::{self, BufRead},
    cmp::min,
    path::{Path, PathBuf},
    time::{Instant, Duration},
    process::{Command, Stdio},
};
use htslib::bam;
use crate::{
    ext,
    seq::{self, Interval, ContigNames},
    algo::HashSet,
    err::{Error, add_path},
};
use super::NamedSeq;

// ------------------------------------------- Helper functions -------------------------------------------

/// Write a single sequence to the FASTA file.
/// Use this function instead of `bio::fasta::Writer` as the latter
/// writes the sequence into a single line, without splitting.
///
/// `name` may include description after a space, if needed.
pub fn write_fasta(
    writer: &mut impl io::Write,
    name: &[u8],
    seq: &[u8]
) -> io::Result<()> {
    writer.write_all(b">")?;
    writer.write_all(name)?;
    writer.write_all(b"\n")?;

    const WIDTH: usize = 120;
    let n = seq.len();
    for i in (0..n).step_by(WIDTH) {
        writer.write_all(&seq[i..min(i + WIDTH, n)])?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

/// Write a single sequence to the FASTA file.
/// Use this function instead of `bio::fasta::Writer` as the latter
/// writes the sequence into a single line, without splitting.
///
/// `name` may include description after a space, if needed.
pub fn write_fastq(
    writer: &mut impl io::Write,
    name: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> io::Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(name)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")
}

/// Reads the next line from a buffered reader, but removes endline characters from the end.
/// If there were no errors, returns the number of written characters.
fn read_line(stream: &mut impl BufRead, buffer: &mut Vec<u8>) -> io::Result<usize> {
    match stream.read_until(b'\n', buffer) {
        Ok(0) => Ok(0),
        Ok(n) => {
            let l = buffer.len();
            if n >= 1 && buffer[l - 1] == b'\n' {
                buffer.pop();
                if n >= 2 && buffer[l - 2] == b'\r' {
                    buffer.pop();
                    Ok(n - 2)
                } else {
                    Ok(n - 1)
                }
            } else {
                Ok(n)
            }
        }
        e => e,
    }
}

/// Qiuckly compare read names, will produce many false positive, but with many reads that is fine.
#[inline]
fn equal_names_fast(name1: &[u8], name2: &[u8]) -> bool {
    let n = name1.len();
    let m = name2.len();
    // Sometimes, last two symbols are `/1`, `/2`, so they will not match.
    // Here, we just compare last symbol before that.
    n == m && (n <= 3 || name1[n - 3] == name2[n - 3])
}

// ------------------------------------------- Traits -------------------------------------------

/// Single-end or paired-end record that can be written to a file.
pub trait WritableRecord: Clone + Default {
    fn write_to(&self, f: &mut impl io::Write) -> io::Result<()>;

    /// Returns sum length and the number of fragments (1 or 2).
    fn sum_len(&self) -> (u32, u32);
}

/// Single read sequence.
pub trait SingleRecord: WritableRecord {
    /// Returns record name.
    fn name(&self) -> &[u8];

    /// Safely converts record name to UTF-8.
    fn name_str(&self) -> Cow<'_, str> {
        String::from_utf8_lossy(self.name())
    }

    /// Read sequence.
    fn seq(&self) -> &[u8];
}

impl<T: SingleRecord + WritableRecord> WritableRecord for [T; 2] {
    /// Writes two FASTQ/FASTA records one after another.
    fn write_to(&self, f: &mut impl io::Write) -> io::Result<()> {
        self.iter().map(|rec| rec.write_to(f)).collect()
    }

    fn sum_len(&self) -> (u32, u32) {
        let (s1, t1) = self[0].sum_len();
        let (s2, t2) = self[1].sum_len();
        (s1 + s2, t1 + t2)
    }
}

/// Trait for reading and writing FASTA/Q single-end and paired-end records.
pub trait FastxRead: Send {
    type Record: WritableRecord;

    /// All filenames, associated with the reader.
    fn filenames(&self) -> &[PathBuf];

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self, record: &mut Self::Record) -> Result<bool, Error>;

    /// Writes the next `count` records to the output stream.
    /// Returns the number of records that were written.
    fn write_next_n(&mut self, writer: &mut impl io::Write, count: usize) -> Result<usize, Error> {
        let mut record = Self::Record::default();
        for i in 0..count {
            match self.read_next(&mut record)? {
                true => record.write_to(writer).map_err(add_path!(!))?,
                false => return Ok(i),
            }
        }
        Ok(count)
    }

    /// Subsamples the reader with the `rate` and optional `seed`.
    fn subsample(&mut self, writer: &mut impl io::Write, rate: f64, rng: &mut impl rand::Rng) -> Result<u64, Error> {
        assert!(rate > 0.0 && rate < 1.0, "Subsampling rate ({}) must be within (0, 1).", rate);
        let mut record = Self::Record::default();
        let mut logger = ProcLogger::new();
        while self.read_next(&mut record)? {
            logger.inc();
            if rng.gen::<f64>() <= rate {
                record.write_to(writer).map_err(add_path!(!))?;
            }
        }
        Ok(logger.reads)
    }

    /// Writes input stream to output.
    fn copy(&mut self, writer: &mut impl io::Write) -> Result<u64, Error> {
        let mut record = Self::Record::default();
        let mut logger = ProcLogger::new();
        while self.read_next(&mut record)? {
            logger.inc();
            record.write_to(writer).map_err(add_path!(!))?;
        }
        Ok(logger.reads)
    }
}

// ------------------------------------------ Single read and various readers ------------------------------------------

/// FASTA/Q record.
#[derive(Default, Clone)]
pub struct FastxRecord {
    /// Name, including description after the space.
    name: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl FastxRecord {
    /// Returns true if the record is empty.
    pub fn is_empty(&self) -> bool {
        self.name.is_empty()
    }
}

impl SingleRecord for FastxRecord {
    /// Returns read name.
    fn name(&self) -> &[u8] {
        &self.name
    }

    // Returns record sequence.
    fn seq(&self) -> &[u8] {
        &self.seq
    }
}

impl WritableRecord for FastxRecord {
    /// Write single record to FASTA/FASTQ file.
    fn write_to(&self, f: &mut impl io::Write) -> io::Result<()> {
        if self.qual.is_empty() {
            write_fasta(f, &self.name, &self.seq)
        } else {
            write_fastq(f, &self.name, &self.seq, &self.qual)
        }
    }

    fn sum_len(&self) -> (u32, u32) {
        (self.seq.len() as u32, 1)
    }
}

impl fmt::Debug for FastxRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut v: Vec<u8> = Vec::new();
        self.write_to(&mut v).unwrap();
        f.write_str(&String::from_utf8_lossy(&v))
    }
}

/// Dynamic FASTA/Q reader.
pub struct Reader<R: BufRead> {
    filenames: Vec<PathBuf>,
    stream: R,
    buffer: Vec<u8>,
}

impl Reader<Box<dyn BufRead + Send>> {
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        Self::new(ext::sys::open(path)?, vec![path.to_path_buf()])
    }
}

impl Reader<ext::sys::FileChain> {
    /// Concatenates multiple FASTA/FASTQ readers.
    pub fn from_paths(paths: &[PathBuf]) -> Result<Self, Error> {
        Self::new(ext::sys::FileChain::new(paths)?, paths.to_vec())
    }
}

impl<R: BufRead> Reader<R> {
    pub fn new(mut stream: R, filenames: Vec<PathBuf>) -> Result<Self, Error> {
        let mut buffer = Vec::new();
        read_line(&mut stream, &mut buffer).map_err(|e| Error::Io(e, filenames.clone()))?;
        Ok(Self { stream, buffer, filenames })
    }

    fn fill_fasta_record(&mut self, record: &mut FastxRecord) -> Result<bool, Error> {
        self.buffer.clear();
        let mut seq_len = 0;
        loop {
            let n = read_line(&mut self.stream, &mut record.seq).map_err(|e| Error::Io(e, self.filenames.clone()))?;
            if n == 0 {
                // File ended.
                if seq_len == 0 {
                    return Err(Error::InvalidData(format!("Fasta record {} has an empty sequence.",
                        record.name_str())));
                }
                return Ok(true);
            } else if record.seq[seq_len] == b'>' || record.seq[seq_len] == b'@' {
                self.buffer.extend_from_slice(&record.seq[seq_len..]);
                record.seq.truncate(seq_len);
                return Ok(true);
            } else {
                seq_len += n;
            }
        }
    }

    fn fill_fastq_record(&mut self, record: &mut FastxRecord) -> Result<bool, Error> {
        // Sequence
        let seq_len = read_line(&mut self.stream, &mut record.seq).map_err(|e| Error::Io(e, self.filenames.clone()))?;
        let prev_buf_len = self.buffer.len();

        // +
        let n = read_line(&mut self.stream, &mut self.buffer).map_err(|e| Error::Io(e, self.filenames.clone()))?;
        if n == 0 {
            return Err(Error::InvalidData(
                format!("Fastq record {} is incomplete", record.name_str())));
        } else if self.buffer[prev_buf_len] != b'+' {
            return Err(Error::InvalidData(
                format!("Fastq record {} has incorrect format", record.name_str())));
        }

        // Qualities
        let qual_len = read_line(&mut self.stream, &mut record.qual).map_err(|e| Error::Io(e, self.filenames.clone()))?;
        if seq_len != qual_len {
            return Err(Error::InvalidData(
                format!("Fastq record {} has non-matching sequence and qualities", record.name_str())));
        }

        // Next record name.
        self.buffer.clear();
        read_line(&mut self.stream, &mut self.buffer).map_err(|e| Error::Io(e, self.filenames.clone()))?;
        Ok(true)
    }
}

impl<R: BufRead + Send> Reader<R> {
    /// Reads the next record and standardizes its sequence.
    pub fn read_next_standardized(&mut self, record: &mut FastxRecord) -> Result<bool, Error> {
        if self.read_next(record)? {
            crate::seq::standardize(&mut record.seq)
                .map_err(|nt| Error::InvalidData(format!("Invalid nucleotide `{}` ({}) for sequence {}",
                    char::from(nt), nt, String::from_utf8_lossy(record.name()))))?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Reads all sequences into memory.
    /// Each sequence is standardized and checked for invalid nucleotides.
    pub fn read_all(&mut self) -> Result<Vec<NamedSeq>, Error> {
        let mut record = FastxRecord::default();
        let mut records = Vec::new();
        while self.read_next_standardized(&mut record)? {
            let name = String::from_utf8(record.name().to_vec())
                .map_err(|_| Error::Utf8("read name", record.name().to_vec()))?;
            records.push(NamedSeq::new(name, record.seq().to_owned()));
        }
        Ok(records)
    }
}

impl<R: BufRead + Send> FastxRead for Reader<R> {
    type Record = FastxRecord;

    /// Reads the next record, and returns true if the read was successful (false if no more reads available).
    fn read_next(&mut self, record: &mut FastxRecord) -> Result<bool, Error> {
        record.name.clear();
        record.seq.clear();
        record.qual.clear();

        // Buffer already contains the first line of the next record.
        // If it is empty, the file has ended.
        if self.buffer.is_empty() {
            return Ok(false);
        }

        // Read only name before any description.
        let end_ix = self.buffer.iter().position(|&ch| ch == b' ').unwrap_or(self.buffer.len());
        record.name.extend_from_slice(&self.buffer[1..end_ix]);
        if self.buffer[0] == b'>' {
            self.fill_fasta_record(record)
        } else {
            self.fill_fastq_record(record)
        }
    }

    fn filenames(&self) -> &[PathBuf] {
        &self.filenames
    }
}

/// Interleaved paired-end FASTA/Q reader, that stores two buffer records to reduce memory allocations.
pub struct PairedEndInterleaved<T: SingleRecord, R: FastxRead<Record = T>> {
    reader: R,
}

impl<T: SingleRecord, R: FastxRead<Record = T>> PairedEndInterleaved<T, R> {
    pub fn new(reader: R) -> Self {
        Self { reader }
    }
}

impl<T: SingleRecord, R: FastxRead<Record = T>> FastxRead for PairedEndInterleaved<T, R> {
    type Record = [T; 2];

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self, [record1, record2]: &mut [T; 2]) -> Result<bool, Error> {
        if !self.reader.read_next(record1)? {
            Ok(false)
        } else if !self.reader.read_next(record2)? {
            Err(Error::InvalidData(format!(
                "Odd number of records in an interleaved input file(s) {}",
                ext::fmt::paths(self.reader.filenames()))))
        } else {
            if !equal_names_fast(record1.name(), record2.name()) {
                Err(Error::InvalidData(format!(
                    "Interleaved input file(s) {} contains non matching first and second mate ({} and {})",
                    ext::fmt::paths(self.reader.filenames()), record1.name_str(), record2.name_str())))
            } else {
                Ok(true)
            }
        }
    }

    fn filenames(&self) -> &[PathBuf] {
        self.reader.filenames()
    }
}

/// Two paired-end FASTA/Q readers, that stores two buffer records to reduce memory allocations.
pub struct PairedEndReaders<T: SingleRecord, R: FastxRead<Record = T>, S: FastxRead<Record = T>> {
    filenames: Vec<PathBuf>,
    reader1: R,
    reader2: S,
}

impl<T: SingleRecord, R: FastxRead<Record = T>, S: FastxRead<Record = T>> PairedEndReaders<T, R, S> {
    pub fn new(reader1: R, reader2: S) -> Self {
        Self {
            filenames: [reader1.filenames(), reader2.filenames()].concat(),
            reader1, reader2,
        }
    }
}

impl<T: SingleRecord, R: FastxRead<Record = T>, S: FastxRead<Record = T>> FastxRead for PairedEndReaders<T, R, S> {
    type Record = [T; 2];

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self, [record1, record2]: &mut [T; 2]) -> Result<bool, Error> {
        let could_read1 = self.reader1.read_next(record1)?;
        let could_read2 = self.reader2.read_next(record2)?;
        match (could_read1, could_read2) {
            (false, false) => Ok(false),
            (true, true) => {
                if !equal_names_fast(record1.name(), record2.name()) {
                    Err(Error::InvalidData(format!(
                        "Paired-end input files {} and {} have non matching first and second mates ({} and {})",
                        ext::fmt::paths(self.reader1.filenames()), ext::fmt::paths(self.reader2.filenames()),
                        record1.name_str(), record2.name_str())))
                } else {
                    Ok(true)
                }
            }
            _ => Err(Error::InvalidData(format!("Different number of records in paired-end input files {} and {}",
                ext::fmt::paths(self.reader1.filenames()), ext::fmt::paths(self.reader2.filenames())))),
        }
    }

    /// Returns the first filename.
    fn filenames(&self) -> &[PathBuf] {
        &self.filenames
    }
}

// Write log messages at most every ten seconds.
pub const UPDATE_SECS: u64 = 10;

struct ProcLogger {
    timer: Instant,
    last_msg: Duration,
    reads: u64,
}

impl ProcLogger {
    fn new() -> Self {
        Self {
            timer: Instant::now(),
            last_msg: Duration::default(),
            reads: 0,
        }
    }

    #[inline]
    fn inc(&mut self) {
        self.reads += 1;
        if self.reads % 100_000 == 0 {
            let elapsed = self.timer.elapsed();
            if (elapsed - self.last_msg).as_secs() >= UPDATE_SECS {
                self.last_msg = elapsed;
                let reads_f64 = self.reads as f64;
                let speed = 1e-3 * f64::from(reads_f64) / elapsed.as_secs_f64();
                log::debug!("    Processed {:7.1}M reads, {:4.0}k reads/s", 1e-6 * reads_f64, speed);
            }
        }
    }
}

// /// Length of a sequence without last \n and possible \r before that.
// #[inline(always)]
// fn seq_length(line: &[u8]) -> u64 {
//     let n = line.len();
//     let carriage_rtn = line.get(n - 2).expect("Sequence is too short") == '\r';
//     (n as u64).saturating_sub()
//     n as u64 - (if *line.get(n - 2).expect("Sequence is too short") }) == b'\r' { 2 } else { 1 })
// }

fn count_reads_fasta(mut stream: impl BufRead) -> io::Result<(u64, u64)> {
    let mut line = Vec::with_capacity(4096);
    // Check first line for carriage return.
    let has_carriage_return = if stream.read_until(b'\n', &mut line)? > 0 {
        line.len() >= 2 && line[line.len() - 2] == b'\r'
    } else {
        return Ok((0, 0));
    };

    // Already read one line.
    let mut count: u64 = 1;
    let mut sum_size: u64 = 0;
    let mut n_seq_lines: u64 = 0;

    line.clear();
    while stream.read_until(b'\n', &mut line)? > 0 {
        if line[0] == b'>' {
            count += 1;
        } else {
            sum_size += line.len() as u64;
            n_seq_lines += 1;
        }
        line.clear();
    }
    // Do not count \n.
    sum_size -= n_seq_lines;
    if has_carriage_return {
        // Do not count \r.
        sum_size -= n_seq_lines;
    }
    Ok((count, sum_size))
}

fn count_reads_fastq(mut stream: impl BufRead) -> io::Result<(u64, u64)> {
    let mut count: u64 = 0;
    let mut sum_size: u64 = 0;

    let mut name = Vec::with_capacity(4096);
    let mut seq = Vec::with_capacity(4096);
    let mut buffer = Vec::with_capacity(4096);
    while stream.read_until(b'\n', &mut name)? > 0 {
        name.clear();
        seq.clear();
        buffer.clear();

        count += 1;
        stream.read_until(b'\n', &mut seq)?;
        sum_size += seq.len() as u64;
        stream.read_until(b'\n', &mut buffer)?;
        stream.read_until(b'\n', &mut buffer)?;
    }

    // Remove \n count times.
    sum_size -= count;
    if seq.len() >= 2 && seq[seq.len() - 2] == b'\r' {
        // Assume that all sequences end with \r, remove corresponding number of characters.
        sum_size -= count;
    }
    Ok((count, sum_size))
}

/// Count the number of FASTA/FASTQ reads in the stream, as well as the total number of basepairs.
/// NOTE: Fasta and Fastq reads should not be present in the same file.
pub fn count_reads_fastx(path: &Path) -> Result<(u64, u64), Error> {
    let mut stream = ext::sys::open(path)?;
    let mut first_byte = [0_u8; 1];
    if stream.read(&mut first_byte).map_err(add_path!(path))? == 1 {
        match first_byte[0] {
            b'>' => count_reads_fasta(stream),
            b'@' => count_reads_fastq(stream),
            _ => {
                log::error!("Unexpected first symbol '{}', assuming FASTQ file", char::from(first_byte[0]));
                count_reads_fastq(stream)
            }
        }.map_err(add_path!(path))
    } else {
        Ok((0, 0))
    }
}

/// Counts the number of reads in a BAM/CRAM file.
pub fn count_reads_bam(path: &Path, samtools: &Path, reference: &Option<PathBuf>, threads: u16) -> Result<u64, Error> {
    let mut command = Command::new(samtools);
    command.args(&["view",
        "-c", // Count reads,
        "-F", "0x980", // Discard secondary, supplemenary alignments + discard second read ends,
        "-@", &threads.to_string()]);
    if let Some(filename) = reference {
        command.arg("-T").arg(filename);
    }
    command.arg(path).stdout(Stdio::piped()).stderr(Stdio::piped());
    let child = command.spawn().map_err(add_path!(!))?;
    let pipe_guard = ext::sys::PipeGuard::new(samtools.to_path_buf(), child);
    let output = pipe_guard.wait()?;
    // Use match instead of `map_err` to suppress brrow checked error.
    let s = match std::str::from_utf8(&output.stdout) {
        Ok(val) => val,
        Err(_) => return Err(Error::Utf8("samtools output", output.stdout)),
    };
    s.trim().parse().map_err(|_| Error::RuntimeError(format!(
        "Cannot parse samtools output: {:?} is not a number", s.trim())))
}

/// Wrapper over bam record, which decodes sequence into vector.
#[derive(Default, Clone)]
pub struct BamRecord {
    inner: bam::Record,
    seq: Vec<u8>,
}

impl BamRecord {
    fn name_str(&self) -> Cow<'_, str> {
        String::from_utf8_lossy(&self.inner.qname())
    }

    /// Read until the next primary alignment is encountered, which is then saved and its sequence decoded.
    /// Returns false if no alignments are available.
    fn read_from(&mut self, reader: &mut impl bam::Read, min_start: i64) -> Result<bool, Error> {
        self.seq.clear();
        loop {
            match reader.read(&mut self.inner) {
                Some(Err(e)) => return Err(e.into()),
                Some(Ok(())) => {}
                None => return Ok(false),
            }
            // Primary alignment.
            if (self.flag() & 0x900) == 0 && self.inner.pos() >= min_start {
                let seq = self.inner.seq();
                let l = seq.len();
                if l == 0 {
                    return Err(Error::InvalidData(format!(
                        "Primary alignment for read {} has no sequence", self.name_str())));
                }
                if self.inner.is_reverse() {
                    self.seq.extend((0..l).rev().map(|i| seq::complement_nt(seq[i])));
                } else {
                    self.seq.extend((0..l).map(|i| seq[i]));
                }
                return Ok(true);
            }
        }
    }

    #[inline]
    fn flag(&self) -> u16 {
        self.inner.flags()
    }
}

impl SingleRecord for BamRecord {
    fn name(&self) -> &[u8] {
        self.inner.qname()
    }

    fn seq(&self) -> &[u8] {
        &self.seq
    }
}

impl WritableRecord for BamRecord {
    /// Write single record to FASTA/FASTQ file.
    fn write_to(&self, f: &mut impl io::Write) -> io::Result<()> {
        let mut qual = self.inner.qual();
        if qual.is_empty() || qual[0] == 255 {
            return write_fasta(f, &self.inner.qname(), &self.seq);
        }

        f.write_all(b"@")?;
        f.write_all(self.inner.qname())?;
        f.write_all(b"\n")?;
        f.write_all(&self.seq)?;
        f.write_all(b"\n+\n")?;

        // Will be enough for short reads.
        const BUF_LEN: usize = 512;
        let mut buffer = [0_u8; BUF_LEN];
        if !self.inner.is_reverse() {
            while !qual.is_empty() {
                for (val, &q) in buffer.iter_mut().zip(qual.iter()) {
                    *val = q + 33;
                }
                let m = min(BUF_LEN, qual.len());
                f.write_all(&buffer[..m])?;
                qual = &qual[m..];
            }
        } else {
            while !qual.is_empty() {
                for (val, &q) in buffer.iter_mut().zip(qual.iter().rev()) {
                    *val = q + 33;
                }
                let m = min(BUF_LEN, qual.len());
                f.write_all(&buffer[..m])?;
                qual = &qual[..qual.len() - m];
            }
        }
        f.write_all(b"\n")
    }

    fn sum_len(&self) -> (u32, u32) {
        (self.seq.len() as u32, 1)
    }
}

/// Hash reads based on their names.
impl std::hash::Hash for BamRecord {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        hasher.write(self.inner.qname())
    }
}

/// Check only read name because we intentionally want to get other read end from a map.
impl PartialEq for BamRecord {
    fn eq(&self, oth: &Self) -> bool {
        self.inner.qname() == oth.inner.qname()
    }
}

impl Eq for BamRecord {}

/// Read BAM/CRAM file directly, without any fetching.
pub struct DirectBamReader {
    filenames: Vec<PathBuf>,
    current_reader: bam::Reader,
    /// Following BAM readers, in reverse order.
    future_readers: Vec<bam::Reader>,
}

impl DirectBamReader {
    pub fn new(
        filenames: &[PathBuf],
        mut build_reader: impl FnMut(&Path) -> Result<bam::Reader, Error>,
    ) -> Result<Self, Error>
    {
        let filenames = filenames.to_vec();
        let current_reader = build_reader(&filenames[0])?;
        let future_readers = filenames[1..].iter().rev().map(|f| build_reader(f))
            .collect::<Result<Vec<_>, Error>>()?;
        Ok(Self { filenames, current_reader, future_readers })
    }
}

impl FastxRead for DirectBamReader {
    type Record = BamRecord;

    fn filenames(&self) -> &[PathBuf] {
        &self.filenames
    }

    /// Reads the next record, and returns true if the read was successful (false if no more reads available).
    fn read_next(&mut self, record: &mut BamRecord) -> Result<bool, Error> {
        while !record.read_from(&mut self.current_reader, i64::MIN)? {
            if let Some(reader) = self.future_readers.pop() {
                self.current_reader = reader;
            } else {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

/// Read single-end BAM/CRAM file across multiple fetches.
pub struct IndexedBamReader {
    filename: [PathBuf; 1],
    reader: bam::IndexedReader,
    /// Fetch regions. Contig IDs must be taken from exactly this bam/cram reader.
    /// All unmapped reads are then fetch in any case.
    /// NOTE: Must go in increasing order. Otherwise, some reads may be missed.
    regions: Vec<Interval>,
    /// Region, which will be fetched next. If `next_ix == regions.len()`, fetch unmapped reads.
    next_ix: usize,
    /// Recruit reads with at least this starting position.
    /// This is needed to discard reads that are fetched twice in nearby regions.
    min_start: i64,
}

impl IndexedBamReader {
    /// Creates bam reader from indexed reader and an iterator over regions.
    /// WARN: Regions should go in increasing order, otherwise some code may fail in paired-end case.
    pub fn new(filename: PathBuf, reader: bam::IndexedReader, regions: Vec<Interval>) -> Result<Self, Error> {
        let mut res = Self {
            filename: [filename],
            reader, regions,
            next_ix: 0,
            min_start: i64::MIN,
        };
        res.fetch_next()?;
        Ok(res)
    }

    /// Fetch the next region.
    /// Returns false if no more fetches available.
    fn fetch_next(&mut self) -> Result<bool, Error> {
        self.min_start = i64::MIN;

        if let Some(region) = self.regions.get(self.next_ix) {
            self.reader.fetch(bam::FetchDefinition::Region(
                region.contig_id().ix() as i32, i64::from(region.start()), i64::from(region.end())))?;
            if self.next_ix > 0 {
                let prev = &self.regions[self.next_ix - 1];
                assert!(prev < region, "Cannot fetch reads from BAM file: regions must be sorted");
                if prev.contig_id() == region.contig_id() {
                    self.min_start = i64::from(prev.end());
                }
            }

        } else if self.next_ix == self.regions.len() {
            self.reader.fetch(bam::FetchDefinition::Unmapped)?;
        } else {
            return Ok(false);
        }
        self.next_ix += 1;
        Ok(true)
    }
}

impl FastxRead for IndexedBamReader {
    type Record = BamRecord;

    /// Reads the next record, and returns true if the read was successful (false if no more reads available).
    fn read_next(&mut self, record: &mut BamRecord) -> Result<bool, Error> {
        loop {
            if record.read_from(&mut self.reader, self.min_start)? {
                return Ok(true);
            } else if !self.fetch_next()? {
                return Ok(false);
            }
        }
    }

    fn filenames(&self) -> &[PathBuf] {
        &self.filename
    }
}

pub struct PairedBamReader {
    reader: IndexedBamReader,
    /// Hash set with records, for which we did not find their mate yet.
    /// Only read name is checked, so we will get a match when we search for the read mate.
    read_pairs: HashSet<BamRecord>,
    discarded: u64,
}

impl PairedBamReader {
    pub fn new(reader: IndexedBamReader) -> Self {
        Self {
            reader,
            read_pairs: Default::default(),
            discarded: 0,
        }
    }
}

impl FastxRead for PairedBamReader {
    type Record = [BamRecord; 2];

    /// Reads records until both read ends are found for some reads, and returns it.
    /// All unpaired reads are saved in a set until their mate is also found.
    /// In the end, all unpaired reads are discarded.
    fn read_next(&mut self, [rec1, rec2]: &mut Self::Record) -> Result<bool, Error> {
        loop {
            if !self.reader.read_next(rec1)? {
                let discarded = self.discarded + self.read_pairs.len() as u64;
                if discarded > 0 {
                    log::warn!("    Discarded {} reads without a pair", discarded);
                }
                return Ok(false);
            }

            if (rec1.flag() & 0xC0) == 0 {
                return Err(Error::InvalidData(format!(
                    "Expected paired-end reads, but read {} is unpaired", rec1.name_str())));
            }

            if let Some(mate) = self.read_pairs.take(rec1) {
                // Read mate was found.
                let rec_first = rec1.inner.is_first_in_template();
                if rec_first == mate.inner.is_first_in_template() {
                    return Err(Error::InvalidData(format!(
                        "Found two primary alignments for {} for the same read end", rec1.name_str())));
                }
                *rec2 = mate;
                if !rec_first {
                    std::mem::swap(rec1, rec2);
                }
                return Ok(true);
            } else {
                // Read mate was not found.
                // Cast to unsigned integer so that -1 becomes the largest value.
                let curr_coord = (rec1.inner.tid() as u32, rec1.inner.pos() as u64);
                let mate_coord = (rec1.inner.mtid() as u32, rec1.inner.mpos() as u64);
                // Do not add record if its mate has smaller coordinates.
                if curr_coord <= mate_coord {
                    self.read_pairs.insert(rec1.clone());
                } else {
                    self.discarded += 1;
                }
            }
        }
    }

    fn filenames(&self) -> &[PathBuf] {
        self.reader.filenames()
    }
}

/// Calculates average read length across the first `n_records`.
pub fn mean_read_len<T, R>(reader: &mut R, n_records: usize) -> Result<f64, Error>
where T: WritableRecord,
      R: FastxRead<Record = T>,
{
    let mut record = Default::default();
    let mut s: u64 = 0;
    let mut n: u64 = 0;
    while reader.read_next(&mut record)? {
        let (l, t) = record.sum_len();
        s += u64::from(l);
        n += u64::from(t);
        if n >= n_records as u64 {
            break;
        }
    }
    if n == 0 {
        Err(Error::InvalidData(format!("Cannot calculate average read length: empty input file(s) {}",
            ext::fmt::paths(reader.filenames()))))
    } else {
        Ok(s as f64 / n as f64)
    }
}

pub trait ExtendedHtslibReader: bam::Read {
    fn set_reference(&mut self, ref_filename: &Path) -> Result<(), htslib::errors::Error>;
}

impl ExtendedHtslibReader for bam::Reader {
    fn set_reference(&mut self, ref_filename: &Path) -> Result<(), htslib::errors::Error> {
        self.set_reference(ref_filename)
    }
}

impl ExtendedHtslibReader for bam::IndexedReader {
    fn set_reference(&mut self, ref_filename: &Path) -> Result<(), htslib::errors::Error> {
        self.set_reference(ref_filename)
    }
}

/// If input file is CRAM, sets its reference. All contigs from the CRAM must match reference.
pub fn set_reference(
    aln_filename: &Path,
    aln_reader: &mut impl ExtendedHtslibReader,
    ref_filename: &Option<PathBuf>,
    ref_contigs: Option<&ContigNames>,
) -> Result<(), Error>
{
    // Possible, that there are no contigs in the BAM/CRAM file: all reads are just stored for compression purposes (?).
    if aln_reader.header().target_count() == 0 {
        return Ok(());
    }

    let ext = aln_filename.extension();
    if ext == Some(OsStr::new("bam")) || ext == Some(OsStr::new("BAM"))
            || ext == Some(OsStr::new("sam")) || ext == Some(OsStr::new("SAM")) {
        return Ok(());
    } else if ext != Some(OsStr::new("cram")) || ext == Some(OsStr::new("CRAM")) {
        log::warn!("Cannot determine alignment file format {}, assuming CRAM", ext::fmt::path(aln_filename));
    }

    let Some(ref_filename) = ref_filename else {
        return Err(Error::InvalidInput("Cannot analyze CRAM file without a reference FASTA".to_owned()))
    };
    let ref_contigs = match ref_contigs {
        Some(val) => Cow::Borrowed(val),
        None => Cow::Owned(ContigNames::load_indexed_fasta("REF", ref_filename)?.0),
    };
    aln_reader.set_reference(ref_filename)?;

    let mut missing = 0;
    let mut incorrect_len = 0;
    let header = aln_reader.header();
    for (tid, contig) in header.target_names().into_iter().enumerate() {
        let contig = std::str::from_utf8(contig)
            .map_err(|_| Error::Utf8("contig name", contig.to_vec()))?;
        if let Some(contig_id) = ref_contigs.try_get_id(&contig) {
            let in_ref_len = u64::from(ref_contigs.get_len(contig_id));
            let in_bam_len = header.target_len(tid as u32).expect("Unknown contig length");
            if in_ref_len != in_bam_len {
                if incorrect_len == 0 {
                    log::error!("Reference {} does not match alignments {}: contig {:?} lengths differ ({} and {})",
                        ext::fmt::path(ref_filename), ext::fmt::path(aln_filename), contig, in_ref_len, in_bam_len);
                }
                incorrect_len += 1;
            }
        } else {
            if missing == 0 {
                log::error!("Reference {} does not match alignments {}: missing contig {:?}",
                    ext::fmt::path(ref_filename), ext::fmt::path(aln_filename), contig);
            }
            missing += 1;
        }
    }

    if missing > 0 || incorrect_len > 0 {
        Err(Error::InvalidInput(format!(
            "Reference does not match alignments: {} missing contigs and {} contigs with incorrect length",
            missing, incorrect_len)))
    } else {
        Ok(())
    }
}

/// Based on the input arguments, selects appropriate reader and performs `action(reader)`.
/// Input arguments must contain:
/// - `args.reads1: Vec<PathBuf>` - first or single-end FASTA/FASTQ filenames,
/// - `args.reads2: Vec<PathBuf>` - second mate FASTA/FASTQ filenames,
/// - `args.alns: Vec<PathBuf>` - zero, one or several BAM/CRAM filename,
/// - `args.reference: Option<PathBuf>` - zero or one FASTA reference file,
/// - `args.interleaved: bool` - are the input reads interleaved.
/// Next: `contigs: Option<&ContigNames>`.
///
/// After semicolon argument: either `let {} reader` or `let {mut} reader`, depending on mutability,
/// Last argument: action on the reader.
macro_rules! process_readers {
    ($args:expr, $contigs:expr; let {$( $mut_:tt )?} $reader:ident; $action:tt) => {
        {
            use crate::seq::fastx;
            let n1 = $args.reads1.len();
            let n2 = $args.reads2.len();
            assert!(n2 == 0 || (n1 == n2));
            let m = $args.alns.len();
            assert!(n1 == 0 || m == 0);

            if n1 >= 2 {
                // Multiple FASTA/Q files.
                let reader1 = fastx::Reader::from_paths(&$args.reads1)?;
                if n2 > 0 {
                    let reader2 = fastx::Reader::from_paths(&$args.reads2)?;
                    let $($mut_)? $reader = fastx::PairedEndReaders::new(reader1, reader2);
                    #[allow(unused_braces)] $action
                } else if $args.interleaved {
                    let $($mut_)? $reader = fastx::PairedEndInterleaved::new(reader1);
                    #[allow(unused_braces)] $action
                } else {
                    let $($mut_)? $reader = reader1;
                    #[allow(unused_braces)] $action
                }
            } else if n1 == 1 {
                // One FASTA/Q files.
                let reader1 = fastx::Reader::from_path(&$args.reads1[0])?;
                if n2 > 0 {
                    let reader2 = fastx::Reader::from_path(&$args.reads2[0])?;
                    let $($mut_)? $reader = fastx::PairedEndReaders::new(reader1, reader2);
                    #[allow(unused_braces)] $action
                } else if $args.interleaved {
                    let $($mut_)? $reader = fastx::PairedEndInterleaved::new(reader1);
                    #[allow(unused_braces)] $action
                } else {
                    let $($mut_)? $reader = reader1;
                    #[allow(unused_braces)] $action
                }
            } else {
                let bam_reader = fastx::DirectBamReader::new(&$args.alns, |path| {
                    let mut inner_reader = bam::Reader::from_path(path)?;
                    fastx::set_reference(path, &mut inner_reader, &$args.reference, $contigs)?;
                    Ok(inner_reader)
                })?;
                if $args.interleaved {
                    let $($mut_)? $reader = fastx::PairedEndInterleaved::new(bam_reader);
                    #[allow(unused_braces)] $action
                } else {
                    let $($mut_)? $reader = bam_reader;
                    #[allow(unused_braces)] $action
                }
            }
        }
    }
}

pub(crate) use process_readers;
