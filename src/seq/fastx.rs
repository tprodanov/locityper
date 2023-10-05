use std::{
    fmt,
    io::{self, BufRead},
    cmp::{min, Ordering},
    path::{Path, PathBuf},
    time::{Instant, Duration},
};
use htslib::bam::{self, Read as BamRead};
use crate::{
    ext,
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
    full_name: &[u8],
    seq: &[u8]
) -> io::Result<()> {
    writer.write_all(b">")?;
    writer.write_all(full_name)?;
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
    full_name: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> io::Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(full_name)?;
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

// ------------------------------------------- Traits -------------------------------------------

/// Single read with name, sequence and qualities.
pub trait SingleRecord: Clone + Default {
    fn name(&self) -> &[u8];

    fn seq(&self) -> &[u8];

    /// If qualities are not available, this slice will be empty.
    fn qual(&self) -> &[u8];
}

/// Single-end or paired-end record that can be written to a file.
pub trait WritableRecord: Clone + Default {
    fn write_to(&self, f: &mut impl io::Write) -> io::Result<()>;
}

impl<T: SingleRecord> WritableRecord for T {
    /// Write single record to FASTA/FASTQ file.
    fn write_to(&self, f: &mut impl io::Write) -> io::Result<()> {
        if self.qual().is_empty() {
            write_fasta(f, self.name(), self.seq())
        } else {
            write_fastq(f, self.name(), self.seq(), self.qual())
        }
    }
}

impl<T: SingleRecord> WritableRecord for [T; 2] {
    /// Writes two FASTQ/FASTA records one after another.
    fn write_to(&self, f: &mut impl io::Write) -> io::Result<()> {
        self.iter().map(|rec| rec.write_to(f)).collect()
    }
}

// ------------------------------------------ Single read and various readers ------------------------------------------

/// FASTA/Q record.
#[derive(Default, Clone)]
pub struct FastxRecord {
    /// Name, including description after the space.
    full_name: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl FastxRecord {
    /// Returns true if the record is empty.
    pub fn is_empty(&self) -> bool {
        self.full_name.is_empty()
    }

    pub fn full_name_str(&self) -> std::borrow::Cow<'_, str> {
        String::from_utf8_lossy(&self.full_name)
    }

    /// Returns only the name of the record, without the optional description.
    pub fn name_only(&self) -> std::borrow::Cow<'_, str> {
        match self.full_name.iter().position(|&c| c == b' ') {
            Some(i) => String::from_utf8_lossy(&self.full_name[..i]),
            None => String::from_utf8_lossy(&self.full_name),
        }
    }
}

impl SingleRecord for FastxRecord {
    /// Returns the name of the record, including description after space, if any.
    fn name(&self) -> &[u8] {
        &self.full_name
    }

    // Returns record sequence.
    fn seq(&self) -> &[u8] {
        &self.seq
    }

    /// Returns record qualities. If unavailable, returns empty slice.
    fn qual(&self) -> &[u8] {
        &self.qual
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
    filename: PathBuf,
    stream: R,
    buffer: Vec<u8>,
}

impl Reader<Box<dyn BufRead + Send>> {
    pub fn from_path(path: &Path) -> Result<Self, Error> {
        Self::new(ext::sys::open(path)?, path.to_path_buf())
    }
}

impl<R: BufRead> Reader<R> {
    pub fn new(mut stream: R, filename: PathBuf) -> Result<Self, Error> {
        let mut buffer = Vec::new();
        read_line(&mut stream, &mut buffer).map_err(add_path!(filename))?;
        Ok(Self { stream, buffer, filename })
    }

    fn fill_fasta_record(&mut self, record: &mut FastxRecord) -> Result<bool, Error> {
        self.buffer.clear();
        let mut seq_len = 0;
        loop {
            let n = read_line(&mut self.stream, &mut record.seq).map_err(add_path!(self.filename))?;
            if n == 0 {
                // File ended.
                if seq_len == 0 {
                    return Err(Error::InvalidData(format!("Fasta record {} has an empty sequence.",
                        record.name_only())));
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
        let seq_len = read_line(&mut self.stream, &mut record.seq).map_err(add_path!(self.filename))?;
        let prev_buf_len = self.buffer.len();

        // +
        let n = read_line(&mut self.stream, &mut self.buffer).map_err(add_path!(self.filename))?;
        if n == 0 {
            return Err(Error::InvalidData(
                format!("Fastq record {} is incomplete.", record.name_only())));
        } else if self.buffer[prev_buf_len] != b'+' {
            return Err(Error::InvalidData(
                format!("Fastq record {} has incorrect format.", record.name_only())));
        }

        // Qualities
        let qual_len = read_line(&mut self.stream, &mut record.qual).map_err(add_path!(self.filename))?;
        if seq_len != qual_len {
            return Err(Error::InvalidData(
                format!("Fastq record {} has non-matching sequence and qualities.", record.name_only())));
        }

        // Next record name.
        self.buffer.clear();
        read_line(&mut self.stream, &mut self.buffer).map_err(add_path!(self.filename))?;
        Ok(true)
    }
}

impl<R: BufRead + Send> Reader<R> {
    /// Calculates mean read length across at most `max_records` entries.
    pub fn mean_read_len(&mut self, max_records: usize) -> Result<f64, Error> {
        let mut count: u64 = 0;
        let mut sum_len: u64 = 0;
        let mut record = FastxRecord::default();
        for _ in 0..max_records {
            match self.read_next(&mut record)? {
                true => {
                    count += 1;
                    sum_len += record.seq().len() as u64;
                }
                false => break,
            }
        }
        Ok(sum_len as f64 / count as f64)
    }
}

impl<R: BufRead + Send> Reader<R> {
    /// Reads all sequences into memory.
    /// Each sequence is standartized and checked for invalid nucleotides.
    pub fn read_all(&mut self) -> Result<Vec<NamedSeq>, Error> {
        let mut record = FastxRecord::default();
        let mut records = Vec::new();
        while self.read_next(&mut record)? {
            let name = record.name_only().into_owned();
            let mut seq = record.seq().to_owned();
            crate::seq::standardize(&mut seq)
                .map_err(|nt| Error::InvalidData(format!("Invalid nucleotide `{}` ({}) for sequence {}",
                    char::from(nt), nt, name)))?;
            records.push(NamedSeq::new(name, seq));
        }
        Ok(records)
    }
}

/// Trait for reading and writing FASTA/Q single-end and paired-end records.
pub trait FastxRead: Send {
    type Record: WritableRecord;

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
        let mut logger = SubsampleLogger::new();
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
        let mut logger = SubsampleLogger::new();
        while self.read_next(&mut record)? {
            logger.inc();
            record.write_to(writer).map_err(add_path!(!))?;
        }
        Ok(logger.reads)
    }
}

impl<R: BufRead + Send> FastxRead for Reader<R> {
    type Record = FastxRecord;

    /// Reads the next record, and returns true if the read was successful (false if no more reads available).
    fn read_next(&mut self, record: &mut FastxRecord) -> Result<bool, Error> {
        record.full_name.clear();
        record.seq.clear();
        record.qual.clear();

        // Buffer already contains the first line of the next record.
        // If it is empty, the file has ended.
        if self.buffer.is_empty() {
            return Ok(false);
        }

        record.full_name.extend_from_slice(&self.buffer[1..]);
        if self.buffer[0] == b'>' {
            self.fill_fasta_record(record)
        } else {
            self.fill_fastq_record(record)
        }
    }
}

/// Interleaved paired-end FASTA/Q reader, that stores two buffer records to reduce memory allocations.
pub struct PairedEndInterleaved<R: BufRead> {
    reader: Reader<R>,
}

impl<R: BufRead> PairedEndInterleaved<R> {
    pub fn new(reader: Reader<R>) -> Self {
        Self { reader }
    }
}

impl<R: BufRead + Send> FastxRead for PairedEndInterleaved<R> {
    type Record = [FastxRecord; 2];

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self, paired_record: &mut Self::Record) -> Result<bool, Error> {
        if !self.reader.read_next(&mut paired_record[0])? {
            Ok(false)
        } else if !self.reader.read_next(&mut paired_record[1])? {
            Err(Error::InvalidData(format!(
                "Odd number of records in an interleaved input file {}", ext::fmt::path(&self.reader.filename))))
        } else {
            Ok(true)
        }
    }
}

/// Two paired-end FASTA/Q readers, that stores two buffer records to reduce memory allocations.
pub struct PairedEndReaders<R: BufRead, S: BufRead> {
    reader1: Reader<R>,
    reader2: Reader<S>,
}

impl<R: BufRead, S: BufRead> PairedEndReaders<R, S> {
    pub fn new(reader1: Reader<R>, reader2: Reader<S>) -> Self {
        Self { reader1, reader2 }
    }
}

impl<R: BufRead + Send, S: BufRead + Send> FastxRead for PairedEndReaders<R, S> {
    type Record = [FastxRecord; 2];

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self, paired_record: &mut Self::Record) -> Result<bool, Error> {
        let could_read1 = self.reader1.read_next(&mut paired_record[0])?;
        let could_read2 = self.reader2.read_next(&mut paired_record[1])?;
        match (could_read1, could_read2) {
            (false, false) => Ok(false),
            (true, true) => Ok(true),
            _ => Err(Error::InvalidData(format!("Different number of records in input files {} and {}",
                ext::fmt::path(&self.reader1.filename), ext::fmt::path(&self.reader2.filename)))),
        }
    }
}

// Write log messages at most every ten seconds.
pub const UPDATE_SECS: u64 = 10;

struct SubsampleLogger {
    timer: Instant,
    last_msg: Duration,
    reads: u64,
}

impl SubsampleLogger {
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

fn count_reads_fasta(mut stream: impl BufRead) -> io::Result<u64> {
    let mut count = 0;
    let mut line = Vec::with_capacity(4096);
    while stream.read_until(b'>', &mut line)? > 0 {
        count += 1;
        line.clear();
    }
    Ok(count)
}

fn count_reads_fastq(mut stream: impl BufRead) -> io::Result<u64> {
    let mut count = 0;
    let mut line = Vec::with_capacity(4096);
    while stream.read_until(b'\n', &mut line)? > 0 {
        count += 1;
        line.clear();
    }
    if count % 4 != 0 {
        log::warn!("Number of lines in FASTQ file does not divide 4 (mod = {})", count % 4);
    }
    Ok(count / 4)
}

/// Count the number of FASTA/FASTQ reads in the stream.
/// NOTE: Fasta and Fastq reads should not be present in the same file.
pub fn count_reads(path: &Path) -> Result<u64, Error> {
    let mut stream = ext::sys::open(path)?;
    let mut first_byte = [0_u8; 1];
    if stream.read(&mut first_byte).map_err(add_path!(path))? == 1 {
        match first_byte[0] {
            b'>' => count_reads_fasta(stream),
            b'@' => count_reads_fastq(stream),
            _ => {
                log::warn!("Unexpected first symbol '{}', assuming FASTQ file", char::from(first_byte[0]));
                count_reads_fastq(stream)
            }
        }.map_err(add_path!(path))
    } else {
        Ok(0)
    }
}

/// Wrapper over bam record, which decodes sequence into vector.
#[derive(Default, Clone)]
pub struct BamRecord {
    inner: bam::Record,
    seq: Vec<u8>,
}

impl BamRecord {
    fn name_str(&self) -> std::borrow::Cow<'_, str> {
        String::from_utf8_lossy(&self.inner.qname())
    }

    /// Read until the next primary alignment is encountered, which is then saved and its sequence decoded.
    /// Returns false if no alignments are available.
    fn read_from(&mut self, reader: &mut impl bam::Read) -> Result<bool, Error> {
        self.seq.clear();
        loop {
            match reader.read(&mut self.inner) {
                Some(Err(e)) => return Err(e.into()),
                Some(Ok(())) => {}
                None => return Ok(false),
            }
            // Primary alignment.
            if (self.inner.flags() & 0x900) == 0 {
                let seq = self.inner.seq();
                let l = seq.len();
                if l == 0 {
                    return Err(Error::InvalidData(format!(
                        "Primary alignment for read {} has no sequence", self.name_str())));
                }
                self.seq.extend((0..l).map(|i| seq[i]));
                return Ok(true);
            }
        }
    }
}

impl SingleRecord for BamRecord {
    #[inline]
    fn name(&self) -> &[u8] {
        self.inner.qname()
    }

    fn seq(&self) -> &[u8] {
        &self.seq
    }

    fn qual(&self) -> &[u8] {
        &self.inner.qual()
    }
}

/// Read single-end BAM/CRAM file starting with some contig and continuing until the end.
pub struct BamReader {
    reader: bam::IndexedReader,
    /// Current contig index.
    /// If `curr_tid == end_tid`, we are currently reading unmapped reads.
    /// If `curr_tid > end_tid`, we completely finished reading.
    curr_tid: i32,
    end_tid: i32,
}

impl BamReader {
    pub fn new(reader: bam::IndexedReader, start_tid: u32) -> Result<Self, Error> {
        let mut res = Self {
            end_tid: reader.header().target_count() as i32,
            curr_tid: i32::try_from(start_tid).unwrap(),
            reader,
        };
        // NOTE: Do not remove this asser.t
        assert!(res.fetch_current_contig()?, "Contig id is too high");
        Ok(res)
    }

    /// Fetches current contig or unmapped reads, depending on `curr_tid`.
    /// Returns false if the reader is completely traversed.
    fn fetch_current_contig(&mut self) -> Result<bool, Error> {
        match self.curr_tid.cmp(&self.end_tid) {
            Ordering::Less => self.reader.fetch(bam::FetchDefinition::CompleteTid(self.curr_tid))?,
            Ordering::Equal => self.reader.fetch(bam::FetchDefinition::Unmapped)?,
            Ordering::Greater => return Ok(false)
        }
        Ok(true)
    }
}

impl FastxRead for BamReader {
    type Record = BamRecord;

    /// Reads the next record, and returns true if the read was successful (false if no more reads available).
    fn read_next(&mut self, record: &mut BamRecord) -> Result<bool, Error> {
        while !record.read_from(&mut self.reader)? {
            // Need to increment contig by one.
            self.curr_tid += 1;
            if !self.fetch_current_contig()? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}
