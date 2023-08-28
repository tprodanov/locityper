use std::{
    fmt,
    io::{self, BufRead},
    cmp::min,
    path::Path,
};
use crate::{
    ext,
    err::{Error, add_path},
};
use super::{
    NamedSeq,
    kmers::{self, Kmer},
};

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

/// FASTA/Q record.
#[derive(Default, Clone)]
pub struct Record {
    /// Name, including description after the space.
    full_name: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl Record {
    /// Returns true if the record is empty.
    pub fn is_empty(&self) -> bool {
        self.full_name.is_empty()
    }

    /// Returns the name of the record, including description after space.
    pub fn full_name(&self) -> &[u8] {
        &self.full_name
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

    /// Returns record sequence.
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    /// Returns record qualities, if available.
    pub fn qual(&self) -> Option<&[u8]> {
        if self.qual.is_empty() {
            None
        } else {
            Some(&self.qual)
        }
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut v: Vec<u8> = Vec::new();
        self.write_to(&mut v).unwrap();
        f.write_str(&String::from_utf8_lossy(&v))
    }
}

/// Extension over a single- or paired-records.
pub trait RecordExt : Default + Clone + Send + 'static {
    /// Writes the single- or paired-read to the stream.
    fn write_to(&self, writer: &mut impl io::Write) -> io::Result<()>;

    /// Extracts all minimizers from the record (see `kmers::minimizers`).
    /// The function does not clear the buffer in advance.
    fn minimizers<K: Kmer>(&self, k: u8, n: u8, buffer: &mut Vec<K>);
}

impl RecordExt for Record {
    fn write_to(&self, writer: &mut impl io::Write) -> io::Result<()> {
        if self.qual.is_empty() {
            write_fasta(writer, &self.full_name, &self.seq)
        } else {
            write_fastq(writer, &self.full_name, &self.seq, &self.qual)
        }
    }

    fn minimizers<K: Kmer>(&self, k: u8, n: u8, buffer: &mut Vec<K>) {
        kmers::minimizers(&self.seq, k, n, buffer);
    }
}

pub type PairedRecord = [Record; 2];

impl RecordExt for PairedRecord {
    fn write_to(&self, writer: &mut impl io::Write) -> io::Result<()> {
        self[0].write_to(writer)?;
        self[1].write_to(writer)
    }

    fn minimizers<K: Kmer>(&self, k: u8, n: u8, buffer: &mut Vec<K>) {
        kmers::minimizers(&self[0].seq, k, n, buffer);
        kmers::minimizers(&self[1].seq, k, n, buffer);
    }
}

/// Dynamic FASTA/Q reader.
pub struct Reader<R: BufRead> {
    stream: R,
    buffer: Vec<u8>,
    total_reads: u64,
}

impl Reader<Box<dyn BufRead + Send>> {
    pub fn from_path(path: &Path) -> Result<Self, Error> {
        Self::new(ext::sys::open(path)?).map_err(add_path!(path))
    }
}

impl<R: BufRead> Reader<R> {
    pub fn new(mut stream: R) -> io::Result<Self> {
        let mut buffer = Vec::new();
        read_line(&mut stream, &mut buffer)?;
        Ok(Self {
            stream, buffer,
            total_reads: 0,
        })
    }

    fn fill_fasta_record(&mut self, record: &mut Record) -> io::Result<bool> {
        self.buffer.clear();
        let mut seq_len = 0;
        loop {
            let n = read_line(&mut self.stream, &mut record.seq)?;
            if n == 0 {
                // File ended.
                if seq_len == 0 {
                    return Err(io::Error::new(io::ErrorKind::InvalidData,
                        format!("Fasta record {} has an empty sequence.", record.name_only())));
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

    fn fill_fastq_record(&mut self, record: &mut Record) -> io::Result<bool> {
        // Sequence
        let seq_len = read_line(&mut self.stream, &mut record.seq)?;
        let prev_buf_len = self.buffer.len();

        // +
        let n = read_line(&mut self.stream, &mut self.buffer)?;
        if n == 0 {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Fastq record {} is incomplete.", record.name_only())));
        } else if self.buffer[prev_buf_len] != b'+' {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Fastq record {} has incorrect format.", record.name_only())));
        }

        // Qualities
        let qual_len = read_line(&mut self.stream, &mut record.qual)?;
        if seq_len != qual_len {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Fastq record {} has non-matching sequence and qualities.", record.name_only())));
        }

        // Next record name.
        self.buffer.clear();
        read_line(&mut self.stream, &mut self.buffer)?;
        Ok(true)
    }
}

impl<R: BufRead + Send> Reader<R> {
    /// Calculates mean read length across at most `max_records` entries.
    pub fn mean_read_len(&mut self, max_records: usize) -> io::Result<f64> {
        let mut count: u64 = 0;
        let mut sum_len: u64 = 0;
        let mut record = Record::default();
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
    pub fn read_all(&mut self) -> io::Result<Vec<NamedSeq>> {
        let mut record = Record::default();
        let mut records = Vec::new();
        while self.read_next(&mut record)? {
            records.push(NamedSeq::new(record.name_only().into_owned(), record.seq().to_owned()));
        }
        Ok(records)
    }
}

/// Trait for reading and writing FASTA/Q single-end and paired-end records.
pub trait FastxRead: Send {
    type Record: RecordExt;

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self, record: &mut Self::Record) -> io::Result<bool>;

    /// Writes the next `count` records to the output stream.
    /// Returns the number of records that were written.
    fn write_next_n(&mut self, writer: &mut impl io::Write, count: usize) -> io::Result<usize> {
        let mut record = Self::Record::default();
        for i in 0..count {
            match self.read_next(&mut record)? {
                true => record.write_to(writer)?,
                false => return Ok(i),
            }
        }
        Ok(count)
    }

    /// Subsamples the reader with the `rate` and optional `seed`.
    fn subsample(&mut self, writer: &mut impl io::Write, rate: f64, rng: &mut impl rand::Rng) -> io::Result<()> {
        assert!(rate > 0.0 && rate < 1.0, "Subsampling rate ({}) must be within (0, 1).", rate);
        let mut record = Self::Record::default();
        while self.read_next(&mut record)? {
            if rng.gen::<f64>() <= rate {
                record.write_to(writer)?;
            }
        }
        Ok(())
    }

    /// Writes input stream to output.
    fn copy(&mut self, writer: &mut impl io::Write) -> io::Result<()> {
        let mut record = Self::Record::default();
        while self.read_next(&mut record)? {
            record.write_to(writer)?;
        }
        Ok(())
    }

    /// Returns the total number of consumed reads/read pairs.
    fn total_reads(&self) -> u64;
}

impl<R: BufRead + Send> FastxRead for Reader<R> {
    type Record = Record;

    /// Reads the next record, and returns true if the read was successful (false if no more reads available).
    fn read_next(&mut self, record: &mut Record) -> io::Result<bool> {
        record.full_name.clear();
        record.seq.clear();
        record.qual.clear();

        // Buffer already contains the first line of the next record.
        // If it is empty, the file has ended.
        if self.buffer.is_empty() {
            return Ok(false);
        }
        self.total_reads += 1;

        record.full_name.extend_from_slice(&self.buffer[1..]);
        if self.buffer[0] == b'>' {
            self.fill_fasta_record(record)
        } else {
            self.fill_fastq_record(record)
        }
    }

    /// Returns the total number of single end reads.
    fn total_reads(&self) -> u64 {
        self.total_reads
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
    type Record = PairedRecord;

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self, paired_record: &mut PairedRecord) -> io::Result<bool> {
        if !self.reader.read_next(&mut paired_record[0])? {
            Ok(false)
        } else if !self.reader.read_next(&mut paired_record[1])? {
            Err(io::Error::new(io::ErrorKind::InvalidData,
                "Odd number of records in an interleaved input file."))
        } else {
            Ok(true)
        }
    }

    /// Returns the total number of consumed paired reads.
    fn total_reads(&self) -> u64 {
        self.reader.total_reads() / 2
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
    type Record = PairedRecord;

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self, paired_record: &mut PairedRecord) -> io::Result<bool> {
        let could_read1 = self.reader1.read_next(&mut paired_record[0])?;
        let could_read2 = self.reader2.read_next(&mut paired_record[1])?;
        match (could_read1, could_read2) {
            (false, false) => Ok(false),
            (true, true) => Ok(true),
            _ => Err(io::Error::new(io::ErrorKind::InvalidData,
                "Different number of records in two input files.")),
        }
    }

    /// Returns the total number of consumed paired reads.
    fn total_reads(&self) -> u64 {
        self.reader1.total_reads()
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
