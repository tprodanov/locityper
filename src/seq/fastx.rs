use std::{
    io::{self, BufRead},
    fmt,
    cmp::min,
};
use rand::{
    SeedableRng,
    rngs::SmallRng,
    distributions::{Bernoulli, Distribution},
};

/// Write a single sequence to the FASTA file.
/// Use this function instead of `bio::fasta::Writer` as the latter
/// writes the sequence into a single line, without splitting.
pub fn write_fasta<W: io::Write>(mut writer: W, name: &[u8], descr: Option<&[u8]>, seq: &[u8]) -> io::Result<()> {
    writer.write_all(b">")?;
    writer.write_all(name)?;
    if let Some(descr) = descr {
        writer.write_all(b" ")?;
        writer.write_all(descr)?;
    }
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
pub fn write_fastq<W: io::Write>(
    mut writer: W,
    name: &[u8],
    descr: Option<&[u8]>,
    seq: &[u8],
    qual: &[u8],
) -> io::Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(name)?;
    if let Some(descr) = descr {
        writer.write_all(b" ")?;
        writer.write_all(descr)?;
    }
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
    name: Vec<u8>,
    descr: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl Record {
    // /// Returns true if the record is empty.
    // pub fn is_empty(&self) -> bool {
    //     self.name.is_empty()
    // }

    /// Returns the name of the record.
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    pub fn name_str(&self) -> std::borrow::Cow<'_, str> {
        String::from_utf8_lossy(&self.name)
    }

    /// Returns record description, if available.
    pub fn descr(&self) -> Option<&[u8]> {
        if self.descr.is_empty() {
            None
        } else {
            Some(&self.descr)
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

    pub fn write_to(&self, writer: impl io::Write) -> io::Result<()> {
        if self.qual.is_empty() {
            write_fasta(writer, &self.name, self.descr(), &self.seq)
        } else {
            write_fastq(writer, &self.name, self.descr(), &self.seq, &self.qual)
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

/// Dynamic FASTA/Q reader.
pub struct Reader<R: BufRead> {
    stream: R,
    buffer: Vec<u8>,
}

impl<R: BufRead> Reader<R> {
    pub fn new(mut stream: R) -> io::Result<Self> {
        let mut buffer = Vec::new();
        read_line(&mut stream, &mut buffer)?;
        Ok(Self { stream, buffer })
    }

    /// Reads the next record, and returns true if the read was successful (false if no more reads available).
    pub fn read_next(&mut self, record: &mut Record) -> io::Result<bool> {
        record.name.clear();
        record.descr.clear();
        record.seq.clear();
        record.qual.clear();

        // Buffer already contains the first line of the next record.
        // If it is empty, the file has ended.
        if self.buffer.is_empty() {
            return Ok(false);
        }

        match self.buffer[1..].iter().position(|&c| c == b' ') {
            Some(i) => {
                record.name.extend(&self.buffer[1..i]);
                record.descr.extend(&self.buffer[i + 1..]);
            }
            None => record.name.extend(&self.buffer[1..]),
        }

        if self.buffer[0] == b'>' {
            self.fill_fasta_record(record)
        } else {
            self.fill_fastq_record(record)
        }
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
                        format!("Fasta record {} has an empty sequence.", record.name_str())));
                }
                return Ok(true);
            } else if record.seq[seq_len] == b'@' || record.seq[seq_len] == b'>' {
                self.buffer.extend(&record.seq[seq_len..]);
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
                format!("Fastq record {} is incomplete.", record.name_str())));
        } else if self.buffer[prev_buf_len] != b'+' {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Fastq record {} has incorrect format.", record.name_str())));
        }

        // Qualities
        let qual_len = read_line(&mut self.stream, &mut record.qual)?;
        if seq_len != qual_len {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Fastq record {} has non-matching sequence and qualities.", record.name_str())));
        }

        // Next record name.
        self.buffer.clear();
        read_line(&mut self.stream, &mut self.buffer)?;
        Ok(true)
    }
}

/// Trait for reading and writing FASTA/Q single-end and paired-end records.
pub trait FastxRead: Send {
    type Record;

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self) -> io::Result<bool>;

    /// Returns the current record.
    fn current(&self) -> &Self::Record;

    /// Writes the current record to the output stream.
    fn write_current(&mut self, writer: impl io::Write) -> io::Result<()>;

    /// Writes the next `count` records to the output stream.
    /// Returns the number of records that were written.
    fn write_next_n(&mut self, mut writer: impl io::Write, count: usize) -> io::Result<usize> {
        for i in 0..count {
            match self.read_next()? {
                true => self.write_current(&mut writer)?,
                false => return Ok(i),
            }
        }
        Ok(count)
    }

    /// Subsamples the reader with the `rate` and optional `seed`.
    fn subsample(&mut self, mut writer: impl io::Write, rate: f64, seed: Option<u64>) -> io::Result<()> {
        assert!(rate > 0.0 && rate < 1.0, "Subsampling rate must be within (0, 1).");
        let rng = if let Some(seed) = seed {
            if seed.count_ones() < 5 {
                log::warn!("Seed ({}) is too simple, consider using a more random number.", seed);
            }
            SmallRng::seed_from_u64(seed)
        } else {
            SmallRng::from_entropy()
        };
        let mut successes = Bernoulli::new(rate).unwrap().sample_iter(rng);

        while self.read_next()? {
            if successes.next().unwrap() {
                self.write_current(&mut writer)?;
            }
        }
        Ok(())
    }
}

/// Single-end FASTA/Q reader, that stores one buffer record to reduce memory allocations.
pub struct SingleEndReader<R: BufRead> {
    reader: Reader<R>,
    record: Record,
}

impl<R: BufRead> SingleEndReader<R> {
    pub fn new(reader: Reader<R>) -> Self {
        Self {
            reader,
            record: Record::default(),
        }
    }
}

impl<R: BufRead + Send> FastxRead for SingleEndReader<R> {
    type Record = Record;

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self) -> io::Result<bool> {
        self.reader.read_next(&mut self.record)
    }

    /// Returns the current record.
    fn current(&self) -> &Self::Record {
        &self.record
    }

    /// Writes the current record to the output stream.
    fn write_current(&mut self, writer: impl io::Write) -> io::Result<()> {
        self.record.write_to(writer)
    }
}

pub type PairedRecord = [Record; 2];

/// Interleaved paired-end FASTA/Q reader, that stores two buffer records to reduce memory allocations.
pub struct PairedEndInterleaved<R: BufRead> {
    reader: Reader<R>,
    records: PairedRecord,
    had_warning: bool,
}

impl<R: BufRead> PairedEndInterleaved<R> {
    pub fn new(reader: Reader<R>) -> Self {
        Self {
            reader,
            records: Default::default(),
            had_warning: false,
        }
    }
}

impl<R: BufRead + Send> FastxRead for PairedEndInterleaved<R> {
    type Record = PairedRecord;

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self) -> io::Result<bool> {
        if !self.reader.read_next(&mut self.records[0])? {
            return Ok(false);
        }
        if !self.reader.read_next(&mut self.records[1])? {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                "Odd number of records in an interleaved input file."))
        }
        if !self.had_warning && self.records[0].name() != self.records[1].name() {
            self.had_warning = true;
            log::warn!("Interleaved input file contains consecutive records with different names: {} and {}",
                self.records[0].name_str(), self.records[1].name_str());
        }
        Ok(true)
    }

    /// Returns the current record.
    fn current(&self) -> &Self::Record {
        &self.records
    }

    /// Writes the current record to the output stream.
    fn write_current(&mut self, mut writer: impl io::Write) -> io::Result<()> {
        self.records[0].write_to(&mut writer)?;
        self.records[1].write_to(&mut writer)
    }
}

/// Two paired-end FASTA/Q readers, that stores two buffer records to reduce memory allocations.
pub struct PairedEndReaders<R: BufRead, S: BufRead> {
    reader1: Reader<R>,
    reader2: Reader<S>,
    records: PairedRecord,
    had_warning: bool,
}

impl<R: BufRead, S: BufRead> PairedEndReaders<R, S> {
    pub fn new(reader1: Reader<R>, reader2: Reader<S>) -> Self {
        Self {
            reader1, reader2,
            records: Default::default(),
            had_warning: false,
        }
    }
}

impl<R: BufRead + Send, S: BufRead + Send> FastxRead for PairedEndReaders<R, S> {
    type Record = PairedRecord;

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self) -> io::Result<bool> {
        let could_read1 = self.reader1.read_next(&mut self.records[0])?;
        let could_read2 = self.reader2.read_next(&mut self.records[1])?;
        match (could_read1, could_read2) {
            (false, false) => return Ok(false),
            (true, true) => {}
            _ => return Err(io::Error::new(io::ErrorKind::InvalidData,
                "Different number of records in two input files."))
        }
        if !self.had_warning && self.records[0].name() != self.records[1].name() {
            self.had_warning = true;
            log::warn!("Paired-end records have different names: {} and {}",
                self.records[0].name_str(), self.records[1].name_str());
        }
        Ok(true)
    }

    /// Returns the current record.
    fn current(&self) -> &Self::Record {
        &self.records
    }

    /// Writes the current record to the output stream.
    fn write_current(&mut self, mut writer: impl io::Write) -> io::Result<()> {
        self.records[0].write_to(&mut writer)?;
        self.records[1].write_to(&mut writer)
    }
}
