use std::{
    io,
    cmp::min,
};
use bio::io::{fasta, fastq};
use rand::{
    SeedableRng,
    rngs::SmallRng,
    distributions::{Bernoulli, Distribution},
};

/// Trait, summarizing `fasta::Record` and `fastq::Record`.
pub trait FastxRecord: Sized {
    /// Creates a new, empty record.
    fn new() -> Self;

    /// Checks if the record is empty.
    fn is_empty(&self) -> bool;

    /// Returns the name of the record.
    fn name(&self) -> &str;

    /// Returns the record sequence.
    fn seq(&self) -> &[u8];

    /// Writes one record to the output stream.
    /// For FASTA files, this function splits sequence into multiple lines.
    fn write<W: io::Write>(&self, writer: W) -> io::Result<()>;

    /// Writes one record to the output stream.
    /// Do not output any description, and for FASTA files, write sequence into one line.
    fn write_simple<W: io::Write>(&self, writer: W) -> io::Result<()>;
}

/// Trait, summarizing `fasta::Reader` and `fastq::Reader`.
pub trait FastxReader {
    type Record: FastxRecord;

    /// Reads the next record.
    /// If there is no record left, `record` will be empty.
    fn read_next(&mut self, record: &mut Self::Record) -> io::Result<()>;
}

/// Implement `FastxRecord` for FASTA record.
impl FastxRecord for fasta::Record {
    fn new() -> Self {
        Self::new()
    }

    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    /// Returns the name of the record.
    fn name(&self) -> &str {
        self.id()
    }

    /// Returns the record sequence.
    fn seq(&self) -> &[u8] {
        self.seq()
    }

    fn write<W: io::Write>(&self, writer: W) -> io::Result<()> {
        write_fasta(writer, self.id(), self.desc(), self.seq())
    }

    fn write_simple<W: io::Write>(&self, mut writer: W) -> io::Result<()> {
        write!(writer, ">{}\n", self.id())?;
        writer.write_all(self.seq())?;
        writer.write_all(b"\n")
    }
}


/// Write a single sequence to the FASTA file.
/// Use this function instead of `bio::fasta::Writer` as the latter
/// writes the sequence into a single line, without splitting.
pub fn write_fasta<W: io::Write>(mut writer: W, name: &str, desc: Option<&str>, seq: &[u8]) -> io::Result<()> {
    write!(writer, ">{}", name)?;
    if let Some(desc) = desc {
        write!(writer, " {}", desc)?;
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

/// Implement `FastxReader` for FASTA reader.
impl<R: io::BufRead> FastxReader for fasta::Reader<R> {
    type Record = fasta::Record;

    fn read_next(&mut self, record: &mut Self::Record) -> io::Result<()> {
        fasta::FastaRead::read(self, record)
    }
}

/// Implement `FastxRecord` for FASTQ record.
impl FastxRecord for fastq::Record {
    fn new() -> Self {
        Self::new()
    }

    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    /// Returns the name of the record.
    fn name(&self) -> &str {
        self.id()
    }

    /// Returns the record sequence.
    fn seq(&self) -> &[u8] {
        self.seq()
    }

    fn write<W: io::Write>(&self, mut writer: W) -> io::Result<()> {
        write!(writer, "@{}", self.id())?;
        if let Some(desc) = self.desc() {
            write!(writer, " {}", desc)?;
        }
        writer.write_all(b"\n")?;
        writer.write_all(self.seq())?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(self.qual())?;
        writer.write_all(b"\n")?;
        Ok(())
    }

    fn write_simple<W: io::Write>(&self, mut writer: W) -> io::Result<()> {
        write!(writer, "@{}", self.id())?;
        writer.write_all(b"\n")?;
        writer.write_all(self.seq())?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(self.qual())?;
        writer.write_all(b"\n")?;
        Ok(())
    }
}

/// Implement `FastxReader` for FASTQ reader.
impl<R: io::BufRead> FastxReader for fastq::Reader<R> {
    type Record = fastq::Record;

    fn read_next(&mut self, record: &mut Self::Record) -> io::Result<()> {
        fastq::FastqRead::read(self, record).map_err(|e| match e {
            fastq::Error::FileOpen { source, .. } => source,
            fastq::Error::ReadError(err) => err,
            _ => io::Error::new(io::ErrorKind::InvalidData, "Failed to process FASTQ file"),
        })
    }
}

/// Trait for reading and writing FASTA/Q single-end and paired-end records.
pub trait ReaderWriter: Sized {
    type Record;

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self) -> io::Result<bool>;

    /// Returns the current record.
    fn current(&self) -> &Self::Record;

    /// Writes the current record to the output stream.
    fn write_current(&mut self) -> io::Result<()>;

    /// Writes the next `count` records to the output stream.
    /// Returns the number of records that were written.
    fn write_next_n(&mut self, count: usize) -> io::Result<usize> {
        for i in 0..count {
            match self.read_next()? {
                true => self.write_current()?,
                false => return Ok(i),
            }
        }
        Ok(count)
    }

    /// Subsamples the reader with the `rate` and optional `seed`.
    fn subsample(mut self, rate: f64, seed: Option<u64>) -> io::Result<()> {
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
                self.write_current()?;
            }
        }
        Ok(())
    }
}

/// Single-end FASTA/Q reader, that stores one buffer record to reduce memory allocations.
pub struct SingleEndReader<R: FastxReader, W> {
    reader: R,
    writer: W,
    record: R::Record,
}

impl<R: FastxReader, W> SingleEndReader<R, W> {
    pub fn new(reader: R, writer: W) -> Self {
        Self {
            reader, writer,
            record: R::Record::new(),
        }
    }
}

impl<R: FastxReader, W: io::Write> ReaderWriter for SingleEndReader<R, W> {
    type Record = R::Record;

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self) -> io::Result<bool> {
        self.reader.read_next(&mut self.record)?;
        Ok(!self.record.is_empty())
    }

    /// Returns the current record.
    fn current(&self) -> &Self::Record {
        &self.record
    }

    /// Writes the current record to the output stream.
    fn write_current(&mut self) -> io::Result<()> {
        self.record.write_simple(&mut self.writer)
    }
}

/// Interleaved paired-end FASTA/Q reader, that stores two buffer records to reduce memory allocations.
pub struct PairedEndInterleaved<R: FastxReader, W> {
    reader: R,
    writer: W,
    records: [R::Record; 2],
    had_warning: bool,
}

impl<R: FastxReader, W> PairedEndInterleaved<R, W> {
    pub fn new(reader: R, writer: W) -> Self {
        Self {
            reader, writer,
            records: [R::Record::new(), R::Record::new()],
            had_warning: false,
        }
    }
}

impl<R: FastxReader, W: io::Write> ReaderWriter for PairedEndInterleaved<R, W> {
    type Record = [R::Record; 2];

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self) -> io::Result<bool> {
        self.reader.read_next(&mut self.records[0])?;
        if self.records[0].is_empty() {
            return Ok(false);
        }
        self.reader.read_next(&mut self.records[1])?;
        if self.records[1].is_empty() {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                "Odd number of records in an interleaved input file."))
        }
        if !self.had_warning && self.records[0].name() != self.records[1].name() {
            self.had_warning = true;
            log::warn!("Interleaved input file contains consecutive records with different names: {} and {}",
                self.records[0].name(), self.records[1].name());
        }
        Ok(true)
    }

    /// Returns the current record.
    fn current(&self) -> &Self::Record {
        &self.records
    }

    /// Writes the current record to the output stream.
    fn write_current(&mut self) -> io::Result<()> {
        self.records[0].write_simple(&mut self.writer)?;
        self.records[1].write_simple(&mut self.writer)
    }
}

/// Two paired-end FASTA/Q readers, that stores two buffer records to reduce memory allocations.
pub struct PairedEndReaders<R: FastxReader, W> {
    readers: [R; 2],
    writer: W,
    records: [R::Record; 2],
    had_warning: bool,
}

impl<R: FastxReader, W> PairedEndReaders<R, W> {
    pub fn new(reader1: R, reader2: R, writer: W) -> Self {
        Self {
            readers: [reader1, reader2],
            writer,
            records: [R::Record::new(), R::Record::new()],
            had_warning: false,
        }
    }
}

impl<R: FastxReader, W: io::Write> ReaderWriter for PairedEndReaders<R, W> {
    type Record = [R::Record; 2];

    /// Read next one/two records, and return true if the read was filled (is not empty).
    fn read_next(&mut self) -> io::Result<bool> {
        self.readers[0].read_next(&mut self.records[0])?;
        self.readers[1].read_next(&mut self.records[1])?;
        match (self.records[0].is_empty(), self.records[1].is_empty()) {
            (true, true) => return Ok(false),
            (false, false) => {}
            _ => return Err(io::Error::new(io::ErrorKind::InvalidData,
                "Different number of records in two input files."))
        }
        if !self.had_warning && self.records[0].name() != self.records[1].name() {
            self.had_warning = true;
            log::warn!("Paired-end records have different names: {} and {}",
            self.records[0].name(), self.records[1].name());
        }
        Ok(true)
    }

    /// Returns the current record.
    fn current(&self) -> &Self::Record {
        &self.records
    }

    /// Writes the current record to the output stream.
    fn write_current(&mut self) -> io::Result<()> {
        self.records[0].write_simple(&mut self.writer)?;
        self.records[1].write_simple(&mut self.writer)
    }
}
