use std::{
    io,
    cmp::min,
};
use bio::io::{fasta, fastq};
use rand::{
    SeedableRng,
    rngs::SmallRng,
    distributions::{Bernoulli, Distribution, DistIter},
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

/// Structure, that subsamples single- and paired-end reads.
pub struct Subsample {
    successes: DistIter<Bernoulli, SmallRng, bool>,
}

impl Subsample {
    /// Creates a new subsampling structure with a given rate and optional seed.
    pub fn new(rate: f64, seed: Option<u64>) -> Self {
        assert!(rate > 0.0 && rate < 1.0, "Subsampling rate must be within (0, 1).");
        let rng = if let Some(seed) = seed {
            if seed.count_ones() < 5 {
                log::warn!("Seed ({}) is too simple, consider using a more random number.", seed);
            }
            SmallRng::seed_from_u64(seed)
        } else {
            SmallRng::from_entropy()
        };
        Self {
            successes: Bernoulli::new(rate).unwrap().sample_iter(rng),
        }
    }

    /// Subsamples one single-end FASTA/Q file.
    pub fn single_end<R: FastxReader>(mut self, mut reader: R, mut writer: impl io::Write) -> io::Result<()> {
        let mut record = R::Record::new();
        loop {
            reader.read_next(&mut record)?;
            if record.is_empty() {
                return writer.flush();
            }
            if self.successes.next().unwrap() {
                record.write_simple(&mut writer)?;
            }
        }
    }

    /// Subsamples an interleaved paired-end FASTA/Q file.
    pub fn interleaved_paired_end<R: FastxReader>(
        mut self,
        mut reader: R,
        mut writer: impl io::Write,
    ) -> io::Result<()>
    {
        let mut record1 = R::Record::new();
        let mut record2 = R::Record::new();
        let mut had_warning = false;
        loop {
            reader.read_next(&mut record1)?;
            if record1.is_empty() {
                return writer.flush();
            }
            reader.read_next(&mut record2)?;
            if record2.is_empty() {
                return Err(io::Error::new(io::ErrorKind::InvalidData,
                    "Odd number of records in an interleaved input file."))
            }
            if !had_warning && record1.name() != record2.name() {
                had_warning = true;
                log::warn!("Interleaved input file contains consecutive records with different names: {} and {}",
                    record1.name(), record2.name());
            }
            if self.successes.next().unwrap() {
                record1.write_simple(&mut writer)?;
                record2.write_simple(&mut writer)?;
            }
        }
    }

    /// Subsamples two paired-end FASTA/Q file into a single interleaved file.
    pub fn paired_end<R: FastxReader, S: FastxReader>(
        mut self,
        mut reader1: R,
        mut reader2: S,
        mut writer: impl io::Write,
    ) -> io::Result<()>
    {
        let mut record1 = R::Record::new();
        let mut record2 = S::Record::new();
        let mut had_warning = false;
        loop {
            reader1.read_next(&mut record1)?;
            reader2.read_next(&mut record2)?;
            match (record1.is_empty(), record2.is_empty()) {
                (true, true) => return writer.flush(),
                (false, false) => {}
                _ => return Err(io::Error::new(io::ErrorKind::InvalidData,
                    "Different number of records in two input files."))
            }
            if !had_warning && record1.name() != record2.name() {
                had_warning = true;
                log::warn!("Paired-end records have different names: {} and {}",
                    record1.name(), record2.name());
            }
            if self.successes.next().unwrap() {
                record1.write_simple(&mut writer)?;
                record2.write_simple(&mut writer)?;
            }
        }
    }
}
