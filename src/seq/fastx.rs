use std::{
    io,
    cmp::min,
};
use bio::io::{fasta, fastq};

/// Trait, summarizing `fasta::Record` and `fastq::Record`.
pub trait FastxRecord: Sized {
    /// Creates a new, empty record.
    fn new() -> Self;

    /// Checks if the record is empty.
    fn is_empty(&self) -> bool;

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