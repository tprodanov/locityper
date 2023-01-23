use std::io::{self, Write};
use std::fs::File;
use crate::{
    // algo::vec_ext::IterExt,
    model::{
        // locs::PairAlignment,
        assgn::ReadAssignment,
    },
};

/// Write debug information about solvers.
pub trait DbgWrite {
    /// Write current read assignments after the corresponding iteration.
    /// If force is true, write in any case.
    fn write<'a>(&mut self, read_assignments: &'a ReadAssignment<'a>, iteration: u32, force: bool) -> io::Result<()>;

    fn flush(&mut self) -> io::Result<()>;
}

/// Do not output any debug information.
pub struct NoDbg;

impl DbgWrite for NoDbg {
    fn write<'a>(&mut self, _: &'a ReadAssignment<'a>, _: u32, _: bool) -> io::Result<()> { Ok(()) }

    fn flush(&mut self) -> io::Result<()> { Ok(()) }
}

/// Debug writer that writes all iterations that are divisible by frequency.
/// Each line starts with prefix. Next, iteration, window index, read depth and likelihood are written.
pub struct DbgWriter<W: Write> {
    writer: W,
    prefix: String,
    frequency: u32,
}

impl<W: Write> DbgWriter<W> {
    pub fn new(writer: W, prefix: String, frequency: u32) -> Self {
        Self { writer, prefix, frequency }
    }
}

impl DbgWriter<io::BufWriter<File>> {
    /// Create debug writer from path, prefix, and frequency.
    pub fn from_path<P: AsRef<std::path::Path>>(path: P, prefix: String, frequency: u32) -> io::Result<Self> {
        let writer = io::BufWriter::new(File::create(path)?);
        Ok(Self::new(writer, prefix, frequency))
    }
}

impl<W: Write> DbgWrite for DbgWriter<W> {
    fn write<'a>(&mut self, read_assignments: &'a ReadAssignment<'a>, iteration: u32, force: bool) -> io::Result<()> {
        if force || iteration % self.frequency == 0 {
            log::debug!("Iteration {:5}:   likelihood {:12.4}", iteration, read_assignments.likelihood());
            writeln!(self.writer, "{}{}\t-1\tNA\t{:.2}", self.prefix, iteration, read_assignments.likelihood())?;
            for (w, (depth, lik)) in read_assignments.depth_lik_iter().enumerate() {
                writeln!(self.writer, "{}{}\t{}\t{}\t{:.2}", self.prefix, iteration, w, depth, lik)?;
            }
        }
        Ok(())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}
