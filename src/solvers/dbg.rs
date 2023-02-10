use std::{
    io::{self, Write},
    fmt,
    fs::File,
};
use crate::{
    // algo::vec_ext::IterExt,
    model::{
        // locs::PairAlignment,
        assgn::ReadAssignment,
    },
};

/// Enum that denotes two-order iteration number.
#[derive(Clone)]
pub enum Iteration {
    /// Initialization was performed. Argument: outer iteration.
    Init(u32),
    /// Single step was performed. Arguments: outer iteration and inner iteration.
    Step(u32, u32),
    /// Last step was performed. Arguments: outer iteration and last inner iteration.
    Last(u32, u32),
    /// Best solution across all iterations.
    Best,
}

impl fmt::Display for Iteration {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Iteration::Init(outer) => write!(f, "{}-init", outer),
            Iteration::Step(outer, inner) => write!(f, "{}-{}", outer, inner),
            Iteration::Last(outer, inner) => write!(f, "{}-{}!", outer, inner),
            Iteration::Best => write!(f, "best"),
        }
    }
}

/// Write debug information about solvers.
pub trait DbgWrite {
    /// Write current read assignments after the corresponding iteration.
    fn write(&mut self, read_assignments: &ReadAssignment, iteration: Iteration) -> io::Result<()>;

    /// Flush cached data to the writer.
    fn flush(&mut self) -> io::Result<()>;
}

/// Do not output any debug information.
pub struct NoDbg;

impl DbgWrite for NoDbg {
    fn write(&mut self, _: &ReadAssignment, _: Iteration) -> io::Result<()> { Ok(()) }

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
    fn write(&mut self, read_assignments: &ReadAssignment, iteration: Iteration) -> io::Result<()> {
        writeln!(self.writer, "{}{}\t-1\tNA\t{:.2}", self.prefix, iteration, read_assignments.likelihood())?;
        match &iteration {
            Iteration::Step(_, inner) if inner % self.frequency != 0 => return Ok(()),
            _ => {},
        };

        log::debug!("{:10}:   likelihood {:12.4}", iteration.to_string(), read_assignments.likelihood());
        for (w, (depth, lik)) in read_assignments.depth_lik_iter().enumerate() {
            writeln!(self.writer, "{}{}\t{}\t{}\t{:.2}", self.prefix, iteration, w, depth, lik)?;
        }
        Ok(())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}
