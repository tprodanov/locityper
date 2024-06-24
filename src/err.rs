use std::{
    io,
    fmt::Write,
    path::PathBuf,
};
use colored::Colorize;
use crate::ext;

/// General enum, representing possible errors.
#[derive(Debug)]
pub enum Error {
    Io(io::Error, Vec<PathBuf>),
    /// Solver finished with an error: `(solver_name, error_description)`.
    Solver(&'static str, String),
    /// Error, produced by an argument parser.
    Lexopt(lexopt::Error),
    /// rust_htslib error.
    Htslib(htslib::errors::Error),
    /// Executable not found.
    NoExec(PathBuf),
    /// Subcommand failed.
    Subprocess(String),
    InvalidInput(String),
    InvalidData(String),
    ParsingError(String),
    RuntimeError(String),
    JsonLoad(String),
    /// Non-UTF8 encountered. First field: what is the object (for example, `read name`), second: failed byte string.
    Utf8(&'static str, Vec<u8>)
}

#[cfg(feature = "gurobi")]
impl From<grb::Error> for Error {
    fn from(e: grb::Error) -> Self {
        Self::Solver("Gurobi", e.to_string())
    }
}

impl From<lexopt::Error> for Error {
    fn from(e: lexopt::Error) -> Self {
        Self::Lexopt(e)
    }
}

impl From<htslib::errors::Error> for Error {
    fn from(e: htslib::errors::Error) -> Self {
        Self::Htslib(e)
    }
}

impl From<json::JsonError> for Error {
    fn from(e: json::JsonError) -> Self {
        Self::JsonLoad(e.to_string())
    }
}

impl Error {
    /// Converts an error, produced by a solver.
    pub fn solver(solver_name: &'static str, s: impl Into<String>) -> Self {
        Self::Solver(solver_name, s.into())
    }

    /// Format error message. TODO: Better formatting.
    pub fn display(&self) -> String {
        let mut s = String::new();
        match self {
            Self::Io(e, files) => {
                write!(s, "{} in relation to ", "Input/Output error".red()).unwrap();
                if files.is_empty() {
                    write!(s, "unnamed streams").unwrap();
                } else {
                    write!(s, "{}", files.iter().map(|f| ext::fmt::path(f).cyan().to_string())
                        .collect::<Vec<_>>().join(", ")).unwrap();
                }
                write!(s, ": {}", e.kind()).unwrap();
                if let Some(e2) = e.get_ref() {
                    write!(s, ", {}", e2).unwrap();
                }
            }
            Self::Solver(solver, e) => write!(s, "ILP solver {} produced error {}", solver.red(), e).unwrap(),
            Self::Lexopt(e) => write!(s, "{} to parse command-line arguments: {}", "Failed".red(), e).unwrap(),
            Self::Htslib(e) => write!(s, "{}: {}", "Htslib error".red(), e).unwrap(),
            Self::NoExec(path) => write!(s, "{} at {}", "Could not find executable".red(),
                ext::fmt::path(path).cyan()).unwrap(),
            Self::Subprocess(e) => write!(s, "{}:\n{}", "Subprocess error".red(), e).unwrap(),
            Self::InvalidInput(e) => write!(s, "{}: {}", "Invalid input".red(), e).unwrap(),
            Self::InvalidData(e) => write!(s, "{}: {}", "Invalid data".red(), e).unwrap(),
            Self::ParsingError(e) => write!(s, "{}: {}", "Parsing error".red(), e).unwrap(),
            Self::RuntimeError(e) => write!(s, "{}: {}", "Runtime error".red(), e).unwrap(),
            Self::JsonLoad(e) => write!(s, "{}: {}", "Could not load JSON".red(), e).unwrap(),
            Self::Utf8(desc, bytes) => write!(s, "{}: {} cannot be decoded ({}, {:?})",
                "UTF-8 error".red(), desc, String::from_utf8_lossy(bytes), bytes).unwrap(),
        };
        s
    }

    /// Converts this error back to io::Error.
    /// If this is already an io::Error, adds information about the filename to the end, otherwise panics.
    pub fn try_into_io_error(self) -> io::Error {
        match self {
            Self::Io(e, files) =>
                io::Error::new(e.kind(), format!("{:?} in relation to {}",
                    e, if files.is_empty() { "unnamed files".to_string() } else { ext::fmt::paths(&files) })),
            e => {
                log::error!("{}", e.display());
                panic!("Cannot convert non-IO error to IO");
            }
        }
    }
}

macro_rules! validate_param {
    ($cond:expr, $($arg:expr),+) => {{
        if !($cond) {
            (
                Err($crate::Error::InvalidInput(format!($($arg),+)))
            ?)
        }
    }};
}
pub(crate) use validate_param;

macro_rules! add_path {
    (!) => {
        |e| $crate::Error::Io(e, Vec::new())
    };
    ($path:expr) => {
        |e| $crate::Error::Io(e, vec![std::convert::AsRef::<std::path::Path>::as_ref(&$path).to_owned()])
    };
    ($($path:expr),+) => {
        |e| {
            let mut v = Vec::new();
            $(
                v.push(std::convert::AsRef::<std::path::Path>::as_ref(&$path).to_owned());
            )*
            $crate::Error::Io(e, v)
        }
    };
}
pub(crate) use add_path;

macro_rules! error {
    ($var:ident, $($arg:expr),+ $(,)?) => {
        $crate::Error::$var(format!($($arg),+))
    }
}
pub(crate) use error;

/// Wrapper around the standard result.
pub type Result<T> = std::result::Result<T, Error>;
