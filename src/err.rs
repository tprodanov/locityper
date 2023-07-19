use std::{
    io,
    path::{Path, PathBuf},
};

/// General enum, representing possible errors.
#[derive(Debug)]
pub enum Error {
    Io(io::Error, Option<PathBuf>),
    /// Solver finished with an error: `(solver_name, error_description)`.
    Solver(&'static str, String),
    /// Error, produced by an argument parser.
    Lexopt(lexopt::Error),
    /// rust_htslib error.
    Htslib(htslib::errors::Error),
    /// Executable not found.
    NoExec(PathBuf),
    /// Subcommand failed.
    SubprocessFail(std::process::Output),
    InvalidInput(String),
    InvalidData(String),
    ParsingError(String),
    RuntimeError(String),
    JsonLoad(String),
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Self::Io(e, None)
    }
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

impl From<std::str::Utf8Error> for Error {
    fn from(_: std::str::Utf8Error) -> Self {
        Self::ParsingError("Failed to parse a string: not a valid UTF-8".to_owned())
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

    pub fn io(e: io::Error, filename: impl AsRef<Path>) -> Self {
        Self::Io(e, Some(filename.as_ref().to_owned()))
    }

    pub fn io_unkn(e: io::Error) -> Self {
        Self::Io(e, None)
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
        |e| $crate::Error::Io(e, None)
    };
    ($path:expr) => {
        |e| $crate::Error::Io(e, Some(std::convert::AsRef::<std::path::Path>::as_ref(&$path).to_owned()))
    };
}
pub(crate) use add_path;
