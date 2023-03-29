use std::io;

/// General enum, representing possible errors.
#[derive(Debug)]
pub enum Error {
    Io(io::Error),
    /// Solver finished with an error: `(solver_name, error_description)`.
    Solver(&'static str, String),
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

#[cfg(feature = "gurobi")]
impl From<grb::Error> for Error {
    fn from(e: grb::Error) -> Self {
        Self::Solver("Gurobi", e.to_string())
    }
}

impl Error {
    pub fn solver(solver_name: &'static str, s: impl Into<String>) -> Self {
        Self::Solver(solver_name, s.into())
    }
}
