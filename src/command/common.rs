use std::{
    io::{BufWriter, Write, BufReader, BufRead, stdin, stdout},
    fs::File,
    path::{Path, PathBuf},
    ffi::OsStr,
};
use colored::Colorize;
use crate::Error;

/// Finds an executable, and returns Error, if executable is not available.
pub(super) fn find_exe(p: PathBuf) -> Result<PathBuf, Error> {
    which::which(&p).map_err(|_| Error::NoExec(p))
}

/// Print tool version and authors.
pub(super) fn print_version() {
    println!("{} {}", env!("CARGO_PKG_NAME").underline(), format!("v{}", env!("CARGO_PKG_VERSION")).green());
    let authors: Vec<_> = env!("CARGO_PKG_AUTHORS").split(':').collect();
    let n = authors.len();
    if n == 0 {
        return;
    }
    print!("Created by ");
    for (i, author) in authors.iter().enumerate() {
        if i == 0 {
            print!("{}", author.bright_blue());
        } else if i < n - 1 {
            print!(", {}", author.bright_blue());
        } else {
            print!(" and {}", author.bright_blue());
        }
    }
    println!();
}

/// Returns a buffered file OR stdin if filename is `"-"`.
pub(super) fn file_or_stdin(filename: &Path) -> Result<Box<dyn BufRead>, Error> {
    if filename == OsStr::new("-") {
        Ok(Box::new(BufReader::new(stdin())))
    } else {
        Ok(Box::new(BufReader::new(File::open(filename)?)))
    }
}

/// Returns a buffered file OR stdout if filename is `"-"`.
pub(super) fn file_or_stdout(filename: &Path) -> Result<Box<dyn Write>, Error> {
    if filename == OsStr::new("-") {
        Ok(Box::new(BufWriter::new(stdout())))
    } else {
        Ok(Box::new(BufWriter::new(File::create(filename)?)))
    }
}