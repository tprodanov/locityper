use std::{
    io::{self, BufWriter, Write, BufReader, BufRead, stdin, stdout},
    fs::{self, File},
    path::{Path, PathBuf},
    ffi::OsStr,
};
use flate2::{
    read::GzDecoder,
    write::GzEncoder,
    Compression,
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

/// Returns
/// - stdin if filename is `-`,
/// - gzip reader if filename ends with `.gz`,
/// - regular text file otherwise.
pub(super) fn open(filename: &Path) -> Result<Box<dyn BufRead>, Error> {
    if filename == OsStr::new("-") {
        Ok(Box::new(BufReader::new(stdin())))
    } else {
        let file = File::open(filename)?;
        if filename.ends_with(".gz") {
            Ok(Box::new(BufReader::new(GzDecoder::new(file))))
        } else {
            Ok(Box::new(BufReader::new(file)))
        }
    }
}

/// Creates a buffered file OR stdout if filename is `-`.
pub(super) fn create_uncompressed(filename: &Path) -> Result<Box<dyn Write>, Error> {
    if filename == OsStr::new("-") {
        Ok(Box::new(BufWriter::new(stdout())))
    } else {
        Ok(Box::new(BufWriter::new(File::create(filename)?)))
    }
}

/// Creates a gzip file.
pub(super) fn create_gzip(filename: &Path) -> Result<BufWriter<GzEncoder<File>>, Error> {
    let file = File::create(filename)?;
    Ok(BufWriter::new(GzEncoder::new(file, Compression::default())))
}

/// Finds all filenames with appropriate extension in the directory.
pub fn find_filenames(dir: &Path, ext: &OsStr) -> io::Result<Vec<PathBuf>> {
    let mut res = Vec::new();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if entry.file_type()?.is_file() && path.extension() == Some(ext) {
            res.push(path);
        }
    }
    Ok(res)
}

// /// Appends current command to the end of the file.
// pub(super) fn update_command_history(filename: &Path) -> io::Result<()> {
//     let args: Vec<_> = std::env::args().collect();
//     let mut f = File::options().append(true).open(filename)?;
//     writeln!(f, "{}\n", args.join(" "))?;
//     f.sync_all()?;
//     Ok(())
// }

/// Returns a path with a new suffix appended to the end.
pub fn append_path(path: &Path, suffix: impl AsRef<OsStr>) -> PathBuf {
    let mut os_string = path.as_os_str().to_owned();
    os_string.push(suffix.as_ref());
    os_string.into()
}
