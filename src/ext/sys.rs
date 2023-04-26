use std::{
    io::{self, Read, BufRead, BufReader, Write, BufWriter, stdin, stdout},
    fs::{self, File},
    path::{Path, PathBuf},
    ffi::OsStr,
    process::Child,
};
use rand::SeedableRng;
use flate2::{
    bufread::MultiGzDecoder,
    write::GzEncoder,
    Compression,
};
use crate::Error;

/// Finds an executable, and returns Error, if executable is not available.
pub fn find_exe(p: impl AsRef<Path>) -> Result<PathBuf, Error> {
    which::which(p.as_ref()).map_err(|_| Error::NoExec(p.as_ref().to_owned()))
}

/// Returns
/// - stdin if filename is `-`,
/// - gzip reader if filename ends with `.gz`,
/// - regular text file otherwise.
pub fn open(filename: &Path) -> io::Result<Box<dyn BufRead + Send>> {
    if filename == OsStr::new("-") || filename == OsStr::new("/dev/stdin") {
        Ok(Box::new(BufReader::new(stdin())))
    } else {
        let mut stream = BufReader::new(File::open(filename)?);
        let mut two_bytes = [0_u8; 2];
        let bytes_read = stream.read(&mut two_bytes)?;
        stream.seek_relative(-(bytes_read as i64))?;
        // Check gzip magic number.
        if two_bytes[0] == 0x1f && two_bytes[1] == 0x8b {
            Ok(Box::new(BufReader::new(MultiGzDecoder::new(stream))))
        } else {
            Ok(Box::new(stream))
        }
    }
}

/// Creates a buffered file OR stdout if filename is `-`.
pub fn create_uncompressed(filename: &Path) -> io::Result<Box<dyn Write>> {
    if filename == OsStr::new("-") {
        Ok(Box::new(BufWriter::new(stdout())))
    } else {
        Ok(Box::new(BufWriter::new(File::create(filename)?)))
    }
}

/// Creates a gzip file.
pub fn create_gzip(filename: &Path) -> io::Result<BufWriter<GzEncoder<File>>> {
    let file = File::create(filename)?;
    Ok(BufWriter::new(GzEncoder::new(file, Compression::default())))
}

/// Finds all filenames with appropriate extension in the directory.
pub fn filenames_with_ext(dir: &Path, ext: impl AsRef<OsStr>) -> io::Result<Vec<PathBuf>> {
    let mut res = Vec::new();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if entry.file_type()?.is_file() && path.extension() == Some(ext.as_ref()) {
            res.push(path);
        }
    }
    Ok(res)
}

/// Returns a path with a new suffix appended to the end.
pub fn append_path(path: &Path, suffix: impl AsRef<OsStr>) -> PathBuf {
    let mut os_string = path.as_os_str().to_owned();
    os_string.push(suffix.as_ref());
    os_string.into()
}

/// Create directory, if it does not exist yet.
pub fn mkdir(path: impl AsRef<Path>) -> io::Result<()> {
    let path = path.as_ref();
    if !path.exists() {
        fs::create_dir(path)
    } else {
        Ok(())
    }
}

/// Counts lines in the stream.
/// Inspired by `linecount` crate.
pub fn count_lines<R: BufRead>(mut stream: R) -> io::Result<u64> {
    const LF: u8 = b'\n';
    let mut count = 0;
    let mut line: Vec<u8> = Vec::new();
    while stream.read_until(LF, &mut line)? > 0 {
        count += 1;
    }
    if line.last() == Some(&LF) {
        count += 1;
    }
    Ok(count)
}


/// RAII child wrapper, that kills the child if it gets dropped.
pub struct ChildGuard {
    child: Child,
    armed: bool,
}

impl ChildGuard {
    pub fn new(child: Child) -> Self {
        Self {
            child,
            armed: true,
        }
    }

    pub fn disarm(&mut self) {
        self.armed = false;
    }
}

impl Drop for ChildGuard {
    fn drop(&mut self) {
        if self.armed {
            match self.child.kill() {
                Err(e) => {
                    // InvalidInput means that the process exited already.
                    if e.kind() != io::ErrorKind::InvalidInput {
                        log::error!("Could not kill child process: {}", e);
                    }
                }
                Ok(_) => log::error!("Successfully killed child process"),
            }
        }
    }
}

/// Inits random number generator from an optional seed.
pub fn init_rng(seed: Option<u64>) -> rand::rngs::SmallRng {
    if let Some(seed) = seed {
        if seed.count_ones() < 5 {
            log::warn!("Seed ({}) is too simple, consider using a more random number.", seed);
        }
        rand::rngs::SmallRng::seed_from_u64(seed)
    } else {
        rand::rngs::SmallRng::from_entropy()
    }
}
