use std::{
    io::{self, BufWriter, Write, BufReader, BufRead, stdin, stdout},
    fmt::{Write as FmtWrite},
    fs::{self, File},
    path::{Path, PathBuf},
    ffi::OsStr,
    borrow::Cow,
    process::Command,
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

/// Replace HOME directory with ~ in the path.
pub(super) fn tilde_home<'a>(path: &'a Path) -> Cow<'a, str> {
    lazy_static::lazy_static!{
        static ref HOME: Option<PathBuf> = std::env::var_os("HOME").map(|s| PathBuf::from(s));
    }
    if let Some(home) = (*HOME).as_ref() {
        match path.strip_prefix(home) {
            Ok(suffix) => Cow::Owned(Path::new("~").join(suffix).to_string_lossy().into_owned()),
            Err(_) => path.to_string_lossy(),
        }
    } else {
        path.to_string_lossy()
    }
}

/// Converts command into a string, removing quotes if argument has no whitespace, and replacing HOME with ~.
pub(super) fn fmt_cmd(cmd: &Command) -> String {
    let prg = tilde_home(cmd.get_program().as_ref());
    let mut s = if prg.as_ref().contains(char::is_whitespace) {
        format!("'{}'", prg)
    } else {
        prg.into_owned()
    };

    for arg in cmd.get_args() {
        let arg = tilde_home(arg.as_ref());
        if arg.contains(char::is_whitespace) {
            write!(s, " '{}'", arg).unwrap();
        } else {
            write!(s, " {}", arg).unwrap();
        }
    }
    s
}

pub(super) fn fmt_duration(duration: std::time::Duration) -> String {
    const IN_HOUR: u64 = 3600;
    const IN_MINUTE: u64 = 60;
    let mut seconds = duration.as_secs();

    let mut res = String::new();
    write!(res, "{}:", seconds / IN_HOUR).unwrap();
    seconds %= IN_HOUR;
    write!(res, "{:02}:", seconds / IN_MINUTE).unwrap();
    seconds %= IN_MINUTE;
    write!(res, "{:02}.{:03}", seconds, duration.subsec_millis()).unwrap();
    res
}
