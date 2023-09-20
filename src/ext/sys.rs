use std::{
    fmt::Write as FmtWrite,
    io::{self, Read, BufRead, BufReader, Write, BufWriter, stdin},
    fs::{self, File},
    path::{Path, PathBuf},
    ffi::OsStr,
    process::{Child, Output},
    mem::ManuallyDrop,
    borrow::Cow,
};
use flate2::{
    bufread::MultiGzDecoder,
    write::GzEncoder,
    Compression,
};
use colored::Colorize;
use crate::err::{Error, add_path};

/// Finds an executable, and returns Error, if executable is not available.
pub fn find_exe(p: impl AsRef<Path>) -> Result<PathBuf, Error> {
    which::which(p.as_ref()).map_err(|_| Error::NoExec(p.as_ref().to_owned()))
}

/// Returns stdin if filename is `-`.
/// Otherwise, tries to guess input format (gzip OR lz4 OR no compression).
pub fn open(filename: &Path) -> Result<Box<dyn BufRead + Send>, Error> {
    if filename == OsStr::new("-") || filename == OsStr::new("/dev/stdin") {
        return Ok(Box::new(BufReader::new(stdin())));
    } else {
        // Guess file format (gz | lz4 | no compression).
        let mut stream = BufReader::new(File::open(filename).map_err(add_path!(filename))?);
        let mut two_bytes = [0_u8; 2];
        let bytes_read = stream.read(&mut two_bytes).map_err(add_path!(filename))?;
        stream.seek_relative(-(bytes_read as i64)).map_err(add_path!(filename))?;
        if two_bytes[0] == 0x1f && two_bytes[1] == 0x8b {
            // gzip magic number
            Ok(Box::new(BufReader::new(MultiGzDecoder::new(stream))))
        } else if two_bytes[0] == 0x04 && two_bytes[1] == 0x22 {
            // lz4 magic number
            Ok(Box::new(BufReader::new(lz4::Decoder::new(stream).map_err(add_path!(filename))?)))
        } else {
            Ok(Box::new(stream))
        }
    }
}

/// Loads full JSON contents from a file.
pub fn load_json(filename: &Path) -> Result<json::JsonValue, Error> {
    let stream = open(filename)?;
    let contents = io::read_to_string(stream).map_err(add_path!(filename))?;
    json::parse(&contents).map_err(|_| Error::InvalidData(
        format!("Failed parsing {}: invalid JSON format", filename.display())))
}

pub type GzFile = BufWriter<GzEncoder<File>>;

/// Creates a gzip compressed file.
pub fn create_gzip(filename: &Path) -> Result<GzFile, Error> {
    let file = File::create(filename).map_err(add_path!(filename))?;
    Ok(BufWriter::new(GzEncoder::new(file, Compression::default())))
}

pub struct AutoFinishLz4<W: Write> {
    encoder: ManuallyDrop<lz4::Encoder<W>>,
}

impl<W: io::Write> Drop for AutoFinishLz4<W> {
    fn drop(&mut self) {
        let encoder = unsafe { ManuallyDrop::take(&mut self.encoder) };
        if let Err(e) = encoder.finish().1 {
            log::error!("Lz4 encoder panicked on drop: {:?}", e);
        }
    }
}

impl<W: Write> Write for AutoFinishLz4<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.encoder.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.encoder.flush()
    }
}

/// Creates an lz4 compressed file.
fn create_lz4(filename: &Path, level: u32) -> Result<AutoFinishLz4<BufWriter<File>>, Error> {
    let file = create_file(filename)?;
    let encoder = lz4::EncoderBuilder::new().level(level).build(file).map_err(add_path!(filename))?;
    Ok(AutoFinishLz4 {
        encoder: ManuallyDrop::new(encoder),
    })
}

pub fn create_lz4_slow(filename: &Path) -> Result<AutoFinishLz4<BufWriter<File>>, Error> {
    create_lz4(filename, 7)
}

/// Creates buffered output file.
pub fn create_file(filename: &Path) -> Result<BufWriter<File>, Error> {
    File::create(filename).map_err(add_path!(filename))
        .map(|w| BufWriter::with_capacity(131_072, w)) // Set buffer capacity to 128 kb.
}

/// Finds all filenames with appropriate extension in the directory.
/// Extension should not include `.`, for example `filenames_with_ext(dir, "csv")`.
pub fn filenames_with_ext(dir: &Path, ext: impl AsRef<OsStr>) -> Result<Vec<PathBuf>, Error> {
    let mut res = Vec::new();
    for entry in fs::read_dir(dir).map_err(add_path!(dir))? {
        let entry = entry.map_err(add_path!(!))?;
        let path = entry.path();
        if entry.file_type().map_err(add_path!(path))?.is_file() && path.extension() == Some(ext.as_ref()) {
            res.push(path);
        }
    }
    Ok(res)
}

/// Returns a path with a new suffix appended to the end.
pub fn path_append(path: &Path, suffix: impl AsRef<OsStr>) -> PathBuf {
    let mut os_string = path.as_os_str().to_owned();
    os_string.push(suffix.as_ref());
    os_string.into()
}

/// Create directory, if it does not exist yet.
pub fn mkdir(path: impl AsRef<Path>) -> Result<(), Error> {
    let path = path.as_ref();
    if !path.exists() {
        fs::create_dir(path).map_err(add_path!(path))
    } else {
        Ok(())
    }
}

/// Directly concantenates files, without trying to decompress them.
/// Therefore, if input files are already gzipped, output writer should be plain, without compression.
pub fn concat_files(filenames: impl Iterator<Item = impl AsRef<Path>>, mut writer: impl Write) -> Result<(), Error> {
    for filename in filenames {
        let mut reader = File::open(&filename).map_err(add_path!(filename))?;
        io::copy(&mut reader, &mut writer).map_err(add_path!(filename))?;
    }
    Ok(())
}

/// Returns path basename, or `???`, if it is empty.
fn safe_basename(path: &Path) -> Cow<'_, str> {
    path.file_name().unwrap_or(path.as_os_str()).to_string_lossy()
}

/// Returns program version.
/// This function does not panic, but can return `version unavailable`.
pub fn get_program_version(exe: &Path) -> String {
    let out = match std::process::Command::new(exe).arg("--version").output() {
        Ok(out) => out,
        Err(_) => return format!("{} version unavailable", safe_basename(exe)),
    };
    let msg = if out.stdout.is_empty() { out.stderr } else { out.stdout };
    if msg.is_empty() {
        return format!("{} version unavailable", safe_basename(exe))
    }
    const MAX_SIZE: usize = 1000;
    let i = msg.iter().take(MAX_SIZE).position(|&ch| ch == b'\n').unwrap_or(MAX_SIZE.min(msg.len()));
    let msg_head = String::from_utf8_lossy(&msg[..i]);
    if msg_head.contains(' ') {
        msg_head.trim().to_string()
    } else {
        format!("{} {}", safe_basename(exe), msg_head.trim())
    }
}

/// Clips errors message at 10000 characters.
fn clip_msg(bytes: &[u8]) -> Cow<'_, str> {
    const MAX_LEN: usize = 10000;
    if bytes.len() > MAX_LEN {
        let mut s = String::from_utf8_lossy(&bytes[..MAX_LEN]).into_owned();
        write!(s, " ...").unwrap();
        Cow::Owned(s)
    } else {
        String::from_utf8_lossy(bytes)
    }
}

fn format_output(s: &mut String, out: &Output) {
    match out.status.code() {
        Some(code) => writeln!(s, "    {}: {}", "Exit code".bold(), code).unwrap(),
        None => writeln!(s, "    {}: {}", "Exit code".bold(), "unknown").unwrap(),
    };
    let stdout = clip_msg(&out.stdout);
    let stdout = stdout.trim();
    if !stdout.is_empty() {
        writeln!(s, "    {}: {}", "Stdout".bold(), stdout).unwrap();
    }
    let stderr = clip_msg(&out.stderr);
    let stderr = stderr.trim();
    if !stderr.is_empty() {
        writeln!(s, "    {}: {}", "Stderr".bold(), stderr).unwrap();
    }
}

pub struct PipeGuard {
    /// Executable, corresponding to each of the children.
    executables: Vec<PathBuf>,
    /// Child processes. Some components may be dropped, then corresponding output is set to self.outputs.
    children: Vec<Option<Child>>,
    /// Outputs, collected from the child processes.
    outputs: Vec<io::Result<Output>>,
    /// Was `wait` or `fail` already called?
    finished: bool,
}

#[inline]
fn init_output() -> io::Result<Output> {
    Err(io::Error::from(io::ErrorKind::NotFound))
}

impl PipeGuard {
    pub fn new(exe: PathBuf, child: Child) -> Self {
        Self {
            executables: vec![exe],
            children: vec![Some(child)],
            outputs: vec![init_output()],
            finished: false,
        }
    }

    pub fn push(&mut self, exe: PathBuf, child: Child) {
        self.executables.push(exe);
        self.children.push(Some(child));
        self.outputs.push(init_output());
    }

    /// We already know that the pipe failed, need to kill all steps and collect available information.
    pub fn fail(&mut self) -> String {
        let mut s = String::new();
        for (exe, (child, output)) in self.executables.iter()
                .zip(self.children.iter_mut().zip(self.outputs.iter_mut())) {
            writeln!(s, "{}", get_program_version(&exe).bold()).unwrap();
            if let Some(mut child) = child.take() {
                if let Err(e) = child.kill() {
                    if e.kind() != io::ErrorKind::InvalidInput {
                        writeln!(s, "    {}: {}", "Could not kill process".red(), e).unwrap();
                    }
                }
                *output = child.wait_with_output();
            }
            match output {
                Ok(out) => format_output(&mut s, out),
                Err(e) => writeln!(s, "    {}: {}", "Could not get child output".red(), e).unwrap(),
            }
        }
        if s.ends_with('\n') {
            s.pop();
        }
        self.finished = true;
        s
    }

    /// Waits for each process from end to start.
    /// If all processed finished successfully, returns output of the last process.
    pub fn wait(mut self) -> Result<Output, Error> {
        for (child, output) in self.children.iter_mut().rev().zip(self.outputs.iter_mut().rev()) {
            if let Some(child) = child.take() {
                *output = child.wait_with_output();
                if !output.as_ref().map(|out| out.status.success()).unwrap_or(false) {
                    return Err(Error::Subprocess(self.fail()));
                }
            }
        }
        self.finished = true;
        Ok(self.outputs.pop().expect("At least one process must be present").expect("Last output must be defined"))
    }
}

impl Drop for PipeGuard {
    fn drop(&mut self) {
        if !self.finished {
            log::error!("Subprocess killed:\n{}", self.fail());
        }
    }
}
