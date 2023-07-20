use std::{
    fmt::{self, Display, Debug},
    path::{Path, PathBuf},
    process::Command,
    ffi::OsStr,
};

/// Pretty path formatting: replace $HOME with ~, put quotes around if needed.
pub fn path(path: &Path) -> String {
    lazy_static::lazy_static!{
        static ref HOME: Option<PathBuf> = std::env::var_os("HOME").map(|s| PathBuf::from(s));
    }
    if let Some(home) = (*HOME).as_ref() {
        if let Ok(suffix) = path.strip_prefix(home) {
            let tilde_path = Path::new("~").join(suffix);
            let s = tilde_path.to_string_lossy();
            return if s.contains(char::is_whitespace) { format!("'{}'", s) } else { s.into_owned() };
        }
    }
    let s = path.to_string_lossy();
    if s.contains(char::is_whitespace) { format!("'{}'", s) } else { s.into_owned() }
}

/// Converts command into a string, removing quotes if argument has no whitespace, and replacing HOME with ~.
pub fn command(cmd: &Command) -> String {
    std::iter::once(cmd.get_program())
        .chain(cmd.get_args())
        .map(OsStr::as_ref)
        .map(path)
        .collect::<Vec<_>>()
        .join(" ")
}

/// Formats duration as `HH:MM:SS.SSS`.
pub struct Duration(pub std::time::Duration);

impl Display for Duration {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        const IN_HOUR: u64 = 3600;
        const IN_MINUTE: u64 = 60;
        let mut seconds = self.0.as_secs();
        write!(f, "{}:", seconds / IN_HOUR)?;
        seconds %= IN_HOUR;
        write!(f, "{:02}:", seconds / IN_MINUTE)?;
        seconds %= IN_MINUTE;
        write!(f, "{:02}.{:03}", seconds, self.0.subsec_millis())?;
        Ok(())
    }
}

impl Debug for Duration {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.0.fmt(f)
    }
}
