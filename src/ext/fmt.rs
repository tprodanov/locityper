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

/// Integer with pretty formatting (k, M, G suffixes, as well as writing `inf` when MAX) and parsing.
macro_rules! impl_pretty_int {
    ($name:ident, $prim:ty) => {
        #[derive(Clone, Copy, Debug, PartialOrd, Ord, PartialEq, Eq)]
        pub struct $name(pub $prim);

        impl $name {
            /// Extract inner value.
            pub fn get(self) -> $prim {
                self.0
            }
        }

        impl Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                // Note: this will break on 16-bit systems.
                const BILLION: $prim = 1_000_000_000;
                const MILLION: $prim = 1_000_000;
                const THOUSAND: $prim = 1000;

                if self.0 == 0 {
                    write!(f, "0")
                } else if self.0 == <$prim>::MAX {
                    write!(f, "inf")
                } else if self.0 % BILLION == 0 {
                    write!(f, "{}G", self.0 / BILLION)
                } else if self.0 % MILLION == 0 {
                    write!(f, "{}M", self.0 / MILLION)
                } else if self.0 % THOUSAND == 0 {
                    write!(f, "{}k", self.0 / THOUSAND)
                } else {
                    write!(f, "{}", self.0)
                }
            }
        }

        impl std::str::FromStr for $name {
            type Err = String;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                if s == "inf" || s == "Inf" || s == "INF" {
                    return Ok(Self(<$prim>::MAX));
                }

                let mut rev_bytes = s.as_bytes().iter().rev();
                let (mut n, mut mult) = match rev_bytes.next().copied() {
                    None => return Err("Cannot parse an empty string into int".to_owned()),
                    Some(b'G') | Some(b'g') => (0, 1_000_000_000),
                    Some(b'M') | Some(b'm') => (0, 1_000_000),
                    Some(b'K') | Some(b'k') => (0, 1000),
                    Some(c @ b'0' ..= b'9') => (<$prim>::from(c - b'0'), 10),
                    Some(c) => return Err(format!("Cannot parse string {:?} to int, unexpected last symbol '{}'",
                        s, char::from(c))),
                };

                let mut was_digit = mult == 10;
                for c in rev_bytes {
                    match *c {
                        c @ b'0' ..= b'9' => {
                            was_digit = true;
                            n += mult * <$prim>::from(c - b'0');
                            mult *= 10;
                        },
                        b',' | b'_' => was_digit = false,
                        c => return Err(format!("Cannot parse string {:?} to int, unexpected symbol '{}'",
                            s, char::from(c))),
                    }
                }

                if was_digit {
                    Ok(Self(n))
                } else {
                    Err(format!("Cannot parse string {:?} to int, unexpected first letter", s))
                }
            }
        }
    }
}

impl_pretty_int!(PrettyU32, u32);
impl_pretty_int!(PrettyU64, u64);
impl_pretty_int!(PrettyUsize, usize);
