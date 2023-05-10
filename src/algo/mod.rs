pub mod loess;
pub mod bisect;

use std::{
    hash::Hasher,
};

/// Parses a string value ending with suffixes G, M or K.
/// Additionally, removes all underscores and commas from the string.
///
/// Allowed strings: `[0-9][0-9_,]*[GgMmKk]?`
pub fn parse_int<T: std::convert::TryFrom<u64>>(s: &str) -> Result<T, String> {
    let mut rev_bytes = s.as_bytes().iter().rev();
    let (mut n, mut mult) = match rev_bytes.next().copied() {
        None => return Err("Cannot parse an empty string into int".to_owned()),
        Some(b'G') | Some(b'g') => (0, 1_000_000_000),
        Some(b'M') | Some(b'm') => (0, 1_000_000),
        Some(b'K') | Some(b'k') => (0, 1000),
        Some(c @ b'0' ..= b'9') => (u64::from(c - b'0'), 10),
        Some(c) => return Err(format!("Cannot parse string {:?} to int, unexpected last symbol '{}'", s, char::from(c))),
    };

    let mut was_digit = mult == 10;
    for c in rev_bytes {
        match *c {
            c @ b'0' ..= b'9' => {
                was_digit = true;
                n += mult * u64::from(c - b'0');
                mult *= 10;
            },
            b',' | b'_' => was_digit = false,
            c => return Err(format!("Cannot parse string {:?} to int, unexpected symbol '{}'", s, char::from(c))),
        }
    }

    if !was_digit {
        return Err(format!("Cannot parse string {:?} to int, unexpected first letter", s));
    }
    n.try_into().map_err(|_| format!("Cannot parse string {:?} to int, value is too large", s))
}

/// Calculates FNV1a function for given bytes.
pub fn fnv1a(bytes: &[u8]) -> u64 {
    let mut hasher = fnv::FnvHasher::default();
    hasher.write(bytes);
    hasher.finish()
}
