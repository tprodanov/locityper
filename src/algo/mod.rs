pub mod loess;
pub mod vec_ext;
pub mod hash;
pub mod bisect;

use std::fmt::Debug;

/// Parses a string value ending with suffixes G, M or K.
/// Additionally, removes all underscores and commas from the string.
///
/// Allowed strings: `[0-9][0-9_,]*[GgMmKk]?`
pub fn parse_int<T: std::convert::TryFrom<u64>>(s: &str) -> Result<T, String> {
    let mut rev_bytes = s.as_bytes().iter().rev();
    let (mut n, mut mult) = match rev_bytes.next().copied() {
        None => return Err("Cannot parse an empty string into int".to_string()),
        Some(b'G') | Some(b'g') => (0, 1_000_000_000),
        Some(b'M') | Some(b'm') => (0, 1_000_000),
        Some(b'K') | Some(b'k') => (0, 1000),
        Some(c @ b'0' ..= b'9') => (u64::from(c - b'0'), 10),
        Some(c) => return Err(format!("Cannot parse string {:?} to int, unexpected last symbol '{}'", s, c as char)),
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
            c => return Err(format!("Cannot parse string {:?} to int, unexpected symbol '{}'", s, c as char)),
        }
    }

    if !was_digit {
        return Err(format!("Cannot parse string {:?} to int, unexpected first letter", s));
    }
    n.try_into().map_err(|_| format!("Cannot parse string {:?} to int, value is too large", s))
}

/// Trait that encodes either `&mut Vec<T>` or `()`, depending on if the information needs to be saved or not.
pub trait VecOrNone<T> : Debug {
    /// True for `()` - meaning that the object does not actually do any work.
    const IS_SINK: bool;

    /// Push a new element.
    fn push(&mut self, val: T);

    /// Returns length, if available.
    fn try_len(&self) -> Option<usize>;
}

impl<T> VecOrNone<T> for () {
    const IS_SINK: bool = true;

    fn push(&mut self, _val: T) {}

    fn try_len(&self) -> Option<usize> { None }
}

impl<T> VecOrNone<T> for &mut Vec<T>
where T: Debug
{
    const IS_SINK: bool = false;

    fn push(&mut self, val: T) {
        Vec::push(*self, val);
    }

    fn try_len(&self) -> Option<usize> {
        Some(self.len())
    }
}
