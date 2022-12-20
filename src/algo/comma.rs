//! Simple trait for writing numbers with comma separator.

use std::fmt::{self, Write};

/// Write
pub trait CommaWrite: Copy {
    fn comma_write_to<W: Write>(self, f: W) -> fmt::Result;

    #[inline]
    fn comma_string(self) -> String {
        let mut s = String::new();
        self.comma_write_to(&mut s).expect("Cannot write to string");
        s
    }
}

impl CommaWrite for u32 {
    fn comma_write_to<W: Write>(mut self, mut f: W) -> fmt::Result {
        const B: u32 = 1_000_000_000;
        const M: u32 = 1_000_000;
        const K: u32 = 1_000;
        let b = self / B;
        if b > 0 {
            write!(f, "{},", b)?;
            self %= B;
        }
        let m = self / M;
        if m > 0 {
            write!(f, "{},", m)?;
            self %= M;
        }
        let k = self / K;
        if k > 0 {
            write!(f, "{},", k)?;
            self %= K;
        }
        write!(f, "{}", self)
    }
}
