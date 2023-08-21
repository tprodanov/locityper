pub mod loess;
pub mod bisect;

use std::{
    hash::Hasher,
};

/// Calculates FNV1a function for given bytes.
pub fn fnv1a(bytes: &[u8]) -> u64 {
    let mut hasher = fnv::FnvHasher::default();
    hasher.write(bytes);
    hasher.finish()
}
