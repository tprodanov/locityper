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

/// Newtype over two u32 values, useful for nohash IntMap.
#[derive(PartialOrd, Ord, PartialEq, Eq, Clone, Copy)]
pub struct TwoU32(pub u32, pub u32);

impl std::hash::Hash for TwoU32 {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        hasher.write_u64((u64::from(self.0) << 32) | u64::from(self.1))
    }
}

impl nohash::IsEnabled for TwoU32 {}
