pub mod loess;
pub mod bisect;

// Use these types instead of the default Hash Map.

pub type Hasher = std::hash::BuildHasherDefault<wyhash::WyHash>;
// Can easily switch between standard and hashbrown implementations, if needed.
pub use std::collections::{hash_map, hash_set};
// pub use hashbrown::{hash_map, hash_set};
pub type HashMap<K, V> = hash_map::HashMap<K, V, Hasher>;
pub type HashSet<T> = hash_set::HashSet<T, Hasher>;
pub type IntMap<K, V> = hash_map::HashMap<K, V, nohash::BuildNoHashHasher<K>>;
pub type IntSet<T> = hash_set::HashSet<T, nohash::BuildNoHashHasher<T>>;

/// Calculates hash function on given bytes using pre-selected hash function.
#[inline]
pub fn get_hash(bytes: &[u8]) -> u64 {
    wyhash::wyhash(bytes, 0)
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
