pub mod loess;
pub mod bisect;

// Use these types instead of the default Hash Map.

pub type BuildHasher = std::hash::BuildHasherDefault<wyhash2::WyHash>;
pub type HashMap<K, V> = std::collections::HashMap<K, V, BuildHasher>;
pub type HashSet<T> = std::collections::HashSet<T, BuildHasher>;

/// Newtype over two u32 values, useful for nohash IntMap.
#[derive(PartialOrd, Ord, PartialEq, Eq, Clone, Copy)]
pub struct TwoU32(pub u32, pub u32);

impl std::hash::Hash for TwoU32 {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        hasher.write_u64((u64::from(self.0) << 32) | u64::from(self.1))
    }
}

impl nohash::IsEnabled for TwoU32 {}
