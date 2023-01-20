
const FNV_PRIME: u64 = 0x100000001b3;
const FNV_OFFSET: u64 = 0xcbf29ce484222325;

/// Calculates FNV-1 hash.
pub fn fnv1(bytes: &[u8]) -> u64 {
    let mut hash = FNV_OFFSET;
    for &b in bytes.iter() {
        hash = hash.wrapping_mul(FNV_PRIME);
        hash ^= u64::from(b);
    }
    hash
}

/// Calculates FNV-1a hash.
pub fn fnv1a(bytes: &[u8]) -> u64 {
    let mut hash = FNV_OFFSET;
    for &b in bytes.iter() {
        hash ^= u64::from(b);
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}
