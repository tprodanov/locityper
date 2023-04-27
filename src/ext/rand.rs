use rand::SeedableRng;

pub type XoshiroRng = rand_xoshiro::Xoshiro256PlusPlus;

/// Inits random number generator from an optional seed.
pub fn init_rng(seed: Option<u64>) -> XoshiroRng {
    if let Some(seed) = seed {
        if seed.count_ones() < 5 {
            log::warn!("Seed ({}) is too simple, consider using a more random number.", seed);
        }
        XoshiroRng::seed_from_u64(seed)
    } else {
        let mut buffer = [0_u8; 8];
        getrandom::getrandom(&mut buffer).unwrap();
        XoshiroRng::seed_from_u64(u64::from_le_bytes(buffer))
    }
}
