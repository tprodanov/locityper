use rand::SeedableRng;

pub type XoshiroRng = rand_xoshiro::Xoshiro256PlusPlus;

/// Inits random number generator from an optional seed.
pub fn init_rng(seed: Option<u64>) -> XoshiroRng {
    if let Some(seed) = seed {
        if seed.count_ones() < 5 {
            log::warn!("Seed ({}) is too simple, consider using a more random number\
                ({}/64 bits set in binary representation)", seed, seed.count_ones());
        }
        XoshiroRng::seed_from_u64(seed)
    } else {
        let mut buffer = [0_u8; 8];
        getrandom::fill(&mut buffer).unwrap();
        let seed = u64::from_le_bytes(buffer);
        log::debug!("Using random seed {}", seed);
        XoshiroRng::seed_from_u64(seed)
    }
}
