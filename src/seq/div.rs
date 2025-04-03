use std::{
    thread, io,
    cmp::Ordering,
    sync::Arc,
    path::Path,
};
use crate::{
    err::{error, add_path},
    math::RoundDiv,
    seq::{NamedSeq, kmers},
    ext::{TriangleMatrix, fmt},
};
use varint_rs::{VarintWriter, VarintReader};

/// Calculates number of non-matching entries and corresponding Jaccard distance given two sorted vectors..
fn jaccard_distance(minimizers1: &[u64], minimizers2: &[u64]) -> (u32, f64) {
    let mut iter1 = minimizers1.iter();
    let mut opt_x = iter1.next();
    let mut iter2 = minimizers2.iter();
    let mut opt_y = iter2.next();

    let mut inters: u32 = 0;
    while let (Some(x), Some(y)) = (opt_x, opt_y) {
        match x.cmp(y) {
            Ordering::Equal => {
                inters += 1;
                opt_x = iter1.next();
                opt_y = iter2.next();
            }
            Ordering::Less => opt_x = iter1.next(),
            Ordering::Greater => opt_y = iter2.next(),
        }
    }
    let n1 = minimizers1.len() as u32;
    let n2 = minimizers2.len() as u32;
    let overlap = n1 + n2 - inters;
    let non_shared = overlap - inters;
    (non_shared, f64::from(non_shared) / f64::from(overlap))
}

/// Calculates all pairwise divergences and returns the distance matrix.
/// Distances are calculated based on the number of shared non-canonical minimizers.
/// Returns two values for each pair: number of non-shared minimizers and Jaccard distance.
pub fn minimizer_divergences(
    entries: &[NamedSeq],
    pairs: &[(u32, u32)],
    k: u8,
    w: u8,
    threads: u16,
) -> Vec<(u32, f64)>
{
    let mut minimizers = Vec::with_capacity(entries.len());
    for entry in entries {
        let seq = entry.seq();
        // Expected num of minimizers = 2L / (w + 1), here we take a bit more to be safe.
        let mut buff = Vec::with_capacity((5 * seq.len()).fast_round_div(2 * usize::from(w) + 2));
        kmers::minimizers::<u64, _, { kmers::NON_CANONICAL }>(seq, k, w, &mut buff);
        buff.sort_unstable();
        minimizers.push(buff);
    }

    if threads == 1 {
        let mut divs = Vec::with_capacity(pairs.len());
        for &(i, j) in pairs.iter() {
            divs.push(jaccard_distance(&minimizers[i as usize], &minimizers[j as usize]));
        }
        divs
    } else {
        divergences_multithread(minimizers, pairs, threads)
    }
}

fn divergences_multithread(
    minimizers: Vec<Vec<u64>>,
    pairs: &[(u32, u32)],
    threads: u16,
) -> Vec<(u32, f64)> {
    let threads = usize::from(threads);
    let pairs = Arc::new(pairs.to_vec());
    let n_pairs = pairs.len();
    let mut handles = Vec::with_capacity(threads);

    let minimizers = Arc::new(minimizers);
    let mut start = 0;
    for worker_ix in 0..threads {
        if start == n_pairs {
            break;
        }
        let rem_workers = threads - worker_ix;
        let end = start + (n_pairs - start).fast_ceil_div(rem_workers);
        assert!(start < end);
        // Closure with cloned data.
        {
            let minimizers = Arc::clone(&minimizers);
            let pairs = Arc::clone(&pairs);
            handles.push(thread::spawn(move ||
                pairs[start..end].iter().map(|&(i, j)|
                    jaccard_distance(&minimizers[i as usize], &minimizers[j as usize]))
                    .collect::<Vec<_>>()
            ));
        }
        start = end;
    }
    assert_eq!(start, n_pairs);

    handles.into_iter()
        .flat_map(|handle| handle.join().expect("Worker process failed"))
        .collect()
}

/// Writes the number of non-shared minimizers to a binary file.
pub fn write_divergences<W>(k: u8, w: u8, divs: &TriangleMatrix<(u32, f64)>, mut f: W) -> io::Result<()>
where W: VarintWriter<Error = io::Error>,
{
    f.write(k)?;
    f.write(w)?;
    f.write_u32_varint(divs.side() as u32)?;
    for &(d, _) in divs.iter() {
        f.write_u32_varint(d)?;
    }
    Ok(())
}

/// Loads minimizer sizes k and w, and the number of non-shared minimizers for each pair of contigs.
pub fn load_divergences<R>(
    mut f: R,
    filename: &Path,
    n: usize,
) -> crate::Result<(u8, u8, TriangleMatrix<u32>)>
where R: VarintReader<Error = io::Error>,
{
    let k = f.read().map_err(add_path!(filename))?;
    let w = f.read().map_err(add_path!(filename))?;
    let m = f.read_u32_varint().map_err(add_path!(filename))?;
    if n != m as usize {
        return Err(error!(InvalidData,
            "Cannot read distances from {}: invalid number of haplotypes (expected {}, found {})",
            fmt::path(filename), n, m));
    }

    let total = TriangleMatrix::<()>::expected_len(n);
    let mut divs = Vec::with_capacity(total);
    for _ in 0..total {
        let d = f.read_u32_varint().map_err(add_path!(filename))?;
        divs.push(d);
    }
    Ok((k, w, TriangleMatrix::from_linear(n, divs)))
}
