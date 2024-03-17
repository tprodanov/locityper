use std::{
    thread,
    cmp::Ordering,
    sync::Arc,
    io::{self, Write, Read},
    path::Path,
};
use crate::{
    err::add_path,
    math::RoundDiv,
    seq::{NamedSeq, kmers},
    ext::{TriangleMatrix, fmt},
};

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
pub fn pairwise_divergences(
    entries: &[NamedSeq],
    k: u8,
    w: u8,
    threads: u16,
) -> TriangleMatrix<(u32, f64)> {
    let mut minimizers = Vec::with_capacity(entries.len());
    for entry in entries {
        let seq = entry.seq();
        // Expected num of minimizers = 2L / (w + 1), here we take a bit more to be safe.
        let mut buff = Vec::with_capacity((5 * seq.len()).fast_round_div(2 * usize::from(w) + 2));
        kmers::minimizers::<u64, u64, { kmers::NON_CANONICAL }>(seq, k, w, &mut buff);
        buff.sort_unstable();
        minimizers.push(buff);
    }

    let n = entries.len();
    let divergences = if threads == 1 {
        let mut divs = Vec::with_capacity(n * (n - 1) / 2);
        for (i, minimizers1) in minimizers.iter().enumerate() {
            for minimizers2 in minimizers[i+1..].iter() {
                divs.push(jaccard_distance(minimizers1, minimizers2));
            }
        }
        divs
    } else {
        divergences_multithread(minimizers, threads)
    };
    TriangleMatrix::from_linear(n, divergences)
}

fn divergences_multithread(
    minimizers: Vec<Vec<u64>>,
    threads: u16,
) -> Vec<(u32, f64)> {
    let threads = usize::from(threads);
    let n = minimizers.len();
    let pairs: Arc<Vec<(u32, u32)>> = Arc::new(TriangleMatrix::indices(n).map(|(i, j)| (i as u32, j as u32)).collect());
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
pub fn write_divergences(divs: &TriangleMatrix<(u32, f64)>, mut f: impl Write) -> io::Result<()> {
    f.write_all(&(divs.side() as u32).to_le_bytes())?;
    for &(n, _) in divs.iter() {
        f.write_all(&n.to_le_bytes())?;
    }
    Ok(())
}

/// Reads the number of non-shared minimizers from a binary file.
pub fn load_divergences(
    mut f: impl Read,
    filename: &Path,
    n: usize,
) -> Result<TriangleMatrix<u32>, crate::Error>
{
    let mut buf = [0_u8; 4];
    f.read_exact(&mut buf).map_err(add_path!(filename))?;
    let m = u32::from_le_bytes(buf);
    if m as usize != n {
        return Err(crate::Error::InvalidData(
            format!("Cannot read distances from {}: invalid number of haplotypes ({} != {})",
            fmt::path(filename), n, m)));
    }

    let total = TriangleMatrix::<()>::expected_len(n);
    let mut divs = Vec::with_capacity(total);
    for _ in 0..total {
        f.read_exact(&mut buf).map_err(add_path!(filename))?;
        divs.push(u32::from_le_bytes(buf));
    }
    Ok(TriangleMatrix::from_linear(n, divs))
}
