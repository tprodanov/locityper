use std::{
    io::{self, Write},
    cmp::min,
    path::Path,
};
use super::depth::ReadDepthParams;
use crate::{
    model::windows::WindowGetter,
    err::{error, add_path},
    seq::{self, Interval, kmers::KmerCounts},
    ext,
};

/// Single window, used for background parameter estimation.
pub struct Window {
    start: u32,
    end: u32,
    /// GC-content.
    gc: f64,
    /// Is this window used or discarded?
    keep: bool,
}

impl Window {
    fn new(start: u32, end: u32) -> Self {
        Self {
            start, end,
            gc: f64::NAN,
            keep: true,
        }
    }

    pub fn gc(&self) -> f64 {
        self.gc
    }

    pub fn keep(&self) -> bool {
        self.keep
    }
}

/// Filter windows based on GC-content and k-mer frequencies,
/// and returns the number of windows that passed filtering.
fn filter_windows(
    windows: &mut [Window],
    interval: &Interval,
    ref_seq: &[u8],
    kmer_coll: &KmerCounts,
    window_size: u32,
    neighb_size: u32,
    uniq_kmer_frac: f64,
    mut dbg_writer: impl Write,
) -> io::Result<usize>
{
    writeln!(dbg_writer, "#chrom\tstart\tend\tgc\tkmer_cdf1\tkeep")?;
    log::debug!("    Total windows:   {:7}", windows.len());
    let seq_len = ref_seq.len();
    let k = kmer_coll.k();
    assert!(neighb_size > k);
    let kmer_counts = kmer_coll.get_first();
    let left_padding = (neighb_size - window_size) / 2;
    let right_padding = neighb_size - window_size - left_padding;

    let chrom_name = interval.contig_name();
    let seq_shift = interval.start();

    let mut selected = 0;
    let mut have_ns = 0;
    let mut have_common_kmers = 0;
    for window in windows.iter_mut() {
        let start_ix = window.start.saturating_sub(left_padding + seq_shift) as usize;
        let end_ix = min((window.end + right_padding - seq_shift) as usize, seq_len);
        let window_seq = &ref_seq[start_ix..end_ix];
        let mut inv_quant1 = f64::NAN;
        window.keep = if seq::has_n(window_seq) {
            have_ns += 1;
            false
        } else {
            let end_ix2 = (end_ix + 1).checked_sub(k as usize).unwrap();
            assert!(end_ix2 > start_ix);
            // Inverse quantile (1) - what percentage of k-mer counts is <= 1.
            // This calculation is actually faster than quantile, as it does not require sorting.
            inv_quant1 = kmer_counts[start_ix..end_ix2].iter().filter(|&&x| x <= 1).count() as f64
                / (end_ix2 - start_ix) as f64;
            window.gc = seq::gc_content(window_seq);
            if inv_quant1 < uniq_kmer_frac {
                have_common_kmers += 1;
                false
            } else {
                true
            }
        };
        selected += usize::from(window.keep);
        writeln!(dbg_writer, "{}\t{}\t{}\t{:.1}\t{:.1}\t{}", chrom_name, window.start, window.end,
            window.gc, inv_quant1, if window.keep { 'T' } else { 'F' })?;
    }
    log::debug!("    Remove {} windows with Ns, {} windows with common k-mers", have_ns, have_common_kmers);
    log::debug!("    After filtering: {:7}", selected);
    Ok(selected)
}

/// Calculate GC-content and k-mer frequencies based on the window neighbourhood of at least this size.
/// Neighbourhood size = max(MIN_NEIGHBOURHOOD, window_size).
const MIN_NEIGHBOURHOOD: u32 = 300;
/// Automatic window size = 2/3 read length.
const AUTO_WINDOW_MULT: f64 = 2.0 / 3.0;
/// Window size must be at least 20 bp.
const AUTO_WINDOW_MIN: u32 = 20;
/// Maximum window size.
const AUTO_WINDOW_MAX: u32 = 5000;

/// Collection of windows.
pub struct Windows {
    windows: Vec<Window>,
    window_getter: WindowGetter,
    window_size: u32,
    neighb_size: u32,
}

impl Windows {
    pub fn create(
        interval: &Interval,
        ref_seq: &[u8],
        kmer_coll: &KmerCounts,
        seq_info: &super::SequencingInfo,
        params: &ReadDepthParams,
        out_dir: Option<&Path>,
    ) -> crate::Result<Self>
    {
        assert_eq!(interval.len() as usize, ref_seq.len(),
            "ReadDepth: interval and reference sequence have different lengths!");
        let window_size = params.window_size
            .unwrap_or_else(|| ((seq_info.mean_read_len() * AUTO_WINDOW_MULT).round() as u32)
                .clamp(AUTO_WINDOW_MIN, AUTO_WINDOW_MAX));
        if window_size < AUTO_WINDOW_MIN || window_size > AUTO_WINDOW_MAX {
            log::warn!("    Window size {} may be too small or too big", window_size);
        }
        let neighb_size = window_size.max(MIN_NEIGHBOURHOOD);
        log::debug!("    Using {} bp windows, (window neighbourhood size {}, boundary {})",
            window_size, neighb_size, params.boundary_size);

        assert!(interval.len() >= window_size + 2 * params.boundary_size, "Input interval is too short!");
        let n_windows = (interval.len() - 2 * params.boundary_size) / window_size;
        let sum_len = n_windows * window_size;
        let interval_start = interval.start();
        let start = interval_start + (interval.len() - sum_len) / 2;
        let mut windows: Vec<_> = (0..n_windows)
            .map(|i| Window::new(start + i * window_size, start + (i + 1) * window_size))
            .collect();
        let window_getter = WindowGetter::new(start, start + sum_len, window_size);

        let selected = if let Some(dir) = out_dir {
            let dbg_filename = dir.join("windows.bed.gz");
            let dbg_writer = ext::sys::create_gzip(&dbg_filename)?;
            filter_windows(&mut windows, interval, ref_seq, kmer_coll, window_size, neighb_size,
                0.01 * params.uniq_kmer_perc, dbg_writer).map_err(add_path!(dbg_filename))?
        } else {
            filter_windows(&mut windows, interval, ref_seq, kmer_coll, window_size, neighb_size,
                0.01 * params.uniq_kmer_perc, io::sink()).map_err(add_path!(!))?
        };
        if selected == 0 {
            Err(error!(RuntimeError, "Retained 0 windows after filtering"))
        } else {
            Ok(Self { windows, window_getter, window_size, neighb_size })
        }
    }

    /// Returns number of windows.
    pub fn len(&self) -> usize {
        self.windows.len()
    }

    /// Returns iterator over windows.
    pub fn iter(&self) -> std::slice::Iter<'_, Window> {
        self.windows.iter()
    }

    pub fn window_getter(&self) -> &WindowGetter {
        &self.window_getter
    }

    /// Returns true if the middle of the read lies within a selected window.
    pub fn keep_window(&self, middle: u32) -> bool {
        if let Some(i) = self.window_getter.middle_window(middle) {
            self.windows[i as usize].keep
        } else {
            false
        }
    }

    pub fn window_size(&self) -> u32 {
        self.window_size
    }

    pub fn neighb_size(&self) -> u32 {
        self.neighb_size
    }
}
