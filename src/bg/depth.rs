use std::{
    io::{self, Write},
    path::Path,
};
use htslib::bam::{
    Record,
    ext::BamRecordExtensions,
};
use crate::{
    seq::{
        self, Interval,
        cigar::Cigar,
        kmers::KmerCounts,
        aln::LightAlignment,
    },
    algo::{bisect, loess::Loess},
    ext::{
        self,
        vec::{VecExt, F64Ext},
    },
    math::{
        Ln,
        distr::{NBinom, WithMoments},
    },
    err::{Error, validate_param},
};
use super::{
    err_prof::ErrorProfile,
    ser::{JsonSer, parse_f64_arr, json_get},
};

pub trait DepthDistr {
    fn window_size(&self) -> u32;

    fn ln_prob(&self, gc_content: u8, depth: u32) -> f64;
}

/// Number of various read types per window.
#[derive(Clone, Debug)]
struct WindowCounts {
    start: u32,
    end: u32,
    depth: [u32; 2],
}

impl WindowCounts {
    fn new(start: u32, end: u32) -> Self {
        WindowCounts {
            start, end,
            depth: [0, 0],
        }
    }

    fn len(&self) -> u32 {
        self.end - self.start
    }
}

/// Count reads in various windows of length `params.window_size` between ~ `interval.start() + params.boundary_size`
/// and ~ `interval.end() - params.boundary_size`.
fn count_reads<'a>(
    alignments: impl Iterator<Item = &'a LightAlignment>,
    interval: &Interval,
    params: &ReadDepthParams,
) -> Vec<WindowCounts> {
    assert!(interval.len() >= params.window_size + 2 * params.boundary_size, "Input interval is too short!");
    let n_windows = (f64::from(interval.len() - 2 * params.boundary_size) / f64::from(params.window_size))
        .floor() as u32;
    let sum_len = n_windows * params.window_size;
    let interval_start = interval.start();
    let start = interval_start + (interval.len() - sum_len) / 2;
    let end = start + sum_len;
    let mut windows: Vec<_> = (0..n_windows)
        .map(|i| WindowCounts::new(start + i * params.window_size, start + (i + 1) * params.window_size))
        .collect();

    for aln in alignments {
        let middle = aln.interval().middle();
        if start <= middle && middle < end {
            let ix = (middle - start) / params.window_size;
            windows[ix as usize].depth[aln.read_end().ix()] += 1;
        }
    }
    windows
}

/// Discard windows, where the region <left_padding><window><right padding> contains
/// - unknown nucleotides (Ns),
/// - frequent k-mers (average k-mer frequency > max_kmer_freq).
///
/// All windows must have the same length! (is not checked)
///
/// Consumes windows, and returns a new vector of selected windows + a vector of GC-content values for each window.
fn filter_windows(
    windows: Vec<WindowCounts>,
    seq_shift: u32,
    ref_seq: &[u8],
    kmer_counts_coll: &KmerCounts,
    params: &ReadDepthParams,
    mut dbg_writer: impl Write,
) -> io::Result<(Vec<WindowCounts>, Vec<f64>)>
{
    writeln!(dbg_writer, "start\tgc\tkmer_freq\tkeep\tdepth1\tdepth2")?;
    log::debug!("    Total windows:   {:7}", windows.len());
    let mut have_ns = 0;
    let mut have_common_kmers = 0;

    let k = kmer_counts_coll.k();
    let kmer_counts = kmer_counts_coll.get_first();
    let window_padding = params.window_padding;
    let neigh_len = (windows[0].len() + 2 * window_padding) as usize;
    let neigh_kmers = neigh_len + 1 - k as usize;
    let norm_fct = 1.0 / (neigh_kmers as f64);

    let mut sel_windows = Vec::with_capacity(windows.len());
    let mut gc_contents = Vec::with_capacity(windows.len());
    for window in windows.into_iter() {
        let start_ix = (window.start - window_padding - seq_shift) as usize;
        let window_seq = &ref_seq[start_ix..start_ix + neigh_len];
        let mut keep = true;
        let mut mean_kmer_freq = f64::NAN;
        let mut gc_content = f64::NAN;
        if seq::has_n(window_seq) {
            have_ns += 1;
            keep = false
        } else {
            mean_kmer_freq = norm_fct * kmer_counts[start_ix..start_ix + neigh_kmers]
                .iter().copied().map(f64::from).sum::<f64>();
            gc_content = seq::gc_content(window_seq);
            if mean_kmer_freq > params.max_kmer_freq {
                have_common_kmers += 1;
                keep = false;
            }
        }

        writeln!(dbg_writer, "{}\t{:.1}\t{:.1}\t{}\t{}\t{}", window.start, gc_content, mean_kmer_freq,
            if keep { 'T' } else { 'F' }, window.depth[0], window.depth[1])?;
        if keep {
            sel_windows.push(window);
            gc_contents.push(gc_content);
        }
    }
    log::debug!("    Remove {} windows with Ns, {} windows with common k-mers", have_ns, have_common_kmers);
    log::debug!("    After filtering: {:7}", sel_windows.len());
    Ok((sel_windows, gc_contents))
}

/// In total, GC-content falls into 101 bins (0..=100).
pub const GC_BINS: usize = 101;

/// For each GC-content in 0..=100, returns start & end indices,
/// where gc_contents[start..end] lie in [GC - 0.5, GC + 0.5).
/// gc_contents must be sorted in advance.
fn find_gc_bins(gc_contents: &[f64]) -> Vec<(usize, usize)> {
    let n = gc_contents.len();
    let mut res = Vec::with_capacity(GC_BINS);
    let mut i = 0;
    for gc in 0..GC_BINS {
        let gc = gc as f64;
        let j = bisect::right_at(gc_contents, &(gc + 0.5), i, n);
        res.push((i, j));
        debug_assert!(i == j || (gc - 0.5 <= gc_contents[i] && gc_contents[j - 1] < gc + 0.5));
        i = j;
    }
    debug_assert_eq!(i, n);
    res
}

/// Predicts read depth means and variance for each GC-content from 0 to 100.
/// GC-contents must be sorted. gc_contents and depth must have the same length and have the go in the same order.
///
/// Mean values: use LOESS with all observations, without weights, and using `frac_windows`.
///
/// Variance: use LOESS with frac = 1, where observations are all GC-content values with at least 10 observations.
/// Y values are read depth variance in each particular GC-content bin.
fn predict_mean_var(gc_contents: &[f64], gc_bins: &[(usize, usize)], depth: &[f64], frac_windows: f64)
    -> (Vec<f64>, Vec<f64>)
{
    const VAR_MIN_WINDOWS: usize = 10;
    let n = depth.len();
    let m = gc_bins.len();
    let xout: Vec<f64> = (0..u32::try_from(m).unwrap()).map(Into::into).collect();
    let means = Loess::new(frac_windows, 1).set_xout(xout.clone()).calculate(&gc_contents, &depth);

    let mut x = Vec::with_capacity(m);
    let mut y = Vec::with_capacity(m);
    let mut w = Vec::with_capacity(m);

    for (gc, &(i, j)) in gc_bins.iter().enumerate() {
        if j - i >= VAR_MIN_WINDOWS {
            x.push(gc as f64);
            y.push(F64Ext::variance(&depth[i..j], None));
            w.push(((j - i) as f64 / n as f64).sqrt());
        }
    }
    let vars = Loess::new(1.0, 1).set_xout(xout).calculate_weighted(&x, &y, &w);
    (means, vars)
}

/// At very small and very large GC-content values, there may be not enough observations
/// to well estimate read depth mean and value.
///
/// To counter this, find GC-content values such that to the left of it there are < min_obs observations.
/// Fill unavailable values according to read depth parameter value `tail_var_mult`.
fn blur_boundary_values(means: &[f64], vars: &[f64], gc_bins: &[(usize, usize)], params: &ReadDepthParams
) -> (Vec<f64>, Vec<f64>)
{
    let min_obs = params.min_tail_obs;
    let n = gc_bins.len();
    let m = gc_bins[gc_bins.len() - 1].1;

    // Index of the first GC-content value, such that there are at least `min_obs` windows to the left.
    let left_ix = gc_bins.iter().position(|(_, end)| *end >= min_obs).unwrap_or(n);
    // Index of the last GC-content value, such that there are at least `min_obs` windows to the right.
    let right_ix = n - gc_bins.iter().rev().position(|(start, _)| m - *start >= min_obs).unwrap_or(n);
    assert!(left_ix < right_ix, "Too few windows to calculate read depth!");
    log::debug!("    Few windows (< {}) with GC-content < {} or > {}, bluring distributions",
        min_obs, left_ix, right_ix);

    let mut blurred_means = means.to_vec();
    let mut blurred_vars = vars.to_vec();
    for i in 0..left_ix {
        blurred_means[i] = means[left_ix];
        let mult = 1.0 + (left_ix - i) as f64 * params.tail_var_mult;
        blurred_vars[i] = (mult * vars[left_ix]).max(vars[i]);
    }
    for i in right_ix + 1..n {
        blurred_means[i] = means[right_ix];
        let mult = 1.0 + (i - right_ix) as f64 * params.tail_var_mult;
        blurred_vars[i] = (mult * vars[right_ix]).max(vars[i]);
    }
    (blurred_means, blurred_vars)
}

/// Read depth parameters.
#[derive(Clone, Debug)]
pub struct ReadDepthParams {
    /// Specie ploidy (2 in most cases).
    pub ploidy: u8,

    /// Calculate background per windows of this size.
    /// Default: 100.
    pub window_size: u32,

    /// Calculate GC-content and k-mer frequencies based on the window of size `window_size + 2 * window_padding`.
    /// Default: 100.
    pub window_padding: u32,

    /// Ignore left-most and right-most `boundary_size` base-pairs. Must not be smaller than `window_padding`.
    /// Default: 1000.
    pub boundary_size: u32,

    /// Filter windows with too many frequent k-mers (average frequency must be at most this value).
    /// Default: 1.2.
    pub max_kmer_freq: f64,

    /// When calculating read depth averages for various GC-content values, use `frac_windows` fraction
    /// across all windows. For example, if frac_windows is 0.1 (default),
    /// use 10% of all observations with the most similar GC-content values.
    pub frac_windows: f64,

    /// Calculate read depth parameters differently for the highest and lowest GC-content values.
    /// First, select such GC-content values L & R such that there `< min_tail_obs` windows with GC-content `< L`
    /// and `< min_tail_obs` windows with GC-content `> R`. Default: 100.
    ///
    /// Next, for all GC-content < L use mean as for L, and multiple variance by
    /// `1 + GC_diff * tail_var_mult` where `GC_diff` is the difference between GC-content of the window and L.
    /// For example, window with GC-content `L-2` will have `mean[L-2] = mean[L]` and
    /// `var[L-2] = var[L] * (1 + 2 * tail_var_mult)`.
    ///
    /// Repeat the same procedure for GC-content values > R.
    pub min_tail_obs: usize,

    /// See `min_tail_obs`. Default: 0.02.
    pub tail_var_mult: f64,
}

impl Default for ReadDepthParams {
    fn default() -> Self {
        Self {
            ploidy: 2,
            window_size: 100,
            window_padding: 100,
            boundary_size: 1000,
            max_kmer_freq: 1.2,

            frac_windows: 0.5,
            min_tail_obs: 100,
            tail_var_mult: 0.02,
        }
    }
}

impl ReadDepthParams {
    /// Validate all parameter values.
    pub fn validate(&self) -> Result<(), Error> {
        validate_param!(self.ploidy > 0, "Ploidy cannot be zero");
        validate_param!(self.boundary_size >= self.window_padding,
            "Edge padding ({}) must not be smaller than window padding ({})", self.boundary_size, self.window_padding);
        if self.boundary_size < self.window_size {
            log::warn!("Edge padding ({}) is smaller than the window size ({}), consider using a larger value",
                self.boundary_size, self.window_size);
        }

        validate_param!(self.max_kmer_freq >= 1.0,
            "Maximum k-mer frequency ({}) must be at least 1", self.max_kmer_freq);
        validate_param!(0.0 < self.frac_windows && self.frac_windows <= 1.0,
            "Fraction of windows ({}) must be within (0, 1]", self.frac_windows);

        if self.min_tail_obs < 100 {
            log::warn!("Number of windows with extreme GC-content ({}) is too small, consider using a larger value",
                self.min_tail_obs);
        }
        validate_param!(self.tail_var_mult >= 0.0, "Extreme GC-content variance factor ({}) must be non-negative",
            self.tail_var_mult);
        if self.tail_var_mult >= 0.5 {
            log::warn!("Extreme GC-content variance factor ({}) is too large, consider using a smaller value",
                self.tail_var_mult);
        }
        Ok(())
    }
}

fn write_summary(
    gc_bins: &[(usize, usize)],
    depth: &[f64],
    loess_means: &[f64],
    loess_vars: &[f64],
    blurred_means: &[f64],
    blurred_vars: &[f64],
    distrs: &[NBinom],
    mut dbg_writer: impl Write,
) -> io::Result<()>
{
    writeln!(dbg_writer, "gc\tnwindows\tmean\tvar\tloess_mean\tloess_var\
        \tblurred_mean\tblurred_var\tnbinom_n\tnbinom_p")?;
    for (gc, &(i, j)) in gc_bins.iter().enumerate() {
        let n = j - i;
        let mean = if n > 0 { F64Ext::mean(&depth[i..j]) } else { f64::NAN };
        let var = if n > 1 { F64Ext::variance(&depth[i..j], Some(mean)) } else { f64::NAN };
        write!(dbg_writer, "{gc}\t{n}\t{mean:.5}\t{var:.5}\t")?;
        write!(dbg_writer, "{:.5}\t{:.5}\t", loess_means[gc], loess_vars[gc])?;
        write!(dbg_writer, "{:.5}\t{:.5}\t", blurred_means[gc], blurred_vars[gc])?;
        writeln!(dbg_writer, "{}\t{}", distrs[gc].n(), distrs[gc].p())?;
    }
    Ok(())
}

/// Background read depth distributions for each GC-content value.
#[derive(Clone)]
pub struct ReadDepth {
    /// Read depth is calculated per windows with this size.
    window_size: u32,
    /// Specie ploidy (2 in most cases).
    ploidy: u8,
    /// For each window, add `window_padding` to the left and right before calculating GC-content.
    window_padding: u32,
    // Read depth distribution for each GC-content in 0..=100.
    distributions: Vec<NBinom>,
}

impl ReadDepth {
    /// Estimates read depth from primary alignments, mapped to the `interval` with sequence `ref_seq`.
    /// Ignore reads with alignment probability < `params.a` and with insert size > `max_insert_size`.
    ///
    /// Write debug information if `out_dir` is Some.
    pub fn estimate<'a>(
        alignments: impl Iterator<Item = &'a LightAlignment>,
        interval: &Interval,
        ref_seq: &[u8],
        kmer_counts: &KmerCounts,
        params: &ReadDepthParams,
        subsampling_rate: f64,
        is_paired_end: bool,
        out_dir: Option<&Path>,
    ) -> io::Result<Self>
    {
        log::info!("Estimating read depth (ploidy {}, subsampling rate {})", params.ploidy, subsampling_rate);
        assert_eq!(interval.len() as usize, ref_seq.len(),
            "ReadDepth: interval and reference sequence have different lengths!");

        let windows = count_reads(alignments, interval, params);
        let (windows, gc_contents) = if let Some(dir) = out_dir {
            let dbg_writer = ext::sys::create_gzip(&dir.join("counts.csv.gz"))?;
            filter_windows(windows, interval.start(), &ref_seq, &kmer_counts, params, dbg_writer)?
        } else {
            filter_windows(windows, interval.start(), &ref_seq, &kmer_counts, params, io::sink())?
        };
        assert!(windows.len() > 0, "ReadDepth: no applicable windows!");

        let ixs = VecExt::argsort(&gc_contents);
        let gc_contents = VecExt::reorder(&gc_contents, &ixs);
        let gc_bins = find_gc_bins(&gc_contents);

        let depth: Vec<f64> = ixs.iter().map(|&i| f64::from(windows[i].depth[0])).collect();
        let (loess_means, loess_vars) = predict_mean_var(&gc_contents, &gc_bins, &depth, params.frac_windows);
        let (blurred_means, blurred_vars) = blur_boundary_values(&loess_means, &loess_vars, &gc_bins, params);

        let nbinom_mul = 1.0 / (subsampling_rate * f64::from(params.ploidy));
        let distributions: Vec<_> = blurred_means.iter().zip(blurred_vars.iter())
            .map(|(&m, &v)| NBinom::estimate(m, v.max(m * 1.00001)).mul(nbinom_mul))
            .collect();

        const GC_VAL: usize = 40;
        let logging_distr = distributions[GC_VAL].mul(f64::from(params.ploidy) * if is_paired_end { 2.0 } else { 1.0 });
        log::info!("    Read depth mean = {:.2},  variance: {:.2}  (per {} bp window at GC-content {})",
            logging_distr.mean(), logging_distr.variance(), params.window_size, GC_VAL);

        if let Some(dir) = out_dir {
            let dbg_writer = ext::sys::create_gzip(&dir.join("depth.csv.gz"))?;
            write_summary(&gc_bins, &depth, &loess_means, &loess_vars, &blurred_means, &blurred_vars, &distributions,
                dbg_writer)?;
        }
        Ok(Self {
            ploidy: params.ploidy,
            window_size: params.window_size,
            window_padding: params.window_padding,
            distributions,
        })
    }

    /// Returns window size.
    pub fn window_size(&self) -> u32 {
        self.window_size
    }

    /// Returns window padding - added to both sides of the window to calculate GC-content and k-mer frequencies.
    pub fn window_padding(&self) -> u32 {
        self.window_padding
    }

    /// Returns read depth distribution at GC-content `gc_content` (between 0 and 100).
    pub fn depth_distribution(&self, gc_content: u8) -> &NBinom {
        &self.distributions[usize::from(gc_content)]
    }

    /// Returns iterator over read-depth distributions at GC-contents between 0 and 100.
    pub fn distributions(&self) -> std::slice::Iter<'_, NBinom> {
        self.distributions.iter()
    }
}

impl JsonSer for ReadDepth {
    fn save(&self) -> json::JsonValue {
        let n_params: Vec<f64> = self.distributions.iter().map(|distr| distr.n()).collect();
        let p_params: Vec<f64> = self.distributions.iter().map(|distr| distr.p()).collect();
        json::object!{
            ploidy: self.ploidy,
            window: self.window_size,
            window_padding: self.window_padding,
            n: &n_params as &[f64],
            p: &p_params as &[f64],
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> ploidy (as_u8), window (as_usize), window_padding (as_u32));
        let mut n_params = vec![0.0; window + 1];
        parse_f64_arr(obj, "n", &mut n_params)?;
        let mut p_params = vec![0.0; window + 1];
        parse_f64_arr(obj, "p", &mut p_params)?;
        Ok(Self {
            ploidy, window_padding,
            window_size: window as u32,
            distributions: n_params.into_iter().zip(p_params).map(|(n, p)| NBinom::new(n, p)).collect(),
        })
    }
}
