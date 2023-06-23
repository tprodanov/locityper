use std::{
    io::{self, Write},
    path::Path,
};
use crate::{
    seq::{
        self, Interval,
        kmers::KmerCounts,
        aln::LightAlignment,
    },
    algo::{bisect, loess::Loess},
    ext::{
        self,
        vec::{VecExt, F64Ext},
    },
    math::{
        distr::{
            WithMoments,
            nbinom::{NBinom, RegularizedEstimator},
        },
    },
    err::{Error, validate_param},
    model::windows::WindowGetter,
};
use super::{
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
}

/// Count reads in various windows of length `params.window_size` between ~ `interval.start() + params.boundary_size`
/// and ~ `interval.end() - params.boundary_size`.
fn count_reads<'a>(
    alignments: impl Iterator<Item = &'a LightAlignment>,
    interval: &Interval,
    params: &ReadDepthParams,
) -> Vec<WindowCounts> {
    assert!(interval.len() >= params.window_size + 2 * params.boundary_size, "Input interval is too short!");
    let n_windows = (interval.len() - 2 * params.boundary_size) / params.window_size;
    let sum_len = n_windows * params.window_size;
    let interval_start = interval.start();
    let start = interval_start + (interval.len() - sum_len) / 2;
    let mut windows: Vec<_> = (0..n_windows)
        .map(|i| WindowCounts::new(start + i * params.window_size, start + (i + 1) * params.window_size))
        .collect();
    let window_getter = WindowGetter::new(start, start + sum_len, params.window_size);

    for aln in alignments {
        let (aln_start, aln_end) = aln.interval().range();
        let (start_ix, end_ix) = window_getter.covered_middle(aln_start, aln_end);
        let read_end = aln.read_end().ix();
        for ix in start_ix..end_ix {
            windows[ix as usize].depth[read_end] += 1;
        }
    }
    windows
}

/// Discard windows, where the region <left_padding><window><right padding> (sum length = neighb_size) contains
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
    let window_padding = (params.neighb_size - params.window_size) / 2;
    let neighb_size = params.neighb_size as usize;
    let neighb_kmers = neighb_size + 1 - k as usize;
    let norm_fct = 1.0 / (neighb_kmers as f64);

    let mut sel_windows = Vec::with_capacity(windows.len());
    let mut gc_contents = Vec::with_capacity(windows.len());
    for window in windows.into_iter() {
        let start_ix = (window.start - window_padding - seq_shift) as usize;
        let window_seq = &ref_seq[start_ix..start_ix + neighb_size];
        let mut keep = true;
        let mut mean_kmer_freq = f64::NAN;
        let mut gc_content = f64::NAN;
        if seq::has_n(window_seq) {
            have_ns += 1;
            keep = false
        } else {
            mean_kmer_freq = norm_fct * kmer_counts[start_ix..start_ix + neighb_kmers]
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

    /// Calculate background read depth per windows of this size.
    /// Default: 100.
    pub window_size: u32,
    /// Calculate GC-content and k-mer frequencies based on the window neighbourhood of this size (includes window).
    /// Default: 300.
    pub neighb_size: u32,
    /// Ignore left-most and right-most `boundary_size` base-pairs. Must be >= `neighb_size`.
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
            neighb_size: 300,
            boundary_size: 1000,
            max_kmer_freq: 1.2,

            frac_windows: 0.5,
            min_tail_obs: 100,
            tail_var_mult: 0.02,
        }
    }
}

impl ReadDepthParams {
    /// Compare window, neighbourhood and boundary sizes.
    /// This function is static, as it is used elsewhere.
    pub fn validate_sizes(window: u32, neighb: u32, boundary: u32) -> Result<(), Error> {
        validate_param!(0 < window,
            "Window size ({}) must be positive", window);
        validate_param!(window <= neighb,
            "Neighbourhood size ({}) must be at least window size ({})", neighb, window);
        validate_param!((neighb - window + 1) / 2 <= boundary,
            "Region boundary ({}) is too small compared to window ({}) and neighbourhood ({}) sizes",
            boundary, window, neighb);
        Ok(())
    }

    /// Validate all parameter values.
    pub fn validate(&self) -> Result<(), Error> {
        validate_param!(self.ploidy > 0, "Ploidy cannot be zero");
        Self::validate_sizes(self.window_size, self.neighb_size, self.boundary_size)?;

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
    /// Specie ploidy (2 in most cases).
    ploidy: u8,
    /// Read depth is calculated per windows with this size.
    window_size: u32,
    /// Calculate GC-content per window neighbourhood of this size (includes window).
    neighb_size: u32,
    // Read depth distribution for each GC-content in 0..=100.
    distributions: Vec<NBinom>,
}

/// Estimate read depth distributions based on observed mean and variance,
/// as well on the subsampling rate and ploidy.
///
/// Output distributions are calculated for ploidy 1 and first read ends.
///
/// Although it is possible to estimate n & p directly, sometimes $p$ becomes close to 1 and $n$ grows very large.
/// Due to this, we perform L1-regularization on $n$, and estimate parameters numerically using Nelder-Mead algorithm.
fn estimate_nbinoms(means: &[f64], vars: &[f64], rate: f64, ploidy: f64) -> Vec<NBinom> {
    let mut estimator = RegularizedEstimator::default();
    estimator.set_subsampling_rate(rate).set_lambda(1e-5);

    means.iter().zip(vars).map(|(&m, &v)| {
        estimator.estimate(m, v).mul(1.0 / ploidy)
    }).collect()
}

impl ReadDepth {
    /// Estimates read depth from primary alignments, mapped to the `interval` with sequence `ref_seq`.
    /// Ignore reads with alignment probability < `params.a` and with insert size > `max_insert_size`.
    ///
    /// Write debug information if `out_dir` is Some.
    pub fn estimate<'a>(
        alignments: &[&'a LightAlignment],
        interval: &Interval,
        ref_seq: &[u8],
        kmer_counts: &KmerCounts,
        params: &ReadDepthParams,
        subsampling_rate: f64,
        is_paired_end: bool,
        out_dir: Option<&Path>,
    ) -> io::Result<Self>
    {
        log::info!("Estimating read depth from {} reads", alignments.len());
        log::info!("    ploidy {}, subsampling rate {}", params.ploidy, subsampling_rate);
        assert_eq!(interval.len() as usize, ref_seq.len(),
            "ReadDepth: interval and reference sequence have different lengths!");

        let windows = count_reads(alignments.iter().copied(), interval, params);
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
        let ploidy = f64::from(params.ploidy);
        let distributions = estimate_nbinoms(&blurred_means, &blurred_vars, subsampling_rate, ploidy);

        const GC_VAL: usize = 40;
        let logging_distr = distributions[GC_VAL].mul(ploidy * if is_paired_end { 2.0 } else { 1.0 });
        log::info!("    Read depth mean = {:.2},  variance: {:.2}  (at GC-content {})",
            logging_distr.mean(), logging_distr.variance(), GC_VAL);

        if let Some(dir) = out_dir {
            let dbg_writer = ext::sys::create_gzip(&dir.join("depth.csv.gz"))?;
            write_summary(&gc_bins, &depth, &loess_means, &loess_vars, &blurred_means, &blurred_vars, &distributions,
                dbg_writer)?;
        }
        Ok(Self {
            ploidy: params.ploidy,
            window_size: params.window_size,
            neighb_size: params.neighb_size,
            distributions,
        })
    }

    /// Returns window size.
    pub fn window_size(&self) -> u32 {
        self.window_size
    }

    /// Returns neighbourhood size.
    pub fn neighb_size(&self) -> u32 {
        self.neighb_size
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
            neighb: self.neighb_size,
            n: &n_params as &[f64],
            p: &p_params as &[f64],
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> ploidy (as_u8), window (as_u32), neighb (as_u32));
        let mut n_params = vec![0.0; GC_BINS];
        let mut p_params = vec![0.0; GC_BINS];
        parse_f64_arr(obj, "n", &mut n_params)?;
        parse_f64_arr(obj, "p", &mut p_params)?;
        Ok(Self {
            ploidy,
            window_size: window,
            neighb_size: neighb,
            distributions: n_params.into_iter().zip(p_params).map(|(n, p)| NBinom::new(n, p)).collect(),
        })
    }
}
