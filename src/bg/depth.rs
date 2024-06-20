use std::{
    io::{self, Write},
    path::Path,
};
use super::{
    Windows,
    ser::{JsonSer, parse_f64_arr, json_get},
};
use crate::{
    seq::aln::Alignment,
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
    err::{Error, validate_param, add_path},
    model::windows::WindowGetter,
};

pub trait DepthDistr {
    fn window_size(&self) -> u32;

    fn ln_prob(&self, gc_content: u8, depth: u32) -> f64;
}

/// Count reads in various windows of length `window_size`.
fn count_reads<'a>(
    alignments: impl Iterator<Item = &'a Alignment>,
    n_windows: usize,
    window_getter: &WindowGetter,
) -> Vec<[u32; 2]> {
    let mut depth = vec![[0, 0]; n_windows];
    for aln in alignments {
        if let Some(window) = window_getter.middle_window(aln.interval().middle()) {
            depth[window as usize][aln.read_end().ix()] += 1;
        }
    }
    depth
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
            y.push(F64Ext::variance(&depth[i..j]));
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
    /// None: automatic.
    pub window_size: Option<u32>,
    /// Ignore left-most and right-most `boundary_size` bp.
    /// Default: 1000.
    pub boundary_size: u32,

    /// Filter windows with too many frequent k-mers:
    /// where less than this percetage of k-mers are unique.
    /// Default: 90.
    pub uniq_kmer_perc: f64,

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
            window_size: None,
            boundary_size: 1000,
            uniq_kmer_perc: 90.0,

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
        validate_param!(1.0 < self.uniq_kmer_perc && self.uniq_kmer_perc <= 100.0,
            "Unique k-mer percentile ({}) must be within (1, 100].", self.uniq_kmer_perc);
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
        let (mean, var) = F64Ext::mean_variance_or_nan(&depth[i..j]);
        write!(dbg_writer, "{gc}\t{n}\t{mean:.5}\t{var:.5}\t")?;
        write!(dbg_writer, "{:.5}\t{:.5}\t", loess_means[gc], loess_vars[gc])?;
        write!(dbg_writer, "{:.5}\t{:.5}\t", blurred_means[gc], blurred_vars[gc])?;
        writeln!(dbg_writer, "{}\t{}", distrs[gc].n(), distrs[gc].p())?;
    }
    Ok(())
}

fn write_summary_without_gc(
    count: usize,
    mean: f64,
    var: f64,
    distr: &NBinom,
    mut dbg_writer: impl Write,
) -> io::Result<()>
{
    writeln!(dbg_writer, "nwindows\tmean\tvar\tnbinom_n\tnbinom_p")?;
    writeln!(dbg_writer, "{}\t{:.5}\t{:.5}\t{}\t{}", count, mean, var, distr.n(), distr.p())
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

/// Estimate read depth distributions based on observed mean, variance, and ploidy.
///
/// Output distributions are calculated for ploidy 1 and first read ends.
///
/// Although it is possible to estimate n & p directly, sometimes $p$ becomes close to 1 and $n$ grows very large.
/// Due to this, we perform L1-regularization on $n$, and estimate parameters numerically using Nelder-Mead algorithm.
fn estimate_nbinoms<'a>(
    means: &'a [f64],
    vars: &'a [f64],
    ploidy: f64,
) -> impl Iterator<Item = NBinom> + 'a
{
    let mut estimator = RegularizedEstimator::default();
    estimator.set_lambda(1e-5);

    means.iter().zip(vars).map(move |(&m, &v)| {
        estimator.estimate(m, v).mul(1.0 / ploidy)
    })
}

/// Returns two vectors (first with read1 depth values, and second with GC-content values).
/// Both vectors are sorted by GC-content.
fn get_depth_and_gc(
    depth: &[[u32; 2]],
    windows: &Windows,
    mut dbg_writer: impl Write,
) -> io::Result<(Vec<f64>, Vec<f64>)>
{
    let mut depth1 = Vec::with_capacity(depth.len());
    let mut gc_contents = Vec::with_capacity(depth.len());
    for (i, (&[d1, d2], window)) in depth.iter().zip(windows.iter()).enumerate() {
        writeln!(dbg_writer, "{}\t{}\t{}", i, d1, d2)?;
        if window.keep() {
            depth1.push(f64::from(d1));
            gc_contents.push(window.gc());
        }
    }
    let ixs = VecExt::argsort(&gc_contents);
    Ok((VecExt::reorder(&depth1, &ixs), VecExt::reorder(&gc_contents, &ixs)))
}

impl ReadDepth {
    /// Estimates read depth from primary alignments.
    ///
    /// Write debug information if `out_dir` is Some.
    pub fn estimate<'a>(
        alignments: &[&'a Alignment],
        windows: &Windows,
        params: &ReadDepthParams,
        is_paired_end: bool,
        seq_info: &super::SequencingInfo,
        out_dir: Option<&Path>,
    ) -> Result<Self, Error>
    {
        let depth = count_reads(alignments.iter().copied(), windows.len(), windows.window_getter());
        let (depth, gc_contents) = if let Some(dir) = out_dir {
            let dbg_filename1 = dir.join("window_depth.csv.gz");
            let mut dbg_writer1 = ext::sys::create_gzip(&dbg_filename1)?;
            writeln!(dbg_writer1, "window_ix\tdepth1\tdepth2").map_err(add_path!(dbg_filename1))?;
            get_depth_and_gc(&depth, windows, dbg_writer1).map_err(add_path!(dbg_filename1))?
        } else {
            get_depth_and_gc(&depth, windows, io::sink()).map_err(add_path!(!))?
        };

        let gc_bins = find_gc_bins(&gc_contents);
        let ploidy = f64::from(params.ploidy);
        let dbg_filename2 = out_dir.map(|dirname| dirname.join("depth.csv.gz"));
        let dbg_writer2 = dbg_filename2.as_ref().map(|filename| ext::sys::create_gzip(filename)).transpose()?;
        let distributions: Vec<NBinom>;
        if seq_info.technology().has_gc_bias() {
            let (loess_means, loess_vars) = predict_mean_var(&gc_contents, &gc_bins, &depth, params.frac_windows);
            let (blurred_means, blurred_vars) = blur_boundary_values(&loess_means, &loess_vars, &gc_bins, params);
            distributions = estimate_nbinoms(&blurred_means, &blurred_vars, ploidy).collect();
            if let Some(writer) = dbg_writer2 {
                write_summary(&gc_bins, &depth, &loess_means, &loess_vars, &blurred_means, &blurred_vars,
                    &distributions, writer).map_err(add_path!(dbg_filename2.as_ref().unwrap()))?;
            }
        } else {
            let (mean, var) = F64Ext::mean_variance(&depth);
            let distr = estimate_nbinoms(&[mean], &[var], ploidy).next().unwrap();
            if let Some(writer) = dbg_writer2 {
                write_summary_without_gc(depth.len(), mean, var, &distr, writer)
                    .map_err(add_path!(dbg_filename2.unwrap()))?;
            }
            distributions = vec![distr; GC_BINS];
        }

        const GC_VAL: usize = 40;
        let logging_distr = distributions[GC_VAL].mul(ploidy * if is_paired_end { 2.0 } else { 1.0 });
        log::info!("    Read depth mean = {:.2},  variance: {:.2}  (at GC-content {})",
            logging_distr.mean(), logging_distr.variance(), GC_VAL);

        Ok(Self {
            ploidy: params.ploidy,
            window_size: windows.window_size(),
            neighb_size: windows.neighb_size(),
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

    /// Assume that read depth changes by a factor of rate.
    pub fn mul_depth(&mut self, rate: f64) {
        for distr in self.distributions.iter_mut() {
            *distr = distr.mul(rate);
        }
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
        json_get!(obj => ploidy (as_u8), window (as_u32), neighb (as_u32));
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
