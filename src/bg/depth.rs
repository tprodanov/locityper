use std::cmp::max;
use htslib::bam::{
    Record,
    ext::BamRecordExtensions,
};
use crate::{
    seq::{
        cigar::{Operation, CigarItem},
        interv::Interval,
        seq::gc_content,
    },
    algo::{
        nbinom::{NBinom, CachedDistr},
        vec_ext::*,
        loess::Loess,
    },
    bg::ser::{JsonSer, LoadError, parse_f64_arr},
};

pub trait DepthDistr {
    fn window_size(&self) -> u32;

    fn ln_prob(&self, gc_content: u8, depth: u32) -> f64;
}

/// Returns maximal of the left and right clipping.
fn soft_clipping(record: &Record) -> u32 {
    let raw_cigar = record.raw_cigar();
    let first_item = CigarItem::from_u32(raw_cigar[0]);
    let left = match first_item.operation() {
        Operation::Soft => first_item.len(),
        _ => 0,
    };

    let n = raw_cigar.len();
    if n == 1 {
        return left;
    }
    let last_item = CigarItem::from_u32(raw_cigar[n - 1]);
    match last_item.operation() {
        Operation::Soft => max(left, last_item.len()),
        _ => left,
    }
}

/// Number of various read types per window.
#[derive(Clone, Debug)]
struct WindowCounts {
    start: u32,
    end: u32,

    depth1: u32,
    depth2: u32,

    low_mapq_reads: u32,
    clipped_reads: u32,
    unpaired_reads: u32,
}

impl WindowCounts {
    #[inline]
    fn new(start: u32, end: u32) -> Self {
        WindowCounts {
            start, end,
            depth1: 0,
            depth2: 0,
            low_mapq_reads: 0,
            clipped_reads: 0,
            unpaired_reads: 0,
        }
    }

    fn add_read(&mut self, record: &Record, max_insert_size: u32) {
        const MAX_CLIPPING: u32 = 4;
        const MIN_MAPQ: u8 = 30;

        debug_assert!(!record.is_unmapped(), "WindowCounts: Cannot add unmapped read!");
        if record.is_last_in_template() {
            self.depth2 += 1;
        } else {
            self.depth1 += 1;
        }
        if record.mapq() < MIN_MAPQ {
            self.low_mapq_reads += 1;
        }
        if record.is_paired() {
            if !record.is_proper_pair() || record.is_mate_unmapped() || record.tid() != record.mtid()
                    || record.insert_size().abs() as u32 > max_insert_size {
                self.unpaired_reads += 1;
            }
        }
        if soft_clipping(record) > MAX_CLIPPING {
            self.clipped_reads += 1;
        }
    }

    /// Returns the total read depth across first and second read mates.
    #[inline]
    fn total_depth(&self) -> u32 {
        self.depth1 + self.depth2
    }

    /// Returns the rate of low-mapq reads (low_mapq_reads / total_depth).
    fn low_mapq_rate(&self) -> f64 {
        self.low_mapq_reads as f64 / max(self.total_depth(), 1) as f64
    }

    /// Returns the rate of clipped reads (clipped_reads / total_depth).
    fn clipped_rate(&self) -> f64 {
        self.clipped_reads as f64 / max(self.total_depth(), 1) as f64
    }

    /// Returns the rate of unpaired reads (unpaired_reads / total_depth).
    fn unpaired_rate(&self) -> f64 {
        self.unpaired_reads as f64 / max(self.total_depth(), 1) as f64
    }

    /// Does the read satisfies set constraints?
    fn satisfies(&self, lim: &LimitingValues) -> bool {
        let dp = self.total_depth() as f64;
        dp == 0.0 ||
            (dp <= lim.depth
            && self.low_mapq_reads as f64 / dp <= lim.low_mapq_rate
            && self.clipped_reads as f64 / dp <= lim.clipped_rate
            && self.unpaired_reads as f64 / dp <= lim.unpaired_rate)
    }
}

/// Count reads in various windows of length `params.window_size` between ~ `interval.start() + params.edge_padding`
/// and ~ `interval.end() - params.edge_padding`.
fn count_reads<'a>(
        interval: &Interval,
        records: impl Iterator<Item = &'a Record>,
        params: &ReadDepthParams,
        max_insert_size: u32,
    ) -> Vec<WindowCounts>
{
    assert!(interval.len() >= params.window_size + 2 * params.edge_padding, "Input interval is too short!");
    let n_windows = ((interval.len() - 2 * params.edge_padding) as f64 / params.window_size as f64).floor() as u32;
    let sum_len = n_windows * params.window_size;
    let shift = interval.start() + (interval.len() - sum_len) / 2;
    let end = shift + sum_len;
    let mut windows: Vec<_> = (0..n_windows)
        .map(|i| WindowCounts::new(shift + i * params.window_size, shift + (i + 1) * params.window_size))
        .collect();

    for record in records {
        if (record.flags() & 3844) == 0 {
            let middle = (record.reference_start() + record.reference_end()) as u32 >> 1;
            if shift <= middle && middle < end {
                let ix = (middle - shift) / params.window_size;
                windows[ix as usize].add_read(record, max_insert_size);
            }
        }
    }
    windows
}

/// Filter windows that have too large read depth, too many low-mapq reads, etc.
#[derive(Debug, Clone)]
struct LimitingValues {
    depth: f64,
    low_mapq_rate: f64,
    clipped_rate: f64,
    unpaired_rate: f64,
}

impl LimitingValues {
    /// Set max_depth to 2.0 * <0.999 quant. of read depth>,
    /// Max low_mapq_rate and others to `filt_quantile` of the corresponding rate across input windows.
    fn new(windows: &[WindowCounts], filt_quantile: f64) -> Self {
        const READ_DP_QUANTILE: f64 = 0.999;
        const READ_DP_MULT: f64 = 2.0;
        /// Maximum depth should be at least 100.
        const DEF_MAX_DEPTH: f64 = 100.0;

        Self {
            depth: (READ_DP_MULT * IterExt::quantile(
                windows.iter().map(|window| window.total_depth() as f64), READ_DP_QUANTILE)).max(DEF_MAX_DEPTH),
            low_mapq_rate: IterExt::quantile(windows.iter().map(WindowCounts::low_mapq_rate), filt_quantile),
            clipped_rate: IterExt::quantile(windows.iter().map(WindowCounts::clipped_rate), filt_quantile),
            unpaired_rate: IterExt::quantile(windows.iter().map(WindowCounts::unpaired_rate), filt_quantile),
        }
    }

    /// Discard windows with extreme values.
    fn filter_windows(&self, windows: Vec<WindowCounts>) -> Vec<WindowCounts> {
        log::debug!(
            "        Filter by:  depth <= {:.0},  <= {:.1}% low mapq,  <= {:.1}% clipped,  <= {:.1}% unpaired reads",
            self.depth, 100.0 * self.low_mapq_rate, 100.0 * self.clipped_rate, 100.0 * self.unpaired_rate);
        windows.into_iter()
            .filter(|window| window.satisfies(self))
            .collect()
    }
}

impl JsonSer for LimitingValues {
    fn save(&self) -> json::JsonValue {
        json::object!{
            depth: self.depth,
            low_mapq: self.low_mapq_rate,
            clipped: self.clipped_rate,
            unpaired: self.unpaired_rate,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        let depth = obj["depth"].as_f64().ok_or_else(|| LoadError(format!(
            "ReadDepth: Failed to parse '{}': missing or incorrect 'depth' field!", obj)))?;
        let low_mapq_rate = obj["low_mapq"].as_f64().ok_or_else(|| LoadError(format!(
            "ReadDepth: Failed to parse '{}': missing or incorrect 'depth' field!", obj)))?;
        let clipped_rate = obj["clipped"].as_f64().ok_or_else(|| LoadError(format!(
            "ReadDepth: Failed to parse '{}': missing or incorrect 'depth' field!", obj)))?;
        let unpaired_rate = obj["unpaired"].as_f64().ok_or_else(|| LoadError(format!(
            "ReadDepth: Failed to parse '{}': missing or incorrect 'depth' field!", obj)))?;
        Ok(Self { depth, low_mapq_rate, clipped_rate, unpaired_rate })
    }
}

/// Return `f64` GC-content for each window in `WindowCounts`.
fn get_window_gc_contents(interval: &Interval, ref_seq: &[u8], windows: &[WindowCounts], padd: u32) -> Vec<f64> {
    let shift = interval.start();
    windows.iter().map(|window|
            gc_content(&ref_seq[(window.start - padd - shift) as usize..(window.end + padd - shift) as usize]))
        .collect()
}

/// In total, GC-content falls into 101 bins (0..=100).
const GC_BINS: usize = 101;

/// For each GC-content in 0..=100, returns start & end indices,
/// where gc_contents[start..end] lie in [GC - 0.5, GC + 0.5).
/// gc_contents must be sorted in advance.
fn find_gc_bins(gc_contents: &[f64]) -> Vec<(usize, usize)> {
    let n = gc_contents.len();
    let mut res = Vec::with_capacity(GC_BINS);
    let mut i = 0;
    for gc in 0..GC_BINS {
        let gc = gc as f64;
        let j = gc_contents.binary_search_right_at(&(gc + 0.5), i, n);
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
/// Mean values: use LOESS with all observations, without weights, and using `mean_loess_frac`.
///
/// Variance: Useuse LOESS, where observations are all GC-content values with at least 10 observations.
/// Y values are read depth variance in each particular GC-content bin.
/// Use loess_frac = 1.
fn predict_mean_var(gc_contents: &[f64], gc_bins: &[(usize, usize)], depth: &[u32], mean_loess_frac: f64)
        -> (Vec<f64>, Vec<f64>)
{
    const VAR_MIN_WINDOWS: usize = 10;
    let n = depth.len();
    let m = gc_bins.len();
    let depth_f: Vec<f64> = depth.iter().cloned().map(Into::into).collect();
    let xout: Vec<f64> = (0..m as u32).map(Into::into).collect();
    let means = Loess::new(mean_loess_frac, 1).set_xout(xout.clone()).calculate(&gc_contents, &depth_f);

    let mut x = Vec::with_capacity(m);
    let mut y = Vec::with_capacity(m);
    let mut w = Vec::with_capacity(m);

    for (gc, &(i, j)) in gc_bins.iter().enumerate() {
        if j - i >= VAR_MIN_WINDOWS {
            x.push(gc as f64);
            y.push(depth_f[i..j].variance(None));
            w.push((j - i) as f64 / n as f64);
        }
    }
    let vars = Loess::new(1.0, 1).set_xout(xout).calculate_weighted(&x, &y, &w);
    (means, vars)
}

/// At very small and very large GC-content values, there may be not enough observations
/// to well estimate read depth mean and value.
///
/// To counter this, find GC-content values such that to the left of it there are < min_obs observations.
/// Fill unavailable values with the last available mean, and with the double of the last available variance.
fn blur_boundary_values(means: &mut [f64], vars: &mut [f64], gc_bins: &[(usize, usize)], params: &ReadDepthParams) {
    let min_obs = params.min_tail_obs;
    let n = gc_bins.len();
    let m = gc_bins[gc_bins.len() - 1].1;

    // Index of the first GC-content value, where mean and variance are available.
    let left_ix = gc_bins
        .iter().enumerate()
        .filter(|(_, (_, end))| *end >= min_obs).next()
        .map(|t| t.0)
        .unwrap_or(0);
    // Index of the first GC-content value, where mean and variance are available.
    let right_ix = gc_bins
        .iter().enumerate().rev()
        .filter(|(_, (start, _))| m - start >= min_obs).next()
        .map(|t| t.0)
        .unwrap_or(n);
    assert!(left_ix < right_ix, "Too few windows to calculate read depth!");
    log::debug!("        Few windows (< {}) with GC-content < {} OR > {}, bluring distribution",
        min_obs, left_ix, right_ix);

    for i in 0..left_ix {
        means[i] = means[left_ix];
        vars[i] = (vars[left_ix] * (1.0 + (left_ix - i) as f64 * params.tail_var_mult)).max(vars[i]);
    }
    for i in right_ix + 1..n {
        means[i] = means[right_ix];
        vars[i] = (vars[right_ix].max(vars[i]) * (1.0 + (i - right_ix) as f64 * params.tail_var_mult)).max(vars[i]);
    }
}

/// Read depth parameters.
pub struct ReadDepthParams {
    /// Calculate background per windows of this size.
    /// Default: 100.
    pub window_size: u32,

    /// Calculate GC-content based on the window of size `window_size + 2 * gc_padding`.
    /// Default: 100.
    pub gc_padding: u32,

    /// Ignore left-most and right-most `edge_padding` base-pairs. Must not be smaller than `gc_padding`.
    /// Default: 1000.
    pub edge_padding: u32,

    /// Filter windows with too many low-MAPQ reads, clipped reads or unpaired reads.
    /// "too many" is defined on the `filter_quantile` across all analyzed windows. Default: 0.99.
    pub filter_quantile: f64,

    /// When calculating read depth averages for various GC-content values, use `mean_loess_frac` fraction
    /// across all windows. For example, if mean_loess_frac is 0.1 (default),
    /// use 10% of all observations with the most similar GC-content values.
    pub mean_loess_frac: f64,

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

    /// See `min_tail_obs`. Default: 0.05.
    pub tail_var_mult: f64,
}

impl Default for ReadDepthParams {
    fn default() -> Self {
        Self {
            window_size: 100,
            gc_padding: 100,
            edge_padding: 1000,
            filter_quantile: 0.99,
            mean_loess_frac: 0.1,
            min_tail_obs: 100,
            tail_var_mult: 0.05,
        }
    }
}

pub struct ReadDepth {
    window_size: u32,
    limits: LimitingValues,
    // Read depth distribution for each GC-content in 0..=100.
    distributions: Vec<CachedDistr<NBinom>>,
}

impl ReadDepth {
    pub fn estimate<'a>(
            interval: &Interval,
            ref_seq: &[u8],
            records: impl Iterator<Item = &'a Record>,
            params: &ReadDepthParams,
            max_insert_size: u32,
        ) -> Self
    {
        log::info!("    Estimating read depth");
        assert_eq!(interval.len() as usize, ref_seq.len(),
            "ReadDepth: interval and reference sequence have different lengths!");
        assert!(params.edge_padding >= params.gc_padding, "Edge padding must not be smaller than GC padding!");
        let windows = count_reads(interval, records, params, max_insert_size);
        log::debug!("        Count reads in  {:7} windows", windows.len());
        let limits = LimitingValues::new(&windows, params.filter_quantile);
        let filt_windows = limits.filter_windows(windows);
        log::debug!("        After filtering {:7} windows", filt_windows.len());
        assert!(filt_windows.len() > 0, "ReadDepth: no applicable windows!");

        let gc_contents = get_window_gc_contents(interval, ref_seq, &filt_windows, params.gc_padding);
        let ixs = gc_contents.argsort();
        let gc_contents = gc_contents.reorder(&ixs);
        let gc_bins = find_gc_bins(&gc_contents);
        let depth1: Vec<u32> = ixs.iter().map(|&i| filt_windows[i].depth1).collect();

        let (mut means, mut variances) = predict_mean_var(&gc_contents, &gc_bins, &depth1, params.mean_loess_frac);
        log::info!("        Read depth:  mean = {:.2},  st.dev. = {:.2}   (GC-content 40, {} bp windows)",
            means[40], variances[40].sqrt(), params.window_size);
        blur_boundary_values(&mut means, &mut variances, &gc_bins, params);

        /// Cache up to 256 values for each GC-content.
        const CACHE_SIZE: usize = 256;
        let distributions: Vec<_> = means.into_iter().zip(variances)
            .map(|(m, v)| NBinom::estimate(m, v).cached(CACHE_SIZE))
            .collect();
        Self {
            window_size: params.window_size,
            limits,
            distributions,
        }
    }
}

impl JsonSer for ReadDepth {
    fn save(&self) -> json::JsonValue {
        let n_params: Vec<f64> = self.distributions.iter().map(|distr| distr.distr().n()).collect();
        let p_params: Vec<f64> = self.distributions.iter().map(|distr| distr.distr().p()).collect();
        json::object!{
            window: self.window_size,
            limits: self.limits.save(),
            n: &n_params as &[f64],
            p: &p_params as &[f64],
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        let window_size = obj["window"].as_usize().ok_or_else(|| LoadError(format!(
            "ReadDepth: Failed to parse '{}': missing or incorrect 'window' field!", obj)))?;
        if !obj.has_key("limits") {
            return Err(LoadError(format!("BgDistr: Failed to parse '{}': missing 'limits' key!", obj)))
        }
        let limits = LimitingValues::load(&obj["limits"])?;
        let mut n_params = vec![0.0; window_size + 1];
        parse_f64_arr(obj, "n", &mut n_params)?;
        let mut p_params = vec![0.0; window_size + 1];
        parse_f64_arr(obj, "p", &mut p_params)?;
        Ok(Self {
            window_size: window_size as u32,
            limits,
            distributions: n_params.into_iter().zip(p_params).map(|(n, p)| NBinom::new(n, p).cached_q999()).collect(),
        })
    }
}
