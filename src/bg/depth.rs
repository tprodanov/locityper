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
            if (!record.is_proper_pair() || record.is_mate_unmapped() || record.tid() != record.mtid())
                    && record.insert_size().abs() as u32 > max_insert_size {
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
    fn satisfies(&self, max_depth: f64, max_low_mapq: f64, max_clipped: f64, max_unpaired: f64) -> bool {
        let dp = self.total_depth() as f64;
        dp > max_depth
            || self.low_mapq_reads as f64 / dp > max_low_mapq
            || self.clipped_reads as f64 / dp > max_clipped
            || self.unpaired_reads as f64 / dp > max_unpaired
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
            if middle < shift || middle >= end {
                let ix = (middle - shift) / params.window_size;
                windows[ix as usize].add_read(record, max_insert_size);
            }
        }
    }
    windows
}

/// Filters windows by removing windows with extreme values
/// (extreme read depth, number of low-mapq reads, clipped reads or unpaired reads).
fn filter_windows(windows: Vec<WindowCounts>, filt_quantile: f64) -> Vec<WindowCounts> {
    const READ_DP_QUANTILE: f64 = 0.999;
    const READ_DP_MULT: f64 = 2.0;

    // Filter windows with surprisingly high read depth.
    let max_depth = READ_DP_MULT * IterExt::quantile(
        windows.iter().map(|window| window.total_depth() as f64), READ_DP_QUANTILE);
    // Filter windows with too many low-mapq reads, reads with soft clipping or unpaired reads.
    let max_low_mapq = IterExt::quantile(windows.iter().map(WindowCounts::low_mapq_rate), filt_quantile);
    let max_clipped = IterExt::quantile(windows.iter().map(WindowCounts::clipped_rate), filt_quantile);
    let max_unpaired = IterExt::quantile(windows.iter().map(WindowCounts::unpaired_rate), filt_quantile);

    windows.into_iter()
        .filter(|window| window.satisfies(max_depth, max_low_mapq, max_clipped, max_unpaired))
        .collect()
}

/// Return `f64` GC-content for each window in `WindowCounts`.
fn get_window_gc_contents(interval: &Interval, ref_seq: &[u8], windows: &[WindowCounts], padd: u32) -> Vec<f64> {
    let shift = interval.start();
    windows.iter().map(|window|
            gc_content(&ref_seq[(window.start - padd - shift) as usize..(window.end + padd - shift) as usize]))
        .collect()
}

/// In total, GC-content falls in 101 bins (0..=100).
const GC_BINS: usize = 101;

/// For each GC-content in 0..=100, returns start & end indices,
/// where gc_contents[start..end] lie in [GC - 0.5, GC + 0.5).
/// gc_contents must be sorted in advance.
fn find_gc_bins(gc_contents: &[f64]) -> Vec<(usize, usize)> {
    let n = gc_contents.len();
    let mut res = Vec::with_capacity(GC_BINS);
    let mut i = 0;
    for gc in 0..=100 {
        let j = gc_contents.binary_search_right_at(&(gc as f64 + 0.5), i, n);
        res.push((i, j));
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
    if let Some((ix, (_, _))) = gc_bins.iter().enumerate().filter(|(_, (_, end))| *end >= min_obs).next() {
        // ix -- index of the first GC-content value, where mean and variance are available.
        for i in 0..ix {
            means[i] = means[ix];
            vars[i] = params.tail_var_mult * vars[ix];
        }
    }

    let n = gc_bins[gc_bins.len() - 1].1;
    if let Some((ix, (_, _))) = gc_bins.iter().enumerate().rev()
            .filter(|(_, (start, _))| n - start >= min_obs).next() {
        // ix -- index of the last GC-content value, where mean and variance are available.
        for i in ix+1..n {
            means[i] = means[ix];
            vars[i] = params.tail_var_mult * vars[ix];
        }
    }
}

/// Read depth parameters.
pub struct ReadDepthParams {
    /// Calculate background per windows of this size.
    pub window_size: u32,

    /// Calculate GC-content based on the window of size `window_size + 2 * gc_padding`.
    pub gc_padding: u32,

    /// Ignore left-most and right-most `edge_padding` base-pairs. Must not be smaller than `gc_padding`.
    pub edge_padding: u32,

    /// Filter windows with too many low-MAPQ reads, clipped reads or unpaired reads.
    /// "too many" is defined on the `filter_quantile` across all analyzed windows.
    pub filter_quantile: f64,

    /// When calculating read depth averages for various GC-content values, use `mean_loess_frac` fraction
    /// across all windows. For example, if mean_loess_frac is 0.1, use 10% of all observations with the most
    /// similar GC-content values.
    pub mean_loess_frac: f64,

    /// Calculate read depth parameters differently for the highest and lowest GC-content values.
    /// First, select such GC-content values L & R such that there `< min_tail_obs` windows with GC-content `< L`.
    /// Next, for all GC-content < L use mean as for L, and multiple variance by `tail_var_mult`.
    /// Repeat the same procedure for GC-content values > R.
    pub min_tail_obs: usize,

    /// See `min_tail_obs`.
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
            min_tail_obs: 500,
            tail_var_mult: 2.0,
        }
    }
}

pub struct ReadDepth {
    window_size: u32,
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
        assert_eq!(interval.len() as usize, ref_seq.len(),
            "ReadDepth: interval and reference sequence have different lengths!");
        assert!(params.edge_padding >= params.gc_padding, "Edge padding must not be smaller than GC padding!");
        let windows = count_reads(interval, records, params, max_insert_size);
        let filt_windows = filter_windows(windows, params.filter_quantile);
        assert!(filt_windows.len() > 0, "ReadDepth: no applicable windows!");

        let gc_contents = get_window_gc_contents(interval, ref_seq, &filt_windows, params.gc_padding);
        let ixs = gc_contents.argsort();
        let gc_contents = gc_contents.reorder(&ixs);
        let gc_bins = find_gc_bins(&gc_contents);
        let depth1: Vec<u32> = ixs.iter().map(|&i| filt_windows[i].depth1).collect();

        let (mut means, mut variances) = predict_mean_var(&gc_contents, &gc_bins, &depth1, params.mean_loess_frac);
        blur_boundary_values(&mut means, &mut variances, &gc_bins, params);

        /// Cache up to 256 values for each GC-content.
        const CACHE_SIZE: usize = 256;
        let distributions: Vec<_> = means.into_iter().zip(variances)
            .map(|(m, v)| NBinom::estimate(m, v).cached(CACHE_SIZE))
            .collect();
        Self {
            window_size: params.window_size,
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
            n: &n_params as &[f64],
            p: &p_params as &[f64],
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        let window_size = obj["window"].as_usize().ok_or_else(|| LoadError(format!(
            "ReadDepth: Failed to parse '{}': missing or incorrect 'window' field!", obj)))?;
        let mut n_params = vec![0.0; window_size + 1];
        parse_f64_arr(obj, "n", &mut n_params)?;
        let mut p_params = vec![0.0; window_size + 1];
        parse_f64_arr(obj, "p", &mut p_params)?;
        Ok(Self {
            window_size: window_size as u32,
            distributions: n_params.into_iter().zip(p_params).map(|(n, p)| NBinom::new(n, p).cached_q999()).collect(),
        })
    }
}
