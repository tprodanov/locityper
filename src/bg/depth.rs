use htslib::bam::{
    Record,
    ext::BamRecordExtensions,
};
use crate::{
    seq::{
        self, Interval,
        cigar::Cigar,
        kmers::KmerCounts,
    },
    algo::{bisect, loess::Loess},
    ext::vec::{VecExt, F64Ext},
    math::distr::{NBinom, WithMoments},
    err::{Error, validate_param},
};
use super::{
    err_prof::ErrorProfile,
    ser::{JsonSer, LoadError, parse_f64_arr},
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
    depth1: u32,
    depth2: u32,
}

impl WindowCounts {
    fn new(start: u32, end: u32) -> Self {
        WindowCounts {
            start, end,
            depth1: 0,
            depth2: 0,
        }
    }

    fn len(&self) -> u32 {
        self.end - self.start
    }

    fn add_read(&mut self, record: &Record) {
        if record.is_last_in_template() {
            self.depth2 += 1;
        } else {
            self.depth1 += 1;
        }
    }
}

/// Count reads in various windows of length `params.window_size` between ~ `interval.start() + params.edge_padding`
/// and ~ `interval.end() - params.edge_padding`.
fn count_reads<'a>(
    records: impl Iterator<Item = &'a Record>,
    interval: &Interval,
    ref_seq: &[u8],
    params: &ReadDepthParams,
    err_prof: &ErrorProfile,
    max_insert_size: i64,
    min_aln_ln_prob: f64,
) -> Vec<WindowCounts> {
    assert!(interval.len() >= params.window_size + 2 * params.edge_padding, "Input interval is too short!");
    let n_windows = (f64::from(interval.len() - 2 * params.edge_padding) / f64::from(params.window_size))
        .floor() as u32;
    let sum_len = n_windows * params.window_size;
    let interval_start = interval.start();
    let start = interval_start + (interval.len() - sum_len) / 2;
    let end = start + sum_len;
    let mut windows: Vec<_> = (0..n_windows)
        .map(|i| WindowCounts::new(start + i * params.window_size, start + (i + 1) * params.window_size))
        .collect();

    for record in records {
        if super::read_unpaired_or_proper_pair(record, max_insert_size) {
            let aln_start = record.reference_start();
            let aln_end = record.reference_end();
            let middle = (aln_start + aln_end) as u32 / 2;
            if start <= middle && middle < end {
                let cigar = Cigar::infer_ext_cigar(record, ref_seq, interval_start);
                if err_prof.ln_prob(&cigar) >= min_aln_ln_prob {
                    let ix = (middle - start) / params.window_size;
                    windows[ix as usize].add_read(record);
                }
            }
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
/// Modifies input `windows` variable, and returns GC-content for each window.
fn filter_windows(
    windows: &mut Vec<WindowCounts>,
    seq_shift: u32,
    ref_seq: &[u8],
    kmer_counts_coll: &KmerCounts,
    params: &ReadDepthParams,
) -> Vec<f64> {
    log::debug!("    Total windows:   {:7}", windows.len());
    let mut have_ns = 0;
    let mut have_common_kmers = 0;

    let k = kmer_counts_coll.k();
    let kmer_counts = kmer_counts_coll.get_first();
    let window_padding = params.window_padding;
    let neigh_len = (windows[0].len() + 2 * window_padding) as usize;
    let neigh_kmers = neigh_len + 1 - k as usize;
    let thresh_sum = params.max_kmer_freq * neigh_kmers as f64;

    let mut gc_contents = Vec::with_capacity(windows.len());
    windows.retain(|window| {
        let start_ix = (window.start - window_padding - seq_shift) as usize;
        let window_seq = &ref_seq[start_ix..start_ix + neigh_len];
        if seq::has_n(window_seq) {
            have_ns += 1;
            return false;
        }
        let sum_kmer_freq = kmer_counts[start_ix..start_ix + neigh_kmers].iter().copied().map(f64::from).sum::<f64>();
        if sum_kmer_freq > thresh_sum {
            have_common_kmers += 1;
            false
        } else {
            gc_contents.push(seq::gc_content(window_seq));
            true
        }
    });
    assert_eq!(gc_contents.len(), windows.len());
    log::debug!("    Remove {} windows with Ns, {} windows with common k-mers", have_ns, have_common_kmers);
    log::debug!("    After filtering: {:7}", windows.len());
    gc_contents
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
fn predict_mean_var(gc_contents: &[f64], gc_bins: &[(usize, usize)], depth: &[u32], frac_windows: f64)
    -> (Vec<f64>, Vec<f64>)
{
    const VAR_MIN_WINDOWS: usize = 10;
    let n = depth.len();
    let m = gc_bins.len();
    let depth_f: Vec<f64> = depth.iter().copied().map(Into::into).collect();
    let xout: Vec<f64> = (0..u32::try_from(m).unwrap()).map(Into::into).collect();
    let means = Loess::new(frac_windows, 1).set_xout(xout.clone()).calculate(&gc_contents, &depth_f);

    let mut x = Vec::with_capacity(m);
    let mut y = Vec::with_capacity(m);
    let mut w = Vec::with_capacity(m);

    for (gc, &(i, j)) in gc_bins.iter().enumerate() {
        if j - i >= VAR_MIN_WINDOWS {
            x.push(gc as f64);
            y.push(F64Ext::variance(&depth_f[i..j], None));
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

    // Index of the first GC-content value, such that there are at least `min_obs` windows to the left.
    let left_ix = gc_bins.iter().position(|(_, end)| *end >= min_obs).unwrap_or(n);
    // Index of the last GC-content value, such that there are at least `min_obs` windows to the right.
    let right_ix = n - gc_bins.iter().rev().position(|(start, _)| m - *start >= min_obs).unwrap_or(n);
    assert!(left_ix < right_ix, "Too few windows to calculate read depth!");
    log::debug!("    Few windows (< {}) with GC-content < {} or > {}, bluring distributions",
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

    /// Ignore left-most and right-most `edge_padding` base-pairs. Must not be smaller than `window_padding`.
    /// Default: 1000.
    pub edge_padding: u32,

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

    /// See `min_tail_obs`. Default: 0.05.
    pub tail_var_mult: f64,
}

impl Default for ReadDepthParams {
    fn default() -> Self {
        Self {
            ploidy: 2,
            window_size: 100,
            window_padding: 100,
            edge_padding: 1000,
            max_kmer_freq: 1.2,

            frac_windows: 0.1,
            min_tail_obs: 100,
            tail_var_mult: 0.05,
        }
    }
}

impl ReadDepthParams {
    /// Validate all parameter values.
    pub fn validate(&self) -> Result<(), Error> {
        validate_param!(self.ploidy > 0, "Ploidy cannot be zero");
        validate_param!(self.edge_padding >= self.window_padding,
            "Edge padding ({}) must not be smaller than window padding ({})", self.edge_padding, self.window_padding);
        if self.edge_padding < self.window_size {
            log::warn!("Edge padding ({}) is smaller than the window size ({}), consider using a larger value",
                self.edge_padding, self.window_size);
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

/// Background read depth distributions for each GC-content value.
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
    /// Ignore reads with alignment probability < `min_aln_ln_prob` and with insert size > `max_insert_size`.
    pub fn estimate<'a>(
        records: impl Iterator<Item = &'a Record>,
        interval: &Interval,
        ref_seq: &[u8],
        kmer_counts: &KmerCounts,
        err_prof: &ErrorProfile,
        max_insert_size: i64,
        params: &ReadDepthParams,
        subsampling_rate: f64,
    ) -> Self {
        log::info!("Estimating read depth (ploidy {}, subsampling rate {})", params.ploidy, subsampling_rate);
        assert_eq!(interval.len() as usize, ref_seq.len(),
            "ReadDepth: interval and reference sequence have different lengths!");

        let min_aln_ln_prob = err_prof.min_aln_ln_prob();
        let mut windows = count_reads(records, interval, ref_seq, params, err_prof, max_insert_size, min_aln_ln_prob);
        let gc_contents = filter_windows(&mut windows, interval.start(), &ref_seq, &kmer_counts, params);
        assert!(windows.len() > 0, "ReadDepth: no applicable windows!");

        let ixs = VecExt::argsort(&gc_contents);
        let gc_contents = VecExt::reorder(&gc_contents, &ixs);
        let gc_bins = find_gc_bins(&gc_contents);
        let depth1: Vec<u32> = ixs.iter().map(|&i| windows[i].depth1).collect();

        let (mut means, mut variances) = predict_mean_var(&gc_contents, &gc_bins, &depth1, params.frac_windows);
        blur_boundary_values(&mut means, &mut variances, &gc_bins, params);

        let nbinom_mul = 1.0 / (subsampling_rate * f64::from(params.ploidy));
        let distributions: Vec<_> = means.into_iter().zip(variances)
            .map(|(m, v)| NBinom::estimate(m, v.max(m * 1.00001)).mul(nbinom_mul))
            .collect();
        const GC_VAL: usize = 40;
        log::info!("    Read depth:  mean = {:.2},  st.dev. = {:.2}   \
            (Ploidy 1, first mates, GC-content {}, {} bp windows)",
            distributions[GC_VAL].mean(), distributions[GC_VAL].variance().sqrt(), GC_VAL, params.window_size);
        Self {
            ploidy: params.ploidy,
            window_size: params.window_size,
            window_padding: params.window_padding,
            distributions,
        }
    }

    /// Returns window size.
    pub fn window_size(&self) -> u32 {
        self.window_size
    }

    pub fn window_padding(&self) -> u32 {
        self.window_padding
    }

    /// Returns read depth distribution at GC-content `gc_content` (between 0 and 100).
    pub fn depth_distribution(&self, gc_content: usize) -> &NBinom {
        &self.distributions[gc_content]
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

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        let ploidy = obj["ploidy"].as_usize().ok_or_else(|| LoadError(format!(
            "ReadDepth: Failed to parse '{}': missing or incorrect 'ploidy' field!", obj)))?;
        let window_size = obj["window"].as_usize().ok_or_else(|| LoadError(format!(
            "ReadDepth: Failed to parse '{}': missing or incorrect 'window' field!", obj)))?;
        let window_padding = obj["window_padding"].as_usize().ok_or_else(|| LoadError(format!(
            "ReadDepth: Failed to parse '{}': missing or incorrect 'window_padding' field!", obj)))?;
        let mut n_params = vec![0.0; window_size + 1];
        parse_f64_arr(obj, "n", &mut n_params)?;
        let mut p_params = vec![0.0; window_size + 1];
        parse_f64_arr(obj, "p", &mut p_params)?;
        Ok(Self {
            ploidy: u8::try_from(ploidy).unwrap(),
            window_size: u32::try_from(window_size).unwrap(),
            window_padding: u32::try_from(window_padding).unwrap(),
            distributions: n_params.into_iter().zip(p_params).map(|(n, p)| NBinom::new(n, p)).collect(),
        })
    }
}
