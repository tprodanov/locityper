use std::cmp::max;
use htslib::bam::{
    Record,
    ext::BamRecordExtensions,
};
use crate::{
    seq::{
        cigar::{Operation, CigarItem},
        interv::Interval,
        seq::gc_count,
    },
    algo::{
        nbinom::{NBinom, CachedDistr},
        vec_ext::*,
        loess::Loess,
    },
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
        const MIN_MAPQ: u8 = 10;

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

/// Count reads in various windows of length `params.window_size` between ~ `interval.start() + params.padding`
/// and ~ `interval.end() - params.padding`.
fn count_reads<'a>(
        interval: &Interval,
        records: impl Iterator<Item = &'a Record>,
        params: &ReadDepthParams,
    ) -> Vec<WindowCounts>
{
    assert!(interval.len() >= params.window_size + 2 * params.padding);
    let n_windows = ((interval.len() - 2 * params.padding) as f64 / params.window_size as f64).floor() as u32;
    let sum_len = n_windows * params.window_size;
    let shift = interval.start() + (interval.len() - sum_len) / 2;
    let end = shift + sum_len;
    let mut counts: Vec<_> = (0..n_windows)
        .map(|i| WindowCounts::new(shift + i * params.window_size, shift + (i + 1) * params.window_size))
        .collect();

    for record in records {
        if (record.flags() & 3844) == 0 {
            let middle = (record.reference_start() + record.reference_end()) as u32 >> 1;
            if middle < shift || middle >= end {
                let ix = (middle - shift) / params.window_size;
                counts[ix as usize].add_read(record, params.max_insert_size);
            }
        }
    }
    counts
}

/// Filters windows by removing windows with extreme values
/// (extreme read depth, number of low-mapq reads, clipped reads or unpaired reads).
fn filter_windows(counts: Vec<WindowCounts>, filt_quantile: f64) -> Vec<WindowCounts> {
    const READ_DP_QUANTILE: f64 = 0.999;
    const READ_DP_MULT: f64 = 2.0;

    // Filter windows with surprisingly high read depth.
    let max_depth = READ_DP_MULT * IterExt::quantile(
        counts.iter().map(|counts| counts.total_depth() as f64), READ_DP_QUANTILE);
    // Filter windows with too many low-mapq reads, reads with soft clipping or unpaired reads.
    let max_low_mapq = IterExt::quantile(counts.iter().map(WindowCounts::low_mapq_rate), filt_quantile);
    let max_clipped = IterExt::quantile(counts.iter().map(WindowCounts::clipped_rate), filt_quantile);
    let max_unpaired = IterExt::quantile(counts.iter().map(WindowCounts::unpaired_rate), filt_quantile);

    counts.into_iter()
        .filter(|counts| counts.satisfies(max_depth, max_low_mapq, max_clipped, max_unpaired))
        .collect()
}

/// Predicts read depth variance for each GC-count from 0 to window_size.
/// Uses LOESS, with input points taken for all GC-count values with at least ten observations.
fn predict_variance(depth: &[u32], gc_counts: &[u32], window_size: u32) -> Vec<f64> {
    const MIN_WINDOWS: usize = 10;
    let n = depth.len();
    let ixs = gc_counts.argsort();
    let gc_counts = gc_counts.reorder(&ixs);
    let depth_f64: Vec<f64> = (0..n).map(|i| depth[ixs[i]] as f64).collect();

    let m = window_size as usize + 1;
    let mut x = Vec::with_capacity(m);
    let mut y = Vec::with_capacity(m);
    let mut w = Vec::with_capacity(m);
    let mut i = 0;
    for gc in 0..=window_size {
        let j = gc_counts.binary_search_right_at(&gc, i, n);
        if j - i >= MIN_WINDOWS {
            x.push(gc as f64);
            y.push(x[i..j].variance(None));
            w.push((j - i) as f64 / n as f64);
        }
        i = j;
    }
    debug_assert_eq!(i, n);

    let xout: Vec<f64> = (0..=window_size).map(Into::into).collect();
    Loess::new().set_frac(1.0).set_xout(xout.clone()).calculate_weighted(&x, &y, &w)
}

/// Read depth parameters.
pub struct ReadDepthParams {
    pub window_size: u32,
    pub padding: u32,
    pub max_insert_size: u32,
    pub filter_quantile: f64,
}

pub struct ReadDepth {
    window_size: u32,
    // Read depth distribution for each GC-count in 0..=window_size.
    distributions: Vec<CachedDistr<NBinom>>,
}

impl ReadDepth {
    pub fn estimate<'a>(
            interval: &Interval,
            ref_seq: &[u8],
            records: impl Iterator<Item = &'a Record>,
            params: &ReadDepthParams,
        ) -> Self
    {
        assert_eq!(interval.len() as usize, ref_seq.len(),
            "ReadDepth: interval and reference sequence have different lengths!");
        let window_counts = count_reads(interval, records, params);
        let filt_window_counts = filter_windows(window_counts, params.filter_quantile);

        let shift = interval.start();
        let gc_counts: Vec<u32> = filt_window_counts.iter()
            .map(|window| gc_count(&ref_seq[(window.start - shift) as usize .. (window.end - shift) as usize]))
            .collect();
        let depth1: Vec<u32> = filt_window_counts.iter().map(|counts| counts.depth1).collect();

        unimplemented!()
    }
}
