use std::{
    cmp::{Ordering, min},
    io::{self, Write},
};
use rand::Rng;
use crate::{
    ext::vec::F64Ext,
    bg::depth::ReadDepth,
    seq::{
        self, ContigId, ContigNames, ContigSet,
        kmers::KmerCounts,
    },
};
use super::{
    locs::{ReadPairAlignments, PairAlignment},
    dp_cache::{CachedDepthDistrs, DistrBox},
};

/// Structure that, for predefined region boundaries and window size,
/// returns windows for each alignment boundaries.
/// This function is used both in read depth calculation and during the locus genotyping.
#[derive(Clone)]
pub struct WindowGetter {
    start: u32,
    end: u32,
    window: u32,
    // Half window: size that needs to be covered in `get_significant` - 1.
    halfw: u32,
}

impl WindowGetter {
    pub fn new(start: u32, end: u32, window: u32) -> Self {
        Self {
            start, end, window,
            halfw: window.checked_sub(1).unwrap() / 2,
        }
    }

    /// Returns range of windows (start_ix..end_ix),
    /// covered by the alignment (at least by 1 bp).
    pub fn covered_any(&self, aln_start: u32, aln_end: u32) -> (u32, u32) {
        (
            min(aln_start, self.end).saturating_sub(self.start) / self.window,
            (min(aln_end, self.end) + self.window - 1).saturating_sub(self.start) / self.window,
        )
    }

    /// Returns range of windows (start_ix..end_ix),
    /// middle of which are covered by the alignment.
    pub fn covered_middle(&self, aln_start: u32, aln_end: u32) -> (u32, u32) {
        (
            (min(aln_start, self.end) + self.halfw).saturating_sub(self.start) / self.window,
            (min(aln_end, self.end) + self.halfw).saturating_sub(self.start) / self.window,
        )
    }
}

/// First window represents unmapped reads (or windows on the boundary).
pub(crate) const UNMAPPED_WINDOW: u32 = 0;
/// Second window represents reads alignments within the boundary.
pub(crate) const BOUNDARY_WINDOW: u32 = 1;
/// Regular windows have indices starting with 2.
pub(crate) const REG_WINDOW_SHIFT: u32 = 2;

/// Alignment of a read pair to a specific multi-contig windows.
pub struct ReadWindows {
    /// Index in the list of read-pair alignments.
    aln_ix: u32,
    /// ln-probability of this location.
    ln_prob: f64,
    /// Index of the contig (within multi-contig windows).
    contig_ix: Option<u8>,
    /// Start and end of the first read end alignment.
    range1: Option<(u32, u32)>,
    /// Start and end of the second read end alignment.
    range2: Option<(u32, u32)>,
    /// Windows, affected by the read.
    windows: AffectedWindows,
}

impl ReadWindows {
    fn new(aln_ix: u32, contig_ix: u8, paln: &PairAlignment) -> Self {
        Self {
            aln_ix,
            contig_ix: Some(contig_ix),
            range1: paln.intervals().range1(),
            range2: paln.intervals().range2(),
            ln_prob: paln.ln_prob(),
            windows: AffectedWindows::new(),
        }
    }

    fn both_unmapped(ln_prob: f64) -> Self {
        Self {
            ln_prob,
            aln_ix: u32::MAX,
            contig_ix: None,
            range1: None,
            range2: None,
            windows: AffectedWindows::new(),
        }
    }

    pub fn define_windows_determ(&mut self, mcontigs: &MultiContigWindows) {
        if let Some(i) = self.contig_ix {
            let contig = &mcontigs.by_contig[i as usize];
            let shift = mcontigs.wshifts[i as usize];
            let (w1s, w1e) = match self.range1 {
                Some((start, end)) => contig.get_windows(start, end, shift),
                None => (UNMAPPED_WINDOW, UNMAPPED_WINDOW + 1),
            };
            let (w2s, w2e) = match self.range2 {
                Some((start, end)) => contig.get_windows(start, end, shift),
                None => (UNMAPPED_WINDOW, UNMAPPED_WINDOW + 1),
            };
            self.windows.init_from_ranges(w1s, w1e, w2s, w2e);
        }
    }

    pub fn define_windows_random(&mut self, mcontigs: &MultiContigWindows, tweak: u32, rng: &mut impl Rng) {
        if let Some(i) = self.contig_ix {
            let contig = &mcontigs.by_contig[i as usize];
            let shift = mcontigs.wshifts[i as usize];
            let r = rng.next_u64();
            let tweak1 = (r >> 32) as u32 % (2 * tweak + 1);
            let tweak2 = r as u32 % (2 * tweak + 1);

            let (w1s, w1e) = match self.range1 {
                Some((start, end)) => contig.get_windows(
                    (start + tweak1).saturating_sub(tweak), (end + tweak1).saturating_sub(tweak), shift),
                None => (UNMAPPED_WINDOW, UNMAPPED_WINDOW + 1),
            };
            let (w2s, w2e) = match self.range2 {
                Some((start, end)) => contig.get_windows(
                    (start + tweak2).saturating_sub(tweak), (end + tweak2).saturating_sub(tweak), shift),
                None => (UNMAPPED_WINDOW, UNMAPPED_WINDOW + 1),
            };
            self.windows.init_from_ranges(w1s, w1e, w2s, w2e);
        }
    }

    /// Returns range of windows, to which the first and the second read ends are aligned.
    pub fn windows(&self) -> &AffectedWindows {
        &self.windows
    }

    /// Returns ln-probability of the alignment.
    pub fn ln_prob(&self) -> f64 {
        self.ln_prob
    }

    /// Index of the read-pair alignment across all alignments for the read pair.
    pub(crate) fn aln_ix(&self) -> u32 {
        self.aln_ix
    }
}

/// Calculates weight of a positive value x.
/// If x <= c1, returns 1.
/// If x = c2, returns 0.5 (returns ln 0.5).
/// For values between c1 and infinity, the weight is distributed according to a sech function.
/// ```
/// f(x) = { 1                              if x in [0, c1].
///         /      x - c1
///        { sech ------- * ln(2 + sqrt3)   if x > c1.
///         \     c2 - c1
/// ```
fn sech_weight(x: f64, c1: f64, c2: f64) -> f64 {
    if x <= c1 {
        return 1.0;
    }
    // ln(2 + sqrt(3)).
    const LN2_SQRT3: f64 = 1.31695789692481;
    let t = (x - c1) / (c2 - c1) * LN2_SQRT3;
    let expt = t.exp();
    2.0 * expt / (expt * expt + 1.0)
}

/// Set of windows for one contig.
#[derive(Clone)]
pub struct ContigWindows {
    contig_id: ContigId,
    /// Structure, that stores start, end and window.
    window_getter: WindowGetter,
    /// GC-content for each window within contig.
    window_gcs: Vec<u8>,
    /// k-mer-based window weights.
    window_weights: Vec<f64>,
}

impl ContigWindows {
    pub fn new(
        contig_id: ContigId,
        contigs: &ContigNames,
        seq: &[u8],
        kmer_counts: &KmerCounts,
        depth: &ReadDepth,
        params: &super::Params,
    ) -> Self
    {
        let contig_len = seq.len() as u32;
        let window = depth.window_size();
        let neighb_size = depth.neighb_size();
        assert!(contig_len > window + 2 * params.boundary_size,
            "Contig {} is too short (len = {})", contigs.get_name(contig_id), contig_len);
        debug_assert_eq!(contig_len, contigs.get_len(contig_id));
        let n_windows = (contig_len - 2 * params.boundary_size) / window;
        let sum_len = n_windows * window;
        let start = (contig_len - sum_len) / 2;
        let end = start + sum_len;

        let k = kmer_counts.k();
        let halfk = k / 2;
        let contig_kmer_counts = kmer_counts.get(contig_id);

        let mut window_gcs = Vec::with_capacity(n_windows as usize);
        let mut window_weights = Vec::with_capacity(n_windows as usize);
        let window_padding = (neighb_size - window) / 2;
        for i in 0..n_windows {
            let padded_start = start + i * window - window_padding;
            let padded_end = padded_start + neighb_size;
            window_gcs.push(seq::gc_content(&seq[padded_start as usize..padded_end as usize]).round() as u8);

            let mean_kmer_freq = F64Ext::mean(&contig_kmer_counts[padded_start.saturating_sub(halfk) as usize
                ..min(padded_end - halfk, contig_len - k + 1) as usize]);
            window_weights.push(sech_weight(mean_kmer_freq, params.rare_kmer, params.semicommon_kmer));
        }
        Self {
            window_getter: WindowGetter::new(start, end, window),
            contig_id, window_gcs, window_weights
        }
    }

    /// Creates a set of contig windows for each contig.
    pub fn new_all(set: &ContigSet, depth: &ReadDepth, params: &super::Params) -> Vec<Self> {
        let contigs = set.contigs();
        contigs.ids().zip(set.seqs())
            .map(|(id, seq)| Self::new(id, contigs, seq, set.kmer_counts(), depth, params))
            .collect()
    }

    pub fn n_windows(&self) -> u32 {
        self.window_gcs.len() as u32
    }

    pub fn window_size(&self) -> u32 {
        self.window_getter.window
    }

    /// Returns window range within the contig based on the read alignment range.
    fn get_windows(&self, aln_start: u32, aln_end: u32, shift: u32) -> (u32, u32) {
        let (start_ix, end_ix) = self.window_getter.covered_middle(aln_start, aln_end);
        if start_ix < end_ix {
            (start_ix + shift, end_ix + shift)
        } else {
            (BOUNDARY_WINDOW, BOUNDARY_WINDOW + 1)
        }
    }

    pub const BED_HEADER: &'static str = "contig\tstart\tend\tgc\tweight";

    /// Writes windows for this contig in a BED format (see `BED_HEADER`).
    pub fn write_to(&self, f: &mut impl Write, contigs: &ContigNames) -> io::Result<()> {
        let name = contigs.get_name(self.contig_id);
        for (i, (&gc, &weight)) in self.window_gcs.iter().zip(&self.window_weights).enumerate() {
            let start = self.window_getter.start + i as u32 * self.window_getter.window;
            writeln!(f, "{}\t{}\t{}\t{}\t{:.5}", name, start, start + self.window_getter.window, gc, weight)?;
        }
        Ok(())
    }
}

/// Stores the contigs and windows corresponding to the windows.
pub struct MultiContigWindows {
    by_contig: Vec<ContigWindows>,
    /// `ln(sum(cns))`.
    ln_ploidy: f64,
    /// Start index for each contig id (length: n-contigs + 1).
    /// Contig with index `i` will have windows between `wshifts[i]..wshifts[i + 1]`.
    wshifts: Vec<u32>,
    window: u32,
}

impl MultiContigWindows {
    pub fn new(contig_ids: &[ContigId], contig_windows: &[ContigWindows]) -> Self {
        let n = contig_ids.len();
        let mut by_contig = Vec::<ContigWindows>::with_capacity(n);
        let mut wshifts = Vec::with_capacity(n + 1);
        let mut curr_wshift = REG_WINDOW_SHIFT;
        wshifts.push(curr_wshift);

        for &id in contig_ids.iter() {
            let curr_contig = contig_windows[id.ix()].clone();
            curr_wshift += curr_contig.n_windows();
            wshifts.push(curr_wshift);
            by_contig.push(curr_contig);
        }
        assert!(by_contig.len() < 256, "Multi-contig collection cannot contain more than 256 entries");
        Self {
            ln_ploidy: (n as f64).ln(),
            window: by_contig[0].window_size(),
            by_contig, wshifts,
        }
    }

    /// Total number of contigs. Can be less than ploidy, as contigs can be repeated.
    pub fn n_contigs(&self) -> usize {
        self.by_contig.len()
    }

    /// Total number of windows in the contig group.
    pub fn total_windows(&self) -> u32 {
        *self.wshifts.last().unwrap()
    }

    pub fn window_size(&self) -> u32 {
        self.window
    }

    pub fn ids(&self) -> impl Iterator<Item = ContigId> + '_ {
        self.by_contig.iter().map(|el| el.contig_id)
    }

    /// Returns string with all contig names through a comma.
    pub fn ids_str(&self, contigs: &ContigNames) -> String {
        contigs.get_names(self.ids())
    }

    /// Returns window shift for the `i`-th contig.
    pub(crate) fn get_wshift(&self, i: usize) -> u32 {
        self.wshifts[i]
    }

    /// Given all pair-alignments for a single read pair,
    /// finds all alignments corresponding to `self`, and appends them to `out_alns`.
    /// Returns the number of added alignments.
    ///
    /// Output pair-alignments with ln-probability worse than `best_prob - prob_diff` are discarded.
    /// Remaining alignments have random order, probabilities are not normalized.
    ///
    /// Any alignments that were in the vector before, stay as they are and in the same order.
    pub fn read_windows(&self,
        pair_alns: &ReadPairAlignments,
        out_alns: &mut Vec<ReadWindows>,
        prob_diff: f64,
    ) -> usize {
        let start_len = out_alns.len();
        // Probability of being unmapped to any of the contigs.
        let unmapped_prob = pair_alns.unmapped_prob() + self.ln_ploidy;
        // Current threshold, is updated during the for-loop.
        let mut thresh_prob = unmapped_prob - prob_diff;
        for (i, contig_id) in self.ids().enumerate() {
            let contig_ix = u8::try_from(i).unwrap();
            let (start_ix, palns) = pair_alns.contig_alns(contig_id);
            for (aln_ix, paln) in (start_ix..).zip(palns) {
                let ln_prob = paln.ln_prob();
                if ln_prob >= thresh_prob {
                    thresh_prob = thresh_prob.max(ln_prob - prob_diff);
                    out_alns.push(ReadWindows::new(aln_ix as u32, contig_ix, paln));
                }
            }
        }

        // As threshold was updated during the for-loop, some alignments in the beginning may need to removed.
        let mut i = start_len;
        while i < out_alns.len() {
            if out_alns[i].ln_prob < thresh_prob {
                out_alns.swap_remove(i);
            } else {
                i += 1;
            }
        }
        if unmapped_prob >= thresh_prob {
            out_alns.push(ReadWindows::both_unmapped(unmapped_prob));
        }
        out_alns.len() - start_len
    }

    pub(crate) fn get_distributions(&self, cached_distrs: &CachedDepthDistrs) -> Vec<DistrBox> {
        const _: () = assert!(UNMAPPED_WINDOW == 0 && BOUNDARY_WINDOW == 1 && REG_WINDOW_SHIFT == 2,
            "Constants were changed!");
        let mut distrs: Vec<DistrBox> = Vec::with_capacity(self.total_windows() as usize);
        distrs.push(Box::new(cached_distrs.unmapped_distr()));
        distrs.push(Box::new(cached_distrs.boundary_distr()));

        for curr_contig in self.by_contig.iter() {
            for (&gc_content, &weight) in curr_contig.window_gcs.iter().zip(&curr_contig.window_weights) {
                distrs.push(cached_distrs.get_distribution(gc_content, weight));
            }
        }
        distrs
    }
}

/// Windows, affected by a single read/read pair.
/// Contains pairs (window index, depth change), sorted by window index.
#[derive(Debug)]
pub struct AffectedWindows(Vec<(u32, i8)>);

impl AffectedWindows {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// Creates list of affected windows from two ranges `w1s..w1e` and `w2s..w2e`.
    fn init_from_ranges(&mut self, w1s: u32, w1e: u32, w2s: u32, w2e: u32) {
        let (s1, s2) = if w1s <= w2s { (w1s, w2s) } else { (w2s, w1s) };
        let (e1, e2) = if w1e <= w2e { (w1e, w2e) } else { (w2e, w1e) };

        self.0.clear();
        if e1 <= s2 {
            // Ranges do not overlap.
            // s1 < e1 <= s2 < e2
            self.0.extend((s1..e1).map(|w| (w, 1)));
            self.0.extend((s2..e2).map(|w| (w, 1)));
        } else {
            // Ranges overlap.
            // s1 <= s2 < e1 <= e2
            self.0.extend((s1..s2).map(|w| (w, 1)));
            self.0.extend((s2..e1).map(|w| (w, 2)));
            self.0.extend((e1..e2).map(|w| (w, 1)));
        }
    }

    /// Calculates difference (w, a1 - a2), and fills buffer `self` with it.
    pub fn fill_diff(&mut self, windows1: &Self, windows2: &Self) {
        self.0.clear();
        let mut it1 = windows1.iter();
        let mut it2 = windows2.iter();
        const UNDEF: (u32, i8) = (u32::MAX, 0);
        let (mut w1, mut a1) = it1.next().unwrap_or(UNDEF);
        let (mut w2, mut a2) = it2.next().unwrap_or(UNDEF);
        while w1 < u32::MAX || w2 < u32::MAX {
            match w1.cmp(&w2) {
                Ordering::Less => {
                    self.0.push((w1, a1));
                    (w1, a1) = it1.next().unwrap_or(UNDEF);
                }
                Ordering::Equal => {
                    if a1 != a2 {
                        self.0.push((w1, a1 - a2));
                    }
                    (w1, a1) = it1.next().unwrap_or(UNDEF);
                    (w2, a2) = it2.next().unwrap_or(UNDEF);
                }
                Ordering::Greater => {
                    self.0.push((w2, -a2));
                    (w2, a2) = it2.next().unwrap_or(UNDEF);
                }
            }
        }
    }

    /// Returns iterator over pairs (window index, depth).
    pub fn iter(&self) -> impl Iterator<Item = (u32, i8)> + std::iter::ExactSizeIterator + '_ {
        self.0.iter().copied()
    }
}
