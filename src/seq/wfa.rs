use std::{
    ffi::c_char,
    cmp::max,
};
use crate::{
    seq::cigar::{Cigar, Operation},
    err::validate_param,
};

#[allow(non_upper_case_globals)]
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused, unnecessary_transmutes)]
mod cwfa {
    include!(concat!(env!("OUT_DIR"), "/bindings_wfa.rs"));
}

/// Maximum accuracy level (need to change `add` help message, if changed this).
pub const MAX_ACCURACY: u8 = 9;

#[derive(Clone, Copy)]
pub struct Penalties {
    // NOTE: Adding match may lead to problems in other places (such as `panvcf` divergence calculation).
    // Similarly, adding separate parameters for insert/deletion penalties may lead to problems elsewhere.
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for Penalties {
    fn default() -> Self {
        Self {
            mismatch: 4,
            gap_open: 6,
            gap_extend: 1,
        }
    }
}

impl Penalties {
    pub fn validate(&self) -> crate::Result<()> {
        validate_param!(self.mismatch >= 0, "Mismatch penalty ({}) must be non-negative", self.mismatch);
        validate_param!(self.gap_open >= 0, "Gap opening penalty ({}) must be non-negative", self.gap_open);
        validate_param!(self.gap_extend >= 0, "Gap extension penalty ({}) must be non-negative", self.gap_extend);
        Ok(())
    }

    /// Simple alignment between two sequences: insert/deletion at the start,
    /// followed by a sequence of =/X without any gaps.
    pub fn align_simple(&self, seq1: &[u8], seq2: &[u8], cigar: &mut Cigar) -> i32 {
        let n = seq1.len();
        let m = seq2.len();
        debug_assert!(n > 0 && m > 0, "Cannot align empty sequences");
        let diff = n as i32 - m as i32;
        let mut score;
        let (i, j) = if diff < 0 {
            cigar.push_unchecked(Operation::Ins, (-diff) as u32);
            score = -self.gap_open + diff * self.gap_extend;
            (0, (-diff) as usize)
        } else if diff > 0 {
            cigar.push_unchecked(Operation::Del, diff as u32);
            score = -self.gap_open - diff * self.gap_extend;
            (diff as usize, 0)
        } else {
            score = 0;
            (0, 0)
        };

        let mut curr_match = seq1[i] == seq2[j];
        let mut curr_len = 1;
        for (&nt1, &nt2) in seq1[i + 1..].iter().zip(&seq2[j + 1..]) {
            if (nt1 == nt2) != curr_match {
                cigar.push_unchecked(if curr_match { Operation::Equal } else { Operation::Diff }, curr_len);
                score -= if curr_match { 0 } else { self.mismatch * curr_len as i32 };
                curr_match = !curr_match;
                curr_len = 1;
            } else {
                curr_len += 1;
            }
        }
        cigar.push_unchecked(if curr_match { Operation::Equal } else { Operation::Diff }, curr_len);
        score -= if curr_match { 0 } else { self.mismatch * curr_len as i32 };
        score
    }
}

/// Number of alignment steps based on the accuracy level (0-9).
fn alignment_steps(accuracy: u8) -> i32 {
    match accuracy {
        0 =>       1,
        1 =>     100,
        2 =>     500,
        3 =>   1_000,
        4 =>   3_000,
        5 =>   5_000,
        6 =>  10_000,
        7 =>  50_000,
        8 => 100_000,
        9 => i32::MAX,
        _ => unreachable!(),
    }
}

/// Used for check for max jump. This way we can turn on/off max jump check during compile time.
pub trait Threshold : Copy {
    /// Returns true if `threshold < val`.
    fn under(self, val: u32) -> bool;
}

impl Threshold for u32 {
    #[inline(always)]
    fn under(self, val: u32) -> bool { self < val }
}

/// Larger than any value, always return false
impl Threshold for () {
    #[inline(always)]
    fn under(self, _val: u32) -> bool { false }
}

pub struct Aligner {
    /// WFA aligner for global alignments.
    global_aligner: *mut cwfa::wavefront_aligner_t,
    /// WFA aligner for left clipping (sequence right side should be aligned, left side can be unaligned)
    /// and right clipping.
    semiglobal_aligner: Option<*mut cwfa::wavefront_aligner_t>,
    penalties: Penalties,
    /// At what size nX could not theoretically be replaced with 1I(n-1)=1D.
    /// For example, with default penalties (4,6,1)
    /// ACGT                          ACGT
    /// |||| will have score -16, but  |||  will have score -14, so safe_mismatch_size < 4.
    /// CGTA                           CGTA
    safe_mismatch_size: u32,
}

impl Drop for Aligner {
    fn drop(&mut self) {
        unsafe { cwfa::wavefront_aligner_delete(self.global_aligner) };
        if let Some(ptr) = self.semiglobal_aligner {
            unsafe { cwfa::wavefront_aligner_delete(ptr) };
        }
    }
}

impl Aligner {
    /// Accuracy level must be in [1, 9], 0 - very fast and inaccurate, 9 - very slow and accurate.
    pub fn new(
        penalties: Penalties,
        accuracy: u8,
        bound: Option<i32>,
        enable_semiglobal: bool,
    ) -> Self {
        assert!(1 <= accuracy && accuracy <= MAX_ACCURACY, "Cannot construct WFA aligner for accuracy {}", accuracy);
        let mut attributes = unsafe { cwfa::wavefront_aligner_attr_default }.clone();
        // Limit the number of alignment steps.
        attributes.system.max_alignment_steps = alignment_steps(accuracy);

        // High memory for some reason produces alignments with M both for X and =.
        attributes.memory_mode = if accuracy < 7 {
            cwfa::wavefront_memory_t_wavefront_memory_med
        } else {
            cwfa::wavefront_memory_t_wavefront_memory_low
        };
        // Compute score and CIGAR as well.
        attributes.alignment_scope = cwfa::alignment_scope_t_compute_alignment;

        if let Some(k) = bound {
            attributes.heuristic.strategy = cwfa::wf_heuristic_strategy_wf_heuristic_banded_adaptive;
            attributes.heuristic.min_k = -k;
            attributes.heuristic.max_k = k;
            attributes.heuristic.steps_between_cutoffs = 3;
        }

        // Set the cost model and parameters.
        attributes.distance_metric = cwfa::distance_metric_t_gap_affine;
        attributes.affine_penalties.mismatch = penalties.mismatch;
        attributes.affine_penalties.gap_opening = penalties.gap_open;
        attributes.affine_penalties.gap_extension = penalties.gap_extend;

        let semiglobal_aligner = if enable_semiglobal {
            // Need positive match score for alignment to work.
            attributes.affine_penalties.match_ = -max(1, attributes.affine_penalties.mismatch / 2);
            attributes.alignment_form.span = cwfa::alignment_span_t_alignment_endsfree;
            Some(unsafe { cwfa::wavefront_aligner_new(&mut attributes.clone()) })
        } else { None };

        attributes.affine_penalties.match_ = 0;
        attributes.alignment_form.span = cwfa::alignment_span_t_alignment_end2end;
        let global_aligner = unsafe { cwfa::wavefront_aligner_new(&mut attributes) };

        // At this value it could not be beneficial to replace nX with 1I(n-1)=1D.
        let safe_mismatch_size = ((2 * penalties.gap_open + 2 * penalties.gap_extend) / penalties.mismatch) as u32;
        Self { global_aligner, semiglobal_aligner, safe_mismatch_size, penalties }
    }

    pub fn penalties(&self) -> &Penalties {
        &self.penalties
    }

    /// Aligns two sequences (first: ref, second: query), extends `cigar`, and returns alignment score.
    /// If the alignment is dropped, returns VERY approximate alignment.
    ///
    /// If `LEFT_CLIPPING`, replace starting operations until "=" with soft clipping.
    /// !!! Resulting score will be incorrect.
    fn align<const LEFT_CLIPPING: bool>(
        &self,
        aligner: *mut cwfa::wavefront_aligner_t,
        seq1: &[u8],
        seq2: &[u8],
        cigar: &mut Cigar,
    ) -> i32 {
        let status = unsafe { cwfa::wavefront_align(aligner,
            seq1.as_ptr() as *const c_char, seq1.len() as i32,
            seq2.as_ptr() as *const c_char, seq2.len() as i32,
        ) };
        if status != 0 {
            // Alignment was dropped, create approximate alignment.
            return self.penalties.align_simple(seq1, seq2, cigar);
        }

        let c_cigar = unsafe { (*aligner).cigar };
        let begin_offset = usize::try_from(unsafe { (*c_cigar).begin_offset }).unwrap();
        let end_offset = usize::try_from(unsafe { (*c_cigar).end_offset }).unwrap();
        let mut no_matches_yet = true;
        if begin_offset < end_offset {
            let cigar_slice: &[u8] = unsafe {
                std::slice::from_raw_parts((*c_cigar).operations as *const u8, end_offset - begin_offset)
            };
            for &ch in cigar_slice {
                let op = op_from_char(ch);
                if LEFT_CLIPPING && no_matches_yet && op == Operation::Equal {
                    no_matches_yet = false;
                    let soft_clipping = cigar.query_len();
                    cigar.clear();
                    if soft_clipping > 0 {
                        cigar.push_unchecked(Operation::Ins, soft_clipping);
                    }
                }
                cigar.push_checked(op, 1);
            }
        }
        if LEFT_CLIPPING && no_matches_yet {
            let soft_clipping = cigar.query_len();
            cigar.clear();
            if soft_clipping > 0 {
                cigar.push_unchecked(Operation::Ins, soft_clipping);
            }
        }
        let score = unsafe { (*c_cigar).score };
        if score <= -0x60000000_i32 {
            log::warn!("WFA produced very small score ({}). Sequences: {} and {}",
                score, String::from_utf8_lossy(seq1), String::from_utf8_lossy(seq2));
        }
        score
    }

    /// If `i1 == i2` or `j1 == j2` inserts simple INS/DEL,
    /// If gap is bigger than `max_gap` inserts INS/DEL followed by (mis)matches,
    /// If gap can be quickly replaced with mismatches, do that.
    /// Otherwise performs proper alignment of the two subsequences.
    #[inline(always)]
    pub fn smart_align(
        &self,
        seq1: &[u8], // ref sequence
        seq2: &[u8], // query sequence
        i1: u32,
        i2: u32,
        j1: u32,
        j2: u32,
        max_gap: impl Threshold,
        cigar: &mut Cigar,
    ) -> i32
    {
        debug_assert!(i1 <= i2 && j1 <= j2);
        let jump1 = i2 - i1;
        let jump2 = j2 - j1;
        match (jump1 > 0, jump2 > 0) {
            (true, true) => {
                let subseq1 = &seq1[i1 as usize..i2 as usize];
                let subseq2 = &seq2[j1 as usize..j2 as usize];
                if max_gap.under(jump1) || max_gap.under(jump2) {
                    self.penalties.align_simple(subseq1, subseq2, cigar)
                } else if jump1 == jump2 && jump1 <= self.safe_mismatch_size {
                    let mut ndiff = 0;
                    for (&c1, &c2) in itertools::izip!(subseq1, subseq2) {
                        cigar.push_checked(if c1 == c2 { Operation::Equal } else { Operation::Diff }, 1);
                        ndiff -= i32::from(c1 != c2);
                    }
                    ndiff * self.penalties.mismatch
                } else {
                    self.align::<false>(self.global_aligner, subseq1, subseq2, cigar)
                }
            }
            (true, false) => {
                cigar.push_unchecked(Operation::Del, jump1);
                -self.penalties.gap_open - jump1 as i32 * self.penalties.gap_extend
            }
            (false, true) => {
                cigar.push_unchecked(Operation::Ins, jump2);
                -self.penalties.gap_open - jump2 as i32 * self.penalties.gap_extend
            }
            (false, false) => 0,
        }
    }

    pub fn align_clipping<const LEFT: bool>(
        &self,
        seq1: &[u8], // ref sequence
        seq2: &[u8], // query sequence
        i1: u32,
        i2: u32,
        j1: u32,
        j2: u32,
        cigar: &mut Cigar,
    ) {
        assert!(j1 != j2);
        if i1 == i2 {
            // Will be replaced with Operation::Soft later in the analysis.
            cigar.push_unchecked(Operation::Ins, j2 - j1);
            return
        }
        let subseq1 = &seq1[i1 as usize..i2 as usize];
        let subseq2 = &seq2[j1 as usize..j2 as usize];

        let aligner = self.semiglobal_aligner.expect("Semi-global aligner undefined");
        if LEFT {
            unsafe { cwfa::wavefront_aligner_set_alignment_free_ends(
                aligner, subseq1.len() as i32, 0, subseq2.len() as i32, 0) };
        } else {
            unsafe { cwfa::wavefront_aligner_set_alignment_free_ends(
                aligner, 0, subseq1.len() as i32, 0, subseq2.len() as i32) };
        }

        let subseq1 = &seq1[i1 as usize..i2 as usize];
        let subseq2 = &seq2[j1 as usize..j2 as usize];
        log::debug!("Align clipping (left? {}) between {} and {}. Current CIGAR = {:?}",
            LEFT, std::str::from_utf8(subseq1).unwrap(), std::str::from_utf8(subseq2).unwrap(), cigar);
        self.align::<LEFT>(aligner, subseq1, subseq2, cigar);
        if !LEFT {
            let mut soft_clipping = 0;
            loop {
                match cigar.last() {
                    Some(item) if item.operation() != Operation::Equal
                        => soft_clipping += u32::from(item.operation().consumes_query()) * item.len(),
                    _ => break,
                }
                cigar.pop();
            }
            if soft_clipping > 0 {
                cigar.push_unchecked(Operation::Ins, soft_clipping);
            }
        }
        log::debug!("    -> {:?}", cigar);
    }
}

/// Convert char into operation, replacing M with X.
#[inline(always)]
fn op_from_char(ch: u8) -> Operation {
    match ch {
        b'M' | b'=' => Operation::Equal,
        b'X' => Operation::Diff,
        b'I' => Operation::Ins,
        b'D' => Operation::Del,
        b'S' => Operation::Soft,
        _ => panic!("Unexpected CIGAR operation {}", ch as char),
    }
}
