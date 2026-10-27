use std::{
    ffi::c_char,
    cmp::min,
};
use crate::{
    seq::cigar::{Cigar, Operation},
    err::validate_param,
    ext::opt::{Optional, Flag, Off, With},
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
            cigar.push_checked(Operation::Ins, (-diff) as u32);
            score = -self.gap_open + diff * self.gap_extend;
            (0, (-diff) as usize)
        } else if diff > 0 {
            cigar.push_checked(Operation::Del, diff as u32);
            score = -self.gap_open - diff * self.gap_extend;
            (diff as usize, 0)
        } else {
            score = 0;
            (0, 0)
        };
        for (&nt1, &nt2) in seq1[i + 1..].iter().zip(&seq2[j + 1..]) {
            cigar.push_checked(if nt1 == nt2 { Operation::Equal } else { Operation::Diff }, 1);
            score -= if nt1 != nt2 { self.mismatch } else { 0 };
        }
        score
    }

    /// Calculates score of an extended CIGAR (does not include M symbol).
    pub fn calculate_score(&self, cigar: &Cigar) -> i32 {
        let mut score = 0;
        for item in cigar.iter() {
            let len = item.len() as i32;
            score -= match item.operation() {
                Operation::Equal => 0,
                Operation::Diff => self.mismatch * len,
                Operation::Soft | Operation::Ins | Operation::Del => self.gap_open + self.gap_extend * len,
                _ => panic!("`calculate_score`: Not supported operation {}", item.operation()),
            }
        }
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
        band: Option<i32>,
        add_semiglobal: bool,
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

        if let Some(k) = band {
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

        let semiglobal_aligner = if add_semiglobal {
            // Need positive match score for alignment to work.
            attributes.affine_penalties.match_ = min(-1, -attributes.affine_penalties.mismatch / 2);
            attributes.alignment_form.span = cwfa::alignment_span_t_alignment_endsfree;
            Some(unsafe { cwfa::wavefront_aligner_new(&mut attributes.clone()) })
        } else {
            None
        };

        attributes.affine_penalties.match_ = 0;
        attributes.alignment_form.span = cwfa::alignment_span_t_alignment_end2end;
        let global_aligner = unsafe { cwfa::wavefront_aligner_new(&mut attributes) };

        // At this value it could not be beneficial to replace nX with 1I(n-1)=1D.
        let safe_mismatch_size = ((2 * penalties.gap_open + 2 * penalties.gap_extend) / penalties.mismatch) as u32;
        Self { global_aligner, semiglobal_aligner, safe_mismatch_size, penalties }
    }

    #[inline(always)]
    pub fn penalties(&self) -> &Penalties {
        &self.penalties
    }

    /// Aligns two sequences (first: ref, second: query), extending `cigar`.
    /// If the alignment is dropped, returns VERY approximate alignment.
    /// Returns offset in the reference sequence and alignment score.
    /// Offset = 0 unless left clipping was used.
    ///
    /// If `CLIP_LEFT`, replace starting operations until "=" with soft clipping.
    /// !!! Resulting score will be incorrect.
    fn align(
        &self,
        aligner: *mut cwfa::wavefront_aligner_t,
        ref_seq: &[u8],
        query_seq: &[u8],
        cigar: &mut Cigar,
        clip_left: impl Flag,
    ) -> (u32, i32) {
        let status = unsafe { cwfa::wavefront_align(aligner,
            ref_seq.as_ptr() as *const c_char, ref_seq.len() as i32,
            query_seq.as_ptr() as *const c_char, query_seq.len() as i32,
        ) };
        if status != 0 {
            // Alignment was dropped, create approximate alignment.
            return (0, self.penalties.align_simple(ref_seq, query_seq, cigar));
        }

        let c_cigar = unsafe { (*aligner).cigar };
        let begin_slice_offset = usize::try_from(unsafe { (*c_cigar).begin_offset }).unwrap();
        let end_slice_offset = usize::try_from(unsafe { (*c_cigar).end_offset }).unwrap();
        let mut ref_offset = 0;
        let mut no_matches_yet = true;
        if begin_slice_offset < end_slice_offset {
            let cigar_slice: &[u8] = unsafe {
                std::slice::from_raw_parts((*c_cigar).operations as *const u8, end_slice_offset - begin_slice_offset)
            };
            for &ch in cigar_slice {
                let op = op_from_char(ch);
                if clip_left.enabled() && no_matches_yet {
                    if op == Operation::Equal {
                        no_matches_yet = false;
                        let soft_clipping = cigar.query_len();
                        cigar.clear();
                        if soft_clipping > 0 {
                            // Will be replaced by soft clipping later.
                            cigar.push_checked(Operation::Ins, soft_clipping);
                        }
                    } else {
                        ref_offset += u32::from(op.consumes_ref());
                    }
                }
                cigar.push_checked(op, 1);
            }
        }
        if clip_left.enabled() && no_matches_yet {
            let soft_clipping = cigar.query_len();
            cigar.clear();
            if soft_clipping > 0 {
                cigar.push_checked(Operation::Ins, soft_clipping);
            }
        }
        let score = unsafe { (*c_cigar).score };
        if score <= -0x60000000_i32 {
            #[cold]
            fn warn(score: i32, ref_seq: &[u8], query_seq: &[u8]) {
                log::warn!("WFA produced very small score ({}). Sequences: {} and {}",
                    score, String::from_utf8_lossy(ref_seq), String::from_utf8_lossy(query_seq));
            }
            warn(score, ref_seq, query_seq);
        }
        (ref_offset, score)
    }

    /// If `i1 == i2` or `j1 == j2` inserts simple INS/DEL,
    /// If gap is bigger than `max_gap` inserts INS/DEL followed by (mis)matches,
    /// If gap can be quickly replaced with mismatches, do that.
    /// Otherwise performs proper alignment of the two subsequences.
    #[inline(always)]
    pub fn smart_align(
        &self,
        seq1: &[u8], // ref sequence
        i1: u32,
        i2: u32,
        seq2: &[u8], // query sequence
        j1: u32,
        j2: u32,
        max_gap: impl Optional<u32>,
        cigar: &mut Cigar,
    ) -> i32 {
        debug_assert!(i1 <= i2 && j1 <= j2);
        let jump1 = i2 - i1;
        let jump2 = j2 - j1;
        match (jump1 > 0, jump2 > 0) {
            (true, true) => {
                let subseq1 = &seq1[i1 as usize..i2 as usize];
                let subseq2 = &seq2[j1 as usize..j2 as usize];
                if max_gap.get().map(|m| jump1 > m || jump2 > m).unwrap_or(false) {
                    self.penalties.align_simple(subseq1, subseq2, cigar)
                } else if jump1 == jump2 && jump1 <= self.safe_mismatch_size {
                    let mut ndiff = 0;
                    for (&c1, &c2) in itertools::izip!(subseq1, subseq2) {
                        cigar.push_checked(if c1 == c2 { Operation::Equal } else { Operation::Diff }, 1);
                        ndiff -= i32::from(c1 != c2);
                    }
                    ndiff * self.penalties.mismatch
                } else {
                    self.align(self.global_aligner, subseq1, subseq2, cigar, Off).1
                }
            }
            (true, false) => {
                cigar.push_checked(Operation::Del, jump1);
                -self.penalties.gap_open - jump1 as i32 * self.penalties.gap_extend
            }
            (false, true) => {
                cigar.push_checked(Operation::Ins, jump2);
                -self.penalties.gap_open - jump2 as i32 * self.penalties.gap_extend
            }
            (false, false) => 0,
        }
    }

    /// Use semi-global alignment to align sequences in either left-free or right-free mode.
    /// Returns start offset.
    /// `clip_left` and `clip_right` are either `Off`, or `With(ref_free_end, read_free_end)`.
    /// If `clip_left` is not `Off`, replace all relevant starting non-matches with soft clipping.
    pub fn align_free_ends(
        &self,
        ref_seq: &[u8],
        query_seq: &[u8],
        cigar: &mut Cigar,
        clip_left: impl Optional<(i32, i32)>,
        clip_right: impl Optional<(i32, i32)>,
    ) -> u32 {
        assert!(!query_seq.is_empty());
        if ref_seq.is_empty() {
            // Will be replaced by Operation::Soft later in the analysis.
            cigar.push_checked(Operation::Ins, query_seq.len() as u32);
            return 0;
        }
        let aligner = self.semiglobal_aligner.expect("Semi-global aligner undefined");

        let (free1_left, free2_left) = clip_left.get().unwrap_or((0, 0));
        let (free1_right, free2_right) = clip_right.get().unwrap_or((0, 0));
        unsafe {
            cwfa::wavefront_aligner_set_alignment_free_ends(aligner, free1_left, free1_right, free2_left, free2_right);
        };
        let (start_offset, _) = self.align(aligner, ref_seq, query_seq, cigar, clip_left.flag());
        if clip_right.enabled() {
            let mut soft_clipping = 0;
            while let Some(item) = cigar.pop_if(|item| item.operation() != Operation::Equal) {
                soft_clipping += u32::from(item.operation().consumes_query()) * item.len();
            }
            if soft_clipping > 0 {
                // Will be replaced by soft clipping later.
                cigar.push_checked(Operation::Ins, soft_clipping);
            }
        }
        start_offset
    }

    /// Aligns two sequences with free left end and fixed right end.
    /// Returns start offset in the reference sequence.
    #[inline(always)]
    pub fn align_left_end(&self, ref_seq: &[u8], query_seq: &[u8], cigar: &mut Cigar) -> u32 {
        self.align_free_ends(ref_seq, query_seq, cigar, With((ref_seq.len() as i32, query_seq.len() as i32)), Off)
    }

    /// Aligns two sequences with fixed left end and free right end.
    #[inline(always)]
    pub fn align_right_end(&self, ref_seq: &[u8], query_seq: &[u8], cigar: &mut Cigar) {
        self.align_free_ends(ref_seq, query_seq, cigar, Off, With((ref_seq.len() as i32, query_seq.len() as i32)));
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
