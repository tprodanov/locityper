#[cfg(feature = "align")]
use std::ffi::c_char;
use crate::{
    seq::cigar::{Cigar, Operation},
    err::{Error, validate_param},
};

#[cfg(feature = "align")]
#[allow(non_upper_case_globals)]
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused)]
mod cwfa {
    include!(concat!(env!("OUT_DIR"), "/bindings_wfa.rs"));
}

/// Maximum accuracy level (need to change `add` help message, if changed this).
pub const MAX_ACCURACY: u8 = 9;

#[derive(Clone)]
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
    pub fn validate(&self) -> Result<(), Error> {
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
#[cfg(feature = "align")]
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

#[cfg(feature = "align")]
pub struct Aligner {
    inner: *mut cwfa::wavefront_aligner_t,
    penalties: Penalties,
}

#[cfg(feature = "align")]
impl Drop for Aligner {
    fn drop(&mut self) {
        unsafe { cwfa::wavefront_aligner_delete(self.inner) }
    }
}

#[cfg(feature = "align")]
impl Aligner {
    /// Accuracy level must be in [1, 9], 0 - very fast and inaccurate, 9 - very slow and accurate.
    pub fn new(penalties: Penalties, accuracy: u8) -> Self {
        assert!(1 <= accuracy && accuracy <= MAX_ACCURACY, "Cannot construct WFA aligner for accuracy {}", accuracy);
        let mut attributes = unsafe { cwfa::wavefront_aligner_attr_default }.clone();
        // Limit the number of alignment steps.
        attributes.system.max_alignment_steps = alignment_steps(accuracy);

        // High memory for some reason produces alignments with M both for X and =.
        attributes.memory_mode = if accuracy < 7 {
            cwfa::wavefront_memory_t_wavefront_memory_med
        } else {
            cwfa::wavefront_memory_t_wavefront_memory_ultralow
        };
        // Compute score and CIGAR as well.
        attributes.alignment_scope = cwfa::alignment_scope_t_compute_alignment;
        // Compute global alignment.
        attributes.alignment_form.span = cwfa::alignment_span_t_alignment_end2end;

        // Set the cost model and parameters.
        attributes.distance_metric = cwfa::distance_metric_t_gap_affine;
        attributes.affine_penalties.mismatch = penalties.mismatch;
        attributes.affine_penalties.gap_opening = penalties.gap_open;
        attributes.affine_penalties.gap_extension = penalties.gap_extend;

        Self {
            inner: unsafe { cwfa::wavefront_aligner_new(&mut attributes) },
            penalties,
        }
    }

    pub fn penalties(&self) -> &Penalties {
        &self.penalties
    }

    /// Aligns two sequences (first: ref, second: query), extends `cigar`, and returns alignment score.
    /// If the alignment is dropped, returns VERY approximate alignment.
    pub fn align(&self, seq1: &[u8], seq2: &[u8], cigar: &mut Cigar) -> Result<i32, Error> {
        let status = unsafe { cwfa::wavefront_align(
            self.inner,
            seq1.as_ptr() as *const c_char,
            seq1.len() as i32,
            seq2.as_ptr() as *const c_char,
            seq2.len() as i32,
        ) };
        if status != 0 {
            // Alignment was dropped, create approximate alignment.
            return Ok(self.penalties.align_simple(seq1, seq2, cigar));
        }

        let c_cigar = unsafe { (*self.inner).cigar };
        let begin_offset = usize::try_from(unsafe { (*c_cigar).begin_offset }).unwrap();
        let end_offset = usize::try_from(unsafe { (*c_cigar).end_offset }).unwrap();
        if begin_offset < end_offset {
            let cigar_slice: &[u8] = unsafe {
                std::slice::from_raw_parts((*c_cigar).operations as *const u8, end_offset - begin_offset)
            };
            for &ch in cigar_slice {
                let op = op_from_char(ch).map_err(|_| Error::RuntimeError(format!(
                    "Could not align two sequences: violating CIGAR character `{}` ({}) in {:?}. Sequences: {} and {}",
                    char::from(ch), ch, cigar_slice, String::from_utf8_lossy(seq1), String::from_utf8_lossy(seq2))))?;
                cigar.push_checked(op, 1);
            }
        }
        Ok(unsafe { (*c_cigar).score })
    }
}

/// Convert char into operation, replacing M with X.
#[inline]
#[cfg(feature = "align")]
fn op_from_char(ch: u8) -> Result<Operation, ()> {
    match ch {
        b'M' | b'=' => Ok(Operation::Equal),
        b'X' => Ok(Operation::Diff),
        b'I' => Ok(Operation::Ins),
        b'D' => Ok(Operation::Del),
        b'S' => Ok(Operation::Soft),
        _ => Err(()),
    }
}
