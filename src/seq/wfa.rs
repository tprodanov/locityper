use std::{
    ffi::c_char,
};
use crate::{
    seq::cigar::{Cigar, Operation},
    err::{Error, validate_param},
};

#[allow(non_upper_case_globals)]
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused)]
mod cwfa {
    include!(concat!(env!("OUT_DIR"), "/bindings_wfa.rs"));
}

#[derive(Clone)]
pub struct Penalties {
    // NOTE: Adding match may lead to problems in other places (such as `panvcf` divergence calculation).
    // Similarly, adding separate parameters for insert/deletion penalties may lead to problems elsewhere.
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub more_heuristics: bool,
}

impl Default for Penalties {
    fn default() -> Self {
        Self {
            mismatch: 4,
            gap_open: 6,
            gap_extend: 1,
            more_heuristics: true,
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
}

pub struct Aligner {
    inner: *mut cwfa::wavefront_aligner_t,
    penalties: Penalties,
}

impl Drop for Aligner {
    fn drop(&mut self) {
        unsafe { cwfa::wavefront_aligner_delete(self.inner) }
    }
}

impl Aligner {
    pub fn new(penalties: Penalties) -> Self {
        let mut attributes = unsafe { cwfa::wavefront_aligner_attr_default }.clone();
        if penalties.more_heuristics {
            // Limit the number of alignment steps.
            attributes.system.max_alignment_steps = 1000;
        }

        // Use X-drop heuristic.
        attributes.heuristic.strategy = cwfa::wf_heuristic_strategy_wf_heuristic_xdrop;
        attributes.heuristic.min_wavefront_length = 100;
        attributes.heuristic.steps_between_cutoffs = 100;
        attributes.heuristic.xdrop = if penalties.more_heuristics { 50 } else { 5000 };

        // High memory for some reason produces alignments with M both for X and =.
        attributes.memory_mode = cwfa::wavefront_memory_t_wavefront_memory_med;
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
            let n = seq1.len() as u32;
            let m = seq2.len() as u32;
            let score = if n < m {
                cigar.push_checked(Operation::Diff, n);
                cigar.push_unchecked(Operation::Ins, m - n);
                -self.penalties.mismatch * n as i32
                    - self.penalties.gap_open - self.penalties.gap_extend * (m - n) as i32
            } else if m < n {
                cigar.push_checked(Operation::Diff, m);
                cigar.push_unchecked(Operation::Del, n - m);
                -self.penalties.mismatch * m as i32
                    - self.penalties.gap_open - self.penalties.gap_extend * (n - m) as i32
            } else {
                cigar.push_checked(Operation::Diff, n);
                -self.penalties.mismatch * n as i32
            };
            return Ok(score);
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
