use std::{
    ffi::c_char,
};
use crate::{
    seq::cigar::{Cigar, CigarItem, Operation},
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
    pub gap_opening: i32,
    pub gap_extension: i32,
}

impl Default for Penalties {
    fn default() -> Self {
        Self {
            mismatch: 4,
            gap_opening: 6,
            gap_extension: 1,
        }
    }
}

impl Penalties {
    pub fn validate(&self) -> Result<(), Error> {
        validate_param!(self.mismatch >= 0, "Mismatch penalty ({}) must be non-negative", self.mismatch);
        validate_param!(self.gap_opening >= 0, "Gap opening penalty ({}) must be non-negative", self.gap_opening);
        validate_param!(self.gap_extension >= 0, "Gap extension penalty ({}) must be non-negative", self.gap_extension);
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
        // Use Adaptive heuristic.
        attributes.heuristic.strategy = cwfa::wf_heuristic_strategy_wf_heuristic_wfadaptive;
        attributes.heuristic.min_wavefront_length = 50;
        attributes.heuristic.max_distance_threshold = 1000;
        attributes.heuristic.steps_between_cutoffs = 1;
        // Use less memory at the expense of running time.
        attributes.memory_mode = cwfa::wavefront_memory_t_wavefront_memory_low;

        // Compute score and CIGAR as well.
        attributes.alignment_scope = cwfa::alignment_scope_t_compute_alignment;
        // Compute global alignment.
        attributes.alignment_form.span = cwfa::alignment_span_t_alignment_end2end;

        // Set the cost model and parameters.
        attributes.distance_metric = cwfa::distance_metric_t_gap_affine;
        attributes.affine_penalties.mismatch = penalties.mismatch as i32;
        attributes.affine_penalties.gap_opening = penalties.gap_opening as i32;
        attributes.affine_penalties.gap_extension = penalties.gap_extension as i32;

        Self {
            inner: unsafe { cwfa::wavefront_aligner_new(&mut attributes) },
            penalties,
        }
    }

    pub fn penalties(&self) -> &Penalties {
        &self.penalties
    }

    /// Aligns two sequences and returns pair (Cigar, alignment score).
    /// If Cigar creation fails, returns violating character.
    pub fn align(&self, seq1: &[u8], seq2: &[u8]) -> Result<(Cigar, i32), Error> {
        let status = unsafe { cwfa::wavefront_align(
            self.inner,
            seq1.as_ptr() as *const c_char,
            seq1.len() as i32,
            seq2.as_ptr() as *const c_char,
            seq2.len() as i32,
        ) };
        assert_eq!(status, 0, "WFA alignment failed");
        let c_cigar = unsafe { (*self.inner).cigar };
        let cigar = convert_cigar(c_cigar)
            .map_err(|(ch, raw_cigar)| Error::RuntimeError(format!(
                "Could not align {} and {}. Violating CIGAR character '{}' ({}) in {:?}",
                String::from_utf8_lossy(seq1), String::from_utf8_lossy(seq2), char::from(ch), ch, raw_cigar)))?;
        let score = unsafe { (*c_cigar).score };
        Ok((cigar, score))
    }
}

fn convert_cigar(c_cigar: *const cwfa::cigar_t) -> Result<Cigar, (u8, Vec<u8>)> {
    let mut cigar = Cigar::new();
    let begin_offset = usize::try_from(unsafe { (*c_cigar).begin_offset }).unwrap();
    let end_offset = usize::try_from(unsafe { (*c_cigar).end_offset }).unwrap();
    if begin_offset >= end_offset {
        return Ok(cigar);
    }

    let cigar_slice: &[u8] = unsafe {
        std::slice::from_raw_parts((*c_cigar).operations as *const u8, end_offset - begin_offset)
    };
    let mut last_op = cigar_slice[0];
    let mut last_len = 1;
    for &curr_op in &cigar_slice[1..] {
        if last_op == curr_op {
            last_len += 1;
        } else {
            let op = op_from_char(last_op).map_err(|_| (last_op, cigar_slice.to_vec()))?;
            cigar.push(CigarItem::new(op, last_len));
            last_op = curr_op;
            last_len = 1;
        }
    }
    let op = op_from_char(last_op).map_err(|_| (last_op, cigar_slice.to_vec()))?;
    cigar.push(CigarItem::new(op, last_len));
    Ok(cigar)
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
