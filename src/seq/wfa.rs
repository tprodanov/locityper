use std::ffi::c_char;
use crate::{
    err::{Error, validate_param},
    seq::cigar::{Cigar, CigarItem, Operation},
};

#[allow(non_upper_case_globals)]
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused)]
mod cwfa {
    include!(concat!(env!("OUT_DIR"), "/bindings_wfa.rs"));
}

pub struct Penalties {
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
        validate_param!(self.mismatch > 0, "Mismatch penalty ({}) must be positive", self.mismatch);
        validate_param!(self.gap_opening > 0, "Open gap penalty ({}) must be positive", self.gap_opening);
        validate_param!(self.gap_extension > 0, "Extend gap penalty ({}) must be positive", self.gap_extension);
        Ok(())
    }
}

pub struct Aligner(*mut cwfa::wavefront_aligner_t);

impl Drop for Aligner {
    fn drop(&mut self) {
        unsafe { cwfa::wavefront_aligner_delete(self.0) }
    }
}

impl Aligner {
    pub fn new(penalties: &Penalties) -> Self {
        let mut attributes = unsafe { cwfa::wavefront_aligner_attr_default }.clone();
        // Do not use a heuristic (enabled by default).
        attributes.heuristic.strategy = cwfa::wf_heuristic_strategy_wf_heuristic_none;
        // Compute score and CIGAR as well.
        attributes.alignment_scope = cwfa::alignment_scope_t_compute_alignment;
        // Compute global alignment.
        attributes.alignment_form.span = cwfa::alignment_span_t_alignment_end2end;

        // Set the cost model and parameters.
        attributes.distance_metric = cwfa::distance_metric_t_gap_affine;
        attributes.affine_penalties.mismatch = penalties.mismatch;
        attributes.affine_penalties.gap_opening = penalties.gap_opening;
        attributes.affine_penalties.gap_extension = penalties.gap_extension;

        Self(unsafe { cwfa::wavefront_aligner_new(&mut attributes) })
    }

    pub fn align(&self, seq1: &[u8], seq2: &[u8]) -> (Cigar, i32) {
        let status = unsafe { cwfa::wavefront_align(
            self.0,
            seq1.as_ptr() as *const c_char,
            seq1.len() as i32,
            seq2.as_ptr() as *const c_char,
            seq2.len() as i32,
        ) };
        assert_eq!(status, 0, "WFA alignment failed");
        let c_cigar = unsafe { (*self.0).cigar };
        let cigar = unsafe { convert_cigar(c_cigar) };
        let score = unsafe { (*c_cigar).score };
        (cigar, score)
    }
}

unsafe fn convert_cigar(cigar: *mut cwfa::cigar_t) -> Cigar {
    let mut res = Cigar::new();
    let begin_offset = (*cigar).begin_offset as isize;
    let end_offset = (*cigar).end_offset as isize;

    if begin_offset >= end_offset {
        return res;
    }

    let operations = (*cigar).operations;
    // Index into c-array.
    let mut last_op = *operations.offset(begin_offset) as u8;
    let mut last_len = 1;
    for i in begin_offset + 1..end_offset {
        let curr_op = *operations.offset(i) as u8;
        if last_op == curr_op {
            last_len += 1;
        } else {
            res.push(CigarItem::new(op_from_char(last_op), last_len));
            last_op = curr_op;
            last_len = 1;
        }
    }
    res.push(CigarItem::new(op_from_char(last_op), last_len));
    res
}


/// Convert char into operation, replacing M with X.
fn op_from_char(ch: u8) -> Operation {
    match ch {
        b'M' => Operation::Equal,
        b'=' => Operation::Equal,
        b'X' => Operation::Diff,
        b'I' => Operation::Ins,
        b'D' => Operation::Del,
        b'S' => Operation::Soft,
        _ => panic!("Unexpected CIGAR operation '{}'", char::from(ch)),
    }
}