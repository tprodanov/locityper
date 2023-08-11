use std::{
    cmp::max,
    ffi::c_char,
};
use crate::{
    seq::cigar::{Cigar, CigarItem, Operation},
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
    pub mismatch: u32,
    pub gap_opening: u32,
    pub gap_extension: u32,
    // Similarly, adding separate parameters for insert/deletion penalties may lead to problems elsewhere.
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
    /// Returns worst possible score (everything is not aligned).
    pub fn worst_score(&self, mut len1: u32, mut len2: u32) -> u32 {
        if len1 > len2 {
            (len1, len2) = (len2, len1);
        }
        max(self.mismatch * len1 + if len1 == len2 { 0 } else { self.gap_opening + self.gap_extension * (len2 - len1) },
            2 * self.gap_opening + self.gap_extension * (len1 + len2))
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
        let cigar = convert_cigar(c_cigar);
        let score = unsafe { (*c_cigar).score };
        (cigar, score)
    }
}

fn convert_cigar(cigar: *const cwfa::cigar_t) -> Cigar {
    let mut res = Cigar::new();
    let begin_offset = usize::try_from(unsafe { (*cigar).begin_offset }).unwrap();
    let end_offset = usize::try_from(unsafe { (*cigar).end_offset }).unwrap();

    if begin_offset >= end_offset {
        return res;
    }

    let operations = unsafe { (*cigar).operations };
    // Index into c-array.
    let mut last_op = unsafe { *operations.add(begin_offset) } as u8;
    let mut last_len = 1;
    for i in begin_offset + 1..end_offset {
        let curr_op = unsafe { *operations.add(i) } as u8;
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
        b'M' | b'=' => Operation::Equal,
        b'X' => Operation::Diff,
        b'I' => Operation::Ins,
        b'D' => Operation::Del,
        b'S' => Operation::Soft,
        _ => panic!("Unexpected CIGAR operation '{}'", char::from(ch)),
    }
}