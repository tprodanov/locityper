//! Variant discovery.

use smallvec::SmallVec;
use crate::{
    solvers::solve,
    seq::{
        aln::Alignment,
        contigs::{ContigId, ContigNames, Genotype},
        cigar::Operation,
    },
    model::{
        windows::ReadGtAlns,
        assgn::GenotypeAlignments,
    },
    algo::{HashMap, IntSet},
};

/// Go through (potentially paired) alignments and output indices of single-end alignments (from `groupped_alns`)
/// that were observed at least once.
fn find_observed_alignments(
    read_alns: &[ReadGtAlns],
    assgn_counts: &mut impl Iterator<Item = u16>,
    use_all_alns: bool,
    indices: &mut IntSet<u32>,
) {
    indices.clear();
    for aln in read_alns {
        let count = assgn_counts.next().expect("Not enough assignment counts");
        if !use_all_alns && count == 0 { continue };
        let Some((_, aln_pair)) = aln.parent() else { continue }; // Not intertested in unmapped alignments.
        if let Some(i) = aln_pair.ix1() {
            indices.insert(i);
        }
        if let Some(j) = aln_pair.ix2() {
            indices.insert(j);
        }
    }
}

/// Convert contigs into consecutive indices.
/// Returns indices, corresponding contig lengths, the total number of used indices.
fn contigs_to_ixs(gt: &Genotype, contigs: &ContigNames) -> (Vec<Option<usize>>, Vec<u32>) {
    // [TODO] Make a structure, store ploidy of each contig.
    let mut contig_ixs = vec![None; contigs.len()];
    let mut contig_lens = Vec::with_capacity(gt.ploidy());
    let mut curr_ix = 0;
    for &id in gt.ids() {
        if !contig_ixs[id.ix()].is_some() {
            contig_ixs[id.ix()] = Some(curr_ix);
            contig_lens.push(contigs.get_len(id));
            curr_ix += 1;
        }
    }
    (contig_ixs, contig_lens)
}

#[inline(always)]
fn is_boundary_clipping(
    operation: Operation,
    hap_pos: u32, next_hap_pos: u32, hap_len: u32,
) -> bool {
    const BOUNDARY_PADDING: u32 = 10;
    operation == Operation::Soft && (hap_pos < BOUNDARY_PADDING || next_hap_pos + BOUNDARY_PADDING > hap_len)
}

// type Allele = SmallVec<[u8; 8]>;

// /// Key = (start, end, alt sequence), value = count.
// type VarCounter = HashMap<(u32, u32, Allele), u16>;

/// Go through read alignment and record start and end positions of all variations into `breakpoints`.
fn enumerate_aln_differences(
    aln: &Alignment,
    contig_len: u32,
    breakpoints: &mut Vec<(u32, bool)>,
) {
    // let read_len = seq.len() as u32;
    let mut hap_pos = aln.interval().start();
    // let mut read_pos = 0;
    log::debug!("    Cigar {:?}  interval {}", aln.cigar(), aln.interval());
    for entry in aln.cigar().iter() {
        let op = entry.operation();
        let next_hap_pos = hap_pos + u32::from(op.consumes_ref()) * entry.len();
        // let next_read_pos = read_pos + u32::from(op.consumes_query()) * entry.len();
        if op != Operation::Equal && !is_boundary_clipping(op, hap_pos, next_hap_pos, contig_len) {
            let var_start = if op.consumes_both() { hap_pos } else { hap_pos.saturating_sub(1) };
            breakpoints.push((var_start, true));
            breakpoints.push((next_hap_pos, false));
            log::debug!("        Variant {}-{}", var_start+1, next_hap_pos);
        } else {
            log::debug!("        Skip {}-{}", hap_pos+1, next_hap_pos);
        }
        hap_pos = next_hap_pos;
        // read_pos = next_read_pos;
    }
}

/// Identify variable regions with at least `min_support` supporting reads.
fn combine_breakpoints(breakpoints: &mut [Vec<(u32, bool)>], min_support: u32) -> Vec<Vec<(u32, u32)>> {
    let mut var_regions: Vec<Vec<(u32, u32)>> = vec![Vec::new(); breakpoints.len()];
    for (hap_breakpoints, hap_regions) in itertools::izip!(breakpoints, &mut var_regions) {
        log::debug!("Combine breakpoints");
        // let add_region = |start, end| {
        //     match hap_regions.last_mut() {
        //         Some(last) if last.1 == start => last.1 = end,
        //         _ => hap_regions.push((start, end)),
        //     }
        // };

        hap_breakpoints.sort();
        let mut last_pos = 0;
        let mut support = 0;
        for &(pos, is_start) in hap_breakpoints.iter() {
            log::debug!("    {:5} {}, [{}]", pos, if is_start { "start" } else { " end " }, support);
            if pos > last_pos && support >= min_support {
                match hap_regions.last_mut() {
                    Some(last) if last.1 == last_pos => last.1 = pos,
                    _ => hap_regions.push((last_pos, pos)),
                }
            }
            support = support.strict_add_signed(if is_start { 1 } else { -1 });
            last_pos = pos;
        }
    }
    var_regions
}

/// Identify potential variants across reads.
pub fn discover_variants(
    gt: &Genotype,
    contigs: &ContigNames,
    data: &solve::Data,
    assgn_counts: &[u16],
    assgn_attempts: u16,
) {
    let (contig_ixs, contig_lens) = contigs_to_ixs(gt, contigs);
    let n_haps = contig_lens.len();

    const MIN_ATTEMPTS: u16 = 10;
    let use_all_alns = assgn_attempts < MIN_ATTEMPTS;
    if use_all_alns {
        log::warn!("Too few attempts ({}) on the last stage, using all alignments for variant discovery",
            assgn_attempts);
    }
    let gt_alns = GenotypeAlignments::new(gt.clone(), &data.contig_infos, &data.all_alns, &data.assgn_params);
    let mut counts_iter = assgn_counts.iter().copied();
    let mut ixs_buffer = IntSet::default();
    let mut breakpoints = vec![Default::default(); n_haps];

    log::debug!("Processing reads");
    for (rp, groupped_alns) in data.all_alns.reads().iter().enumerate() {
        // let read_data = groupped_alns.read_data();
        let read_alns = gt_alns.possible_read_alns(rp);
        find_observed_alignments(read_alns, &mut counts_iter, use_all_alns, &mut ixs_buffer);
        for &i in &ixs_buffer {
            let aln = groupped_alns.ith_aln(i);
            // let seq = read_data.mate_data(aln.read_end()).get_seq(aln.strand());
            let contig_ix = contig_ixs[aln.contig_id().ix()].expect("Contig must be in the list");
            enumerate_aln_differences(aln, contig_lens[contig_ix], &mut breakpoints[contig_ix]);
        }
    }

    // [TODO] Pass min_support from arguments.
    let var_regions = combine_breakpoints(&mut breakpoints, 3);
    for (i, hap_regions) in var_regions.iter().enumerate() {
        log::debug!("Haplotype {}:", i);
        for &(start, end) in hap_regions {
            log::debug!("    {}-{}", start+1, end);
        }
    }
}
