//! Transfer alignments from one haplotype to another.

use std::{
    sync::Arc,
    cmp::{min, max},
};
use vers_vecs::rmq::fast_rmq::FastRmq;
use crate::{
    seq::{
        Interval,
        contigs::{ContigId, ContigNames, ContigSet},
        cigar::{Cigar, CigarIndex},
        paf::{PafEntry},
        aln::Alignment,
        wfa::Aligner,
        compl::linguistic_complexity,
    },
    bg::ErrorProfile,
    model::locs::{ReadData, PrelimAlignments},
    ext::{
        TriangleMatrix,
        opt::With,
    }
};

/// Store ambiguity for each 32 bp.
const STEP_PWR: u8 = 5;
/// 64 bp.
const WINDOW_SIZE: usize = 1 << (STEP_PWR + 1);

/// For one haplotype answers queries where alignments are ambiguous.
struct AmbiguousRegions {
    /// For each haplotype, store a range-minimum-query structure on the number of unique k-mers in
    /// the next 64 bp. So, we have counts for 64 bp windows stored for each mode 32 bp position.
    haplotype_rmqs: Vec<FastRmq>,
}

impl AmbiguousRegions {
    fn new(contig_set: &ContigSet) -> Self {
        const KMER_SIZE: u8 = 5;
        // const ALL_UNIQUE: u16 = (WINDOW_SIZE + 1 - KMER_SIZE as usize) as u16;

        let mut haplotype_rmqs = Vec::with_capacity(contig_set.len());
        let mut buffer = Vec::<u32>::new();
        let mut n_unique = Vec::<u16>::new();

        for contig_id in contig_set.contigs().ids() {
            let seq = contig_set.get_seq(contig_id);
            linguistic_complexity(seq, KMER_SIZE, WINDOW_SIZE, &mut buffer, &mut n_unique);
            let n = (seq.len() - 1) >> STEP_PWR;
            let mut rmq = Vec::<u64>::with_capacity(n);
            let m = min(n, ((n_unique.len() - 1) >> STEP_PWR) + 1);
            for i in 0..m {
                rmq.push(u64::from(n_unique[i << STEP_PWR]));
            }
            let last = u64::from(*n_unique.last().unwrap());
            rmq.resize(n, last);
            haplotype_rmqs.push(FastRmq::from_vec(rmq));
        }
        Self { haplotype_rmqs }
    }

    fn min_unique_kmers(&self, contig_id: ContigId, start: u32, end: u32) -> u64 {
        let rmq = &self.haplotype_rmqs[contig_id.ix()];
        let a = (start >> STEP_PWR) as usize;
        // Because windows overlap, we want to take a previous window for the end index.
        let b = max(a, ((end - 1) >> STEP_PWR).saturating_sub(1) as usize);
        let i = self.haplotype_rmqs[contig_id.ix()].range_min(a, b);
        rmq[i]
    }

    #[inline]
    fn ambiguous(&self, contig_id: ContigId, start: u32, end: u32) -> bool {
        /// Expect at least 2/3 of unique k-mers.
        const THRESHOLD: u64 = WINDOW_SIZE as u64 * 2 / 3;
        self.min_unique_kmers(contig_id, start, end) < THRESHOLD
    }
}

/// Haplotype-haplotype alignments.
pub struct HapAlns {
    /// Alignments from each contig to other contigs.
    /// matrix[i, j] contains alignment where contig i is query, j is reference.
    aln_matrix: TriangleMatrix<Option<(Cigar, CigarIndex)>>,
    ambiguities: AmbiguousRegions,
    transfer_fails: u32,

    local_step: u32,
    /// Haplotype -> Window -> offset in `locally_similar_ixs`. Inner vector has size of (n_windows + 1).
    contig_offsets: Vec<Vec<usize>>,
    /// Slice [contig_offset[contig][window] .. contig_offset[contig][window + 1]] will contain a sorted list of
    /// other locally similar contigs.
    locally_similar_ixs: Vec<ContigId>,
}

impl HapAlns {
    pub fn new(contig_set: &ContigSet, transfer_fails: u32) -> Self {
        Self {
            aln_matrix: TriangleMatrix::new(contig_set.len(), None),
            ambiguities: AmbiguousRegions::new(contig_set),
            transfer_fails,
            local_step: 0,
            contig_offsets: Vec::new(),
            locally_similar_ixs: Vec::new(),
        }
    }

    pub fn add(&mut self, entry: &mut PafEntry, contigs: &ContigNames, max_div: f64) {
        let id1 = entry.query_id();
        let id2 = entry.target_id();
        let matrix_cell = self.aln_matrix.get_symmetric_mut(id1.ix(), id2.ix());
        if matrix_cell.is_some() { return };

        if !entry.full_positive_alignment() {
            log::warn!("Alignment between {} and {} is on the reverse strand or does not fully cover both sequences",
                contigs.get_name(id1), contigs.get_name(id2));
            return;
        }
        if entry.divergence().unwrap_or(f64::INFINITY) > max_div { return };
        let Some(mut cigar) = entry.take_cigar() else { return };
        if id1 > id2 {
            cigar = cigar.invert();
        }
        let cigar_index = CigarIndex::new(&cigar);
        *matrix_cell = Some((cigar, cigar_index));
    }

    /// Identify closest local haplotypes.
    pub fn create_locally_similar_lists(
        &mut self,
        contigs: &ContigNames,
        basis_contigs: &[ContigId],
        window: u32,
        step: u32,
        max_div: f64,
    ) {
        let max_edit = (f64::from(window) * max_div).floor() as u32;
        let n_contigs = contigs.len();
        let mut contig_offsets = vec![Vec::new(); n_contigs];
        let mut locally_similar_ixs = Vec::new();

        let mut buffer = Vec::new();
        let max_len = contigs.lengths().iter().max().expect("At least one contig must be present");
        // Temporary structure: window index -> [(haplotype index, edit size)]
        let mut contig_locally_similar = vec![Vec::new(); max_len.saturating_sub(window).div_ceil(step) as usize + 1];
        for &contig_id in basis_contigs {
            let contig_len = contigs.get_len(contig_id);
            let n_windows = contig_len.saturating_sub(window).div_ceil(step) as usize + 1;
            // This way we will for certain not get out of bounds.
            let curr_locally_similar = &mut contig_locally_similar[..n_windows];

            let i = contig_id.ix();
            for j in 0..n_contigs {
                if i == j { continue };
                let Some((cigar, _)) = self.aln_matrix.get_symmetric(i, j).as_ref() else { continue };
                let in_query = i < j;
                cigar.locally_similar::<(u32, u32)>(in_query, window, step, max_edit, &mut buffer);
                for &(window, edit) in &buffer {
                    curr_locally_similar[window as usize].push((j, edit));
                }
            }

            let curr_offsets = &mut contig_offsets[i];
            curr_offsets.reserve_exact(n_windows + 1);
            curr_offsets.push(locally_similar_ixs.len());
            for similar_at_window in curr_locally_similar {
                similar_at_window.sort_unstable_by_key(|tup| tup.1);
                locally_similar_ixs.extend(similar_at_window.iter().map(|tup| ContigId::new(tup.0)));
                curr_offsets.push(locally_similar_ixs.len());
                similar_at_window.clear();
            }
            assert!(curr_offsets.len() == n_windows + 1);
        }

        self.local_step = step;
        self.contig_offsets = contig_offsets;
        self.locally_similar_ixs = locally_similar_ixs;
    }
}

/// Constructs new alignment to the target region.
/// Returns alignment start and CIGAR.
fn align_anew(
    aligner: &Aligner,
    target_seq: &[u8],
    approx_target_start: u32,
    approx_target_end: u32,
    read_seq: &[u8],
    read_clipping: (u32, u32),
) -> (u32, Cigar) {
    // How much to add to the sides of the reference sequence
    const PADDING: u32 = 20;
    // Inside the region corresponding to the reference alignment,
    // allow this size of free indels on the sides.
    const FREE_MIDDLE: u32 = 50;
    let target_len = target_seq.len() as u32;
    let padding_left = min(PADDING, approx_target_start);
    let padded_target_start = approx_target_start - padding_left;
    let padding_right = min(PADDING, target_len - approx_target_end);
    let padded_target_end = approx_target_end + padding_right;
    let target_subseq = &target_seq[padded_target_start as usize..padded_target_end as usize];

    let read_free_left = min((read_clipping.0 + FREE_MIDDLE) as i32, read_seq.len() as i32);
    let read_free_right = min((read_clipping.1 + FREE_MIDDLE) as i32, read_seq.len() as i32);
    let padded_len = padded_target_end - padded_target_start;
    let target_free_left = min(padding_left + FREE_MIDDLE, padded_len) as i32;
    let target_free_right = min(padding_right + FREE_MIDDLE, padded_len) as i32;

    let mut cigar = Cigar::new();
    let aln_start = padded_target_start + aligner.align_free_ends(target_subseq, read_seq, &mut cigar,
        With((target_free_left, read_free_left)), With((target_free_right, read_free_right)));
    cigar.boundary_ins_to_soft();
    (aln_start, cigar)
}

impl HapAlns {
    /// Try to identify any new alignments for a given read/read pair.
    /// Returns the number of new alignments and a total number of attempted transfers.
    pub(crate) fn transfer_alignments(
        &self,
        prelim_alignments: &mut PrelimAlignments,
        read_data: &ReadData,
        contig_set: &ContigSet,
        aligner: &Aligner,
        err_prof: &ErrorProfile,
    ) -> (usize, usize) {
        let contigs = contig_set.contigs();
        let n = prelim_alignments.len();
        let mut seen = vec![false; n];
        let mut total_transfers = 0;

        for i in 0..n {
            if std::mem::replace(&mut seen[i], true) {
                continue;
            }
            let source_aln = &prelim_alignments[i];
            // Copy various variables before `source_aln` is dropped to make `prelim_alignments` mutable.
            let source_interval = source_aln.interval().clone();
            let source_cigar = source_aln.cigar().clone();
            let read_clipping = source_cigar.soft_clipping();
            let source_strand = source_aln.strand();
            let read_end = source_aln.read_end();
            let read_seq = read_data.mate_data(read_end).get_seq(source_strand);
            let mut n_fails = 0..self.transfer_fails;

            let curr_offsets = &self.contig_offsets[source_interval.contig_id().ix()];
            assert!(!curr_offsets.is_empty(), "Trying to reconstruct alignments not from a basis haplotype");
            let n_windows = curr_offsets.len().strict_sub(1);
            // If we are in the last window, start / step can put us beyond the size of the vector.
            let k = min(n_windows, (source_interval.start() / self.local_step + 1) as usize);
            // Safety: we know that k is at most n_windows, and length = n_windows + 1.
            let start_offset = unsafe { *curr_offsets.get_unchecked(k.strict_sub(1)) };
            let end_offset = unsafe { *curr_offsets.get_unchecked(k) };
            for &target_contig_id in &self.locally_similar_ixs[start_offset..end_offset] {
                let target_seq = contig_set.get_seq(target_contig_id);
                let (haps_cigar, cigar_index) = self.aln_matrix
                    .get_symmetric(source_interval.contig_id().ix(), target_contig_id.ix())
                    .as_ref().expect("Alignment between haplotypes must exist");
                let query_to_ref = source_interval.contig_id() < target_contig_id;

                let offset = cigar_index.find_cigar_offset(source_interval.start(), query_to_ref);
                if let Some(prev_aln_entry) = prelim_alignments.pos_collection()
                    .get(read_end, target_contig_id, offset.rpos())
                {
                    if let Some(seen_val) = seen.get_mut(prev_aln_entry.index() as usize) {
                        *seen_val = true;
                    }
                    if !prev_aln_entry.allowed_to_improve() {
                        continue;
                    }
                }

                let approx_target_start = offset.rpos();
                let approx_target_end = cigar_index.find_cigar_offset(
                    source_interval.end() - 1, query_to_ref).rpos() + 1;
                let ambiguous = self.ambiguities.ambiguous(target_contig_id, approx_target_start, approx_target_end);
                // log::debug!("Transferring {} R{} to {}, ambiguous? {}, n unique {}",
                //     read_data.name(), read_end, contigs.get_name(target_contig_id), ambiguous,
                //     self.ambiguities.min_unique_kmers(target_contig_id, approx_target_start, approx_target_end),
                // );
                // let (new_start, new_cigar) = if ambiguous {
                //     align_anew(aligner, target_seq, approx_target_start, approx_target_end, read_seq, read_clipping)
                // } else {
                //     haps_cigar.transfer_read_alignment(
                //         offset, &source_cigar, read_seq, target_seq, aligner, query_to_ref)
                // };
                // log::debug!("   -> {}", new_cigar);


                let (new_start, new_cigar) = haps_cigar.transfer_read_alignment(
                    offset, &source_cigar, read_seq, target_seq, aligner, query_to_ref);
                let (alt_start, alt_cigar) = align_anew(aligner, target_seq, approx_target_start, approx_target_end,
                    read_seq, read_clipping);



                total_transfers += 1;
                let cigar_ref_len = new_cigar.ref_len();
                const MIN_ALN_SIZE: u32 = 50;
                // Either too high edit distance or too short alignment.
                if cigar_ref_len < MIN_ALN_SIZE {
                    prelim_alignments.record_bad_alignment(read_end, target_contig_id, new_start, false);
                    // Break if produced too many fails for this specific alignment.
                    if n_fails.next().is_none() { break } else { continue };
                }

                let interval = Interval::new(Arc::clone(contigs), target_contig_id,
                    new_start, new_start + cigar_ref_len);
                let mut new_aln = Alignment::new(interval, new_cigar, source_strand, read_end);
                new_aln.calculate_edit_and_prob(contigs, err_prof);

                let alt_interval = Interval::new(Arc::clone(contigs), target_contig_id,
                    alt_start, alt_start + alt_cigar.ref_len());
                let mut alt_aln = Alignment::new(alt_interval, alt_cigar, source_strand, read_end);
                alt_aln.calculate_edit_and_prob(contigs, err_prof);

                println!("{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{:.4}\t{}",
                    read_data.name(), read_data.name_hash(), read_end,
                    self.ambiguities.min_unique_kmers(target_contig_id, approx_target_start, approx_target_end),
                    new_aln.interval(), new_aln.cigar(), new_aln.ln_prob(), new_aln.distance().unwrap().edit(),
                    alt_aln.interval(), alt_aln.cigar(), alt_aln.ln_prob(), alt_aln.distance().unwrap().edit(),
                );

                if let Err(opt_prev_ix) = prelim_alignments.insert(new_aln, false) {
                    if let Some(j) = opt_prev_ix && let Some(seen_val) = seen.get_mut(j as usize) {
                        *seen_val = true;
                    }
                    if n_fails.next().is_none() {
                        break;
                    }
                }
            }
        }
        (prelim_alignments.len() - n, total_transfers)
    }
}
