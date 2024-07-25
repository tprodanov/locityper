//! Output BAM file with read alignments.

use std::{
    path::Path,
    rc::Rc,
};
use htslib::bam;
use crate::{
    solvers::scheme,
    seq::{
        contigs::{Genotype, ContigNames, ContigId},
        cigar::CigarItem,
        aln::{Alignment, Strand, ReadEnd},
    },
    model::{
        locs::{GrouppedAlignments, ReadData, PairAlignment},
        windows::ReadGtAlns,
        assgn::GenotypeAlignments,
    },
    algo::{TwoU32, IntMap},
};

/// Updates `contig_to_tid` vector, and returns
/// - BAM header,
/// - list of unique contig ids.
fn create_bam_header(
    gt: &Genotype,
    contigs: &ContigNames,
    contig_to_tid: &mut Vec<Option<i32>>,
) -> (bam::header::Header, Vec<ContigId>)
{
    let mut unique_ids = Vec::new();
    contig_to_tid.clear();
    contig_to_tid.resize(contigs.len(), None);
    let mut header = bam::header::Header::new();
    let mut tid = 0;
    for &id in gt.ids().iter() {
        if contig_to_tid[id.ix()].is_some() {
            // This contig was already included in the genotype.
            continue;
        }
        unique_ids.push(id);
        contig_to_tid[id.ix()] = Some(tid);
        tid += 1;
        let mut record = bam::header::HeaderRecord::new(b"SQ");
        record.push_tag(b"SN", contigs.get_name(id));
        record.push_tag(b"LN", contigs.get_len(id));
        header.push_record(&record);
    }
    (header, unique_ids)
}

/// Based on the number of read assignments and the total number of attempts,
/// returns probability and mapping quality.
fn count_to_prob(count: u16, attempts: u16) -> (f32, u8) {
    if count == 0 {
        (0.0, 0)
    } else if count == attempts {
        (1.0, 60)
    } else {
        assert!(count < attempts);
        let prob = f32::from(count) / f32::from(attempts);
        let mapq = (-10.0 * (1.0 - prob).log10()).round().min(60.0) as u8;
        (prob, mapq)
    }
}

/// Returns insert size, which will be set to record1. Record2 will be assigned minus this insert size.
fn calc_insert_size(opt_aln1: Option<&Alignment>, opt_aln2: Option<&Alignment>) -> i64 {
    match (opt_aln1, opt_aln2) {
        (Some(aln1), Some(aln2)) => {
            let (start1, end1) = aln1.interval().range();
            let (start2, end2) = aln2.interval().range();
            if start1 <= start2 {
                i64::from(end2) - i64::from(start1)
            } else {
                // This value will be negative (because first mate is on the right)
                i64::from(start2) - i64::from(end1)
            }
        }
        _ => 0,
    }
}

/// Creates mapped/unmapped record, sets position, sequence and qualities.
fn create_record(
    header: Rc<bam::HeaderView>,
    contig_to_tid: &[Option<i32>],
    opt_aln: Option<&Alignment>,
    read_data: &ReadData,
    read_end: ReadEnd,
    cigar_buffer: &mut bam::record::CigarString,
) -> htslib::errors::Result<bam::Record>
{
    let mut record = bam::Record::new();
    record.set_header(header);
    let (cigar_view, strand) = if let Some(aln) = opt_aln {
        record.set_tid(contig_to_tid[aln.contig_id().ix()].expect("Contig ID undefined"));
        record.set_pos(i64::from(aln.interval().start()));
        let strand = aln.strand();
        if strand == Strand::Reverse {
            record.set_reverse();
        }
        cigar_buffer.0.clear();
        cigar_buffer.0.extend(aln.cigar().iter().copied().map(CigarItem::to_htslib));
        (Some(cigar_buffer as &bam::record::CigarString), strand)
    } else {
        record.set_unmapped();
        record.set_tid(-1);
        record.set_pos(0);
        (None, Strand::Forward)
    };

    let mate_data = read_data.mate_data(read_end);
    let (seq, qual) = mate_data.get_seq_and_qual(strand);
    record.set(read_data.name().as_bytes(), cigar_view, seq, qual);

    if let Some(aln) = opt_aln {
        if let Some(dist) = aln.distance() {
            record.push_aux(EDIT_TAG, bam::record::Aux::U32(dist.edit()))?;
        }
        record.push_aux(ALN_LIK_TAG, bam::record::Aux::Double(aln.ln_prob()))?;
    }
    record.push_aux(UNIQ_KMERS_TAG, bam::record::Aux::U16(mate_data.unique_kmers()))?;
    Ok(record)
}

const UNMAPPED: u32 = u32::MAX;

const EDIT_TAG: &'static [u8; 2] = b"NM";
const PROB_TAG: &'static [u8; 2] = b"pr";
const USED_TAG: &'static [u8; 2] = b"us";
const UNIQ_KMERS_TAG: &'static [u8; 2] = b"uk";
const ALN_LIK_TAG: &'static [u8; 2] = b"al";

/// For each alignment pair, sums corresponding `assign_counts`.
/// Alignment pair may appear several times if the same contig appears multiple times in the genotype.
fn count_alignments<const PAIRED: bool>(
    read_alns: &[ReadGtAlns],
    assgn_counts: &mut impl Iterator<Item = u16>,
    buffer1: &mut IntMap<TwoU32, u16>,
    buffer2: &mut Vec<(TwoU32, u16)>,
) {
    buffer1.clear();
    for aln in read_alns {
        let (i, j) = match aln.parent() {
            None => (UNMAPPED, UNMAPPED),
            Some((_, aln_pair)) => (
                aln_pair.ix1().unwrap_or(UNMAPPED),
                if PAIRED { aln_pair.ix2().unwrap_or(UNMAPPED) } else { UNMAPPED },
            ),
        };
        let count = assgn_counts.next().expect("Not enough assignment counts");
        if count > 0 {
            *buffer1.entry(TwoU32(i, j)).or_default() += count;
        }
    }
    buffer2.clear();
    buffer2.extend(buffer1.iter().map(|(&u, &v)| (u, v)));
    buffer2.sort_unstable_by(|a, b| b.1.cmp(&a.1));
}

/// Connect two mates: sets corresponding flags and mate positions.
fn connect_pair(record1: &mut bam::Record, record2: &mut bam::Record, insert_size: i64) {
    match (!record1.is_unmapped(), !record2.is_unmapped()) {
        (false, false) => {
            record1.set_mate_unmapped();
            record2.set_mate_unmapped();
        }
        (false, true) => {
            record2.set_mate_unmapped();
            // Set the same position as the second mate.
            record1.set_tid(record2.tid());
            record1.set_pos(record2.pos());
        }
        (true, false) => {
            record1.set_mate_unmapped();
            // Set the same position as the first mate.
            record2.set_tid(record1.tid());
            record2.set_pos(record1.pos());
        }
        (true, true) => {
            record1.set_proper_pair();
            record2.set_proper_pair();
            if record1.is_reverse() {
                record2.set_mate_reverse();
            }
            if record2.is_reverse() {
                record1.set_mate_reverse();
            }
        }
    }
    record1.set_paired();
    record1.set_first_in_template();
    record1.set_mtid(record2.tid());
    record1.set_mpos(record2.pos());
    record1.set_insert_size(insert_size);

    record2.set_paired();
    record2.set_last_in_template();
    record2.set_mtid(record1.tid());
    record2.set_mpos(record1.pos());
    record2.set_insert_size(-insert_size);
}

/// Generate BAM record for one paired-end read.
fn generate_paired_end_records(
    groupped_alns: &GrouppedAlignments,
    read_alns: &[ReadGtAlns],
    header: &Rc<bam::HeaderView>,
    contig_to_tid: &[Option<i32>],
    records: &mut Vec<bam::Record>,
    assgn_counts: &mut impl Iterator<Item = u16>,
    attempts: u16,
    buffer1: &mut IntMap<TwoU32, u16>,
    buffer2: &mut Vec<(TwoU32, u16)>,
    buffer3: &mut bam::record::CigarString,
) -> htslib::errors::Result<()>
{
    count_alignments::<true>(read_alns, assgn_counts, buffer1, buffer2);
    let read_data = groupped_alns.read_data();
    let mut secondary = false;
    for &(TwoU32(i, j), count) in buffer2.iter() {
        let opt_aln1 = if i == UNMAPPED { None } else { Some(groupped_alns.ith_aln(i)) };
        let opt_aln2 = if j == UNMAPPED { None } else { Some(groupped_alns.ith_aln(j)) };
        let mut record1 = create_record(Rc::clone(header), contig_to_tid, opt_aln1, read_data,
            ReadEnd::First, buffer3)?;
        let mut record2 = create_record(Rc::clone(header), contig_to_tid, opt_aln2, read_data,
            ReadEnd::Second, buffer3)?;
        let (prob, mapq) = count_to_prob(count, attempts);
        record1.set_mapq(mapq);
        record2.set_mapq(mapq);
        record1.push_aux(PROB_TAG, bam::record::Aux::Float(prob))?;
        record2.push_aux(PROB_TAG, bam::record::Aux::Float(prob))?;
        record1.push_aux(USED_TAG, bam::record::Aux::Char(b'T'))?;
        record2.push_aux(USED_TAG, bam::record::Aux::Char(b'T'))?;
        if secondary {
            record1.set_secondary();
            record2.set_secondary();
        }
        secondary = true;
        connect_pair(&mut record1, &mut record2, calc_insert_size(opt_aln1, opt_aln2));
        records.push(record1);
        records.push(record2);
    }
    Ok(())
}

/// Generate BAM record for one paired-end read that was not used in the read assignment.
fn generate_unused_paired_end_records(
    groupped_alns: &GrouppedAlignments,
    aln_pairs: &[&PairAlignment],
    header: &Rc<bam::HeaderView>,
    contig_to_tid: &[Option<i32>],
    records: &mut Vec<bam::Record>,
    buffer: &mut bam::record::CigarString,
) -> htslib::errors::Result<()>
{
    let read_data = groupped_alns.read_data();
    let mut secondary = false;
    for pair in aln_pairs {
        let opt_aln1 = pair.ix1().map(|i| groupped_alns.ith_aln(i));
        let opt_aln2 = pair.ix2().map(|j| groupped_alns.ith_aln(j));
        let mut record1 = create_record(Rc::clone(header), contig_to_tid, opt_aln1, read_data, ReadEnd::First, buffer)?;
        let mut record2 = create_record(Rc::clone(header), contig_to_tid, opt_aln2, read_data, ReadEnd::Second, buffer)?;
        if secondary {
            record1.set_secondary();
            record2.set_secondary();
        }
        secondary = true;
        connect_pair(&mut record1, &mut record2, calc_insert_size(opt_aln1, opt_aln2));
        record1.push_aux(USED_TAG, bam::record::Aux::Char(b'F'))?;
        record2.push_aux(USED_TAG, bam::record::Aux::Char(b'F'))?;
        records.push(record1);
        records.push(record2);
    }
    Ok(())
}

/// Generate BAM record for one single-end read.
fn generate_single_end_records(
    groupped_alns: &GrouppedAlignments,
    read_alns: &[ReadGtAlns],
    header: &Rc<bam::HeaderView>,
    contig_to_tid: &[Option<i32>],
    records: &mut Vec<bam::Record>,
    assgn_counts: &mut impl Iterator<Item = u16>,
    attempts: u16,
    buffer1: &mut IntMap<TwoU32, u16>,
    buffer2: &mut Vec<(TwoU32, u16)>,
    buffer3: &mut bam::record::CigarString,
) -> htslib::errors::Result<()>
{
    count_alignments::<false>(read_alns, assgn_counts, buffer1, buffer2);
    let read_data = groupped_alns.read_data();
    let mut secondary = false;
    for &(TwoU32(i, _), count) in buffer2.iter() {
        let opt_aln = if i == UNMAPPED { None } else { Some(groupped_alns.ith_aln(i)) };
        let mut record = create_record(Rc::clone(header), contig_to_tid, opt_aln, read_data, ReadEnd::First, buffer3)?;
        let (prob, mapq) = count_to_prob(count, attempts);
        record.set_mapq(mapq);
        record.push_aux(PROB_TAG, bam::record::Aux::Float(prob))?;
        record.push_aux(USED_TAG, bam::record::Aux::Char(b'T'))?;
        if secondary {
            record.set_secondary();
        }
        secondary = true;
        records.push(record);
    }
    Ok(())
}

/// Generate BAM record for one single-end read that was not used in the read assignment.
fn generate_unused_single_end_records(
    groupped_alns: &GrouppedAlignments,
    aln_pairs: &[&PairAlignment],
    header: &Rc<bam::HeaderView>,
    contig_to_tid: &[Option<i32>],
    records: &mut Vec<bam::Record>,
    buffer: &mut bam::record::CigarString,
) -> htslib::errors::Result<()>
{
    let read_data = groupped_alns.read_data();
    let mut secondary = false;
    for pair in aln_pairs {
        let opt_aln = pair.ix1().map(|i| groupped_alns.ith_aln(i));
        let mut record = create_record(Rc::clone(header), contig_to_tid, opt_aln, read_data, ReadEnd::First, buffer)?;
        if secondary {
            record.set_secondary();
        }
        secondary = true;
        record.push_aux(USED_TAG, bam::record::Aux::Char(b'F'))?;
        records.push(record);
    }
    Ok(())
}

/// Writes read alignments to a specific genotype.
/// `contig_to_tid`: converts contig id into BAM tid. Is updated automatically for each genotype.
pub fn write_bam(
    bam_path: &Path,
    gt: &Genotype,
    data: &scheme::Data,
    contig_to_tid: &mut Vec<Option<i32>>,
    assgn_counts: &[u16],
) -> htslib::errors::Result<()>
{
    let (header, unique_ids) = create_bam_header(gt, &data.contigs, contig_to_tid);
    let header_view = Rc::new(bam::HeaderView::from_header(&header));
    let mut records = Vec::new();

    let gt_alns = GenotypeAlignments::new(gt.clone(), &data.all_contig_infos, &data.all_alns, &data.assgn_params);
    let mut counts_iter = assgn_counts.iter().copied();
    let mut buffer1 = Default::default();
    let mut buffer2 = Default::default();
    let mut buffer3 = bam::record::CigarString(Vec::new());
    for (rp, groupped_alns) in data.all_alns.reads().iter().enumerate() {
        let read_alns = gt_alns.possible_read_alns(rp);
        if data.is_paired_end {
            generate_paired_end_records(groupped_alns, read_alns, &header_view, contig_to_tid, &mut records,
                &mut counts_iter, data.assgn_params.attempts, &mut buffer1, &mut buffer2, &mut buffer3)?;
        } else {
            generate_single_end_records(groupped_alns, read_alns, &header_view, contig_to_tid, &mut records,
                &mut counts_iter, data.assgn_params.attempts, &mut buffer1, &mut buffer2, &mut buffer3)?;
        }
    }
    assert!(counts_iter.next().is_none(), "Too many assignment counts");

    let mut aln_pairs: Vec<&PairAlignment> = Vec::new();
    // Reads unused in the analysis.
    for groupped_alns in data.all_alns.unused_reads() {
        aln_pairs.clear();
        for &id in unique_ids.iter() {
            aln_pairs.extend(groupped_alns.contig_alns(id).iter());
        }
        aln_pairs.sort_unstable_by(|a, b| b.ln_prob().total_cmp(&a.ln_prob()));
        if data.is_paired_end {
            generate_unused_paired_end_records(groupped_alns, &aln_pairs, &header_view, contig_to_tid, &mut records,
                &mut buffer3)?;
        } else {
            generate_unused_single_end_records(groupped_alns, &aln_pairs, &header_view, contig_to_tid, &mut records,
                &mut buffer3)?;
        }
    }

    // Stable sort so that we do not separate mates.
    // Convert tid to u32 so that unmapped reads (-1) appear at the end.
    records.sort_by(|rec1, rec2| (rec1.tid() as u32, rec1.pos()).cmp(&(rec2.tid() as u32, rec2.pos())));
    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam)?;
    for record in records.into_iter() {
        writer.write(&record)?;
    }
    std::mem::drop(writer);
    // Use 1 thread.
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1)?;
    Ok(())
}
