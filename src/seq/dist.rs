use std::{
    cmp::min,
    io::{self, Write},
};
use once_cell::sync::OnceCell;
use nohash::IntMap;
use bio::alignment::sparse::lcskpp;
use smallvec::SmallVec;
use super::{
    contigs::{ContigId, ContigNames},
    wfa::{Aligner, Penalties},
    cigar::{Cigar, Operation},
};

/// Calculates all pairwise distances between all sequences, writes alignments to PAF file,
/// and returns condensed distance matrix (see `kodama` crate).
pub fn calculate_all_distances(
    names_seqs: &[(String, Vec<u8>)],
    mut paf_writer: impl Write,
    penalties: &Penalties,
    threads: u16,
) -> io::Result<Vec<f64>> {
    if threads == 1 {
        let n = names_seqs.len();
        let mut distances = Vec::with_capacity(n * (n - 1) / 2);
        let aligner = Aligner::new(penalties);
        for (i, (name1, seq1)) in names_seqs.iter().enumerate() {
            for (name2, seq2) in names_seqs[i + 1..].iter() {
                let (cigar, score) = aligner.align(seq1, seq2);
                assert_eq!(seq1.len() as u32, cigar.ref_len());
                assert_eq!(seq2.len() as u32, cigar.query_len());
                let divergence = write_paf(&mut paf_writer, name2, name1, score, &cigar)?;
                distances.push(divergence);
            }
        }
        Ok(distances)
    } else {
        unimplemented!("Multi-thread seq. distance cannot be calculated yet.")
    }
}

/// Writes the alignment to PAF file and returns sequence divergence.
fn write_paf(
    writer: &mut impl Write,
    qname: &str,
    rname: &str,
    score: i32,
    cigar: &Cigar,
) -> io::Result<f64>
{
    let qlen = cigar.query_len();
    let rlen = cigar.ref_len();
    write!(writer, "{qname}\t{qlen}\t0\t{qlen}\t+\t{rname}\t{rlen}\t0\t{rlen}\t")?;

    let mut nmatches = 0;
    let mut total_size = 0;
    for item in cigar.iter() {
        total_size += item.len();
        if item.operation() == Operation::Equal {
            nmatches += item.len();
        }
    }
    let edit_dist = total_size - nmatches;
    let diverg = f64::from(edit_dist) / f64::from(total_size);
    writeln!(writer, "{nmatches}\t{total_size}\t60\tNM:i:{edit_dist}\tAS:i:{score}\tdv:f:{diverg:.7}\tcg:Z:{cigar}")?;
    Ok(diverg)
}
