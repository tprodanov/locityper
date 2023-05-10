//! Functions related to the pangenome VCF file.

use std::str::from_utf8;
use fnv::FnvHashSet;
use htslib::bcf::{
    header::HeaderView,
    record::{Record, GenotypeAllele},
};
use crate::Error;
use super::NamedSeq;

/// Formats variant as `chrom:pos`.
fn format_var(var: &Record, header: &HeaderView) -> String {
    let chrom = var.rid()
        .and_then(|rid| header.rid2name(rid).ok())
        .and_then(|bytes| from_utf8(bytes).ok())
        .unwrap_or("??");
    format!("{}:{}", chrom, var.pos() + 1)
}

/// Reconstructs sample sequences by adding variants to the reference sequence.
/// Returns Vector of pairs `(name, sequence)`.
pub fn reconstruct_sequences(
    ref_start: u32,
    ref_seq: &[u8],
    ref_name: &str,
    header: &HeaderView,
    recs: &[Record],
    leave_out: &FnvHashSet<String>,
) -> Result<Vec<NamedSeq>, Error>
{
    const PLOIDY: usize = 2;
    let ref_end = ref_start + ref_seq.len() as u32;
    let n_samples = header.sample_count() as usize;
    let mut haplotype_names = FnvHashSet::with_capacity_and_hasher(1 + n_samples * PLOIDY, Default::default());
    haplotype_names.insert(ref_name.to_owned());

    let mut seqs = vec![NamedSeq::new(ref_name.to_owned(), ref_seq.to_vec())];
    // Boolean vector, do we keep the sequence according to `leave_out`?
    let mut keep_seq = vec![leave_out.is_empty() || !leave_out.contains(ref_name)];

    let capacity = ref_seq.len() * 3 / 2;
    let sample_names = header.samples();
    for sample in sample_names.iter() {
        let sample = from_utf8(sample)?;
        for i in 1..=PLOIDY {
            let name = format!("{}.{}", sample, i);
            if !haplotype_names.insert(name.clone()) {
                return Err(Error::InvalidData(format!("Haplotype name {} is not unique", name)));
            }
            keep_seq.push(leave_out.is_empty() || (!leave_out.contains(sample) && !leave_out.contains(&name)));
            seqs.push(NamedSeq::new(name, Vec::with_capacity(capacity)));
        }
    }

    let mut ref_pos = ref_start;
    for var in recs.iter() {
        let alleles = var.alleles();
        let var_start = u32::try_from(var.pos()).unwrap();
        let var_end = var_start + alleles[0].len() as u32;
        if var_start < ref_pos {
            if var_end > ref_pos {
                return Err(Error::InvalidData(format!("Input VCF file contains overlapping variants: see {}",
                    format_var(var, header))));
            }
            continue;
        }
        if ref_end <= var_start {
            break;
        }
        let seq_between_vars = &ref_seq[(ref_pos - ref_start) as usize..(var_start - ref_start) as usize];

        let gts = var.genotypes()?;
        for sample_id in 0..n_samples {
            let gt = gts.get(sample_id);
            if gt.len() != PLOIDY {
                return Err(Error::InvalidData(format!("Variant {} in sample {} has ploidy {} (expected {})",
                    format_var(var, header), from_utf8(&sample_names[sample_id]).unwrap(),
                    gt.len(), PLOIDY)));
            }
            // Check only last allele in the genotype, as only last allele has phasing marking,
            // (see https://docs.rs/rust-htslib/latest/rust_htslib/bcf/record/struct.Genotypes.html).
            match gt[PLOIDY - 1] {
                GenotypeAllele::Unphased(_) => return Err(Error::InvalidData(format!(
                    "Variant {} is unphased in sample {}", format_var(var, header),
                    from_utf8(&sample_names[sample_id]).unwrap()))),
                GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing => return Err(
                    Error::InvalidData(format!(
                    "Variant {} is missing genotype for sample {}", format_var(var, header),
                    from_utf8(&sample_names[sample_id]).unwrap()))),
                _ => {}
            }
            for haplotype in 0..PLOIDY {
                let mut_seq = seqs[PLOIDY * sample_id + haplotype + 1].seq_mut();
                mut_seq.extend_from_slice(seq_between_vars);
                match gt[haplotype] {
                    GenotypeAllele::Phased(allele_ix) | GenotypeAllele::Unphased(allele_ix) =>
                        mut_seq.extend_from_slice(alleles[allele_ix as usize]),
                    _ => return Err(Error::InvalidData(format!(
                        "Variant {} is missing genotype for sample {}", format_var(var, header),
                        from_utf8(&sample_names[sample_id]).unwrap()))),
                }
            }
        }
        ref_pos = var_end;
    }
    let suffix_seq = &ref_seq[(ref_pos - ref_start) as usize..];
    for entry in seqs[1..].iter_mut() {
        entry.seq_mut().extend_from_slice(suffix_seq);
    }

    let n_keep = keep_seq.iter().fold(0, |acc, &keep| acc + usize::from(keep));
    if n_keep < seqs.len() {
        log::warn!("        Leave out {} sequences", seqs.len() - n_keep);
        Ok(
            keep_seq.into_iter().zip(seqs.into_iter())
                .filter_map(|(keep, seq)| if keep { Some(seq) } else { None })
                .collect()
        )
    } else {
        Ok(seqs)
    }
}
