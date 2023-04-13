//! Functions related to the pangenome VCF file.

use std::str::from_utf8;
use fnv::FnvHashSet;
use htslib::bcf::{
    header::HeaderView,
    record::{Record, GenotypeAllele},
};
use crate::Error;

/// Formats variant as `chrom:pos`.
fn format_var(var: &Record, header: &HeaderView) -> String {
    let chrom = var.rid()
        .and_then(|rid| header.rid2name(rid).ok())
        .and_then(|bytes| from_utf8(bytes).ok())
        .unwrap_or("??");
    format!("{}:{}", chrom, var.pos() + 1)
}

/// Reconstructs sample sequences by adding variants to the reference sequence.
/// Returns Vector of pairs `(name, sequence)`. Reference sequence is added if `ref_name` is not None.
/// All names must be unique.
pub fn reconstruct_sequences(
    ref_start: u32, ref_seq: &[u8], ref_name: &Option<String>,
    header: &HeaderView, recs: &[Record],
) -> Result<Vec<(String, Vec<u8>)>, Error>
{
    const PLOIDY: usize = 2;
    let ref_end = ref_start + ref_seq.len() as u32;
    let n_samples = header.sample_count() as usize;
    let mut haplotype_names = FnvHashSet::with_capacity_and_hasher(1 + n_samples * PLOIDY, Default::default());
    let mut res = Vec::new();
    if let Some(name) = ref_name.as_ref() {
        haplotype_names.insert(name.clone());
        res.push((name.clone(), ref_seq.to_vec()));
    }
    let capacity = ref_seq.len() * 3 / 2;
    let sample_names = header.samples();
    for sample in sample_names.iter() {
        let sample = from_utf8(sample)?;
        for i in 1..=PLOIDY {
            let name = format!("{}.{}", sample, i);
            if !haplotype_names.insert(name.clone()) {
                return Err(Error::InvalidData(format!("Haplotype name {} is not unique", name)));
            }
            res.push((name, Vec::with_capacity(capacity)));
        }
    }

    let shift = ref_name.is_some() as usize;
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
                return Err(Error::InvalidData(format!("")))
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
                let curr_pair = &mut res[shift + PLOIDY * sample_id + haplotype];
                curr_pair.1.extend(seq_between_vars);
                match gt[haplotype] {
                    GenotypeAllele::Phased(allele_ix) | GenotypeAllele::Unphased(allele_ix) =>
                        curr_pair.1.extend(alleles[allele_ix as usize]),
                    _ => return Err(Error::InvalidData(format!(
                        "Variant {} is missing genotype for sample {}", format_var(var, header),
                        from_utf8(&sample_names[sample_id]).unwrap()))),
                }
            }
        }
        ref_pos = var_end;
    }
    let suffix_seq = &ref_seq[(ref_pos - ref_start) as usize..];
    for i in shift..res.len() {
        res[i].1.extend(suffix_seq);
    }
    Ok(res)
}