//! Functions related to the pangenome VCF file.

use htslib::bcf::{
    self, Read,
    header::HeaderView,
    record::{Record, GenotypeAllele},
};
use crate::{
    err::{Error, error},
    seq::Interval,
    ext::vec::F64Ext,
    algo::HashSet,
};
use super::NamedSeq;

/// Stores one sample in the VCF file.
pub struct HaplotypeName {
    /// Haplotype index for this sample.
    hap_ix: usize,
    /// Shift in index across all remaining haplotypes.
    shift_ix: usize,
    /// Haplotype name.
    name: String,
}

// impl HaplotypeName {
//     pub fn hap_ix(&self) -> usize {
//         self.hap_ix
//     }
// }

/// One sample.
pub struct Sample {
    /// Sample index.
    sample_id: usize,
    name: String,
    ploidy: usize,
    /// Sample haplotypes (some may be discarded).
    subnames: Vec<HaplotypeName>,
}

// impl Sample {
//     pub fn sample_id(&self) -> usize {
//         self.sample_id
//     }

//     pub fn subnames(&self) -> &[HaplotypeName] {
//         &self.subnames
//     }
// }

/// Collection of all haplotypes.
pub struct HaplotypeNames {
    /// None if reference is skipped.
    ref_name: Option<String>,
    /// All samples.
    samples: Vec<Sample>,
    /// Total number of haplotypes.
    total: usize,
}

impl HaplotypeNames {
    /// Examines VCF header and the first record to extract sample ploidy.
    /// Samples and haplotypes in `leave_out` are discarded.
    pub fn new(
        reader: &mut impl bcf::Read,
        ref_name: &str,
        leave_out: &HashSet<String>,
    ) -> crate::Result<Self>
    {
        let mut discarded = 0;
        let mut left_out = Vec::new();
        let mut total = 0;
        let mut hap_names = HashSet::default();
        let ref_name = if leave_out.contains(ref_name) {
            left_out.push(ref_name.to_owned());
            None
        } else {
            total += 1;
            hap_names.insert(ref_name.to_owned());
            Some(ref_name.to_owned())
        };

        let mut samples = Vec::new();
        let record = reader.records().next()
            .ok_or_else(|| error!(InvalidData, "Input VCF file does not contain any records"))??;
        let genotypes = record.genotypes()?;
        for (sample_id, sample) in reader.header().samples().into_iter().enumerate() {
            let sample = std::str::from_utf8(sample).map_err(|_| Error::Utf8("VCF sample name", sample.to_vec()))?;

            let ploidy = genotypes.get(sample_id).len();
            if leave_out.contains(sample) {
                left_out.push(format!("{} x{}", sample, ploidy));
                discarded += ploidy;
                continue;
            }
            if ploidy == 0 {
                return Err(error!(InvalidData, "Sample {} has zero ploidy", sample));
            } else if ploidy > 255 {
                return Err(error!(InvalidData, "Sample {} has extremely high ploidy", sample));
            }
            let mut subnames = Vec::with_capacity(ploidy);
            for hap_ix in 0..ploidy {
                let haplotype = if ploidy == 1 { sample.to_owned() } else { format!("{}.{}", sample, hap_ix + 1) };
                if leave_out.contains(&haplotype) {
                    left_out.push(haplotype);
                    discarded += 1;
                    continue;
                }
                if !hap_names.insert(haplotype.clone()) {
                    return Err(error!(InvalidData, "Duplicate haplotype name ({})", haplotype));
                }
                subnames.push(HaplotypeName {
                    hap_ix,
                    shift_ix: total,
                    name: haplotype,
                });
                total += 1;
            }
            samples.push(Sample {
                name: sample.to_owned(),
                sample_id, ploidy, subnames,
            });
        }
        log::info!("VCF file contains {} haplotypes", total);
        if discarded > 0 {
            if left_out.len() > 5 {
                left_out.truncate(5);
                left_out.push("...".to_owned());
            }
            log::warn!("    Leave out {} haplotypes ({})", discarded, left_out.join(", "));
        } else if !leave_out.is_empty() {
            log::warn!("Zero matches between leave-out and VCF samples");
        }
        if samples.is_empty() {
            return Err(error!(InvalidData, "Loaded zero haplotypes"));
        }
        Ok(Self { ref_name, samples, total })
    }

    // pub fn samples(&self) -> &[Sample] {
    //     &self.samples
    // }

    // /// Total number of haplotypes (including the reference).
    // pub fn total(&self) -> usize {
    //     self.total
    // }
}

/// Discard variants where there is no known variation.
/// NOTE: Can also filter out poorly known variants, as well as variants with bad quality.
pub fn filter_variants(
    reader: &mut impl bcf::Read,
    hap_names: &HaplotypeNames,
) -> crate::Result<Vec<Record>>
{
    let mut vars = Vec::new();
    for rec in reader.records() {
        let var = rec?;
        let gts = var.genotypes()?;
        let mut has_variation = false;
        for sample in hap_names.samples.iter() {
            let gt = gts.get(sample.sample_id);
            if gt.len() != sample.ploidy {
                return Err(error!(InvalidData, "Variant {} in sample {} has ploidy {} (expected {})",
                    format_var(&var, reader.header()), sample.name, gt.len(), sample.ploidy));
            }
            // Check only last allele in the genotype, as only last allele has phasing marking,
            // (see https://docs.rs/rust-htslib/latest/rust_htslib/bcf/record/struct.Genotypes.html).
            match gt[sample.ploidy - 1] {
                GenotypeAllele::Unphased(_) | GenotypeAllele::UnphasedMissing if sample.ploidy > 1 =>
                    return Err(error!(InvalidData, "Variant {} is unphased in sample {}",
                        format_var(&var, reader.header()), sample.name)),
                _ => {}
            }
            has_variation |= sample.subnames.iter().any(|haplotype|
                match gt[haplotype.hap_ix] {
                    GenotypeAllele::Phased(1..) | GenotypeAllele::Unphased(1..) => true,
                    _ => false,
                });
        }
        if has_variation {
            vars.push(var);
        }
    }
    Ok(vars)
}

/// Formats variant as `chrom:pos`.
fn format_var(var: &Record, header: &HeaderView) -> String {
    let chrom = var.rid()
        .and_then(|rid| header.rid2name(rid).ok())
        .and_then(|bytes| std::str::from_utf8(bytes).ok())
        .unwrap_or("??");
    format!("{}:{}", chrom, var.pos() + 1)
}

/// Discard sequences with too many unknown nucleotides.
fn discard_unknown(
    seqs: &mut Vec<NamedSeq>,
    unknown_nts: &[u32],
    unknown_frac: f64,
) {
    let aver_unknown = F64Ext::mean(unknown_nts);
    if aver_unknown > 0.0 {
        log::debug!("        On average, {:.1} bp unknown per haplotype", aver_unknown);
    }
    let keep_seqs: Vec<_> = seqs.iter().zip(unknown_nts)
        .map(|(seq, &unknown)| f64::from(unknown) <= unknown_frac * f64::from(seq.len()))
        .collect();
    let total = seqs.len();
    let n_remain: usize = keep_seqs.iter().copied().map(usize::from).sum();
    let n_discard = total - n_remain;
    if n_discard == 0 {
        return;
    }

    let unfilt_seqs = std::mem::replace(seqs, Vec::with_capacity(n_remain));
    seqs.extend(unfilt_seqs.into_iter().zip(&keep_seqs)
        .filter_map(|(seq, &keep)| if keep { Some(seq) } else { None }));
    log::warn!("        Reconstructed {} alleles ({} unavailable)", n_remain, n_discard);
}

/// Reconstructs sample sequences by adding variants to the reference sequence.
/// Returns a vector of named sequences.
pub fn reconstruct_sequences(
    interval: &Interval,
    ref_seq: &[u8],
    vcf_file: &mut bcf::IndexedReader,
    hap_names: &HaplotypeNames,
    unknown_frac: f64,
    overlaps_allowed: bool,
) -> crate::Result<Vec<NamedSeq>>
{
    assert_eq!(interval.len(), ref_seq.len() as u32);
    let (ref_start, ref_end) = interval.range();
    let vcf_rid = vcf_file.header().name2rid(interval.contig_name().as_bytes())?;
    vcf_file.fetch(vcf_rid, u64::from(ref_start), Some(u64::from(ref_end)))?;
    let recs = filter_variants(vcf_file, hap_names)?;

    let mut seqs = Vec::with_capacity(hap_names.total);
    if let Some(ref_name) = hap_names.ref_name.as_ref() {
        // Sequence will be extended at the end of the function, during suffix extension.
        seqs.push(NamedSeq::new(ref_name.to_owned(), Vec::with_capacity(ref_seq.len())));
    }
    for sample in hap_names.samples.iter() {
        for haplotype in sample.subnames.iter() {
            seqs.push(NamedSeq::new(haplotype.name.clone(), Vec::with_capacity(ref_seq.len() * 3 / 2)));
        }
    }
    // Number of unknown nucleotides for each sequence.
    let mut unknown_nts = vec![0_u32; hap_names.total];
    let header = vcf_file.header();
    let mut ref_pos = vec![ref_start; hap_names.total];
    let mut total_overlaps = 0;
    const MAX_OVERLAP_MSGS: u32 = 3;

    for var in recs.iter() {
        let alleles = var.alleles();
        // if alleles.iter().copied().any(seq::has_n) {
        //     return Err(error!(InvalidData, "Input VCF file contains Ns in one of the alleles of {}",
        //         format_var(var, header)));
        // }
        let var_start = u32::try_from(var.pos()).unwrap();
        let ref_len = alleles[0].len() as u32;
        let var_end = var_start + ref_len;
        if var_end <= ref_start {
            continue;
        } else if ref_end <= var_start {
            break;
        } else if var_start < ref_start || ref_end < var_end {
            return Err(error!(RuntimeError, "Variant {} overlaps the boundary of the region {}-{}",
                format_var(var, header), ref_start + 1, ref_end));
        }

        let gts = var.genotypes()?;
        for sample in hap_names.samples.iter() {
            let gt = gts.get(sample.sample_id);
            assert_eq!(gt.len(), sample.ploidy);
            for haplotype in sample.subnames.iter() {
                let allele_ix = match gt[haplotype.hap_ix] {
                    GenotypeAllele::Phased(allele_ix) | GenotypeAllele::Unphased(allele_ix) => allele_ix as usize,
                    GenotypeAllele::PhasedMissing | GenotypeAllele::UnphasedMissing =>
                    {
                        unknown_nts[haplotype.shift_ix] += ref_len;
                        0 // Use reference allele
                    }
                };
                // Simply do not update `ref_pos`, nothing more to do.
                if allele_ix == 0 {
                    continue;
                }

                let prev_end = ref_pos[haplotype.shift_ix];
                if var_start < prev_end {
                    if !overlaps_allowed {
                        return Err(error!(InvalidData,
                            "Overlapping variants forbidden ({} for {})", format_var(var, header), haplotype.name));
                    } else if total_overlaps < MAX_OVERLAP_MSGS {
                        log::warn!("One of the overlapping variants ignored for {} ({})",
                            haplotype.name, format_var(var, header));
                    }
                    total_overlaps += 1;
                    continue;
                }
                let mut_seq = seqs[haplotype.shift_ix].seq_mut();
                mut_seq.extend_from_slice(&ref_seq[(prev_end - ref_start) as usize..(var_start - ref_start) as usize]);
                mut_seq.extend_from_slice(alleles[allele_ix]);
                ref_pos[haplotype.shift_ix] = var_end;
            }
        }
    }
    for (&prev_end, entry) in ref_pos.iter().zip(seqs.iter_mut()) {
        if prev_end < ref_end {
            entry.seq_mut().extend_from_slice(&ref_seq[(prev_end - ref_start) as usize..]);
        }
    }
    if total_overlaps > 0 {
        log::warn!("In total, {} overlapping variants discarded", total_overlaps);
    }

    discard_unknown(&mut seqs, &unknown_nts, unknown_frac);
    Ok(seqs)
}
