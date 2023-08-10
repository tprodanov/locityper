//! Functions related to the pangenome VCF file.

use std::{
    str::from_utf8,
    cmp::Ordering,
};
use fnv::FnvHashSet;
use htslib::bcf::{
    self,
    header::HeaderView,
    record::{Record, GenotypeAllele},
};
use crate::{
    Error,
    seq::{
        self,
        cigar::Cigar,
        wfa::{Aligner, Penalties},
    },
    ext::vec::F64Ext,
};
use super::NamedSeq;

/// Stores one sample in the VCF file.
pub struct Haplotype {
    /// Haplotype index for this sample.
    hap_ix: usize,
    /// Shift in index across all remaining haplotypes.
    shift_ix: usize,
    /// Haplotype name.
    name: String,
}

impl Haplotype {
    pub fn hap_ix(&self) -> usize {
        self.hap_ix
    }
}

/// One sample.
pub struct Sample {
    /// Sample index.
    sample_id: usize,
    name: String,
    ploidy: usize,
    /// Sample haplotypes (some may be discarded).
    haplotypes: Vec<Haplotype>,
}

impl Sample {
    pub fn sample_id(&self) -> usize {
        self.sample_id
    }

    pub fn haplotypes(&self) -> &[Haplotype] {
        &self.haplotypes
    }
}

/// Collection of all haplotypes.
pub struct AllHaplotypes {
    /// None if reference is skipped.
    ref_name: Option<String>,
    /// All samples.
    samples: Vec<Sample>,
    /// Total number of haplotypes.
    total: usize,
}

impl AllHaplotypes {
    pub fn samples(&self) -> &[Sample] {
        &self.samples
    }
}

impl AllHaplotypes {
    /// Examines VCF header and the first record to extract sample ploidy.
    /// Samples and haplotypes in `leave_out` are discarded.
    pub fn new(
        reader: &mut impl bcf::Read,
        ref_name: &str,
        leave_out: &FnvHashSet<String>,
    ) -> Result<Self, Error>
    {
        let mut discarded = 0;
        let mut total = 0;
        let mut hap_names = FnvHashSet::default();
        let ref_name = if leave_out.contains(ref_name) {
            discarded += 1;
            None
        } else {
            total += 1;
            hap_names.insert(ref_name.to_owned());
            Some(ref_name.to_owned())
        };

        let mut samples = Vec::new();
        let record = reader.records().next()
            .ok_or_else(|| Error::InvalidData("Input VCF file does not contain any records".to_owned()))??;
        let genotypes = record.genotypes()?;
        for (sample_id, sample) in reader.header().samples().into_iter().enumerate() {
            let sample = std::str::from_utf8(sample).map_err(|_|
                Error::InvalidData(format!("Input VCF file contains non-UTF8 sample names ({:?})",
                String::from_utf8_lossy(sample))))?;

            let ploidy = genotypes.get(sample_id).len();
            if leave_out.contains(sample) {
                discarded += ploidy;
                continue;
            }
            if ploidy == 0 {
                return Err(Error::InvalidData(format!("Sample {} has zero ploidy", sample)));
            } else if ploidy > 255 {
                return Err(Error::InvalidData(format!("Sample {} has extremely high ploidy", sample)));
            }
            let mut haplotypes = Vec::with_capacity(ploidy);
            for hap_ix in 0..ploidy {
                let haplotype = if ploidy == 1 { sample.to_owned() } else { format!("{}.{}", sample, hap_ix + 1) };
                if leave_out.contains(&haplotype) {
                    discarded += 1;
                    continue;
                }
                if !hap_names.insert(haplotype.clone()) {
                    return Err(Error::InvalidData(format!("Duplicate haplotype name ({})", haplotype)));
                }
                haplotypes.push(Haplotype {
                    hap_ix,
                    shift_ix: total,
                    name: haplotype,
                });
                total += 1;
            }
            samples.push(Sample {
                name: sample.to_owned(),
                sample_id, ploidy, haplotypes,
            });
        }
        log::info!("Total {} haplotypes", total);
        if discarded > 0 {
            log::warn!("Leave out {} haplotypes", discarded);
        }
        if samples.is_empty() {
            return Err(Error::InvalidData("Loaded zero haplotypes".to_owned()));
        }
        Ok(Self { ref_name, samples, total })
    }
}

/// Discard variants where there is no known variation.
/// NOTE: Can also filter out poorly known variants, as well as variants with bad quality.
pub fn filter_variants(
    reader: &mut impl bcf::Read,
    haplotypes: &AllHaplotypes,
) -> Result<Vec<Record>, Error>
{
    let mut vars = Vec::new();
    for rec in reader.records() {
        let var = rec?;
        let gts = var.genotypes()?;
        let mut has_variation = false;
        for sample in haplotypes.samples.iter() {
            let gt = gts.get(sample.sample_id);
            if gt.len() != sample.ploidy {
                return Err(Error::InvalidData(format!("Variant {} in sample {} has ploidy {} (expected {})",
                    format_var(&var, reader.header()), sample.name, gt.len(), sample.ploidy)));
            }
            // Check only last allele in the genotype, as only last allele has phasing marking,
            // (see https://docs.rs/rust-htslib/latest/rust_htslib/bcf/record/struct.Genotypes.html).
            match gt[sample.ploidy - 1] {
                GenotypeAllele::Unphased(_) | GenotypeAllele::UnphasedMissing if sample.ploidy > 1 =>
                    return Err(Error::InvalidData(format!("Variant {} is unphased in sample {}",
                        format_var(&var, reader.header()), sample.name))),
                _ => {}
            }
            has_variation |= sample.haplotypes.iter().any(|haplotype|
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
        .and_then(|bytes| from_utf8(bytes).ok())
        .unwrap_or("??");
    format!("{}:{}", chrom, var.pos() + 1)
}

/// Reconstructs sample sequences by adding variants to the reference sequence.
/// Returns
/// - vector of named sequences,
/// - all pairwise alignments in a vector form.
pub fn reconstruct_sequences(
    ref_start: u32,
    ref_seq: &[u8],
    recs: &[Record],
    haplotypes: &AllHaplotypes,
    header: &HeaderView,
    unavail_rate: f64,
    aln_penalties: &Penalties,
) -> Result<(Vec<NamedSeq>, Vec<(Cigar, i32)>), Error>
{
    let ref_end = ref_start + ref_seq.len() as u32;
    let capacity = ref_seq.len() * 3 / 2;

    let has_reference = haplotypes.ref_name.is_some();
    let mut seqs = Vec::with_capacity(haplotypes.total);
    if let Some(ref_name) = haplotypes.ref_name.as_ref() {
        seqs.push(NamedSeq::new(ref_name.to_owned(), ref_seq.to_vec()));
    }
    for sample in haplotypes.samples.iter() {
        for haplotype in sample.haplotypes.iter() {
            seqs.push(NamedSeq::new(haplotype.name.clone(), Vec::with_capacity(capacity)));
        }
    }
    let aligner = Aligner::new(aln_penalties);
    let mism_penalty = -(aln_penalties.mismatch as i32);
    let mut pairwise_alns = vec![(Cigar::new(), 0); haplotypes.total * (haplotypes.total - 1) / 2];
    // Number of unavailable nucleotides for each sequence.
    let mut unavail_size = vec![0_u32; haplotypes.total];

    // Store last allele for each sample in this buffer.
    let mut last_allele = vec![0; haplotypes.total];
    let mut ref_pos = ref_start;
    for var in recs.iter() {
        let alleles = var.alleles();
        if alleles.iter().copied().any(seq::has_n) {
            return Err(Error::InvalidData(format!("Input VCF file contains Ns in one of the alleles of {}",
                format_var(var, header))));
        }
        let var_start = u32::try_from(var.pos()).unwrap();
        let var_end = var_start + alleles[0].len() as u32;
        if var_end <= ref_start {
            continue;
        }
        if ref_end <= var_start {
            break;
        }
        assert!(ref_start <= var_start && var_end <= ref_end,
            "Variant {} overlaps the boundary of the region {}-{}", format_var(var, header), ref_start + 1, ref_end);
        if var_start < ref_pos {
            return Err(Error::InvalidData(format!("Input VCF file contains overlapping variants: see {}. \
                Consider running `vcfbub -l 0 -i VCF | bgzip > VCF2 && tabix -p vcf VCF2`.",
                format_var(var, header))));
        }
        let seq_between_vars = &ref_seq[(ref_pos - ref_start) as usize..(var_start - ref_start) as usize];

        let gts = var.genotypes()?;
        let mut last_allele_iter = last_allele[usize::from(has_reference)..].iter_mut();
        for sample in haplotypes.samples.iter() {
            let gt = gts.get(sample.sample_id);
            assert_eq!(gt.len(), sample.ploidy);
            for haplotype in sample.haplotypes.iter() {
                let mut_seq = seqs[haplotype.shift_ix].seq_mut();
                mut_seq.extend_from_slice(seq_between_vars);
                let allele_ix = match gt[haplotype.hap_ix] {
                    GenotypeAllele::Phased(allele_ix) | GenotypeAllele::Unphased(allele_ix) => allele_ix as usize,
                    GenotypeAllele::PhasedMissing | GenotypeAllele::UnphasedMissing =>
                    {
                        unavail_size[haplotype.shift_ix] += alleles[0].len() as u32;
                        0 // Use reference allele
                    }
                };
                mut_seq.extend_from_slice(alleles[allele_ix]);
                *last_allele_iter.next().unwrap() = allele_ix;
            }
        }
        add_variant_to_alns(var_start - ref_pos, &alleles, &aligner, mism_penalty,
            &last_allele, &mut pairwise_alns);
        ref_pos = var_end;
    }
    let suffix_size = ref_end - ref_pos;
    if suffix_size > 0 {
        let suffix_seq = &ref_seq[(ref_pos - ref_start) as usize..];
        for entry in seqs[haplotypes.samples[0].haplotypes[0].shift_ix..].iter_mut() {
            entry.seq_mut().extend_from_slice(suffix_seq);
        }
        let cigar_suffix = Cigar::new_full_match(suffix_size);
        pairwise_alns.iter_mut().for_each(|(cigar, _)| cigar.extend(&cigar_suffix));
    }

    let aver_unavail = F64Ext::mean(&unavail_size);
    if aver_unavail > 0.0 {
        log::debug!("        On average, {:.1} bp unavailable per haplotype", aver_unavail);
    }
    let avail_seqs: Vec<NamedSeq> = seqs.into_iter().zip(unavail_size.into_iter())
        .filter(|(seq, unavail)| f64::from(*unavail) <= unavail_rate * f64::from(seq.len()))
        .map(|(seq, _)| seq)
        .collect();
    let n_remain = avail_seqs.len();
    if n_remain < haplotypes.total {
        log::warn!("        Reconstructed {} haplotypes ({} unavailable)", n_remain, haplotypes.total - n_remain);
    }
    if n_remain < 2 {
        return Err(Error::InvalidData("Less than two haplotypes reconstructed".to_owned()));
    }
    Ok((avail_seqs, pairwise_alns))
}

/// Finds alignments between all alleles.
/// Returns vector of size `n * n`, where `n` is the number of alleles.
fn add_variant_to_alns(
    between_dist: u32,
    alleles: &[&[u8]],
    aligner: &Aligner,
    mism_penalty: i32,
    last_allele: &[usize],
    pairwise_alns: &mut [(Cigar, i32)],
) {
    let n = alleles.len();
    // Vector with all alignment CIGARs and scores between two alleles.
    let mut allele_alns = Vec::with_capacity(n * n);
    for (i, allele_i) in alleles.iter().enumerate() {
        let i_len = allele_i.len() as u32;
        for (j, allele_j) in alleles.iter().enumerate() {
            match i.cmp(&j) {
                Ordering::Less => {
                    if i_len == 1 && allele_j.len() == 1 {
                        assert_ne!(allele_i, allele_j, "Two variant alleles are identical");
                        allele_alns.push((Cigar::new_full_mismatch(1), mism_penalty));
                    } else {
                        allele_alns.push(aligner.align(allele_i, allele_j));
                    }
                }
                Ordering::Equal => allele_alns.push((Cigar::new_full_match(i_len), 0)),
                Ordering::Greater => {
                    let (cigar, score) = &allele_alns[j * n + i];
                    allele_alns.push((cigar.inverse(), *score));
                }
            }
        }
    }

    let cigar_between = if between_dist > 0 { Cigar::new_full_match(between_dist) } else { Cigar::new() };
    let mut pairwise_alns_iter = pairwise_alns.iter_mut();
    for (i, &ix1) in last_allele.iter().enumerate() {
        for &ix2 in &last_allele[i + 1..] {
            let aln = pairwise_alns_iter.next().expect("Incorrect number of alignments");
            aln.0.extend(&cigar_between);

            let (allele_aln, allele_score) = &allele_alns[ix1 * n + ix2];
            aln.0.extend(allele_aln);
            aln.1 += *allele_score;
        }
    }
}
