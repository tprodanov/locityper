use std::{
    path::Path,
    fmt::{self, Display, Formatter, Write},
};
use smallvec::{smallvec, SmallVec};
use htslib::{
    bcf::{
        self,
        Read, IndexedReader,
        record::{Record, Genotype, GenotypeAllele},
    },
};
use crate::{
    err::error,
};
use super::{
    interv::ReferenceRegion,
};

#[derive(Clone, Debug)]
pub struct Params {
    pub sample: Option<String>,
    pub min_gq: f32,
    pub max_len: u32,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            sample: None,
            min_gq: 10.0,
            max_len: 50,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Variant {
    start: u32,
    end: u32,
    /// Consensus sequences for each haplotype.
    haplotypes: SmallVec<[Vec<u8>; 2]>,
    genotype_qual: f32,
}

impl Variant {
    fn new(rec: &Record, genotype: &Genotype, genotype_qual: f32) -> Self {
        let start = rec.pos() as u32;
        let end = rec.end() as u32;
        let alleles = rec.alleles();
        Self {
            start, end, genotype_qual,
            haplotypes: genotype.iter()
                .map(|allele| alleles[allele.index().expect("All alleles must be defined") as usize].to_vec())
                .collect(),
        }
    }
}

// [DEBUG] Remove later.
/// Writes a sequence, eliding the middle of long ones.
fn write_seq(f: &mut Formatter<'_>, seq: &[u8]) -> fmt::Result {
    const MAX: usize = 40;
    let s = |b: &[u8]| String::from_utf8_lossy(b).into_owned();
    if seq.len() <= MAX {
        f.write_str(&s(seq))
    } else {
        write!(f, "{}..{}({}bp)..{}", s(&seq[..15]), "", seq.len(), s(&seq[seq.len() - 15..]))
    }
}

impl Display for Variant {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}-{} ", self.start + 1, self.end)?;
        for (i, hap) in self.haplotypes.iter().enumerate() {
            if i > 0 {
                f.write_char('|')?;
            }
            write_seq(f, hap)?;
        }
        write!(f, " (GQ {:.1})", self.genotype_qual)
    }
}

/// Block of phased variants.
#[derive(Debug, Clone)]
pub struct VarBlock {
    start: u32,
    end: u32,
    variants: SmallVec<[Variant; 1]>,
}

impl VarBlock {
    fn new() -> Self {
        Self {
            start: u32::MAX,
            end: 0,
            variants: SmallVec::new(),
        }
    }

    fn is_empty(&self) -> bool {
        self.variants.is_empty()
    }

    fn push(&mut self, var: Variant) {
        // if var.start < self.end {
        //     // Add overlapping variant.
        //     if
        // }
        self.start = self.start.min(var.start);
        self.end = self.end.max(var.end);
        self.variants.push(var);
    }
}

impl Display for VarBlock {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "Block {}-{}, {} variant(s)", self.start + 1, self.end, self.variants.len())?;
        for var in &self.variants {
            write!(f, "\n  {}", var)?;
        }
        Ok(())
    }
}

/// Process cluster of overlapping variants.
/// This skips homologous reference variants,
/// and constructs `Variant` with consensus sequences for each non-overlapping region.
fn process_cluster(
    records: &[(Record, Genotype, f32)],
    ref_shift: u32,
    ref_seq: &[u8],
    block: &mut VarBlock,
    min_gq: f32,
    max_len: u32,
) {
    #[inline]
    fn extend_haplotypes(haplotypes: &mut [Vec<u8>], ends: &mut [u32], up_to: u32, ref_seq: &[u8], ref_shift: u32) {
        let rel_end = (up_to - ref_shift) as usize;
        for (hap, end) in itertools::izip!(haplotypes, ends) {
            if *end < up_to {
                hap.extend_from_slice(&ref_seq[(*end - ref_shift) as usize..rel_end]);
                *end = up_to;
            }
        }
    }

    #[inline]
    fn reset_positions(pos: u32, start: &mut u32, ends: &mut [u32], max_end: &mut u32) {
        *start = pos;
        ends.fill(pos);
        *max_end = pos;
    }

    debug_assert!(records.len() > 1);
    let ploidy = records[0].1.len();
    let mut start = records[0].0.pos() as u32;
    let mut ends: SmallVec<[u32; 2]> = smallvec![start; ploidy];
    let mut max_end = start;
    let mut haplotypes: SmallVec<[Vec<u8>; 2]> = smallvec![Vec::new(); ploidy];
    let mut genotype_qual = f32::INFINITY;

    for (record, genotype, rec_gq) in records {
        if genotype.iter().all(|allele| allele.index() == Some(0)) {
            continue;
        }

        let curr_start = record.pos() as u32;
        let curr_end = record.end() as u32;
        if genotype.len() != ploidy {
            log::warn!("hint-vcf contains variants of different ploidy, skipping this variant.");
            reset_positions(curr_end, &mut start, &mut ends, &mut max_end);
            haplotypes.iter_mut().for_each(Vec::clear);
            genotype_qual = f32::INFINITY;
            continue;
        }
        if start == max_end {
            // Between last recorded variant and this position there was no variation, moving start position.
            reset_positions(curr_start, &mut start, &mut ends, &mut max_end);
        } else if max_end <= curr_start {
            // There is a reference gap between haplotypes, splitting into multiple variants.
            if genotype_qual >= min_gq && max_end - start <= max_len {
                extend_haplotypes(&mut haplotypes, &mut ends, max_end, ref_seq, ref_shift);
                block.push(Variant {
                    start, genotype_qual,
                    end: max_end,
                    haplotypes: std::mem::replace(&mut haplotypes, smallvec![Vec::new(); ploidy]),
                });
            } else {
                haplotypes.iter_mut().for_each(Vec::clear);
            }
            reset_positions(curr_start, &mut start, &mut ends, &mut max_end);
            genotype_qual = f32::INFINITY;
        } else {
            // There is a variant in progress, extending some haplotypes up to `curr_start`.
            extend_haplotypes(&mut haplotypes, &mut ends, curr_start, ref_seq, ref_shift);
            max_end = max_end.max(curr_start);
        }

        let alleles = record.alleles();
        genotype_qual = f32::min(genotype_qual, *rec_gq);
        for (&allele, end, hap) in itertools::izip!(genotype.iter(), &mut ends, &mut haplotypes) {
            let i = allele.index().unwrap() as usize;
            if i != 0 {
                if *end != curr_start {
                    log::warn!("Cannot reconcile overlapping variants, skipping some this variant.");
                    reset_positions(curr_end, &mut start, &mut ends, &mut max_end);
                    haplotypes.iter_mut().for_each(Vec::clear);
                    genotype_qual = f32::INFINITY;
                    break;
                }
                hap.extend_from_slice(&alleles[i]);
                *end = curr_end;
                max_end = max_end.max(curr_end);
            }
        }
    }
    if start != max_end && genotype_qual >= min_gq && max_end - start <= max_len {
        extend_haplotypes(&mut haplotypes, &mut ends, max_end, ref_seq, ref_shift);
        block.push(Variant {
            start, genotype_qual,
            end: max_end,
            haplotypes: std::mem::replace(&mut haplotypes, smallvec![Vec::new(); ploidy]),
        });
    }
}

const GQ_TAG: &'static [u8; 2] = b"GQ";

fn has_unkown_allele(gt: &Genotype) -> bool {
    use GenotypeAllele::*;
    gt.iter().copied().any(|allele| match allele {
        UnphasedMissing | PhasedMissing => true,
        _ => false,
    })
}

/// Returns true if all alleles are phased, or if the variant is homozygous.
fn retains_phasing(gt: &Genotype) -> bool {
    use GenotypeAllele::*;
    let allele0 = gt[0].index().unwrap();
    let mut homozygous = true;
    let mut phased = true;
    for &allele in gt.iter() {
        homozygous &= allele.index().unwrap() == allele0;
        if let Unphased(_) = allele {
            phased = false;
        }
    };
    homozygous || phased
}

struct Cluster {
    last_end: u32,
    retains_phasing: bool,
    // VCF record, genotype and genotype quality.
    records: Vec<(Record, Genotype, f32)>,
}

impl Cluster {
    fn new() -> Self {
        Self {
            last_end: 0,
            retains_phasing: true,
            records: Vec::new(),
        }
    }

    fn push(&mut self, record: Record, genotype: Genotype, genotype_qual: f32) {
        self.retains_phasing &= retains_phasing(&genotype);
        self.last_end = self.last_end.max(record.end() as u32);
        self.records.push((record, genotype, genotype_qual));
    }

    fn save(&mut self,
        blocks: &mut Vec<VarBlock>,
        ref_region: &ReferenceRegion,
        ref_seq: &[u8],
        min_gq: f32,
        max_len: u32,
    ) {
        if self.records.is_empty() { return };
        let first_start = self.records[0].0.pos() as u32;
        if first_start < ref_region.start() || self.last_end >= ref_region.end() {
            self.retains_phasing = true;
            self.records.clear();
            return;
        }

        if !blocks.last().unwrap().is_empty() && !self.retains_phasing {
            blocks.push(VarBlock::new());
        }
        let block = blocks.last_mut().unwrap();

        if self.records.len() == 1 {
            let (record, genotype, gq) = &self.records[0];
            block.push(Variant::new(record, genotype, *gq));
        } else {
            process_cluster(&self.records, ref_region.start(), ref_seq, block, min_gq, max_len);
        }

        if !self.retains_phasing {
            blocks.push(VarBlock::new());
        }
        self.retains_phasing = true;
        self.records.clear();
        // Do not reset last_end, because we will anyway set it later.
    }
}

impl Display for Cluster {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "Cluster (end {}, {}phased, {} record(s))",
            self.last_end,
            if self.retains_phasing { "" } else { "not " },
            self.records.len())?;
        for (record, gt, qual) in &self.records {
            let alleles = record.alleles();
            write!(f, "\n  {} ", record.pos() + 1)?;
            for (i, allele) in alleles.iter().enumerate() {
                f.write_str(if i == 0 { "" } else if i == 1 { ">" } else { "," })?;
                write_seq(f, allele)?;
            }
            write!(f, " GT={} GQ={:.1}", gt, qual)?;
        }
        Ok(())
    }
}

/// Obtains genotype quality from the record, if present, otherwise returns 0.0
fn get_gq(record: &Record, sample_id: usize, gq_type: bcf::header::TagType) -> crate::Result<f32> {
    Ok(if gq_type == bcf::header::TagType::Integer {
        record.format(GQ_TAG).integer()?[sample_id].get(0).map(|&x| x as f32)
    } else {
        record.format(GQ_TAG).float()?[sample_id].get(0).copied()
    }.unwrap_or(0.0))
}

/// Previously called VCF file.
pub struct HintVcf {
    vcf: IndexedReader,
    sample_id: usize,
    gq_type: bcf::header::TagType,
    params: Params,
}

impl HintVcf {
    pub fn new(filename: impl AsRef<Path>, params: Params) -> crate::Result<Self> {
        let filename = filename.as_ref();
        let vcf = IndexedReader::from_path(filename)?;
        let header = vcf.header();
        let sample_id = match &params.sample {
            Some(name) => header.sample_id(name.as_bytes())
                .ok_or_else(|| error!(InvalidInput, "Sample {} not found in `{}`", name, filename.display()))?,
            None if header.sample_count() == 1 => 0,
            None => return Err(error!(InvalidInput,
                "VCF `{}` has more than one sample, please specify via --hint-sample", filename.display())),
        };
        let (gq_type, gq_length) = header.format_type(GQ_TAG)
            .map_err(|_| error!(InvalidInput, "Input VCF must have genotype qualities (GQ format field)"))?;
        if gq_type != bcf::header::TagType::Float && gq_type != bcf::header::TagType::Integer {
            return Err(error!(InvalidInput, "Genotype quality (GQ) must be float or integer"));
        } else if gq_length != bcf::header::TagLength::Fixed(1) {
            return Err(error!(InvalidInput, "Genotype quality (GQ) must have Number=1 in the header"));
        }
        Ok(Self { vcf, sample_id, gq_type, params })
    }

    pub fn fetch_blocks(
        &mut self,
        ref_region: &ReferenceRegion,
        ref_seq: &[u8],
    ) -> crate::Result<Vec<VarBlock>> {
        let mut blocks = vec![VarBlock::new()];
        let mut cluster = Cluster::new();

        let chrom = ref_region.chrom();
        let rid = self.vcf.header().name2rid(chrom.as_bytes())
            .map_err(|_| error!(RuntimeError, "Cannot find chromosome {} in the hint VCF file", chrom))?;
        self.vcf.fetch(rid, u64::from(ref_region.start()), Some(u64::from(ref_region.end())))
            .map_err(|_| error!(RuntimeError, "Cannot fetch region {} in the hint VCF file", ref_region))?;
        let sample_id = self.sample_id;
        let gq_type = self.gq_type;
        let params = self.params.clone();
        for item in self.vcf.records() {
            let record = item.map_err(|e| error!(ParsingError, "Error parsing VCF file ({})", e))?;
            let genotype = record.genotypes()?.get(sample_id);
            if has_unkown_allele(&genotype) { continue };
            let gt_qual = get_gq(&record, sample_id, gq_type)?;

            if cluster.last_end <= record.pos() as u32 {
                cluster.save(&mut blocks, ref_region, ref_seq, params.min_gq, params.max_len);
            }
            cluster.push(record, genotype, gt_qual);
        }
        cluster.save(&mut blocks, ref_region, ref_seq, params.min_gq, params.max_len);

        for block in &blocks {
            eprintln!("{}", block);
        }
        Ok(blocks)
    }
}
