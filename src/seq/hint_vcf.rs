use std::{
    path::{Path, PathBuf},
};
use htslib::{
    bcf::{
        self,
        Read, IndexedReader,
        record::{Record as VcfRecord, Genotype},
    },
};
use crate::{
    err::error,
};
use super::{
    Interval,
};

pub struct Params {
    pub vcf: Option<PathBuf>,
    pub already_phased: bool,
    pub ref_hap_or_bed: Option<PathBuf>,
    pub min_gq: f32,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            vcf: None,
            already_phased: false,
            ref_hap_or_bed: None,
            min_gq: 10.0,
        }
    }
}

const GQ_TAG: &'static [u8; 2] = b"GQ";

fn homozygous_or_unknown(gt: &Genotype) -> bool {
    use bcf::record::GenotypeAllele::*;
    gt.iter().copied().all(|allele| match allele {
        Unphased(0) | Phased(0) | UnphasedMissing | PhasedMissing => true,
        _ => false,
    })
}

/// Previously called VCF file.
pub struct HintVcf {
    vcf: IndexedReader,
    sample_id: usize,
    gq_type: bcf::header::TagType,
}

impl HintVcf {
    pub fn new(filename: impl AsRef<Path>, sample: Option<&str>) -> crate::Result<Self> {
        let filename = filename.as_ref();
        let vcf = IndexedReader::from_path(filename)?;
        let header = vcf.header();
        let sample_id = match sample {
            Some(name) => header.sample_id(name.as_bytes())
                .ok_or_else(|| error!(InvalidInput, "Sample {} not found in `{}`", name, filename.display()))?,
            None if header.sample_count() == 1 => 0,
            None => return Err(error!(InvalidInput,
                "VCF `{}` has more than one sample, please specify via --hint-sample", filename.display())),
        };
        let (gq_type, gq_length) = header.format_type(GQ_TAG)
            .map_err(|_| error!(InvalidInput, "Input VCF must have genotype qualities (GQ format field)"))?;
        if gq_type != bcf::header::TagType::Float || gq_type != bcf::header::TagType::Integer {
            return Err(error!(InvalidInput, "Genotype quality (GQ) must be float or integer"));
        } else if gq_length != bcf::header::TagLength::Fixed(1) {
            return Err(error!(InvalidInput, "Genotype quality (GQ) must have Number=1 in the header"));
        }
        Ok(Self { vcf, sample_id, gq_type })
    }

    pub fn fetch(&mut self, interval: &Interval, min_qual: f32) -> crate::Result<HintVariants> {
        let contig_name = interval.contig_name();
        let rid = self.vcf.header().name2rid(contig_name.as_bytes())
            .map_err(|_| error!(RuntimeError, "Cannot find chromosome {} in the hint VCF file", contig_name))?;
        self.vcf.fetch(rid, u64::from(interval.start()), Some(u64::from(interval.end())))
            .map_err(|_| error!(RuntimeError, "Cannot fetch region {} in the hint VCF file", interval))?;
        for item in self.vcf.records() {
            let record = item.map_err(|e| error!(ParsingError, "Error parsing VCF file ({})", e))?;
            let genotype = record.genotypes()?.get(self.sample_id);
            if homozygous_or_unknown(&genotype) {
                continue;
            }
            let gt_qual = if self.gq_type == bcf::header::TagType::Integer {
                record.format(GQ_TAG).integer()?[self.sample_id][0] as f32
            } else {
                record.format(GQ_TAG).float()?[self.sample_id][0]
            };
            if gt_qual < min_qual {
                continue;
            }
        }
        todo!()
    }
}

pub struct HintVariants {
    variants: Vec<(VcfRecord, Genotype)>,
}
