use std::path::Path;
use htslib::{
    bcf::{
        Read, IndexedReader,
        record::{Record as VcfRecord, Genotype},
    },
};
use crate::{
    err::error,
};

/// Previously called VCF file.
pub struct HintVcf {
    vcf: IndexedReader,
    sample_id: usize,
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
        Ok(Self { vcf, sample_id })
    }
}

pub struct HintVariants {
    variants: Vec<(VcfRecord, Genotype)>,
}
