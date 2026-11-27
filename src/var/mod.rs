use smallvec::{smallvec, SmallVec};
use htslib::{
    bcf::{
        self,
        record::{Record, GenotypeAllele},
    },
};

/// Wrapper around a slice of the reference sequence.
#[derive(Clone, Copy)]
pub struct RefSeq<'a> {
    shift: u32,
    seq: &'a [u8],
}

impl<'a> RefSeq<'a> {
    #[inline]
    pub fn new(shift: u32, seq: &'a [u8]) -> Self {
        Self { shift, seq }
    }

    #[inline(always)]
    pub fn get(&self, start: u32, end: u32) -> &[u8] {
        &self.seq[(start - self.shift) as usize..(end - self.shift) as usize]
    }
}

/// Growable vector with 2 elements on stack.
pub type SmallVec2<T> = SmallVec<[T; 2]>;

/// Growable vector with 2 elements on stack.
pub type SmallVec3<T> = SmallVec<[T; 2]>;

const GQ_TAG: &'static [u8; 2] = b"GQ";
const PS_TAG: &'static [u8; 2] = b"PS";

#[derive(Clone, Debug)]
pub struct Variant {
    start: u32,
    end: u32,
    alleles: SmallVec3<Vec<u8>>,
    /// Can actually represent genotype for multiple samples with fixed ploidy, one after another,
    /// the same way it is stored in htslib.
    genotype: Vec<GenotypeAllele>,
    genotype_qual: Option<i32>,
    phased: bool,
    phase_block: Option<i32>,
}

impl Variant {
    pub fn new(start: u32, end: u32, ref_allele: Vec<u8>, phased: bool) -> Self {
        Self {
            start, end,
            alleles: smallvec![ref_allele],
            genotype: Vec::new(),
            genotype_qual: None,
            phased,
            phase_block: None,
        }
    }

    #[inline]
    pub fn start(&self) -> u32 { self.start }

    #[inline]
    pub fn end(&self) -> u32 { self.end }

    #[inline]
    pub fn alleles(&self) -> &[Vec<u8>] {
        &self.alleles
    }

    /// Adds a new allele (if not seen already) and returns its index.
    pub fn add_allele(&mut self, allele: &[u8]) -> usize {
        self.alleles.iter().position(|old_allele| allele == old_allele)
            .unwrap_or_else(|| {
                self.alleles.push(allele.to_vec());
                self.alleles.len() - 1
            })
    }

    #[inline]
    pub fn genotype_mut(&mut self) -> &mut Vec<GenotypeAllele> {
        &mut self.genotype
    }

    /// If current sample has lower than max ploidy, padd genotype with `count` markers.
    pub fn padd_genotype(&mut self, count: usize) {
        const TRAILING_ALLELE: GenotypeAllele = GenotypeAllele::Phased((i32::MIN >> 1) - 1);
        self.genotype.resize(self.genotype.len() + count, TRAILING_ALLELE);
    }

    #[inline]
    pub fn genotype(&self) -> &[GenotypeAllele] {
        &self.genotype
    }

    /// Iterates over alleles.
    #[inline]
    pub fn genotype_iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.genotype.iter().map(|allele_ix| allele_ix.index().expect("Undefined allele") as usize)
    }

    #[inline]
    pub fn genotype_qual(&self) -> i32 {
        self.genotype_qual.unwrap_or(0)
    }

    #[inline]
    pub fn phase_block(&self) -> Option<i32> {
        self.phase_block
    }

    /// Creats record.
    /// WARN: Setting GQ and PS only works correctly if there is only processed sample.
    pub fn create_record(&self, writer: &bcf::Writer) -> crate::Result<Record> {
        let mut record = writer.empty_record();
        record.set_rid(Some(0));
        record.set_pos(i64::from(self.start));
        let alleles: Vec<_> = self.alleles.iter().map(|seq| seq as &[u8]).collect();
        record.set_alleles(&alleles)?;
        record.push_genotypes(&self.genotype)?;
        if let Some(gq) = self.genotype_qual {
            record.push_format_integer(GQ_TAG, &[gq])?;
        }
        if let Some(ps) = self.phase_block {
            record.push_format_integer(PS_TAG, &[ps])?;
        }
        Ok(record)
    }
}
