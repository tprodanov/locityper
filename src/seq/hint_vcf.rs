use std::{
    path::Path,
    fmt::{self, Display, Formatter},
    cmp::max,
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
    err::{error, Error},
    algo::bisect,
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
    genotype_qual: f32,
    /// Consensus sequences for each haplotype.
    haplotypes: SmallVec<[Vec<u8>; 2]>,
    block_ix: usize,
}

impl Variant {
    fn new(start: u32, end: u32, genotype_qual: f32, haplotypes: SmallVec<[Vec<u8>; 2]>) -> Self {
        Self {
            start, end, genotype_qual, haplotypes,
            block_ix: usize::MAX,
        }
    }

    fn from_record(rec: &Record, genotype: &Genotype, genotype_qual: f32) -> Self {
        let start = rec.pos() as u32;
        let end = rec.end() as u32;
        let alleles = rec.alleles();
        Self::new(start, end, genotype_qual,
            genotype.iter()
                .map(|allele| alleles[allele.index().expect("All alleles must be defined") as usize].to_vec())
                .collect()
        )
    }

    fn is_heterozygous(&self) -> bool {
        let allele0 = self.haplotypes.first().unwrap();
        self.haplotypes[1..].iter().any(|allele| allele != allele0)
    }
}

impl Display for Variant {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}-{} ", self.start + 1, self.end)?;
        for (i, hap) in self.haplotypes.iter().enumerate() {
            if i > 0 {
                f.write_str(" | ")?;
            }
            f.write_str(std::str::from_utf8(hap).unwrap())?;
        }
        write!(f, " (GQ {:.1})", self.genotype_qual)
    }
}

#[derive(Clone, Debug)]
struct VarBlock {
    /// This block contains variants start_ix..end_ix
    start_ix: usize,
    end_ix: usize,
    /// Region between these indices is heterozygous (at least first and last variants in this region are).
    heteroz_range: Option<(usize, usize)>,
}

impl VarBlock {
    #[inline]
    fn new(ix: usize) -> Self {
        Self {
            start_ix: ix,
            end_ix: ix,
            heteroz_range: None,
        }
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.start_ix == self.end_ix
    }
}

/// Variants, groupped into blocks.
/// Variants from the same block are phased in respect to each other.
#[derive(Clone, Debug)]
pub struct VarBlocks {
    variants: Vec<Variant>,
    /// A block `i` combines variants `blocks[i].0..blocks[i].1`.
    /// During block construction there should always at least one (potentially empty) block.
    blocks: Vec<VarBlock>,
}

impl VarBlocks {
    #[inline]
    pub fn new_empty() -> Self {
        VarBlocks {
            variants: Vec::new(),
            blocks: Vec::new(),
        }
    }

    #[inline]
    fn new() -> Self {
        Self {
            variants: Vec::new(),
            blocks: vec![VarBlock::new(0)],
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.variants.is_empty()
    }

    #[inline]
    fn last_block_empty(&self) -> bool {
        self.blocks.last().unwrap().is_empty()
    }

    #[inline]
    fn new_block(&mut self) {
        if !self.last_block_empty() {
            let n = self.variants.len();
            self.blocks.push(VarBlock::new(n));
        }
    }

    fn push_variant(&mut self, mut var: Variant) {
        let n = self.variants.len();
        var.block_ix = self.blocks.len() - 1;
        let last_block = self.blocks.last_mut().unwrap();
        if var.is_heterozygous() {
            last_block.heteroz_range.get_or_insert((n - 1, n)).1 = n;
        }
        last_block.end_ix = n;
        self.variants.push(var);
    }

    fn pop_last_if_empty(&mut self) {
        if self.last_block_empty() {
            self.blocks.pop();
        }
    }
}

impl Display for VarBlocks {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        for (i, block) in self.blocks.iter().enumerate() {
            if i > 0 {
                writeln!(f)?;
            }
            write!(f, "Block #{}, {} variant(s)", i, block.end_ix - block.start_ix)?;
            for var in &self.variants[block.start_ix..block.end_ix] {
                write!(f, "\n  {}", var)?;
            }
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
    blocks: &mut VarBlocks,
    min_gq: f32,
    max_len: u32,
    ploidy: usize,
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
    fn reset_positions(pos: u32, start: &mut u32, ends: &mut [u32], max_end: &mut u32, effect_end: &mut u32) {
        *start = pos;
        ends.fill(pos);
        *max_end = pos;
        *effect_end = pos;
    }

    let mut start = records[0].0.pos() as u32;
    let mut ends: SmallVec<[u32; 2]> = smallvec![start; ploidy];
    let mut max_end = start;
    // Reference start + alt length.
    let mut effect_end = start;
    let mut haplotypes: SmallVec<[Vec<u8>; 2]> = smallvec![Vec::new(); ploidy];
    let mut genotype_qual = f32::INFINITY;

    for (record, genotype, rec_gq) in records {
        if genotype.iter().all(|allele| allele.index() == Some(0)) {
            continue;
        }

        let curr_start = record.pos() as u32;
        let curr_end = record.end() as u32;
        if start == max_end {
            // Between last recorded variant and this position there was no variation, moving start position.
            reset_positions(curr_start, &mut start, &mut ends, &mut max_end, &mut effect_end);
        } else if effect_end <= curr_start {
            // There is a reference gap between haplotypes, splitting into multiple variants.
            if genotype_qual >= min_gq && effect_end - start <= max_len {
                extend_haplotypes(&mut haplotypes, &mut ends, max_end, ref_seq, ref_shift);
                blocks.push_variant(Variant::new(start, max_end, genotype_qual,
                    std::mem::replace(&mut haplotypes, smallvec![Vec::new(); ploidy])));
            } else {
                haplotypes.iter_mut().for_each(Vec::clear);
            }
            reset_positions(curr_start, &mut start, &mut ends, &mut max_end, &mut effect_end);
            genotype_qual = f32::INFINITY;
        } else {
            // There is a variant in progress, extending some haplotypes up to `curr_start`.
            extend_haplotypes(&mut haplotypes, &mut ends, curr_start, ref_seq, ref_shift);
            max_end = max_end.max(curr_start);
            effect_end = effect_end.max(curr_start);
        }

        let alleles = record.alleles();
        genotype_qual = f32::min(genotype_qual, *rec_gq);
        for (&allele, end, hap) in itertools::izip!(genotype.iter(), &mut ends, &mut haplotypes) {
            let i = allele.index().unwrap() as usize;
            if i != 0 {
                if *end != curr_start {
                    log::warn!("Cannot reconcile overlapping variants, skipping some this variant.");
                    reset_positions(curr_end, &mut start, &mut ends, &mut max_end, &mut effect_end);
                    haplotypes.iter_mut().for_each(Vec::clear);
                    genotype_qual = f32::INFINITY;
                    break;
                }
                let allele_seq = &alleles[i];
                hap.extend_from_slice(allele_seq);
                *end = curr_end;
                max_end = max_end.max(curr_end);
                effect_end = effect_end.max(max_end).max(curr_start + allele_seq.len() as u32);
            }
        }
    }
    if start != max_end && genotype_qual >= min_gq && effect_end - start <= max_len {
        extend_haplotypes(&mut haplotypes, &mut ends, max_end, ref_seq, ref_shift);
        blocks.push_variant(Variant::new(start, max_end, genotype_qual, haplotypes));
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
    effect_end: u32,
    retains_phasing: bool,
    // VCF record, genotype and genotype quality.
    records: Vec<(Record, Genotype, f32)>,
}

impl Cluster {
    fn new() -> Self {
        Self {
            last_end: 0,
            effect_end: 0,
            retains_phasing: true,
            records: Vec::new(),
        }
    }

    fn push(&mut self, record: Record, genotype: Genotype, genotype_qual: f32) {
        self.retains_phasing &= retains_phasing(&genotype);
        let alleles = record.alleles();
        let max_len = genotype.iter().map(|i| alleles[i.index().unwrap() as usize].len()).max().unwrap();
        self.last_end = self.last_end.max(record.end() as u32).max(record.pos() as u32 + max_len as u32);
        self.records.push((record, genotype, genotype_qual));
    }

    fn save(&mut self,
        blocks: &mut VarBlocks,
        ref_region: &ReferenceRegion,
        ref_seq: &[u8],
        min_gq: f32,
        max_len: u32,
        ploidy: usize,
    ) {
        let Some(first_record) = self.records.first() else { return };
        if (first_record.0.pos() as u32) < ref_region.start() || self.last_end > ref_region.end() {
            self.retains_phasing = true;
            self.records.clear();
            return;
        }

        if !self.retains_phasing {
            blocks.new_block();
        }
        if self.records.len() == 1 {
            let (record, genotype, gq) = first_record;
            if record.rlen() as u32 <= max_len {
                blocks.push_variant(Variant::from_record(record, genotype, *gq));
            }
        } else {
            process_cluster(&self.records, ref_region.start(), ref_seq, blocks, min_gq, max_len, ploidy);
        }
        if !self.retains_phasing {
            blocks.new_block();
        }
        self.retains_phasing = true;
        self.records.clear();
        // Do not reset last_end, because we will anyway set it later.
    }
}

impl Display for Cluster {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "Cluster (end {}, effective end {}, {}phased, {} record(s))",
            self.last_end, self.effect_end,
            if self.retains_phasing { "" } else { "not " },
            self.records.len())?;
        for (record, gt, qual) in &self.records {
            let alleles = record.alleles();
            write!(f, "\n  {} ", record.pos() + 1)?;
            for (i, allele) in alleles.iter().enumerate() {
                f.write_str(if i == 0 { "" } else if i == 1 { ">" } else { "," })?;
                f.write_str(std::str::from_utf8(allele).unwrap())?;
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

/// Pre-loads contig lengths from the VCF file in O(n).
fn load_contig_lengths(header: &bcf::header::HeaderView) -> Vec<Option<u32>> {
    let mut lengths = Vec::new();
    for record in header.header_records() {
        if let bcf::header::HeaderRecord::Contig { key, values } = record {
            let Some(contig_len_str) = values.get("length") else { continue };
            let Ok(contig_len) = contig_len_str.parse() else {
                log::warn!("Cannot parse contig length in the VCF file");
                continue
            };
            let Ok(rid) = header.name2rid(key.as_bytes()) else {
                log::warn!("Cannot fully parse VCF header");
                continue
            };
            lengths.resize(max(lengths.len(), rid as usize + 1), None);
            lengths[rid as usize] = Some(contig_len);
        }
    }
    lengths
}

/// Previously called VCF file.
pub struct HintVcf {
    vcf: IndexedReader,
    sample: String,
    sample_id: usize,
    contig_lengths: Vec<Option<u32>>,
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
        let sample = String::from_utf8(header.id_to_sample(bcf::header::Id(sample_id as u32)))
            .map_err(|_| Error::utf8("sample name", &header.id_to_sample(bcf::header::Id(sample_id as u32))))?;
        let (gq_type, gq_length) = header.format_type(GQ_TAG)
            .map_err(|_| error!(InvalidInput, "Input VCF must have genotype qualities (GQ format field)"))?;
        if gq_type != bcf::header::TagType::Float && gq_type != bcf::header::TagType::Integer {
            return Err(error!(InvalidInput, "Genotype quality (GQ) must be float or integer"));
        } else if gq_length != bcf::header::TagLength::Fixed(1) {
            return Err(error!(InvalidInput, "Genotype quality (GQ) must have Number=1 in the header"));
        }
        let contig_lengths = load_contig_lengths(vcf.header());
        Ok(Self { vcf, sample, sample_id, contig_lengths, gq_type, params })
    }

    pub fn get_chrom_len(&self, name: &str) -> Option<u32> {
        self.vcf.header().name2rid(name.as_bytes())
            .ok()
            .and_then(|rid| self.contig_lengths.get(rid as usize).copied().flatten())
    }

    pub fn sample(&self) -> &str {
        &self.sample
    }

    pub fn fetch_blocks(
        &mut self,
        ref_region: &ReferenceRegion,
        ref_seq: &[u8],
        ploidy: usize,
    ) -> crate::Result<VarBlocks> {
        let mut blocks = VarBlocks::new();
        let mut cluster = Cluster::new();

        let chrom = ref_region.chrom();
        let rid = self.vcf.header().name2rid(chrom.as_bytes())
            .map_err(|_| error!(RuntimeError, "Cannot find chromosome {} in the hint VCF file", chrom))?;
        self.vcf.fetch(rid, u64::from(ref_region.start()), Some(u64::from(ref_region.end())))
            .map_err(|_| error!(RuntimeError, "Cannot fetch region {} in the hint VCF file", ref_region))?;
        let sample_id = self.sample_id;
        let gq_type = self.gq_type;
        let params = self.params.clone();
        let mut warned_ploidy = false;

        for item in self.vcf.records() {
            let record = item.map_err(|e| error!(ParsingError, "Error parsing VCF file ({})", e))?;
            let genotype = record.genotypes()?.get(sample_id);
            if has_unkown_allele(&genotype) { continue };
            if genotype.len() != ploidy {
                if !warned_ploidy {
                    log::warn!("hint-vcf contains variants of different ploidy, skipping this variant.");
                }
                warned_ploidy = true;
                continue;
            }
            let gt_qual = get_gq(&record, sample_id, gq_type)?;

            if cluster.last_end <= record.pos() as u32 {
                cluster.save(&mut blocks, ref_region, ref_seq, params.min_gq, params.max_len, ploidy);
            }
            cluster.push(record, genotype, gt_qual);
        }
        cluster.save(&mut blocks, ref_region, ref_seq, params.min_gq, params.max_len, ploidy);
        blocks.pop_last_if_empty();
        eprintln!("{}", blocks);
        Ok(blocks)
    }
}

/// Combine Locityper variants with variant blocks from the hint VCF.
pub fn combine(
    writer: &bcf::Writer,
    locityper_vars: Vec<Record>,
    blocks: &VarBlocks,
    ref_region: &ReferenceRegion,
    ref_seq: &[u8],
) -> crate::Result<()> {
    let n = locityper_vars.len();
    let variants = &blocks.variants;
    let m = variants.len();

    for block in &blocks.blocks {
        eprintln!("Processing {:?}", block);
        assert!(!block.is_empty());
        let start = variants[block.start_ix].start;
        let end = variants[block.end_ix - 1].end;
        // variants[i - 1].end() <= start
        let i = bisect::right_by(variants, |var| var.end.cmp(&start));
        let j = bisect::position_of(variants, |var| end <= var.start, i, n);
    }
    todo!()
}
