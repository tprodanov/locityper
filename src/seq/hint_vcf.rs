use std::{
    path::Path,
    fmt::{self, Display, Formatter},
    cmp::max,
    sync::Arc,
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
    pub min_gq: i32,
    pub max_len: u32,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            sample: None,
            min_gq: 10,
            max_len: 50,
        }
    }
}

const GQ_TAG: &'static [u8; 2] = b"GQ";
const PS_TAG: &'static [u8; 2] = b"PS";

#[inline]
fn has_unkown_allele(gt: &Genotype) -> bool {
    gt.iter().copied().any(|allele| matches!(allele, GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing))
}

/// Obtains genotype quality from the record, if present, otherwise returns 0.
fn get_gq(record: &Record, sample_id: usize, gq_type: bcf::header::TagType) -> crate::Result<i32> {
    Ok(if gq_type == bcf::header::TagType::Integer {
        record.format(GQ_TAG).integer()?[sample_id].get(0).copied()
    } else {
        record.format(GQ_TAG).float()?[sample_id].get(0).map(|&x| x as i32)
    }.unwrap_or(0))
}

/// Returns the phase set of the sample, if present.
fn get_ps(record: &Record, sample_id: usize) -> Option<i32> {
    record.format(PS_TAG).integer().ok()?
        .get(sample_id)?
        .first()
        .copied()
        .filter(|x| *x >= i32::MIN + 1) // According to some sources, this can happen on missing entries (PS=.)
}

/// Growable vector with 2 elements on stack.
type SmallVec2<T> = SmallVec<[T; 2]>;

/// Wrapper around a record with additional helpful information.
#[derive(Clone, Debug)]
struct HeavyRecord {
    start: u32,
    end: u32,
    record: Option<Record>,
    genotype: Option<Genotype>,
    haplotypes: SmallVec2<Vec<u8>>,
    genotype_qual: i32,
    phase_block: Option<i32>,
    phased: bool,
    homozygous: bool,
    all_reference: bool,
}

impl HeavyRecord {
    fn from_record(
        record: Record,
        sample_id: usize,
        gq_type: bcf::header::TagType,
        has_ps: bool,
        ploidy: usize,
        warn_ploidy: &mut bool,
    ) -> crate::Result<Option<Self>> {
        let genotype = record.genotypes()?.get(sample_id);
        if has_unkown_allele(&genotype) { return Ok(None) };
        if genotype.len() != ploidy {
            if *warn_ploidy {
                log::warn!("hint-vcf contains variants of different ploidy, skipping this variant.");
                *warn_ploidy = false;
            }
            return Ok(None);
        }

        let alleles = record.alleles();
        let mut haplotypes: SmallVec2<Vec<u8>> = SmallVec::new();
        let allele_ix0 = genotype[0].index().unwrap();
        let mut homozygous = true;
        let mut phased = true;
        for &allele_ix in genotype.iter() {
            phased &= matches!(allele_ix, GenotypeAllele::Phased(_));
            let allele_ix = allele_ix.index().unwrap();
            homozygous &= allele_ix == allele_ix0;
            haplotypes.push(alleles[allele_ix as usize].to_vec());
        };
        let all_reference = homozygous && allele_ix0 == 0;

        let genotype_qual = get_gq(&record, sample_id, gq_type)?;
        let phase_block = (has_ps && phased).then(|| get_ps(&record, sample_id)).flatten();
        Ok(Some(Self {
            start: record.pos() as u32,
            end: record.end() as u32,
            record: Some(record),
            genotype: Some(genotype),
            haplotypes, genotype_qual, phase_block, phased, homozygous, all_reference,
        }))
    }

    fn from_locityper_rec(record: Record) -> Self {
        let genotype = record.genotypes().expect("Undefined genotype").get(0);
        let alleles = record.alleles();
        let mut haplotypes: SmallVec2<Vec<u8>> = SmallVec::new();
        let allele_ix0 = genotype[0].index().unwrap();
        let mut homozygous = true;
        for &allele_ix in genotype.iter() {
            let allele_ix = allele_ix.index().unwrap();
            homozygous &= allele_ix == allele_ix0;
            haplotypes.push(alleles[allele_ix as usize].to_vec());
        };
        assert!(!homozygous || allele_ix0 != 0, "Locityper variants should not be all-reference");
        Self {
            start: record.pos() as u32,
            end: record.end() as u32,
            record: Some(record),
            genotype: Some(genotype),
            genotype_qual: 0,
            phase_block: None,
            phased: true,
            all_reference: false,
            haplotypes, homozygous,
        }
    }

    fn from_haplotypes(
        start: u32,
        end: u32,
        genotype_qual: i32,
        haplotypes: SmallVec2<Vec<u8>>,
    ) -> Self {
        let hap0 = haplotypes.first().unwrap();
        let homozygous = haplotypes[1..].iter().all(|hap| hap == hap0);
        Self {
            record: None,
            genotype: None,
            phase_block: None,
            all_reference: false,
            phased: true,
            start, end, haplotypes, genotype_qual, homozygous,
        }
    }

    fn record(&self) -> &Record {
        self.record.as_ref().expect("Record must be defined")
    }

    fn genotype(&self) -> &Genotype {
        self.genotype.as_ref().expect("Genotype must be defined")
    }
}

impl Display for HeavyRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}-{} ", self.start + 1, self.end)?;
        for (i, hap) in self.haplotypes.iter().enumerate() {
            if i > 0 {
                f.write_str(" | ")?;
            }
            f.write_str(std::str::from_utf8(hap).unwrap())?;
        }
        write!(f, " (GQ {})", self.genotype_qual)
    }
}

/// Variants, groupped into blocks.
/// Variants from the same block are phased in respect to each other.
#[derive(Clone, Debug)]
pub struct VarBlocks {
    variants: Vec<HeavyRecord>,
    /// A block `i` combines variants `blocks[i].0..blocks[i].1`.
    /// During block construction there should always at least one (potentially empty) block.
    blocks: Vec<(usize, usize)>,
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
            blocks: vec![(0, 0)],
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.variants.is_empty()
    }

    #[inline]
    fn last_block_empty(&self) -> bool {
        let last_block = *self.blocks.last().unwrap();
        last_block.0 == last_block.1
    }

    #[inline]
    fn new_block(&mut self) {
        if !self.last_block_empty() {
            let n = self.variants.len();
            self.blocks.push((n, n));
        }
    }

    #[inline]
    fn push_variant(&mut self, var: HeavyRecord) {
        self.variants.push(var);
        self.blocks.last_mut().unwrap().1 = self.variants.len();
    }

    #[inline]
    fn pop_last_if_empty(&mut self) {
        if self.last_block_empty() {
            self.blocks.pop();
        }
    }
}

impl Display for VarBlocks {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        for (i, &(start, end)) in self.blocks.iter().enumerate() {
            if i > 0 {
                writeln!(f)?;
            }
            write!(f, "Block #{}, {} variant(s)", i, end - start)?;
            for var in &self.variants[start..end] {
                write!(f, "\n  {}", var)?;
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct RefSeq<'a> {
    shift: u32,
    seq: &'a [u8],
}

impl<'a> RefSeq<'a> {
    #[inline]
    fn new(shift: u32, seq: &'a [u8]) -> Self {
        Self { shift, seq }
    }

    #[inline(always)]
    fn get(&self, start: u32, end: u32) -> &[u8] {
        &self.seq[(start - self.shift) as usize..(end - self.shift) as usize]
    }
}

/// Process cluster of overlapping variants.
/// This skips homozygous reference variants,
/// and constructs `HeavyRecord`s with consensus sequences for each non-overlapping region.
fn process_cluster(
    records: &[HeavyRecord],
    ref_seq: RefSeq,
    blocks: &mut VarBlocks,
    min_gq: i32,
    max_len: u32,
    ploidy: usize,
) {
    #[inline]
    fn extend_haplotypes(haplotypes: &mut [Vec<u8>], ends: &mut [u32], up_to: u32, ref_seq: RefSeq) {
        for (hap, end) in itertools::izip!(haplotypes, ends) {
            if *end < up_to {
                hap.extend_from_slice(ref_seq.get(*end, up_to));
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

    let mut start = records[0].start;
    let mut ends: SmallVec2<u32> = smallvec![start; ploidy];
    let mut max_end = start;
    // Reference start + alt length.
    let mut effect_end = start;
    let mut haplotypes: SmallVec2<Vec<u8>> = smallvec![Vec::new(); ploidy];
    let mut genotype_qual = i32::MAX;

    for record in records {
        if record.all_reference { continue };

        let curr_start = record.start;
        let curr_end = record.end;
        if start == max_end {
            // Between last recorded variant and this position there was no variation, moving start position.
            reset_positions(curr_start, &mut start, &mut ends, &mut max_end, &mut effect_end);
        } else if effect_end <= curr_start {
            // There is a reference gap between haplotypes, splitting into multiple variants.
            if genotype_qual >= min_gq && effect_end - start <= max_len {
                extend_haplotypes(&mut haplotypes, &mut ends, max_end, ref_seq);
                blocks.push_variant(HeavyRecord::from_haplotypes(start, max_end, genotype_qual,
                    std::mem::replace(&mut haplotypes, smallvec![Vec::new(); ploidy])));
            } else {
                haplotypes.iter_mut().for_each(Vec::clear);
            }
            reset_positions(curr_start, &mut start, &mut ends, &mut max_end, &mut effect_end);
            genotype_qual = i32::MAX;
        } else {
            // There is a variant in progress, extending some haplotypes up to `curr_start`.
            extend_haplotypes(&mut haplotypes, &mut ends, curr_start, ref_seq);
            max_end = max_end.max(curr_start);
            effect_end = effect_end.max(curr_start);
        }

        genotype_qual = genotype_qual.min(record.genotype_qual);
        for (&allele, allele_seq, end, hap) in
            itertools::izip!(record.genotype().iter(), &record.haplotypes, &mut ends, &mut haplotypes) {
            let i = allele.index().unwrap() as usize;
            if i != 0 {
                if *end != curr_start {
                    log::warn!("Cannot reconcile overlapping variants, skipping some this variant.");
                    reset_positions(curr_end, &mut start, &mut ends, &mut max_end, &mut effect_end);
                    haplotypes.iter_mut().for_each(Vec::clear);
                    genotype_qual = i32::MAX;
                    break;
                }
                hap.extend_from_slice(allele_seq);
                *end = curr_end;
                max_end = max_end.max(curr_end);
                effect_end = effect_end.max(max_end).max(curr_start + allele_seq.len() as u32);
            }
        }
    }
    if start != max_end && genotype_qual >= min_gq && effect_end - start <= max_len {
        extend_haplotypes(&mut haplotypes, &mut ends, max_end, ref_seq);
        blocks.push_variant(HeavyRecord::from_haplotypes(start, max_end, genotype_qual, haplotypes));
    }
}

struct Cluster {
    last_end: u32,
    /// Does this cluster connect to the previous one via phase blocks?
    connects_back: bool,
    /// Same with the next block.
    connects_forward: bool,
    /// VCF record, genotype and genotype quality.
    records: Vec<HeavyRecord>,
}

impl Cluster {
    fn new() -> Self {
        Self {
            last_end: 0,
            connects_back: true,
            connects_forward: true,
            records: Vec::new(),
        }
    }

    fn push(&mut self, record: HeavyRecord, current_phase_block: &mut Option<i32>) {
        if !record.homozygous {
            match (*current_phase_block, record.phase_block) {
                (Some(curr_ps), Some(rec_ps)) => {
                    if curr_ps != rec_ps {
                        self.connects_back = false;
                        // If there is phase mismatch, connect forward only if this is the first variant of the cluster.
                        self.connects_forward &= self.records.is_empty();
                    }
                }
                (Some(_), None) => {
                    self.connects_back = false;
                    self.connects_forward = false;
                }
                (None, Some(_)) => {
                    self.connects_back = false;
                }
                (None, None) => {
                    self.connects_back &= record.phased;
                    self.connects_forward &= record.phased;
                }
            }
            *current_phase_block = record.phase_block;
        }
        let max_len = record.haplotypes.iter().map(|hap| hap.len()).max().unwrap();
        self.last_end = self.last_end.max(record.end).max(record.start + max_len as u32);
        self.records.push(record);
    }

    #[inline]
    fn reset(&mut self) {
        self.connects_back = true;
        self.connects_forward = true;
        self.records.clear();
        // Do not reset last_end because they will be updated on push.
    }

    fn save(&mut self,
        blocks: &mut VarBlocks,
        ref_region: &ReferenceRegion,
        ref_seq: RefSeq,
        min_gq: i32,
        max_len: u32,
        ploidy: usize,
    ) {
        let Some(first_record) = self.records.first() else { return };
        if (first_record.start) < ref_region.start() || self.last_end > ref_region.end() {
            self.reset();
            return;
        }

        if !self.connects_back {
            blocks.new_block();
        }
        if self.records.len() == 1 {
            let record = self.records.pop().expect("Length is exactly 1");
            if record.record().rlen() as u32 <= max_len && record.genotype_qual >= min_gq {
                blocks.push_variant(record);
            }
        } else {
            process_cluster(&self.records, ref_seq, blocks, min_gq, max_len, ploidy);
        }
        if !self.connects_forward {
            blocks.new_block();
        }
        self.reset();
    }
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
    has_ps: bool,
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
        let has_ps = matches!(header.format_type(PS_TAG), Ok((bcf::header::TagType::Integer, _)));
        let contig_lengths = load_contig_lengths(vcf.header());
        Ok(Self { vcf, sample, sample_id, contig_lengths, gq_type, has_ps, params })
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
        let has_ps = self.has_ps;
        let params = self.params.clone();
        let ref_seq = RefSeq::new(ref_region.start(), ref_seq);
        let mut phase_block = None;
        let mut warn_ploidy = true;

        for item in self.vcf.records() {
            let record = item.map_err(|e| error!(ParsingError, "Error parsing VCF file ({})", e))?;
            let Some(record) = HeavyRecord::from_record(record, sample_id, gq_type, has_ps, ploidy, &mut warn_ploidy)?
                else { continue };
            if cluster.last_end <= record.start {
                cluster.save(&mut blocks, ref_region, ref_seq, params.min_gq, params.max_len, ploidy);
            }
            cluster.push(record, &mut phase_block);
        }
        cluster.save(&mut blocks, ref_region, ref_seq, params.min_gq, params.max_len, ploidy);
        blocks.pop_last_if_empty();
        eprintln!("{}", blocks);
        Ok(blocks)
    }
}

// trait VarOrRec {
//     fn start(&self) -> u32;

//     fn end(&self) -> u32;

//     fn haplotypes(&self) -> impl Iterator<Item = &[u8]>;
// }

// impl VarOrRec for Variant {
//     #[inline]
//     fn start(&self) -> u32 { self.start }

//     #[inline]
//     fn end(&self) -> u32 { self.end }

//     #[inline]
//     fn haplotypes(&self) -> impl Iterator<Item = &[u8]> {
//         self.haplotypes.iter().map(AsRef::as_ref)
//     }
// }

// impl VarOrRec for HeavyRecord {
//     #[inline]
//     fn start(&self) -> u32 { self.start }

//     #[inline]
//     fn end(&self) -> u32 { self.end }

//     #[inline]
//     fn haplotypes(&self) -> impl Iterator<Item = &[u8]> {
//         let alleles = self.record.alleles();
//         self.genotype.iter()
//             .map(move |allele_ix| alleles[allele_ix.index().expect("Undefined allele") as usize])
//     }
// }

// fn consensus_sequence<'a>(
//     ref_seq: RefSeq,
//     start: u32,
//     end: u32,
//     variants: impl Iterator<Item = &'a Variant>,
//     ploidy: usize,
// ) -> SmallVec2<Vec<u8>> {
//     let mut last_end = start;
//     let mut consensuses: SmallVec2<Vec<u8>> = smallvec![Vec::new(); ploidy];
//     for var in variants {
//         let between_seq = ref_seq.get(last_end, var.start);
//         for (consensus, hap) in itertools::izip!(&mut consensuses, &var.haplotypes) {
//             consensus.extend_from_slice(between_seq);
//             consensus.extend_from_slice(hap);
//         }
//         last_end = var.end;
//     }
//     if last_end < end {
//         let between_seq = ref_seq.get(last_end, end);
//         for consensus in &mut consensuses {
//             consensus.extend_from_slice(between_seq);
//         }
//     }
//     consensuses
// }

// fn create_record(
//     header: &Arc<bcf::header::HeaderView>,
//     start: u32,
//     end: u32,
//     ref_seq: RefSeq,
//     hap_alleles: &[impl AsRef<[u8]>],
//     phase_block: Option<i32>,
// ) -> Record {
//     let mut record = Record::new(Arc::clone(header));
//     record.set_rid(Some(0));
//     record.set_pos(i64::from(start));
//     record.set_qual(60.0);

//     let ref_allele = ref_seq.get(start, end);
//     let mut alleles: SmallVec<[&[u8]; 3]> = smallvec![ref_allele];
//     let mut genotype: SmallVec2<GenotypeAllele> = SmallVec::new();
//     let phased = phase_block.is_some();
//     for allele in hap_alleles {
//         let allele = allele.as_ref();
//         let allele_ix = alleles.iter().position(|&existing_allele| existing_allele == allele)
//             .unwrap_or_else(|| {
//                 alleles.push(allele);
//                 alleles.len() - 1
//             }) as i32;
//         genotype.push(if phased { GenotypeAllele::Phased(allele_ix) } else { GenotypeAllele::Unphased(allele_ix) });
//     }
//     record.set_alleles(&alleles).expect("Failed to set alleles");
//     record.push_genotypes(&genotype).expect("Failed to set genotype");
//     if let Some(ps) = phase_block {
//         record.push_format_integer(PS_TAG, &[ps]).expect("Failed to set PS tag");
//     }
//     record
// }

// /// Identify, which of the hint haplotypes is more similar to which of the Locityper haplotypes.
// /// Returns conversion: `conv[hint_haplotype] = locityper_haplotype`, if phasing is possible.
// fn phase_block(
//     block_variants: &[Variant],
//     locityper_recs: &[Record],
//     ref_seq: RefSeq,
// ) -> Option<SmallVec2<usize>> {
//     // Returns None if there were no heterozygous variants, so there is nothing to phase.
//     let het_i = block_variants.iter().position(|var| !var.homozygous)?;
//     let het_j = block_variants.iter().rposition(|var| !var.homozygous).unwrap();
//     let het_reg_variants = &block_variants[het_i..=het_j];
//     let mut het_start = het_reg_variants.first().unwrap().start;
//     let mut het_end = het_reg_variants.last().unwrap().end;
//     // Returns None if no Locityper variants in the heterozygous region.
//     let lt_i = locityper_recs.iter().position(|rec| het_start < rec.end() as u32)?;
//     let lt_j = locityper_recs.iter().rposition(|rec| (rec.pos() as u32) < het_end).unwrap();
//     if lt_i > lt_j {
//         return None;
//     }
//     let het_reg_lt_recs = &locityper_recs[lt_i..=lt_j];
//     het_start = het_start.min(locityper_recs.first().unwrap().pos() as u32);
//     het_end = het_end.max(locityper_recs.last().unwrap().end() as u32);
//     todo!()
// }

// fn process_block(
//     block_ps: Option<i32>,
//     block_variants: &[Variant],
//     locityper_recs: &[Record],
//     header: &Arc<bcf::header::HeaderView>,
//     ref_seq: RefSeq,
//     new_records: &mut Vec<Record>,
// ) {
//     if locityper_recs.is_empty() {
//         // Case 1: No Locityper records.
//         new_records.extend(
//             block_variants.iter()
//                 .filter(|var| !var.all_reference)
//                 .map(|var| create_record(header, var.start, var.end, ref_seq, &var.haplotypes, block_ps))
//         );
//         return;
//     }

// }

// /// Combine Locityper variants with variant blocks from the hint VCF.
// pub fn combine(
//     header: &Arc<bcf::header::HeaderView>,
//     mut locityper_recs: Vec<Record>,
//     blocks: &VarBlocks,
//     ref_region: &ReferenceRegion,
//     ref_seq: &[u8],
// ) -> crate::Result<Vec<Record>> {
//     let n = locityper_recs.len();
//     let ref_seq = RefSeq::new(ref_region.start, ref_seq);
//     let variants = &blocks.variants;
//     let mut new_records = Vec::new();
//     let copied_up_to = 0;

//     for &(start_ix, end_ix) in &blocks.blocks {
//         assert!(start_ix < end_ix);
//         let block_variants = &variants[start_ix..end_ix];
//         let start = block_variants.first().unwrap().start;
//         let end = block_variants.last().unwrap().end;
//         let block_ps = (block_variants.len() > 1).then_some(start);
//         let i = bisect::right_by(&locityper_recs, |rec| (rec.end() as u32).cmp(&start));
//         let j = bisect::position_of(&locityper_recs, |rec| end <= rec.pos() as u32, i, n);
//         new_records.extend_from_slice(&locityper_recs[copied_up_to..i]);
//         // [TODO] Update copied_up_to. (at least i).
//         // [TODO] Add source.
//     }
//     Ok(new_records)
// }
