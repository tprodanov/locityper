use std::{
    path::{Path, PathBuf},
    time::Instant,
    ffi::OsStr,
    sync::Arc,
};
use colored::Colorize;
use regex::Regex;
use htslib::{
    bcf::{self, header::Header, record::GenotypeAllele},
};
use crate::{
    err::{error, add_path, validate_param},
    ext::{
        self,
        fmt::PrettyU32,
    },
    seq::{
        interv::{self, ReferenceRegion},
        paf::{PafParseResult, PafFile},
        cigar::{Cigar, Operation},
        contigs::{ContigId, ContigNames, ContigSet, DiscardedHaplotypes},
    },
    var::Variant,
    algo::{bisect, HashMap, IntMap},
    math::ifelse0,
    command::paths,
};

struct Args {
    input: Option<PathBuf>,
    paf: Option<PathBuf>,
    fasta: Option<PathBuf>,
    disc_filename: String,
    out_merged: Option<PathBuf>,
    out_separate: Option<PathBuf>,
    ref_hap: Option<String>,
    region: String,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: None,
            paf: None,
            fasta: None,
            disc_filename: String::from("auto"),
            out_merged: None,
            out_separate: None,
            ref_hap: None,
            region: String::from("auto"),
        }
    }
}

impl Args {
    fn replace_input(mut self) -> Self {
        if let Some(dir) = &self.input {
            self.paf = self.paf.or_else(|| paths::LOCUS_PAFS.iter()
                .map(|path| dir.join(path)).filter(|path| path.exists()).next());
            self.fasta = self.fasta.or_else(|| Some(dir.join(paths::LOCUS_FASTA)));
            self.out_merged = self.out_merged.or_else(|| Some(dir.join("haplotypes.vcf.gz")));
        }
        self
    }

    fn validate(self) -> crate::Result<Self> {
        validate_param!(self.paf.is_some(), "Input PAF file is not provided (see -p/--paf)");
        validate_param!(self.fasta.is_some(), "Input FASTA file is not provided (see -f/--fasta)");
        validate_param!(self.out_merged.is_some(), "Output VCF file must be provided (see -o/--output)");
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 15;
    const VAL: usize = 4;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    // let defaults = Args::default();
    println!("{}", "Convert PAF file to VCF file(s).".yellow());

    println!("\n{}", "Usage:".bold());
    println!("    {} -i dir [-r ref_hap] [args]", super::PROGRAM);
    println!("    {} -p input.paf -f input.fa -r ref_hap [-d discarded_haplotypes.txt] \\", super::PROGRAM);
    println!("        -o merged.vcf.gz [separate.vcf.gz] [args]");

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Input locityper/loci/{{}} directory. Unless specified\n\
        {EMPTY}  otherwise, sets {}, {}, {}, {} to files in that directory.",
        "-i, --input".green(), "DIR".yellow(), "-p".green(), "-f".green(), "-d".green(), "-o".green());
    println!("    {:KEY$} {:VAL$}  Input PAF[.gz] file with alignments between haplotypes.",
        "-p, --paf".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Input FASTA[.gz] file with haplotypes.",
        "-f, --fasta".green(), "FILE".yellow());
    println!("    {:KEY$} {}\n\
        {EMPTY}  Text file with discarded haplotypes [{}].\n\
        {EMPTY}  {} = <Directory of {}>/discarded_haplotypes.txt.",
        "-d, --discarded".green(), "FILE|auto|none".yellow(), super::fmt_def("auto"),
        "auto".yellow(), "-f".green());
    println!("    {:KEY$} {:VAL$}\n\
        {EMPTY}  Two output VCF[.gz] files, first with merged variants\n\
        {EMPTY}  and second (optional) with unmerged variants.",
        "-o, --output".green(), "FILE [FILE]".yellow());

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Reference haplotypes name. Optional if\n\
        {EMPTY}  <-i directory>/ref.bed or -R region.bed contain haplotype name\n\
        {EMPTY}  in the format \"ref-hap=...\" in one of the columns.",
        "-r, --ref-hap".green(), "STR".yellow());
    println!("    {:KEY$} {:VAL$}  Adjust variant positions relative to this region [{}].\n\
        {EMPTY}  Should be either {} (search for ref.bed relative to haplotypes PAF), {},\n\
        {EMPTY}  `chrom:start-end` or a BED file with one entry.\n\
        {EMPTY}  {}: This will only work if the corresponding haplotype is on the forward strand.",
        "-R, --region".green(), "STR".yellow(), super::fmt_def("auto"), "auto".yellow(), "none".yellow(), "NOTE".red());

    println!("\n{}", "Other arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Show this help message.", "-h, --help".green(), "");
    println!("    {:KEY$} {:VAL$}  Show version.", "-V, --version".green(), "");
}

fn parse_args(argv: &[String]) -> crate::Result<Args> {
    if argv.is_empty() {
        print_help();
        std::process::exit(1);
    }
    use lexopt::prelude::*;
    let mut args = Args::default();
    let mut parser = lexopt::Parser::from_args(argv);

    while let Some(arg) = parser.next()? {
        match arg {
            Short('i') | Long("input") => args.input = Some(parser.value()?.parse()?),
            Short('p') | Long("paf") => args.paf = Some(parser.value()?.parse()?),
            Short('f') | Long("fasta") => args.fasta = Some(parser.value()?.parse()?),
            Short('d') | Long("discarded") => args.disc_filename = parser.value()?.parse()?,
            Short('o') | Long("output") => {
                let mut values = parser.values()?.take(2);
                args.out_merged = Some(values.next().expect("First argument is always present").parse()?);
                if let Some(val) = values.next() {
                    args.out_separate = Some(val.parse()?);
                }
            }
            Short('r') | Long("ref-hap") => args.ref_hap = Some(parser.value()?.parse()?),
            Short('R') | Long("region") => args.region = parser.value()?.parse()?,

            Short('V') | Long("version") => {
                super::print_version();
                std::process::exit(0);
            }
            Short('h') | Long("help") | Long("full-help") | Long("hidden-help") => {
                print_help();
                std::process::exit(0);
            }
            _ => Err(arg.unexpected())?,
        }
    }
    Ok(args)
}

pub fn load_region(
    s: &str,
    fasta_filename: &Path,
) -> crate::Result<Option<ReferenceRegion>> {
    if s == "none" { return Ok(None) };

    let fname = if s == "auto" {
        Path::new(fasta_filename).parent()
            .ok_or_else(|| error!(RuntimeError, "Cannot get parent directory for {}", s))?
            .join(super::paths::LOCUS_BED)
    } else {
        let pos_re = Regex::new(interv::CHROM_POS_PATTERN).unwrap();
        let interval_re = Regex::new(interv::INTERVAL_PATTERN).unwrap();
        if let Some(captures) = pos_re.captures(s).or_else(|| interval_re.captures(s)) {
            let chrom = captures[1].to_string();
            let start = captures[2].parse::<PrettyU32>()
                .map_err(|_| error!(ParsingError, "Cannot parse interval '{}'", s))?.get().strict_sub(1);
            let end = captures[3].parse::<PrettyU32>()
                .map_err(|_| error!(ParsingError, "Cannot parse interval '{}'", s))?.get();
            return Ok(Some(ReferenceRegion::new(chrom, start, end, None)));
        }
        PathBuf::from(s)
    };

    if !fname.exists() {
        log::warn!("Cannot find BED file {}, using relative locus coordinates", fname.display());
        return Ok(None);
    }
    ReferenceRegion::load(fname).map(Some)
}

#[derive(Clone, Debug)]
pub(super) struct VarRange {
    ref_start: u32,
    ref_end: u32,
    hap_start: u32,
    hap_end: u32,
}

impl VarRange {
    #[inline]
    pub(super) fn new(ref_start: u32, ref_end: u32, hap_start: u32, hap_end: u32) -> Self {
        Self { ref_start, ref_end, hap_start, hap_end }
    }

    #[inline]
    pub(super) fn ref_len(&self) -> u32 {
        self.ref_end - self.ref_start
    }

    #[inline]
    pub(super) fn alt_len(&self) -> u32 {
        self.hap_end - self.hap_start
    }
}

/// Canonize gap by moving its position at max to the left.
/// Resulting gap start will be >= min_start.
/// Gap can appear both on reference or alternate haplotypes (does not matter), gap start is reference-based.
/// Between min_start and gap_start sequence must be non-variable.
fn gap_move_left(full_ref_seq: &[u8], mut gap_start: u32, gap_seq: &[u8], min_start: u32) -> u32 {
    let last_ix = gap_seq.len() - 1;
    let mut k = last_ix;
    while gap_start > min_start && gap_seq[k] == full_ref_seq[gap_start as usize - 1] {
        gap_start -= 1;
        k = k.checked_sub(1).unwrap_or(last_ix);
    }
    gap_start
}

/// Move variants to canonical representation.
fn move_all_left(vars: &mut Vec<VarRange>, ref_seq: &[u8], hap_seq: &[u8]) {
    let mut last_end = 0;
    for var in vars {
        // Update here so that we can safely use continue.
        let min_start = last_end;
        last_end = var.ref_end;

        if var.ref_len() == var.alt_len() { continue; }
        let var_ref_seq = &ref_seq[var.ref_start as usize..var.ref_end as usize];
        let var_alt_seq = &hap_seq[var.hap_start as usize..var.hap_end as usize];
        // Prefixes do not match.
        if var_ref_seq.iter().zip(var_alt_seq).any(|(c1, c2)| c1 != c2) { continue; }

        let prefix_size = usize::min(var_ref_seq.len(), var_alt_seq.len());
        let gap_seq = if prefix_size == var_ref_seq.len() {
            &var_alt_seq[prefix_size..]
        } else {
            &var_ref_seq[prefix_size..]
        };
        let gap_start = var.ref_start + prefix_size as u32;
        let new_start = gap_move_left(ref_seq, gap_start, gap_seq, min_start + prefix_size as u32);
        let shift = gap_start - new_start;
        assert!(var.hap_start >= shift);
        assert!(var.ref_start - shift >= min_start);
        var.ref_start -= shift;
        var.ref_end -= shift;
        var.hap_start -= shift;
        var.hap_end -= shift;
    }
}

/// Variant ranges on various haplotypes.
#[derive(Default)]
pub(super) struct VarRanges {
    /// For each haplotype, vector of variant ranges.
    haplotypes: IntMap<ContigId, Vec<VarRange>>,
}

impl VarRanges {
    pub fn add_haplotype_ranges(&mut self,
        ref_seq: &[u8],
        hap_id: ContigId,
        hap_seq: &[u8],
        cigar: &Cigar,
    ) -> crate::Result<()> {
        let mut vars = Vec::<VarRange>::new();
        let mut rpos = 0;
        let mut qpos = 0;
        for item in cigar.iter() {
            let op = item.operation();
            if op == Operation::Equal {
                rpos += item.len();
                qpos += item.len();
                continue;
            }

            if op == Operation::Match || op == Operation::Hard {
                return Err(error!(RuntimeError, "Unexpected operation (M/H) in CIGAR {}", cigar));
            }
            let (cons_query, cons_ref) = op.consumes_query_ref();
            let qdiff = ifelse0(cons_query, item.len());
            let rdiff = ifelse0(cons_ref, item.len());

            let mut need_new = true;
            if let Some(last_var) = vars.last_mut() {
                if rpos <= last_var.ref_end && qpos <= last_var.hap_end {
                    last_var.ref_end = last_var.ref_end.max(rpos + rdiff);
                    last_var.hap_end = last_var.hap_end.max(qpos + qdiff);
                    need_new = false;
                } else {
                    // Need to check that the new variant does not overlap the previous one on both haplotypes.
                    assert!(rpos > last_var.ref_end && qpos > last_var.hap_end);
                }
            }

            if need_new {
                let var = if rdiff == qdiff {
                    VarRange::new(rpos, rpos + rdiff, qpos, qpos + qdiff)
                } else if rpos == 0 || qpos == 0 {
                    VarRange::new(rpos, rpos + rdiff + 1, qpos, qpos + qdiff + 1)
                } else {
                    VarRange::new(rpos - 1, rpos + rdiff, qpos - 1, qpos + qdiff)
                };
                vars.push(var);
            }
            rpos += rdiff;
            qpos += qdiff;
        }
        if let Some(var) = vars.last() {
            if var.ref_end > ref_seq.len() as u32 || var.hap_end > hap_seq.len() as u32 {
                return Err(error!(RuntimeError, "CIGAR operation out of range of the sequence"));
            }
        }
        move_all_left(&mut vars, ref_seq, hap_seq);
        let old_val = self.haplotypes.insert(hap_id, vars);
        if old_val.is_some() {
            log::warn!("Defined variants for the same haplotype twice");
        }
        Ok(())
    }

    /// Constructs variant ranges from PAF file with alignments between haplotypes and the reference haplotype.
    pub fn from_paf(
        paf_filename: &Path,
        contig_set: &ContigSet,
        ref_id: ContigId,
    ) -> crate::Result<Self> {
        let contigs = contig_set.contigs();
        let mut this = Self::default();
        let ref_hap = contigs.get_name(ref_id);
        let ref_seq = contig_set.get_seq(ref_id);
        this.haplotypes.insert(ref_id, Vec::new());

        let mut file = ext::sys::open(paf_filename).map(PafFile::new)?;
        while let Some(entry) = file.next_wo_tags(contigs)? {
            let mut entry = match entry {
                PafParseResult::Entry(entry) => entry,
                PafParseResult::UnknownContig(name) => {
                    log::warn!("Cannot find sequence for haplotype {}", name);
                    continue
                }
            };
            let (hap_id, invert) = if ref_id == entry.query_id() {
                (entry.target_id(), true)
            } else if ref_id == entry.target_id() {
                (entry.query_id(), false)
            } else {
                continue
            };
            if !entry.full_positive_alignment() {
                log::warn!("Alignment between {} and {} is on the reverse strand or does not fully cover both sequences",
                    contigs.get_name(hap_id), ref_hap);
                continue
            }

            let Some(mut cigar) = entry.take_cigar() else {
                log::warn!("CIGAR missing for {} and {}", ref_hap, contigs.get_name(hap_id));
                continue
            };
            if invert {
                cigar = cigar.invert();
            }
            let hap_seq = contig_set.get_seq(hap_id);
            if cigar.query_len() != hap_seq.len() as u32 || cigar.ref_len() != ref_seq.len() as u32 {
                log::error!("Incorrect CIGAR for {} and {} (expected lengths {},{}, got {},{})",
                    ref_hap, contigs.get_name(hap_id), cigar.query_len(), cigar.ref_len(), hap_seq.len(), ref_seq.len());
                continue;
            }
            this.add_haplotype_ranges(ref_seq, hap_id, hap_seq, &cigar)?;
        }
        if this.haplotypes.len() < contigs.len() {
            log::warn!("    {} / {} haplotype-reference alignments are missing",
                contigs.len() - this.haplotypes.len(), contigs.len());
        }
        Ok(this)
    }

    /// Returns vector (length = n_haps) of vectors (inner length = ref_ranges.len())
    /// where hap_ranges[hap_ix][var_ix] = (hap_start, hap_end) corresponding to ref_ranges[var_ix].
    /// Undefined haplotypes will received `None`.
    ///
    /// Also, individual haplotype ranges can be None when reference range starts or stops inside haplotype variant.
    pub fn get_hap_ranges(
        &self,
        ref_ranges: &[(u32, u32)],
    ) -> IntMap<ContigId, Vec<Option<(u32, u32)>>> {
        self.haplotypes.iter()
            .map(|(&id, hap_vars)| (id, get_hap_ranges(&ref_ranges, hap_vars)))
            .collect()
    }

    /// Returns
    /// - sorted and de-duplicated reference ranges (ref_start, ref_end),
    /// - merged reference ranges. Combined ranges can touch (end == next_start) but do not overlap.
    pub fn get_ref_ranges(&self) -> (Vec<(u32, u32)>, Vec<(u32, u32)>) {
        let mut unique_ranges: Vec<_> = self.haplotypes.values()
            .flat_map(|v| v.iter())
            .map(|v| (v.ref_start, v.ref_end)).collect();
        unique_ranges.sort_unstable();
        unique_ranges.dedup();

        let mut merged_ranges = Vec::<(u32, u32)>::new();
        for &(start, end) in &unique_ranges {
            if let Some(last) = merged_ranges.last_mut() {
                if last.1 <= start {
                    merged_ranges.push((start, end));
                } else {
                    last.1 = last.1.max(end);
                }
            } else {
                merged_ranges.push((start, end));
            }
        }
        (unique_ranges, merged_ranges)
    }

    pub fn generate_variants(
        &self,
        shift: u32,
        ref_ranges: &[(u32, u32)],
        ref_seq: &[u8],
        contig_set: &ContigSet,
        groupped_haps: &[(String, Vec<Option<ContigId>>)],
    ) -> crate::Result<Vec<Variant>> {
        let hap_ranges = self.get_hap_ranges(&ref_ranges);
        let ploidy = groupped_haps.iter().map(|(_, ids)| ids.len()).max().unwrap_or(0);
        let mut variants = Vec::with_capacity(ref_ranges.len());

        let mut hap_allele_ixs = IntMap::with_capacity_and_hasher(self.haplotypes.len(), Default::default());
        for (i, &(ref_start, ref_end)) in ref_ranges.iter().enumerate() {
            let mut var = Variant::new(shift + ref_start, shift + ref_end,
                ref_seq[ref_start as usize..ref_end as usize].to_vec(), true);
            hap_allele_ixs.clear();
            for (&id, curr_ranges) in hap_ranges.iter() {
                let Some((start, end)) = curr_ranges[i] else { continue };
                let allele_ix = var.add_allele(&contig_set.get_seq(id)[start as usize..end as usize]);
                hap_allele_ixs.insert(id, GenotypeAllele::Phased(allele_ix as i32));
            }
            if var.alleles().len() == 1 {
                // Can happen in separate.vcf.gz as some alleles become unknown.
                continue;
            }

            let genotypes = var.genotype_mut();
            for (_, hap_ids) in groupped_haps {
                genotypes.extend(hap_ids.iter().map(|opt_id|
                    opt_id.and_then(|id| hap_allele_ixs.get(&id).copied()).unwrap_or(GenotypeAllele::PhasedMissing)));
                const TRAILING_ALLELE: GenotypeAllele = GenotypeAllele::Phased((i32::MIN >> 1) - 1);
                genotypes.resize(genotypes.len() + ploidy - hap_ids.len(), TRAILING_ALLELE);
            }
            variants.push(var);
        }
        Ok(variants)
    }
}

/// Based on new reference ranges and old haplotype variants, get new haplotype ranges.
/// It is known that outside of the haplotype variants, reference and haplotype are identical.
/// Returns None when reference range start or stops inside haplotype variant.
fn get_hap_ranges(ref_ranges: &[(u32, u32)], hap_vars: &[VarRange]) -> Vec<Option<(u32, u32)>> {
    let n = hap_vars.len();
    if n == 0 {
        return ref_ranges.iter().copied().map(Some).collect();
    }

    let mut hap_ranges = Vec::with_capacity(n);
    for &(ref_start, ref_end) in ref_ranges {
        let ref_diff = ref_end - ref_start;
        let i = bisect::right_by(hap_vars, |var| var.ref_end.cmp(&ref_start));
        let j = bisect::left_by_at(hap_vars, |var| var.ref_start.cmp(&ref_end), i, n);

        if i == n {
            let last_var = hap_vars[n - 1].clone();
            assert!(ref_start >= last_var.ref_end);
            // Insert variant after the last haplotype variant.
            let shift = ref_start - last_var.ref_end;
            hap_ranges.push(Some((last_var.hap_end + shift, last_var.hap_end + shift + ref_diff)));
            continue;
        }

        let var1 = hap_vars[i].clone();
        if i == j {
            // No variants overlapped.
            assert!(ref_end <= var1.ref_start);
            let left_shift = var1.ref_start - ref_start;
            hap_ranges.push(Some((var1.hap_start.strict_sub(left_shift), var1.hap_start + ref_diff - left_shift)));
            continue;
        }

        let var2 = hap_vars[j - 1].clone();
        if ref_start <= var1.ref_start && var2.ref_end <= ref_end {
            let left_shift = var1.ref_start - ref_start;
            let right_shift = ref_end - var2.ref_end;
            hap_ranges.push(Some((var1.hap_start.strict_sub(left_shift), var2.hap_end + right_shift)));
        } else {
            hap_ranges.push(None);
        }
    }
    hap_ranges
}

pub(super) fn create_vcf_header<'a>(
    chrom: &str,
    length: Option<u32>,
    samples: impl Iterator<Item = &'a str>,
) -> Header {
    let mut header = Header::new();
    let contig_line = match length {
        Some(l) => format!("##contig=<ID={},length={}>", chrom, l),
        None => format!("##contig=<ID={}>", chrom),
    };
    header.push_record(contig_line.as_bytes());
    for sample in samples {
        header.push_sample(sample.as_bytes());
    }
    header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    header
}

pub(super) fn create_header_view(header: &bcf::Header) -> Arc<bcf::header::HeaderView> {
    // duplicates the header, so the copy is synced.
    let inner = unsafe { htslib::htslib::bcf_hdr_dup(header.inner_ptr()) };
    Arc::new(unsafe { bcf::header::HeaderView::from_ptr(inner) })
}

pub(super) fn write_vcf(
    vcf_filename: impl AsRef<Path>,
    header: &bcf::header::Header,
    variants: &[Variant],
) -> crate::Result<()>
{
    let vcf_filename = vcf_filename.as_ref();
    let tmp_filename = ext::sys::temp_filename(vcf_filename);
    let bgzip = vcf_filename.extension().and_then(OsStr::to_str) == Some("gz");
    let mut writer = bcf::Writer::from_path(&tmp_filename, header, !bgzip, bcf::Format::Vcf)?;
    for var in variants {
        let record = var.create_record(&writer)?;
        writer.write(&record)?;
    }
    std::mem::drop(writer);

    std::fs::rename(&tmp_filename, vcf_filename).map_err(add_path!(tmp_filename, vcf_filename))?;
    if bgzip {
        if let Err(e) = htslib::bcf::index::build(vcf_filename, None, 1, htslib::bcf::index::Type::Tbx) {
            log::error!("Could not construct tabix index for VCF file {}: {}", ext::fmt::path(vcf_filename), e);
        }
    }
    Ok(())
}

/// Groups haplotypes into samples.
/// Returns
/// - vector of samples (name, vector [contig index or None if missing]),
/// - reference contig ID.
fn group_haplotypes(
    contigs: &ContigNames,
    ref_hap: &str,
    disc_haps: &DiscardedHaplotypes,
) -> crate::Result<(Vec<(String, Vec<Option<ContigId>>)>, ContigId)>
{
    // First, regular contig name, with lazy *? specifier.
    // Then, single digit after either . or _
    const PATTERN: &'static str = r"^([0-9A-Za-z][0-9A-Za-z+._|~=@^-]*?)([._][1-9])?$";
    let re = regex::Regex::new(PATTERN).unwrap();

    let mut ref_id = None;
    let mut map = HashMap::<String, Vec<Option<ContigId>>>::default();
    let mut add = |id: ContigId, name: &str| -> crate::Result<()> {
        let c = re.captures(name).ok_or_else(|| error!(ParsingError, "Cannot parse contig name `{}`", name))?;
        let sample = &c[1];
        let opt_hap = c.get(2);
        if name == ref_hap {
            ref_id = Some(id);
            if let Some(hap) = &opt_hap {
                log::warn!("Reference name {} has haplotype suffix ({}); will keep it in the VCF file",
                    ref_hap, hap.as_str());
            } else {
                // Return from the inner function.
                return Ok(());
            }
        }
        let hap = opt_hap.map(|hap| hap.as_str().as_bytes()[1] - b'1').unwrap_or(0) as usize;
        let vec = map.entry(sample.to_string()).or_default();
        // If haplotype has .1 suffix, make vector length at least 2,
        // so that we have indication that the sample is diploid.
        let new_len = vec.len().max(hap + 1).max(if opt_hap.is_none() { 1 } else { 2 });
        vec.resize(new_len, None);
        vec[hap] = Some(id);
        Ok(())
    };

    for (id, contig) in itertools::izip!(contigs.ids(), contigs.names()) {
        add(id, contig)?;
        for (hap, _) in disc_haps.identical_to(id).into_iter().flat_map(|vec| vec.iter()) {
            add(id, hap)?;
        }
    }

    let Some(ref_id) = ref_id else {
        return Err(error!(InvalidInput,
            "Cannot find reference haplotype `{}` in neither the fasta file nor among the discarded haplotypes",
            ref_hap));
    };

    let mut groupped: Vec<_> = map.into_iter().collect();
    groupped.sort_unstable();
    Ok((groupped, ref_id))
}

pub(super) fn convert_to_vcf(
    paf_filename: &Path,
    contig_set: &ContigSet,
    disc_haps: &DiscardedHaplotypes,
    opt_ref_hap: Option<&str>,
    region: Option<ReferenceRegion>,
    out_merged: &Path,
    out_separate: Option<&Path>,
) -> crate::Result<()>
{
    if !disc_haps.all_identical() {
        log::warn!("Haplotypes were previously pruned (~ for some lines), VCF will be inaccurate");
    }
    let ref_hap = opt_ref_hap.or(region.as_ref().and_then(|region| region.ref_hap()))
        .ok_or_else(|| error!(InvalidInput, "Reference haplotype not found in the BED file and not provided with -r"))?;
    let (groupped_haps, ref_id) = group_haplotypes(contig_set.contigs(), ref_hap, &disc_haps)?;
    let ref_len = contig_set.contigs().get_len(ref_id);
    let (chrom, shift) = match region {
        Some(region) if region.len() != ref_len =>
            return Err(error!(InvalidData,
                "Region {}:{}-{} (len = {}) does not match saved reference haplotype (len = {})",
                region.chrom(), region.start() + 1, region.end(), region.len(), ref_len)),
        Some(region) => (region.chrom().to_string(), region.start()),
        None => (ref_hap.to_string(), 0),
    };

    let vars = VarRanges::from_paf(paf_filename, &contig_set, ref_id)?;
    let (unique_ranges, merged_ranges) = vars.get_ref_ranges();

    let header = create_vcf_header(&chrom, None, groupped_haps.iter().map(|(name, _)| name as &str));
    let ref_seq = contig_set.get_seq(ref_id);
    let merged_records = vars.generate_variants(shift, &merged_ranges, ref_seq, contig_set, &groupped_haps)?;
    write_vcf(&out_merged, &header, &merged_records)?;
    if let Some(separate_vcf_fname) = &out_separate {
        let sep_records = vars.generate_variants(shift, &unique_ranges, ref_seq, contig_set, &groupped_haps)?;
        write_vcf(&separate_vcf_fname, &header, &sep_records)?;
    }
    Ok(())
}

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let args = parse_args(argv)?.replace_input().validate()?;
    super::greet();
    let timer = Instant::now();

    let fasta_filename = args.fasta.as_ref().expect("Fasta file must be provided");
    let contig_set = ContigSet::load(String::new(), fasta_filename)?;

    let disc_filename = match &args.disc_filename as &str {
        "auto" => Some(Path::new(fasta_filename).parent()
            .ok_or_else(|| error!(RuntimeError, "Cannot get parent directory for {}", fasta_filename.display()))?
            .join(super::paths::DISCARDED_HAPS)),
        "none" => None,
        s => Some(PathBuf::from(s)),
    };
    let disc_haps = match disc_filename {
        Some(fname) if fname.exists() =>
            DiscardedHaplotypes::load(ext::sys::open(fname)?, contig_set.contigs())?,
        Some(fname) if args.disc_filename != "auto" => {
            log::debug!("Skip discarded haplotypes, file {} does not exist", fname.display());
            Default::default()
        }
        _ => Default::default(),
    };

    let region = load_region(&args.region, fasta_filename)?;
    let paf_filename = args.paf.as_ref().expect("PAF filename must be provided");
    convert_to_vcf(paf_filename, &contig_set, &disc_haps, args.ref_hap.as_ref().map(AsRef::as_ref), region,
        args.out_merged.as_ref().expect("Merged output VCF must be present"),
        args.out_separate.as_ref().map(|fname| fname as &Path))?;

    log::info!("Success! Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
