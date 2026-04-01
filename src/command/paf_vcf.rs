use std::{
    path::{Path, PathBuf},
    time::Instant,
    io::{Write, BufRead},
    ffi::OsStr,
};
use colored::Colorize;
use regex::Regex;
use crate::{
    err::{error, add_path, validate_param},
    ext::{
        self,
        fmt::PrettyU32,
    },
    seq::{
        interv,
        cigar::{Cigar, Operation},
        contigs::{ContigId, ContigNames, ContigSet},
    },
    algo::{
        bisect,
        HashMap, IntMap,
    },
    command::prune,
};

struct Args {
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
    fn validate(self) -> crate::Result<Self> {
        validate_param!(self.paf.is_some(), "Input PAF file is not provided (see -p/--paf)");
        validate_param!(self.fasta.is_some(), "Input FASTA file is not provided (see -f/--fasta)");
        validate_param!(self.out_merged.is_some(), "Output VCF file must be provided (see -o/--output)");
        validate_param!(self.ref_hap.is_some(), "Reference haplotype must be provided");
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 16;
    const VAL: usize = 4;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    // let defaults = Args::default();
    println!("{}", "Convert PAF file to VCF file(s).".yellow());

    print!("\n{}", "Usage:".bold());
    println!(" {} -p input.paf -f input.fa -r ref_hap [-d discarded_haplotypes.txt] \\", super::PROGRAM);
    println!("    -o merged.vcf.gz [separate.vcf.gz] [args]");

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Input PAF[.gz] file with alignments between haplotypes.",
        "-p, --paf".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Input FASTA[.gz] file with haplotypes.",
        "-f, --fasta".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Reference haplotypes name.",
        "-r, --ref-hap".green(), "STR".yellow());
    println!("    {:KEY$} {}  Text file with discarded haplotypes [{}].\n\
        {EMPTY}  {} = (DIRNAME of {})/discarded_haplotypes.txt.",
        "-d, --discarded".green(), "FILE|auto|none".yellow(), super::fmt_def("auto"),
        "auto".yellow(), "-f".green());
    println!("    {:KEY$} {:VAL$}  Two output VCF[.gz] files, first with merged variants\n\
        {EMPTY}  and second (optional) with unmerged variants.",
        "-o, --output".green(), "FILE [FILE]".yellow());

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Adjust variant positions relative to this region [{}].\n\
        {EMPTY}  Should be either {} (take from ref.bed), {},\n\
        {EMPTY}  `chrom:start`, `chrom:start-end` or a BED file with one entry.",
        "-R, --region".green(), "STR".yellow(), super::fmt_def("auto"), "auto".yellow(), "none".yellow());

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

/// Parse region (either chrom:start, chrom:start-end, or BED file with a single entry).
/// Returns chromosome name and 0-based start position.
fn load_region(s: &str, fasta_filename: &Path) -> crate::Result<Option<(String, u32)>> {
    if s == "none" { return Ok(None) };

    let filename = if s == "auto" {
        Path::new(fasta_filename).parent()
            .ok_or_else(|| error!(RuntimeError, "Cannot get parent directory for {}", s))?
            .join(super::paths::LOCUS_BED)
    } else {
        let pos_re = Regex::new(interv::CHROM_POS_PATTERN).unwrap();
        let interval_re = Regex::new(interv::INTERVAL_PATTERN).unwrap();
        if let Some(captures) = pos_re.captures(s).or_else(|| interval_re.captures(s)) {
            let chrom = captures[1].to_string();
            let pos: PrettyU32 = captures[2].parse()
                .map_err(|_| error!(ParsingError, "Cannot parse interval '{}'", s))?;
            return Ok(Some((chrom, pos.get().strict_sub(1))));
        }
        PathBuf::from(s)
    };

    if !filename.exists() {
        log::error!("Cannot find BED file {}, using relative locus coordinates", filename.display());
        return Ok(None);
    }

    let mut lines = ext::sys::open(&filename)?.lines();
    let line = lines.next().ok_or_else(|| error!(InvalidInput, "BED file {} is empty", filename.display()))?
        .map_err(add_path!(filename))?;
    let trimmed_line = line.trim();
    let mut split = trimmed_line.split('\t');
    let chrom = split.next().ok_or_else(|| error!(InvalidInput, "BED file {} is empty", filename.display()))?
        .to_string();
    let pos: u32 = split.next()
        .ok_or_else(|| error!(InvalidInput, "Not enough columns in BED file {} ({})",
            filename.display(), trimmed_line))?
        .parse().map_err(|_| error!(ParsingError, "Cannot parse line `{}` in BED file {}",
            trimmed_line, filename.display()))?;
    Ok(Some((chrom, pos)))
}

#[derive(Clone, Debug)]
struct VarRange {
    ref_start: u32,
    ref_end: u32,
    hap_start: u32,
    hap_end: u32,
}

impl VarRange {
    #[inline]
    fn new(ref_start: u32, ref_end: u32, hap_start: u32, hap_end: u32) -> Self {
        Self { ref_start, ref_end, hap_start, hap_end }
    }

    #[inline]
    fn ref_len(&self) -> u32 {
        self.ref_end - self.ref_start
    }

    #[inline]
    fn alt_len(&self) -> u32 {
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

/// Convert one haplotype into VCF.
/// `start` represent shift of the variants relative to the reference haplotype.
/// Returns a vector of variant ranges.
fn process_haplotype(
    ref_seq: &[u8],
    hap_seq: &[u8],
    cigar: &Cigar,
) -> crate::Result<Vec<VarRange>>
{
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
        let rdiff = u32::from(op.consumes_ref()) * item.len();
        let qdiff = u32::from(op.consumes_query()) * item.len();

        let mut need_new = true;
        if let Some(last_var) = vars.last_mut() {
            if last_var.ref_end == rpos && last_var.hap_end == qpos {
                // Current indel is preceded by a mismatch.
                last_var.ref_end = rpos + rdiff;
                last_var.hap_end = qpos + qdiff;
                need_new = false;
            } else {
                assert!(last_var.ref_end < rpos && last_var.hap_end < qpos);
            }
        }

        if need_new {
            let var = if rdiff == qdiff {
                VarRange::new(rpos,      rpos + rdiff,      qpos,      qpos + qdiff)
            } else if rpos == 0 || qpos == 0 {
                VarRange::new(rpos,      rpos + rdiff + 1,  qpos,      qpos + qdiff + 1)
            } else {
                VarRange::new(rpos - 1,  rpos + rdiff,      qpos - 1,  qpos + qdiff)
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
    Ok(vars)
}

/// Creates either bgzip or plain output file.
fn create_vcf_writer<'a>(
    filename: impl AsRef<Path>,
    samples: impl Iterator<Item = &'a str>,
) -> crate::Result<Box<dyn Write>>
{
    let filename = filename.as_ref();
    let mut out = if filename.extension().and_then(OsStr::to_str) == Some("gz") {
        htslib::bgzf::Writer::from_path(&filename)
            .map_err(|_| error!(RuntimeError, "Cannot create file {}", filename.display()))
            .map(|w| Box::new(w) as Box<dyn Write>)?
    } else {
        ext::sys::create_file(filename)
            .map(|w| Box::new(w) as Box<dyn Write>)?
    };
    const HEADER: &'static [u8] = b"\
        ##fileformat=VCFv4.2\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    out.write_all(HEADER).map_err(add_path!(filename))?;
    for sample in samples {
        write!(out, "\t{}", sample).map_err(add_path!(filename))?;
    }
    writeln!(out).map_err(add_path!(filename))?;
    Ok(out)
}

/// Reads PAF file and returns a vector of variants for each haplotype found in the contig set and in the PAF file.
fn process_paf(
    paf_filename: &Path,
    contig_set: &ContigSet,
    ref_id: ContigId,
) -> crate::Result<Vec<Option<Vec<VarRange>>>>
{
    let contigs = contig_set.contigs();
    let ref_hap = contigs.get_name(ref_id);
    let ref_seq = contig_set.get_seq(ref_id);
    let mut var_ranges = vec![None; contigs.len()];
    var_ranges[ref_id.ix()] = Some(Vec::new());

    let file = ext::sys::open(paf_filename)?;
    for line in file.lines() {
        let line = line.map_err(add_path!(paf_filename))?;
        if line.starts_with('#') { continue };
        let trimmed = line.trim();
        let split: Vec<_> = trimmed.split('\t').collect();
        if split.len() < 12 {
            return Err(error!(InvalidInput, "PAF line ({}) has too few columns", trimmed));
        }

        let hap1: &str = &split[0];
        let hap2: &str = &split[5];
        let (hap, invert) = if ref_hap == hap1 {
            (hap2, true)
        } else if ref_hap == hap2 {
            (hap1, false)
        } else {
            continue
        };
        let Some(hap_id) = contigs.try_get_id(hap) else {
            log::warn!("Cannot find sequence for haplotype {}", hap);
            continue
        };
        let hap_seq = contig_set.get_seq(hap_id);

        let hap_column = 5 * usize::from(invert);
        let ref_column = 5 * usize::from(!invert);
        let hap_start: u32 = split[hap_column + 2].parse()
            .map_err(|_| error!(ParsingError, "Cannot parse PAF line {}", trimmed))?;
        let hap_end: u32 = split[hap_column + 3].parse()
            .map_err(|_| error!(ParsingError, "Cannot parse PAF line {}", trimmed))?;
        let ref_start: u32 = split[ref_column + 2].parse()
            .map_err(|_| error!(ParsingError, "Cannot parse PAF line {}", trimmed))?;
        let ref_end: u32 = split[ref_column + 3].parse()
            .map_err(|_| error!(ParsingError, "Cannot parse PAF line {}", trimmed))?;

        if ref_start != 0 || ref_end != ref_seq.len() as u32 || hap_start != 0 || hap_end != hap_seq.len() as u32 {
            log::warn!("Alignment between {} and {} does not cover the full sequence", ref_hap, hap);
            continue;
        }
        let mut cigar = None;
        for col in &split[12..] {
            if col.starts_with("cg:Z:") {
                cigar = Some(Cigar::from_str(&col.as_bytes()[5..])?);
                break;
            }
        }
        let Some(mut cigar) = cigar else {
            log::warn!("CIGAR missing for {} and {}", ref_hap, hap);
            continue
        };

        if invert {
            cigar = cigar.invert();
        }
        if cigar.query_len() != hap_end || cigar.ref_len() != ref_end {
            log::error!("Incorrect CIGAR for {} and {} (expected lengths {},{}, got {},{})", ref_hap, hap,
                cigar.query_len(), cigar.ref_len(), hap_end, ref_end);
            continue;
        }
        var_ranges[hap_id.ix()] = Some(process_haplotype(ref_seq, hap_seq, &cigar)?);
    }
    Ok(var_ranges)
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

fn write_vcf(
    chrom: &str,
    shift: u32,
    ref_ranges: &[(u32, u32)],
    vars: &[Option<Vec<VarRange>>],
    contig_set: &ContigSet,
    ref_id: ContigId,
    groupped_haps: &[(String, Vec<Option<usize>>)],
    vcf_filename: &Path,
) -> crate::Result<()>
{
    let hap_ranges: Vec<_> = vars.iter()
        .map(|hap_vars| hap_vars.as_ref().map(|vars| get_hap_ranges(&ref_ranges, vars))
            .unwrap_or_else(|| vec![None; ref_ranges.len()]))
        .collect();

    let ref_seq = contig_set.get_seq(ref_id);
    let mut writer = create_vcf_writer(vcf_filename, groupped_haps.iter().map(|(name, _)| name as &str))?;

    for (i, &(ref_start, ref_end)) in ref_ranges.iter().enumerate() {
        let mut alleles = vec![&ref_seq[ref_start as usize..ref_end as usize]];
        let hap_allele_ixs: Vec<Option<usize>> = contig_set.seqs().iter().zip(&hap_ranges)
            .map(|(seq, hap)| hap[i].map(|(start, end)| {
                let allele = &seq[start as usize..end as usize];
                match alleles.iter().position(|&existing_allele| existing_allele == allele) {
                    Some(allele_ix) => allele_ix,
                    None => {
                        alleles.push(allele);
                        alleles.len() - 1
                    }
                }
        })).collect();
        if alleles.len() == 1 {
            // Can happen in separate.vcf.gz as some alleles become unknown.
            continue;
        }

        write!(writer, "{}\t{}\t.", chrom, ref_start + shift + 1).map_err(add_path!(vcf_filename))?;
        for (i, allele) in alleles.iter().enumerate() {
            writer.write_all(if i <= 1 { b"\t" } else { b"," }).map_err(add_path!(vcf_filename))?;
            writer.write_all(allele).map_err(add_path!(vcf_filename))?;
        }
        writer.write_all(b"\t60\t.\t.\tGT").map_err(add_path!(vcf_filename))?;

        for (_, hap_ixs) in groupped_haps {
            for (i, &opt_j) in hap_ixs.iter().enumerate() {
                writer.write_all(if i == 0 { b"\t" } else { b"|" }).map_err(add_path!(vcf_filename))?;
                match opt_j.and_then(|j| hap_allele_ixs[j]) {
                    Some(allele_ix) => write!(writer, "{}", allele_ix),
                    None => writer.write_all(b"."),
                }.map_err(add_path!(vcf_filename))?;
            }
        }
        writeln!(writer).map_err(add_path!(vcf_filename))?;
    }
    Ok(())
}

/// Merge overlapping variants and return:
/// - Reference ranges,
/// - Corresponding haplotype ranges.
/// All haplotype vectors have the same length as the the reference vector.
fn combine_variants(
    chrom: &str,
    shift: u32,
    vars: &[Option<Vec<VarRange>>],
    contig_set: &ContigSet,
    ref_id: ContigId,
    groupped_haps: &[(String, Vec<Option<usize>>)],
    merged_vcf_filename: &Path,
    separate_vcf_filename: Option<&Path>,
) -> crate::Result<()>
{
    let mut unique_ranges: Vec<_> = vars.iter()
        .filter_map(|v| v.as_ref())
        .flat_map(|v| v.iter())
        .map(|v| (v.ref_start, v.ref_end)).collect();
    // Sort increasing by start, decreasing by end.
    unique_ranges.sort_unstable();
    unique_ranges.dedup();

    // Get merged reference ranges. Do not merge variants that touch, but don't overlap (prev_end == next_start).
    let mut merged_ref_ranges = Vec::<(u32, u32)>::new();
    for &(start, end) in &unique_ranges {
        if let Some(last) = merged_ref_ranges.last_mut() {
            if last.1 <= start {
                merged_ref_ranges.push((start, end));
            } else {
                last.1 = last.1.max(end);
            }
        } else {
            merged_ref_ranges.push((start, end));
        }
    }

    write_vcf(chrom, shift, &merged_ref_ranges, &vars, contig_set, ref_id, groupped_haps, merged_vcf_filename)?;

    if let Some(separate_vcf_fname) = separate_vcf_filename {
        write_vcf(chrom, shift, &unique_ranges, &vars, contig_set, ref_id, groupped_haps, separate_vcf_fname)?;
    }
    Ok(())
}

/// Groups haplotypes into samples.
/// Returns vector of samples (name, vector [contig index or None if missing]).
fn group_haplotypes(
    contigs: &ContigNames,
    ref_id: ContigId,
    disc_haps: &IntMap<ContigId, Vec<(String, bool)>>,
) -> crate::Result<Vec<(String, Vec<Option<usize>>)>>
{
    // First, regular contig name, with lazy *? specifier.
    // Then, single digit after either . or _
    const PATTERN: &'static str = r"^([0-9A-Za-z][0-9A-Za-z+._|~=@^-]*?)([._][1-9])?$";
    let re = regex::Regex::new(PATTERN).unwrap();

    let mut map = HashMap::<String, Vec<Option<usize>>>::default();
    let mut add = |i: usize, name: &str| -> crate::Result<()> {
        let c = re.captures(name).ok_or_else(|| error!(ParsingError, "Cannot parse contig name `{}`", name))?;
        let sample = &c[1];
        let opt_hap = c.get(2);
        let hap = opt_hap.map(|hap| hap.as_str().as_bytes()[1] - b'1').unwrap_or(0) as usize;
        let vec = map.entry(sample.to_string()).or_default();
        // If haplotype has .1 suffix, make vector length at least 2,
        // so that we have indication that the sample is diploid.
        let new_len = vec.len().max(hap + 1).max(if opt_hap.is_none() { 1 } else { 2 });
        vec.resize(new_len, None);
        vec[hap] = Some(i);
        Ok(())
    };

    for (i, contig) in contigs.names().iter().enumerate() {
        if i != ref_id.ix() {
            add(i, contig)?;
        }
        for (hap, _) in disc_haps.get(&ContigId::new(i)).into_iter().flat_map(|vec| vec.iter()) {
            add(i, hap)?;
        }
    }
    let mut groupped: Vec<_> = map.into_iter().collect();
    groupped.sort_unstable();
    Ok(groupped)
}

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();

    let fasta_filename = args.fasta.as_ref().expect("Fasta file must be provided");
    let contig_set = ContigSet::load(String::new(), fasta_filename)?;
    let ref_hap = args.ref_hap.clone().expect("Reference haplotype must be provided");
    let ref_id = contig_set.contigs().try_get_id(&ref_hap)
        .ok_or_else(|| error!(InvalidInput, "Reference haplotype {} not found in the fasta file", ref_hap))?;

    let disc_filename = match &args.disc_filename as &str {
        "auto" => Some(Path::new(fasta_filename).parent()
            .ok_or_else(|| error!(RuntimeError, "Cannot get parent directory for {}", fasta_filename.display()))?
            .join(super::paths::DISCARDED_HAPS)),
        "none" => None,
        s => Some(PathBuf::from(s)),
    };
    let disc_haps = match disc_filename {
        Some(fname) if fname.exists() =>
            prune::load_discarded_haplotypes(ext::sys::open(fname)?, contig_set.contigs())?,
        Some(fname) => {
            log::debug!("Skip discarded haplotypes, file {} does not exist", fname.display());
            Default::default()
        }
        None => Default::default()
    };
    if !prune::all_identical(&disc_haps) {
        log::warn!("Haplotypes were previously pruned (~ for some lines), VCF will be inaccurate");
    }
    let groupped_haps = group_haplotypes(contig_set.contigs(), ref_id, &disc_haps)?;

    let (chrom, shift) = load_region(&args.region, fasta_filename)?.unwrap_or((ref_hap.clone(), 0));
    let vars = process_paf(args.paf.as_ref().expect("PAF filename must be provided"), &contig_set, ref_id)?;
    combine_variants(&chrom, shift, &vars, &contig_set, ref_id, &groupped_haps,
        args.out_merged.as_ref().expect("Merged output VCF must be present"),
        args.out_separate.as_ref().map(|fname| fname as &Path))?;

    log::info!("Success! Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
