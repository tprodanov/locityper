use std::{
    path::{Path, PathBuf},
    time::Instant,
    io::{Write, BufRead},
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
        interv, fastx,
        cigar::{Cigar, Operation},
    },
    algo::{HashSet, HashMap},
};

struct Args {
    paf: Option<PathBuf>,
    fasta: Option<PathBuf>,
    out_fmt: Option<String>,
    ref_hap: Option<String>,
    region: Option<String>,
    haplotypes: Option<PathBuf>,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            paf: None,
            fasta: None,
            out_fmt: None,
            ref_hap: None,
            region: None,
            haplotypes: None,
        }
    }
}

impl Args {
    fn validate(self) -> crate::Result<Self> {
        validate_param!(self.paf.is_some(), "Input PAF file is not provided (see -p/--paf)");
        validate_param!(self.fasta.is_some(), "Input FASTA file is not provided (see -f/--fasta)");
        validate_param!(self.out_fmt.is_some(), "Output VCF path is not provided (see -o/--output)");
        validate_param!(self.out_fmt.as_ref().map(|s| s.contains("{}")).unwrap_or(false),
            "Output (-o/--output) must contain {{}}, where haplotype name is inserted");
        validate_param!(self.ref_hap.is_some(), "Reference haplotype must be provided");
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 16;
    const VAL: usize = 5;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    // let defaults = Args::default();
    println!("{}", "Convert PAF file to haplotype-level VCF files.".yellow());

    print!("\n{}", "Usage:".bold());
    println!(" {} -p input.paf -f input.fa -o out.{{}}.vcf.gz -r ref_hap [args]", super::PROGRAM);

    println!("\n{}", "Required arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Input PAF[.gz] file with alignments between haplotypes.",
        "-p, --paf".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Input FASTA[.gz] file with haplotypes.",
        "-f, --fasta".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Format string for output VCF.gz files,\n\
        {EMPTY}  where haplotype name will be inserted in {{}}.\n\
        {EMPTY}  Tries to create parent directory only if it does not contain {{}}.",
        "-o, --output".green(), "STR".yellow());
    println!("    {:KEY$} {:VAL$}  Reference haplotypes name.",
        "-r, --ref-hap".green(), "STR".yellow());

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Adjust variant positions relative to this region.\n\
        {EMPTY}  Either `chrom:start`, `chrom:start-end` or a BED file with one entry.",
        "-R, --region".green(), "STR".yellow());
    println!("    {:KEY$} {:VAL$}  Only create VCF for haplotypes from this file.",
        "-H, --haplotypes".green(), "FILE".yellow());

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
            Short('o') | Long("output") => args.out_fmt = Some(parser.value()?.parse()?),
            Short('r') | Long("ref-hap") => args.ref_hap = Some(parser.value()?.parse()?),

            Short('R') | Long("region") => args.region = Some(parser.value()?.parse()?),
            Short('H') | Long("haplotypes") => args.haplotypes = Some(parser.value()?.parse()?),

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
fn load_region(s: &str) -> crate::Result<(String, u32)> {
    let pos_re = Regex::new(interv::CHROM_POS_PATTERN).unwrap();
    let interval_re = Regex::new(interv::INTERVAL_PATTERN).unwrap();
    if let Some(captures) = pos_re.captures(s).or_else(|| interval_re.captures(s)) {
        let chrom = captures[1].to_string();
        let pos: PrettyU32 = captures[2].parse()
            .map_err(|_| error!(ParsingError, "Cannot parse interval '{}'", s))?;
        return Ok((chrom, pos.get().checked_sub(1).unwrap()));
    }

    let mut lines = ext::sys::open(s)?.lines();
    let line = lines.next().ok_or_else(|| error!(InvalidInput, "BED file {} is empty", s))?
        .map_err(add_path!(s))?;
    let trimmed_line = line.trim();
    let mut split = trimmed_line.split('\t');
    let chrom = split.next().ok_or_else(|| error!(InvalidInput, "BED file {} is empty", s))?.to_string();
    let pos: u32 = split.next()
        .ok_or_else(|| error!(InvalidInput, "Not enough columns in BED file {} ({})", s, trimmed_line))?
        .parse().map_err(|_| error!(ParsingError, "Cannot parse line `{}` in BED file {}", trimmed_line, s))?;
    Ok((chrom, pos))
}

/// Loads lines into HashSet.
fn load_haplotypes(filename: &Path) -> crate::Result<HashSet<String>> {
    ext::sys::open(filename)?
        .lines()
        .map(|line_or_err| line_or_err.map(|line| line.trim().to_string()).map_err(add_path!(filename)))
        .collect()
}

fn load_fasta(filename: &Path) -> crate::Result<HashMap<String, Vec<u8>>> {
    use fastx::SingleRecord;
    let mut fasta_reader = fastx::Reader::from_path(filename)?;
    let mut seqs = HashMap::default();
    let mut record = fastx::FastxRecord::default();
    while fasta_reader.read_next_standardized(&mut record)? {
        let name = String::from_utf8(record.name().to_vec())
            .map_err(|_| crate::Error::Utf8("read name", record.name().to_vec()))?;
        seqs.insert(name, record.seq().to_owned());
    }
    Ok(seqs)
}

/// Convert one haplotype into VCF.
/// `start` represent shift of the variants relative to the reference haplotype.
fn process_haplotype(
    chrom: &str,
    start: u32,
    ref_seq: &[u8],
    hap_seq: &[u8],
    cigar: &Cigar,
    out: &mut impl Write,
) -> crate::Result<()>
{
    fn checked_get(seq: &[u8], start: u32, end: u32) -> crate::Result<&[u8]> {
        if end as usize > seq.len() {
            Err(error!(RuntimeError, "CIGAR operation out of range of the sequence"))
        } else {
            Ok(&seq[start as usize..end as usize])
        }
    }

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

        let (pos, ref_allele, alt_allele) = if rdiff == qdiff {
            (rpos, checked_get(ref_seq, rpos, rpos + rdiff)?, checked_get(hap_seq, qpos, qpos + qdiff)?)
        } else {
            if rpos == 0 || qpos == 0 {
                (rpos, checked_get(ref_seq, rpos, rpos + rdiff + 1)?, checked_get(hap_seq, qpos, qpos + qdiff + 1)?)
            } else {
                (rpos - 1, checked_get(ref_seq, rpos - 1, rpos + rdiff)?,
                    checked_get(hap_seq, qpos - 1, qpos + qdiff)?)
            }
        };

        write!(out, "{}\t{}\t.\t", chrom, start + pos + 1).map_err(add_path!(!))?;
        out.write_all(ref_allele).map_err(add_path!(!))?;
        out.write_all(b"\t").map_err(add_path!(!))?;
        out.write_all(alt_allele).map_err(add_path!(!))?;
        out.write_all(b"\t60\t.\t.\tGT\t1\n").map_err(add_path!(!))?;
        rpos += rdiff;
        qpos += qdiff;
    }
    Ok(())
}

#[inline]
fn write_header(
    sample: &str,
    out: &mut impl Write,
) -> std::io::Result<()>
{
    const HEADER: &'static [u8] = b"\
        ##fileformat=VCFv4.2\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
    out.write_all(HEADER)?;
    writeln!(out, "{}", sample)?;
    Ok(())
}

fn process_paf(
    paf_filename: &Path,
    chrom: &str,
    start: u32,
    ref_hap: &str,
    ref_seq: &[u8],
    seqs: &HashMap<String, Vec<u8>>,
    only_haplotypes: Option<&HashSet<String>>,
    out_fmt: &str,
) -> crate::Result<()>
{
    let mut processed = 0;
    let mut cigar_missing = 0;
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
        if !only_haplotypes.map(|set| set.contains(hap)).unwrap_or(true) { continue };

        let Some(hap_seq) = seqs.get(hap) else {
            log::warn!("Cannot find sequence for haplotype {}", hap);
            continue
        };
        let out_filename = out_fmt.replace("{}", hap);
        // Cannot cast to crate IO error as we get htslib::Error
        let mut out_file = htslib::bgzf::Writer::from_path(&out_filename)
            .map_err(|_| error!(RuntimeError, "Cannot create file {}", out_filename))?;
        write_header(hap, &mut out_file).map_err(add_path!(out_filename))?;

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

        if ref_start != 0 || ref_end as usize != ref_seq.len() || hap_start != 0 || hap_end as usize != hap_seq.len() {
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
            cigar_missing += 1;
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
        process_haplotype(chrom, start, ref_seq, hap_seq, &cigar, &mut out_file)?;
        processed += 1;
    }
    if cigar_missing > 0 {
        log::info!("Processed {} haplotypes, CIGAR missing for {} haplotypes", processed, cigar_missing);
    } else {
        log::info!("Processed {} haplotypes", processed);
    }
    Ok(())
}

// cigar.invert()

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();

    let seqs = load_fasta(args.fasta.as_ref().expect("Fasta file must be provided"))?;
    let ref_hap = args.ref_hap.clone().expect("Reference haplotype must be provided");
    let ref_seq = seqs.get(&ref_hap)
        .ok_or_else(|| error!(InvalidInput, "Reference haplotype {} not found in the fasta file", ref_hap))?
        .clone();

    let region = args.region.as_ref().map(|s| load_region(s)).transpose()?
        .unwrap_or((ref_hap.clone(), 0));
    let only_haplotypes = args.haplotypes.as_ref().map(|filename| load_haplotypes(filename)).transpose()?;

    let out_fmt = args.out_fmt.as_ref().expect("Output path must be provided");
    if let Some(dirname) = Path::new(out_fmt).parent() {
        if !dirname.exists() && !dirname.to_string_lossy().contains("{}") {
            std::fs::create_dir(dirname).map_err(add_path!(dirname))?;
        }
    }
    process_paf(
        args.paf.as_ref().expect("PAF filename must be provided"),
        &region.0, region.1, &ref_hap, &ref_seq, &seqs, only_haplotypes.as_ref(), out_fmt,
    )?;

    log::info!("Success! Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
