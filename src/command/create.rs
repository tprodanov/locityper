//! Create a new database.

use std::{
    fs,
    io::{Read, Seek},
    cmp::max,
    rc::Rc,
    process::Command,
    path::{Path, PathBuf},
};
use bio::io::fasta::IndexedReader;
use colored::Colorize;
use const_format::str_repeat;
use crate::{
    Error,
    seq::{Interval, ContigNames},
};

struct Args {
    reference: Option<PathBuf>,
    database: Option<PathBuf>,
    kmer_size: u8,
    bg_region: Option<String>,
    threads: u16,
    force: bool,
    jellyfish: PathBuf,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            reference: None,
            database: None,
            kmer_size: 25,
            bg_region: None,
            threads: 4,
            force: false,
            jellyfish: PathBuf::from("jellyfish"),
        }
    }
}

fn print_help() {
    const KEY: usize = 16;
    const VAL: usize = 4;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Create an EMPTY database of complex loci.".yellow());

    println!("\n{} {} create -d db -r reference.fa [arguments]",
        "Usage:".bold(), env!("CARGO_PKG_NAME"));

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Output database directory.",
        "-d, --db".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file.",
        "-r, --reference".green(), "FILE".yellow());

    println!("\n{}", "Optional parameters:".bold());
    println!("    {:KEY$} {:VAL$}  k-mer size [{}].",
        "-k, --kmer".green(), "INT".yellow(), defaults.kmer_size);
    println!("    {:KEY$} {:VAL$}  Calculate background distributions based on reads, mapped to this region.\n\
        {EMPTY}  Defaults to: chr17:72062001-76562000 (GRCh38).",
        "-b, --background".green(), "STR".yellow());
    // TODO: Parse default regions from BG_REGIONS.

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads);
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), "");
    println!("    {:KEY$} {:VAL$}  Jellyfish executable [{}].",
        "    --jellyfish".green(), "EXE".yellow(), defaults.jellyfish.display());

    println!("\n{}", "Other parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Show this help message.", "-h, --help".green(), "");
    println!("    {:KEY$} {:VAL$}  Show version.", "-V, --version".green(), "");
}

fn parse_args(argv: &[String]) -> Result<Args, lexopt::Error> {
    use lexopt::prelude::*;
    let mut args = Args::default();
    let mut parser = lexopt::Parser::from_args(argv);

    while let Some(arg) = parser.next()? {
        match arg {
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),

            Short('k') | Long("kmer") => args.kmer_size = parser.value()?.parse()?,
            Short('b') | Long("bg") | Long("background") => args.bg_region = Some(parser.value()?.parse()?),
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,
            Short('j') | Long("jellyfish") => args.jellyfish = parser.value()?.parse()?,

            Short('V') | Long("version") => {
                super::print_version();
                std::process::exit(0);
            }
            Short('h') | Long("help") => {
                print_help();
                std::process::exit(0);
            }
            _ => Err(arg.unexpected())?,
        }
    }
    Ok(args)
}

fn process_args(mut args: Args) -> Result<Args, Error> {
    args.threads = max(args.threads, 1);
    if args.database.is_none() {
        Err(lexopt::Error::from("Database directory is not provided (see -d/--database)"))?;
    }
    if args.reference.is_none() {
        Err(lexopt::Error::from("Reference fasta file is not provided (see -r/--reference)"))?;
    }
    args.jellyfish = super::find_exe(args.jellyfish)?;
    Ok(args)
}

/// Calculate background read distributions based on one of these regions.
const BG_REGIONS: [&'static str; 1] = ["chr17:72062001-76562000"];

/// Returns the first appropriate interval for the contig names.
/// If `bg_region` is set, try to parse it. Otherwise, iterate over `BG_REGIONS`.
///
/// Returns error if no interval is appropriate (chromosome not in the contig set, or interval is out of bounds).
fn select_bg_interval(
    ref_filename: &Path,
    contigs: &Rc<ContigNames>,
    bg_region: &Option<String>
) -> Result<Interval, Error>
{
    if let Some(s) = bg_region {
        let region = Interval::parse(s, contigs).map_err(|_| Error::InvalidInput(
            format!("Reference file '{}' does not contain region {}", ref_filename.display(), s)))?;
        return if contigs.in_bounds(&region) {
            Ok(region)
        } else {
            Err(Error::InvalidInput(format!("Region {} is out of bounds", s)))
        };
    }

    for s in BG_REGIONS.iter() {
        if let Ok(region) = Interval::parse(s, contigs) {
            if contigs.in_bounds(&region) {
                return Ok(region);
            } else {
                log::error!("Chromosome {} is in the reference file '{}', but is shorter than expected",
                    ref_filename.display(), region.contig_name());
            }
        }
    }
    Err(Error::InvalidInput(format!(
        "Reference file '{}' does not contain any of the default background regions. \
        Consider setting region via --region or use a different reference file.", ref_filename.display())))
}

/// Extracts the sequence of a background region, used to estimate the parameters of the sequencing data.
/// The sequence is then written to the `$bg_path/bg.fa.gz`.
fn extract_bg_region<R: Read + Seek>(
    fasta: &mut IndexedReader<R>,
    db_path: &Path,
    region: &Interval,
) -> Result<(), Error>
{
    let seq = region.fetch_seq(fasta)?;
    let mut bg_path = db_path.join("bg");
    fs::create_dir(&bg_path)?;
    bg_path.push("bg.fa.gz");
    log::info!("Writing background region {} to '{}'", region, bg_path.display());
    let mut fasta_writer = super::common::create_gzip(&bg_path)?;
    crate::seq::write_fasta(&mut fasta_writer, "bg", Some(&region.to_string()), &seq)?;
    Ok(())
}

fn run_jellyfish(
    db_path: &Path,
    ref_filename: &Path,
    jellyfish: &Path,
    kmer_size: u8,
    threads: u16,
    genome_size: u64,
) -> Result<(), Error> {
    let mut jf_path = db_path.join("jf");
    fs::create_dir(&jf_path)?;
    jf_path.push(&format!("{}.jf", kmer_size));

    let mut command = Command::new(jellyfish);
    command
        .args(&["count", "--canonical", "--lower-count=2", "--out-counter-len=1",
            &format!("--mer-len={}", kmer_size),
            &format!("--threads={}", threads),
            &format!("--size={}", genome_size)])
        .arg("--output").arg(&jf_path)
        .arg(ref_filename);
    log::info!("Counting {}-mers in {} threads, output: {}", kmer_size, threads, jf_path.display());
    log::trace!("    {:?}", command);

    let output = command.output()?;
    if output.status.success() {
        Ok(())
    } else {
        Err(Error::SubcommandFail(output))
    }
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args: Args = process_args(parse_args(argv)?)?;
    // unwrap as args.database was previously checked to be Some.
    let db_path = args.database.as_ref().unwrap();
    if db_path.exists() {
        if args.force {
            log::warn!("Completely removing output directory '{}'", db_path.display());
            fs::remove_dir_all(db_path)?;
        } else {
            panic!("Output directory '{}' already exists. Remove it or use -F/--force flag.", db_path.display());
        }
    }
    fs::create_dir(db_path)?;

    // unwrap as args.reference was previously checked to be Some.
    let ref_filename = args.reference.as_ref().unwrap();
    let (contigs, mut fasta) = ContigNames::load_indexed_fasta(&ref_filename, "reference".to_string())?;
    let region = select_bg_interval(&ref_filename, &contigs, &args.bg_region)?;

    extract_bg_region(&mut fasta, &db_path, &region)?;
    run_jellyfish(&db_path, &ref_filename, &args.jellyfish, args.kmer_size, args.threads, contigs.genome_size())?;
    fs::create_dir(&db_path.join("loci"))?;
    log::info!("Success!");
    Ok(())
}
