//! Create a new database.

use std::{
    fs,
    io::{self, Read, Seek},
    cmp::max,
    rc::Rc,
    process::Command,
    path::{Path, PathBuf},
};
use bio::io::fasta::IndexedReader;
use colored::Colorize;
use crate::{
    Error,
    seq::{Interval, ContigNames},
};

struct Args {
    fasta: Option<PathBuf>,
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
            fasta: None,
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
    let defaults = Args::default();
    println!("{}", "Create an empty database of complex loci.".yellow());
    println!("Counts k-mers and extracts region for background parameters calculation.");

    println!("\n{} {} create -d db -f genome.fa [arguments]",
        "Usage:".bold(), env!("CARGO_PKG_NAME"));

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:<15} {:<6}  Output database directory.",
        "-d, --db".green(), "DIR".yellow());
    println!("    {:<15} {:<6}  Input FASTA file. Used for k-mer counting and must contain --region.",
        "-f, --fasta".green(), "FILE".yellow());

    println!("\n{}", "Optional parameters:".bold());
    println!("    {:<15} {:<6}  k-mer size [{}].",
        "-k, --kmer".green(), "INT".yellow(), defaults.kmer_size);
    println!("    {:<15} {:<6}  Calculate background distributions based on reads, mapped to this region.\n\
        {:<26}  Defaults to: chr17:72062001-76562000 (GRCh38).",
        "-r, --region".green(), "REGION".yellow(), "");
    // TODO: Write default regions from some array.

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:<15} {:<6}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads);
    println!("    {:<15} {:<6}  Force rewrite output directory.",
        "-F, --force".green(), "");
    println!("    {:<15} {:<6}  Jellyfish executable [{}].",
        "    --jellyfish".green(), "EXE".yellow(), defaults.jellyfish.display());

    println!("\n{}", "Other parameters:".bold());
    println!("    {:<15} {:<6}  Show this help message.", "-h, --help".green(), "");
    println!("    {:<15} {:<6}  Show version.", "-V, --version".green(), "");
}

fn parse_args(argv: &[String]) -> Result<Args, lexopt::Error> {
    use lexopt::prelude::*;
    let mut args = Args::default();
    let mut parser = lexopt::Parser::from_args(argv);

    while let Some(arg) = parser.next()? {
        match arg {
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('f') | Long("fasta") => args.fasta = Some(parser.value()?.parse()?),

            Short('k') | Long("kmer") => args.kmer_size = parser.value()?.parse()?,
            Short('r') | Long("region") => args.bg_region = Some(parser.value()?.parse()?),
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
    if args.fasta.is_none() {
        Err(lexopt::Error::from("Fasta file is not provided (see -f/--fasta)"))?;
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
    fasta_filename: &Path,
    contigs: &Rc<ContigNames>,
    bg_region: &Option<String>
) -> Result<Interval, Error>
{
    if let Some(s) = bg_region {
        let region = Interval::parse(s, contigs).map_err(|_| Error::InvalidInput(
            format!("Input fasta '{}' does not contain region {}", fasta_filename.display(), s)))?;
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
                log::error!("Chromosome {} is in the input fasta file '{}', but is shorter than expected",
                    fasta_filename.display(), region.contig_name());
            }
        }
    }
    Err(Error::InvalidInput(format!(
        "Input fasta '{}' does not contain any of the default background regions. \
        Consider setting region via --region or use a different fasta file.", fasta_filename.display())))
}

/// Extracts the sequence of a background region, used to estimate the parameters of the sequencing data.
/// The sequence is then written to the `$bg_path/bg.fasta`.
fn extract_bg_region<R: Read + Seek>(
    fasta: &mut IndexedReader<R>,
    db_path: &Path,
    region: &Interval,
) -> Result<(), Error>
{
    fasta.fetch(region.contig_name(), u64::from(region.start()), u64::from(region.end()))?;
    let mut seq = Vec::new();
    fasta.read(&mut seq)?;
    crate::seq::standardize(&mut seq);

    let mut bg_path = db_path.join("bg");
    fs::create_dir(&bg_path)?;
    bg_path.push("bg.fasta");
    log::info!("Writing background region {} to '{}'", region, bg_path.display());
    let mut fasta_writer = fs::File::create(&bg_path).map(io::BufWriter::new)?;
    crate::seq::write_fasta(&mut fasta_writer, b"bg", &seq)?;
    Ok(())
}

fn run_jellyfish(
    db_path: &Path,
    fasta_filename: &Path,
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
        .arg(fasta_filename);
    log::info!("Counting {}-mers in {} threads, output: {}", kmer_size, threads, jf_path.display());
    log::debug!("    {:?}", command);

    let output = command.output()?;
    if output.status.success() {
        Ok(())
    } else {
        Err(Error::SubcommandFail(output))
    }
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args: Args = process_args(parse_args(argv)?)?;
    let db_path = args.database.as_ref().expect("Error is impossible");
    if db_path.exists() {
        if args.force {
            log::warn!("Completely removing output directory '{}'", db_path.display());
            fs::remove_dir_all(db_path)?;
        } else {
            panic!("Output directory '{}' already exists. Remove it or use -F/--force flag.", db_path.display());
        }
    }
    fs::create_dir(db_path)?;

    let fasta_filename = args.fasta.as_ref().expect("Error is impossible");
    let mut fasta = IndexedReader::from_file(&fasta_filename)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
    let contigs = Rc::new(ContigNames::from_index("genome".to_string(), &fasta.index));
    let region = select_bg_interval(&fasta_filename, &contigs, &args.bg_region)?;
    extract_bg_region(&mut fasta, &db_path, &region)?;

    run_jellyfish(&db_path, &fasta_filename, &args.jellyfish, args.kmer_size, args.threads, contigs.genome_size())?;
    log::info!("Success!");
    Ok(())
}
