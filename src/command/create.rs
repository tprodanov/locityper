//! Create a new database.

use std::{
    io, fs,
    cmp::max,
    path::{Path, PathBuf},
};
use bio::io::fasta::IndexedReader;
use colored::Colorize;
use crate::seq::{
    contigs::ContigNames,
    interv::Interval,
};

struct Args {
    fasta: Option<PathBuf>,
    database: Option<PathBuf>,
    kmer_size: u8,
    bg_region: Option<String>,
    threads: u8,
    force: bool,
    jellyfish: String,
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
            jellyfish: "jellyfish".to_string(),
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
        "    --jellyfish".green(), "EXE".yellow(), defaults.jellyfish);

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
    process_args(args)
}

fn process_args(mut args: Args) -> Result<Args, lexopt::Error> {
    args.threads = max(args.threads, 1);
    if args.database.is_none() {
        Err("Database directory is not provided (see -d/--database)")?;
    }
    if args.fasta.is_none() {
        Err("Fasta file is not provided (see -f/--fasta)")?;
    }
    Ok(args)
}

const BG_REGIONS: [&'static str; 1] = ["chr17:72062001-76562000"];

fn extract_bg_region(
    fasta_filename: impl AsRef<Path>,
    mut bg_path: PathBuf,
    bg_region: &Option<String>
) -> io::Result<()> {
    // let fasta = fasta::IndexedReader::from_file(fasta_filename)?;
    // let contigs = ContigNames::from_index("genome".to_string(), &fasta.index);
    // if let Some(region) = bg_region {
    //     // Interval::parse(region, &contigs).map_err(|_| Error::new(ErrorKind::))
    // }
    // for region in contigs.iter()

    // fs::create_dir(&bg_path)?;
    // bg_path.push("bg.fasta");
    // log::info!("Writing background regions to '{}'", bg_path.display());
    Ok(())
}

pub(super) fn run(argv: &[String]) -> io::Result<()> {
    let mut args = parse_args(argv).unwrap();
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

    extract_bg_region(args.fasta.as_ref().unwrap(), db_path.join("bg"), &args.bg_region)?;
    Ok(())
}
