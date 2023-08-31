//! Create a new database.

use std::{
    cmp::{min, max},
    process::{Command, Stdio},
    path::{Path, PathBuf},
    time::Instant,
};
use colored::Colorize;
use crate::{
    err::{Error, validate_param, add_path},
    seq::{
        ContigNames,
        kmers::Kmer,
    },
    ext::{self, fmt::PrettyU64},
};
use super::paths;

struct Args {
    reference: Option<PathBuf>,
    database: Option<PathBuf>,
    kmer_size: u8,
    threads: u16,
    jellyfish: PathBuf,
    jf_size: Option<u64>,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            reference: None,
            database: None,
            kmer_size: 25,
            threads: 8,
            jellyfish: PathBuf::from("jellyfish"),
            jf_size: None,
        }
    }
}

impl Args {
    fn validate(mut self) -> Result<Self, Error> {
        self.threads = max(self.threads, 1);
        validate_param!(self.database.is_some(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.reference.is_some(), "Reference fasta file is not provided (see -r/--reference)");
        self.jellyfish = ext::sys::find_exe(self.jellyfish)?;
        validate_param!(self.kmer_size % 2 == 1 && self.kmer_size <= u64::MAX_KMER_SIZE,
            "k-mer size ({}) must be odd, and at most {}", self.kmer_size, u64::MAX_KMER_SIZE);
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 16;
    const VAL: usize = 4;
    // const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", format!("Create an {} database of complex loci.", "empty".bold()).yellow());

    println!("\n{} {} create -d db -r reference.fa [arguments]",
        "Usage:".bold(), super::PROGRAM);

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Output database directory.",
        "-d, --database".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file. Must contain FAI index.",
        "-r, --reference".green(), "FILE".yellow());

    println!("\n{}", "Jellyfish arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Jellyfish executable [{}].",
        "    --jellyfish".green(), "EXE".yellow(), super::fmt_def(defaults.jellyfish.display()));
    println!("    {:KEY$} {:VAL$}  k-mer size [{}].",
        "-k, --kmer".green(), "INT".yellow(), super::fmt_def(defaults.kmer_size));
    println!("    {:KEY$} {:VAL$}  Override Jellyfish cache size.",
        "    --jf-size".green(), "INT".yellow());

    println!("\n{}", "Execution arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), super::fmt_def(defaults.threads));

    println!("\n{}", "Other arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Show this help message.", "-h, --help".green(), "");
    println!("    {:KEY$} {:VAL$}  Show version.", "-V, --version".green(), "");
}

fn parse_args(argv: &[String]) -> Result<Args, lexopt::Error> {
    if argv.is_empty() {
        print_help();
        std::process::exit(1);
    }
    use lexopt::prelude::*;
    let mut args = Args::default();
    let mut parser = lexopt::Parser::from_args(argv);

    while let Some(arg) = parser.next()? {
        match arg {
            Short('d') | Long("db") | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),

            Short('k') | Long("kmer") => args.kmer_size = parser.value()?.parse()?,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Long("jellyfish") => args.jellyfish = parser.value()?.parse()?,
            Long("jf-size") => args.jf_size = Some(parser.value()?.parse::<PrettyU64>()?.get()),

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

fn run_jellyfish(db_path: &Path, ref_filename: &Path, args: &Args, genome_size: u64) -> Result<PathBuf, Error> {
    let mut jf_path = db_path.join(paths::JF_DIR);
    ext::sys::mkdir(&jf_path)?;
    jf_path.push(&format!("{}.jf", args.kmer_size));
    if jf_path.exists() {
        log::warn!("{} already exists, skipping k-mer counting!", ext::fmt::path(&jf_path));
        return Ok(jf_path)
    }

    let timer = Instant::now();
    let jf_size = args.jf_size.unwrap_or(min(genome_size, 10_000_000_000));
    let mut command = Command::new(&args.jellyfish);
    command.args(&["count", "--canonical", "--lower-count=2", "--out-counter-len=1",
            &format!("--mer-len={}", args.kmer_size),
            &format!("--threads={}", args.threads),
            &format!("--size={}", jf_size)])
        .arg("--output").arg(&jf_path)
        .arg(ref_filename)
        .stderr(Stdio::piped());
    log::info!("Counting {}-mers in {} threads", args.kmer_size, args.threads);
    log::debug!("    {}", ext::fmt::command(&command));
    let child = command.spawn().map_err(add_path!(args.jellyfish))?;
    let guard = ext::sys::PipeGuard::new(args.jellyfish.clone(), child);
    guard.wait()?;
    log::debug!("    Finished in {}", ext::fmt::Duration(timer.elapsed()));
    Ok(jf_path)
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args = parse_args(argv)?.validate()?;
    super::greet();
    // unwrap as args.database was previously checked to be Some.
    let db_path = args.database.as_ref().unwrap();
    let success_path = db_path.join(paths::JF_DIR).join(paths::SUCCESS);
    if success_path.exists() && db_path.read_dir().map_err(add_path!(db_path))?.next().is_some() {
        log::error!("Output directory {} is not empty.", ext::fmt::path(db_path));
        log::warn!("Please remove it manually or select a different path.");
        std::process::exit(1);
    }
    ext::sys::mkdir(db_path)?;

    // unwrap as args.reference was previously checked to be Some.
    let ref_filename = args.reference.as_ref().unwrap();
    let (contigs, _) = ContigNames::load_indexed_fasta("ref", &ref_filename)?;
    run_jellyfish(&db_path, &ref_filename, &args, contigs.genome_size())?;

    ext::sys::mkdir(&db_path.join(paths::LOCI_DIR))?;
    super::write_success_file(&success_path)?;
    log::info!("Success!");
    Ok(())
}
