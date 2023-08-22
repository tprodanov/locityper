//! Create a new database.

use std::{
    io::{Read, Write, Seek},
    cmp::{min, max},
    sync::Arc,
    process::Command,
    path::{Path, PathBuf},
    time::Instant,
};
use bio::io::fasta::IndexedReader;
use colored::Colorize;
use crate::{
    err::{Error, validate_param, add_path},
    seq::{
        Interval, ContigNames,
        kmers::{JfKmerGetter, Kmer},
        contigs::GenomeVersion,
    },
    ext::{self, fmt::PrettyU64},
};
use super::paths;

struct Args {
    reference: Option<PathBuf>,
    database: Option<PathBuf>,
    kmer_size: u8,
    bg_region: Option<String>,
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
            bg_region: None,
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
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", format!("Create an {} database of complex loci.", "empty".bold()).yellow());

    println!("\n{} {} create -d db -r reference.fa [arguments]",
        "Usage:".bold(), super::PROGRAM);

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Output database directory.",
        "-d, --db".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file. Must contain FAI index.",
        "-r, --reference".green(), "FILE".yellow());

    println!("\n{}", "Optional parameters:".bold());
    println!("    {:KEY$} {:VAL$}  k-mer size [{}].",
        "-k, --kmer".green(), "INT".yellow(), super::fmt_def(defaults.kmer_size));
    println!("    {:KEY$} {:VAL$}  Preprocess WGS data based on this background region,\n\
        {EMPTY}  preferably >3 Mb and without many duplications.\n\
        {EMPTY}  Default regions are defined for CHM13, GRCh38 and GRCh37.",
        "-b, --bg-region".green(), "STR".yellow());

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), super::fmt_def(defaults.threads));
    println!("    {:KEY$} {:VAL$}  Jellyfish executable [{}].",
        "    --jellyfish".green(), "EXE".yellow(), super::fmt_def(defaults.jellyfish.display()));
    println!("    {:KEY$} {:VAL$}  Override Jellyfish cache size.",
        "    --jf-size".green(), "INT".yellow());

    println!("\n{}", "Other parameters:".bold());
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
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),

            Short('k') | Long("kmer") => args.kmer_size = parser.value()?.parse()?,
            Short('b') | Long("bg") | Long("bg-region") => args.bg_region = Some(parser.value()?.parse()?),
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

/// Returns default background region.
fn default_region(ver: GenomeVersion) -> &'static str {
    match ver {
        GenomeVersion::Chm13 => "chr17:72950001-77450000",
        GenomeVersion::GRCh38 => "chr17:72062001-76562000",
        GenomeVersion::GRCh37 => "chr17:70060001-74560000",
    }
}

/// Returns the first appropriate interval for the contig names.
/// If `bg_region` is set, try to parse it. Otherwise, iterate over `BG_REGIONS`.
///
/// Returns error if no interval is appropriate (chromosome not in the contig set, or interval is out of bounds).
fn select_bg_interval(
    ref_filename: &Path,
    contigs: &Arc<ContigNames>,
    bg_region: &Option<String>
) -> Result<Interval, Error>
{
    if let Some(s) = bg_region {
        let region = Interval::parse(s, contigs).map_err(|_| Error::InvalidInput(
            format!("Reference file {} does not contain region {}", ext::fmt::path(ref_filename), s)))?;
        return if contigs.in_bounds(&region) {
            Ok(region)
        } else {
            Err(Error::InvalidInput(format!("Region {} is out of bounds", s)))
        };
    }

    let genome = GenomeVersion::guess(contigs)
        .ok_or_else(|| Error::RuntimeError("Could not recognize reference genome. \
            Please provide background region (-b) explicitely, preferably >3 Mb long and without many duplications."
            .to_owned()))?;
    let region = default_region(genome);
    log::info!("Recognized {} reference genome, using background region {}", genome, region);
    // Try to crop `chr` if region cannot be found.
    for &crop in &[0, 3] {
        if let Ok(region) = Interval::parse(&region[crop..], contigs) {
            if contigs.in_bounds(&region) {
                return Ok(region);
            } else {
                return Err(Error::InvalidInput(format!("Region {} is too long for the reference genome", region)));
            }
        }
    }
    Err(Error::InvalidInput(format!("Reference file {} does not contain any of the default background regions. \
        Please provide background region (-b) explicitely, preferably >3 Mb long and without many duplications.",
        ext::fmt::path(ref_filename))))
}

/// Extracts the sequence of a background region, used to estimate the parameters of the sequencing data.
/// The sequence is then written to the `$bg_path/bg.fa.gz`.
fn extract_bg_region<R: Read + Seek>(
    region: &Interval,
    fasta: &mut IndexedReader<R>,
    db_path: &Path,
    kmer_getter: &JfKmerGetter,
) -> Result<(), Error>
{
    let seq = region.fetch_seq(fasta).map_err(add_path!(!))?;
    let bg_dir = db_path.join(paths::BG_DIR);
    ext::sys::mkdir(&bg_dir)?;
    let bed_filename = bg_dir.join(paths::BG_BED);
    let mut bed_writer = ext::sys::create_file(&bed_filename)?;
    writeln!(bed_writer, "{}", region.bed_fmt()).map_err(add_path!(bed_filename))?;
    std::mem::drop(bed_writer);

    log::info!("Calculating k-mer counts on the background region.");
    let kmer_counts = kmer_getter.fetch_one(seq)?;
    let kmers_filename = bg_dir.join(paths::KMERS);
    let mut kmers_out = ext::sys::create_lz4_slow(&kmers_filename)?;
    kmer_counts.save(&mut kmers_out).map_err(add_path!(kmers_filename))?;
    Ok(())
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
    let mut command = Command::new(&args.jellyfish);
    let jf_size = args.jf_size.unwrap_or(min(genome_size, 10_000_000_000));
    command.args(&["count", "--canonical", "--lower-count=2", "--out-counter-len=1",
            &format!("--mer-len={}", args.kmer_size),
            &format!("--threads={}", args.threads),
            &format!("--size={}", jf_size)])
        .arg("--output").arg(&jf_path)
        .arg(ref_filename);
    log::info!("Counting {}-mers in {} threads", args.kmer_size, args.threads);
    log::debug!("    {}", ext::fmt::command(&command));

    let output = command.output().map_err(add_path!(!))?;
    log::debug!("    Finished in {}", ext::fmt::Duration(timer.elapsed()));
    if output.status.success() {
        Ok(jf_path)
    } else {
        Err(Error::Subprocess(output))
    }
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args = parse_args(argv)?.validate()?;
    super::greet();
    // unwrap as args.database was previously checked to be Some.
    let db_path = args.database.as_ref().unwrap();
    let success_path = db_path.join(paths::BG_DIR).join(paths::SUCCESS);
    if success_path.exists() && db_path.read_dir().map_err(add_path!(db_path))?.next().is_some() {
        log::error!("Output directory {} is not empty.", ext::fmt::path(db_path));
        log::warn!("Please remove it manually or select a different path.");
        std::process::exit(1);
    }
    ext::sys::mkdir(db_path)?;

    // unwrap as args.reference was previously checked to be Some.
    let ref_filename = args.reference.as_ref().unwrap();
    let (contigs, mut fasta) = ContigNames::load_indexed_fasta("reference", &ref_filename)?;
    let region = select_bg_interval(&ref_filename, &contigs, &args.bg_region)?;
    let jf_path = run_jellyfish(&db_path, &ref_filename, &args, contigs.genome_size())?;
    let kmer_getter = JfKmerGetter::new(args.jellyfish.clone(), jf_path)?;
    extract_bg_region(&region, &mut fasta, &db_path, &kmer_getter)?;

    ext::sys::mkdir(&db_path.join(paths::LOCI_DIR))?;
    super::write_success_file(&success_path)?;
    log::info!("Success!");
    Ok(())
}
