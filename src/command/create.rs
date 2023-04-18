//! Create a new database.

use std::{
    fs,
    io::{Read, Seek},
    cmp::max,
    rc::Rc,
    process::Command,
    path::{Path, PathBuf},
    time::Instant,
};
use bio::io::fasta::IndexedReader;
use colored::Colorize;
use crate::{
    err::{Error, validate_param},
    seq::{
        Interval, ContigNames,
        kmers::{MAX_KMER, JfKmerGetter},
    },
    ext::{sys as sys_ext, fmt as fmt_ext},
};
use super::paths;

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

impl Args {
    fn validate(mut self) -> Result<Self, Error> {
        self.threads = max(self.threads, 1);
        validate_param!(self.database.is_some(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.reference.is_some(), "Reference fasta file is not provided (see -r/--reference)");
        self.jellyfish = sys_ext::find_exe(self.jellyfish)?;
        validate_param!(self.kmer_size % 2 == 1 && self.kmer_size <= MAX_KMER,
            "k-mer size ({}) must be odd, and at most {}", self.kmer_size, MAX_KMER);
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
        "Usage:".bold(), super::PKG_NAME);

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Output database directory.",
        "-d, --db".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file. Must contain FAI index.",
        "-r, --reference".green(), "FILE".yellow());

    println!("\n{}", "Optional parameters:".bold());
    println!("    {:KEY$} {:VAL$}  k-mer size [{}].",
        "-k, --kmer".green(), "INT".yellow(), defaults.kmer_size);
    println!("    {:KEY$} {:VAL$}  Preprocess WGS data based on this background region. Must not be\n\
        {EMPTY}  duplicated in the genome. Defaults to: chr17:72062001-76562000 (GRCh38).",
        "-b, --bg-region".green(), "STR".yellow());
    // TODO: Parse default regions from BG_REGIONS.

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads);
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());
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
            Short('b') | Long("bg") | Long("bg-region") => args.bg_region = Some(parser.value()?.parse()?),
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,
            Long("jellyfish") => args.jellyfish = parser.value()?.parse()?,

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
            format!("Reference file {} does not contain region {}", fmt_ext::path(ref_filename), s)))?;
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
                log::error!("Chromosome {} is in the reference file {}, but is shorter than expected",
                    region.contig_name(), fmt_ext::path(ref_filename));
            }
        }
    }
    Err(Error::InvalidInput(format!(
        "Reference file {} does not contain any of the default background regions. \
        Consider setting region via --region or use a different reference file.", fmt_ext::path(ref_filename))))
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
    let seq = region.fetch_seq(fasta)?;
    let bg_dir = db_path.join(paths::BG_DIR);
    sys_ext::mkdir(&bg_dir)?;
    let bg_fasta_filename = bg_dir.join(paths::BG_FASTA);
    log::info!("Writing background region {} to {}", region, fmt_ext::path(&bg_fasta_filename));
    let mut fasta_writer = sys_ext::create_gzip(&bg_fasta_filename)?;
    crate::seq::write_fasta(&mut fasta_writer, b"bg", Some(region.to_string().as_bytes()), &seq)?;

    log::info!("Calculating k-mer counts on the background region.");
    let kmer_counts = kmer_getter.fetch_one(seq)?;
    let mut kmers_out = sys_ext::create_gzip(&bg_dir.join(paths::KMERS))?;
    kmer_counts.save(&mut kmers_out)?;
    Ok(())
}

fn run_jellyfish(db_path: &Path, ref_filename: &Path, args: &Args, genome_size: u64) -> Result<PathBuf, Error> {
    let mut jf_path = db_path.join(paths::JF_DIR);
    sys_ext::mkdir(&jf_path)?;
    jf_path.push(&format!("{}.jf", args.kmer_size));
    if jf_path.exists() {
        log::warn!("{} already exists, skipping k-mer counting!", fmt_ext::path(&jf_path));
        return Ok(jf_path)
    }

    let mut command = Command::new(&args.jellyfish);
    command.args(&["count", "--canonical", "--lower-count=2", "--out-counter-len=1",
            &format!("--mer-len={}", args.kmer_size),
            &format!("--threads={}", args.threads),
            &format!("--size={}", genome_size)])
        .arg("--output").arg(&jf_path)
        .arg(ref_filename);
    log::info!("Counting {}-mers in {} threads", args.kmer_size, args.threads);
    log::debug!("    {}", fmt_ext::command(&command));

    let start = Instant::now();
    let output = command.output()?;
    log::debug!("    Finished in {:?}", start.elapsed());
    if output.status.success() {
        Ok(jf_path)
    } else {
        Err(Error::SubprocessFail(output))
    }
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args = parse_args(argv)?.validate()?;
    // unwrap as args.database was previously checked to be Some.
    let db_path = args.database.as_ref().unwrap();
    if db_path.exists() {
        if args.force {
            log::warn!("Completely removing output directory {}", fmt_ext::path(db_path));
            fs::remove_dir_all(db_path)?;
        }
    }
    sys_ext::mkdir(db_path)?;

    // unwrap as args.reference was previously checked to be Some.
    let ref_filename = args.reference.as_ref().unwrap();
    let (contigs, mut fasta) = ContigNames::load_indexed_fasta("reference", &ref_filename)?;
    let jf_path = run_jellyfish(&db_path, &ref_filename, &args, contigs.genome_size())?;
    let kmer_getter = JfKmerGetter::new(args.jellyfish.clone(), jf_path)?;

    let region = select_bg_interval(&ref_filename, &contigs, &args.bg_region)?;
    extract_bg_region(&region, &mut fasta, &db_path, &kmer_getter)?;

    sys_ext::mkdir(&db_path.join(paths::LOCI_DIR))?;
    log::info!("Success!");
    Ok(())
}
