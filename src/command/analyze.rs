use std::{
    io, fs,
    cmp::max,
    path::{Path, PathBuf},
};
use colored::Colorize;
use const_format::{str_repeat, concatcp};
use fnv::FnvHashSet;
use crate::{
    err::{Error, validate_param},
    seq::{
        ContigSet,
        recruit, fastx,
    },
    ext::{fmt as fmt_ext, sys as sys_ext},
};
use super::paths;

struct Args {
    input: Vec<PathBuf>,
    database: Option<PathBuf>,
    output: Option<PathBuf>,
    subset_loci: FnvHashSet<String>,

    interleaved: bool,
    threads: u16,
    force: bool,

    params: recruit::Params,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::with_capacity(2),
            database: None,
            output: None,
            subset_loci: FnvHashSet::default(),

            interleaved: false,
            threads: 4,
            force: false,

            params: Default::default(),
        }
    }
}

impl Args {
    /// Validate arguments, modifying some, if needed.
    fn validate(mut self) -> Result<Self, Error> {
        self.threads = max(self.threads, 1);
        let n_input = self.input.len();
        validate_param!(n_input > 0, "Read files are not provided (see -i/--input)");
        validate_param!(n_input != 2 || !self.interleaved,
            "Two read files (-i/--input) are provided, however, --interleaved is specified");
        if !self.interleaved && n_input == 1 {
            log::warn!("Running in single-end mode.");
        }

        validate_param!(self.database.is_some(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.output.is_some(), "Output directory is not provided (see -o/--output)");
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Analyze WGS dataset.".yellow());

    println!("\n{} {} analyze -i reads1.fq [reads2.fq] -d db -o out [arguments]",
        "Usage:".bold(), super::PKG_NAME);

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Database directory (initialized with {} & {}).",
        "-d, --db".green(), "DIR".yellow(), concatcp!(super::PKG_NAME, " create").underline(), "add".underline());
    println!("    {:KEY$} {:VAL$}  Output directory   (initialized with {}).",
        "-o, --output".green(), "DIR".yellow(), concatcp!(super::PKG_NAME, " preproc").underline());
    println!("    {:KEY$} {:VAL$}  Input reads are interleaved.",
        "-^, --interleaved".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Optional: only analyze loci with names from this list.",
        "    --subset-loci".green(), "STR+".yellow());

    println!("\n{}", "Read recruitment:".bold());
    println!("    {:KEY$} {:VAL$}  k-mer size (no larger than 16) [{}].",
        "-k, --recr-kmer".green(), "INT".yellow(), defaults.params.minimizer_k);
    println!("    {:KEY$} {:VAL$}  Take k-mers with smallest hash across {} consecutive k-mers [{}].",
        "-w, --recr-window".green(), "INT".yellow(), "INT".yellow(), defaults.params.minimizer_w);
    println!("    {:KEY$} {:VAL$}  Recruit single-/paired-reads with at least {} k-mers matches [{}].",
        "-m, --min-matches".green(), "INT".yellow(), "INT".yellow(), defaults.params.min_matches);
    println!("    {:KEY$} {:VAL$}  Recruit reads in chunks of this size [{}].",
        "-c, --chunk-size".green(), "INT".yellow(), defaults.params.chunk_size);

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads);
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());

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
            Short('i') | Long("input") => args.input
                .extend(parser.values()?.take(2).map(|s| s.parse()).collect::<Result<Vec<_>, _>>()?),
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),
            Long("subset-loci") => {
                args.subset_loci.extend(parser.values()?.map(|s| s.parse()).collect::<Result<Vec<_>, _>>()?);
                if args.subset_loci.is_empty() {
                    return Err(lexopt::Error::MissingValue { option: Some("subset-loci".to_owned()) });
                }
            }

            Short('k') | Long("recr-kmer") => args.params.minimizer_k = parser.value()?.parse()?,
            Short('w') | Long("recr-window") => args.params.minimizer_w = parser.value()?.parse()?,
            Short('m') | Long("min-matches") => args.params.min_matches = parser.value()?.parse()?,
            Short('c') | Long("chunk") | Long("chunk-size") => args.params.chunk_size = parser.value()?.parse()?,

            Short('^') | Long("interleaved") => args.interleaved = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,

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

fn locus_name_matches<'a>(path: &'a Path, subset_loci: &FnvHashSet<String>) -> Option<&'a str> {
    match path.file_name().unwrap().to_str() {
        None => {
            log::error!("Skipping directory {:?} - filename is not a valid UTF-8", path);
            None
        }
        Some(name) => {
            if !subset_loci.is_empty() && !subset_loci.contains(name) {
                log::trace!("Skipping locus {} (not in the subset loci)", name);
                None
            } else if name.starts_with(".") {
                log::trace!("Skipping hidden directory {}", name);
                None
            } else {
                Some(name)
            }
        }
    }
}

struct LocusData {
    set: ContigSet,
    /// Output directory with locus data.
    locus_dir: PathBuf,
    /// Temporary file with recruited reads.
    tmp_reads_filename: PathBuf,
    /// Final file with recruited reads.
    reads_filename: PathBuf,
    /// Temporary file with read alignments to the haplotypes.
    tmp_aln_filename: PathBuf,
    /// Final file with read alignments to the haplotypes.
    aln_filename: PathBuf,
}

impl LocusData {
    fn new(set: ContigSet, loci_dir: &Path, force: bool) -> io::Result<Self> {
        let locus_dir = loci_dir.join(set.tag());
        if locus_dir.exists() && force {
            fs::remove_dir_all(&locus_dir)?;
        }
        sys_ext::mkdir(&locus_dir)?;
        Ok(Self {
            tmp_reads_filename: locus_dir.join("reads.tmp.fq.gz"),
            reads_filename: locus_dir.join("reads.fq.gz"),
            tmp_aln_filename: locus_dir.join("aln.tmp.bam"),
            aln_filename: locus_dir.join("aln.bam"),
            locus_dir, set,
        })
    }
}

/// Loads all loci from the database. If `subset_loci` is not empty, only loads loci that are contained in it.
/// If `force`, delete old data in the corresponding output directories.
fn load_loci(
    db_path: &Path,
    out_path: &Path,
    subset_loci: &FnvHashSet<String>,
    force: bool
) -> io::Result<Vec<LocusData>>
{
    log::info!("Loading database.");
    let db_loci_dir = db_path.join(paths::LOCI_DIR);
    let out_loci_dir = out_path.join(paths::LOCI_DIR);
    sys_ext::mkdir(&out_loci_dir)?;
    if force {
        log::warn!("Force flag is set: overwriting output directories {}/*", fmt_ext::path(&out_loci_dir));
    }

    let mut loci = Vec::new();
    let mut total_entries = 0;
    for entry in fs::read_dir(&db_loci_dir)? {
        let entry = entry?;
        if !entry.file_type()?.is_dir() {
            continue;
        }

        total_entries += 1;
        let path = entry.path();
        if let Some(name) = locus_name_matches(&path, subset_loci) {
            let fasta_filename = path.join(paths::LOCUS_FASTA);
            let kmers_filename = path.join(paths::KMERS);
            match ContigSet::load(name, &fasta_filename, &kmers_filename, ()) {
                Ok(set) => loci.push(LocusData::new(set, &out_loci_dir, force)?),
                Err(e) => log::error!("Could not load locus information from {}: {:?}", fmt_ext::path(&path), e),
            }
        }
    }
    let n = loci.len();
    if n < total_entries {
        log::info!("Loaded {} loci, skipped {} directories", n, total_entries - n);
    } else {
        log::info!("Loaded {} loci", n);
    }
    Ok(loci)
}

/// Recruits reads to all loci, where neither reads nor alignments are available.
fn recruit_reads(loci: &[LocusData], args: &Args) -> io::Result<()> {
    let filt_loci: Vec<&LocusData> = loci.iter()
        .filter(|locus| !locus.reads_filename.exists() && !locus.aln_filename.exists())
        .collect();
    if filt_loci.is_empty() {
        log::info!("Skipping read recruitment");
        return Ok(());
    }
    if filt_loci.len() < loci.len() {
        log::info!("Skipping read recruitment to {} loci", loci.len() - filt_loci.len());
    }

    let targets = recruit::Targets::new(filt_loci.iter().map(|locus| &locus.set), &args.params);
    let writers: Vec<_> = filt_loci.iter()
        .map(|locus| sys_ext::create_gzip(&locus.tmp_reads_filename))
        .collect::<Result<_, _>>()?;

    // Cannot put reader into a box, because `FastxRead` has a type parameter.
    if args.input.len() == 1 && !args.interleaved {
        let reader = fastx::Reader::new(sys_ext::open(&args.input[0])?)?;
        targets.recruit(reader, writers, args.threads, args.params.chunk_size)?;
    } else if args.interleaved {
        let reader = fastx::PairedEndInterleaved::new(fastx::Reader::new(sys_ext::open(&args.input[0])?)?);
        targets.recruit(reader, writers, args.threads, args.params.chunk_size)?;
    } else {
        let reader = fastx::PairedEndReaders::new(
            fastx::Reader::new(sys_ext::open(&args.input[0])?)?,
            fastx::Reader::new(sys_ext::open(&args.input[1])?)?);
        targets.recruit(reader, writers, args.threads, args.params.chunk_size)?;
    }
    for locus in filt_loci.iter() {
        fs::rename(&locus.tmp_reads_filename, &locus.reads_filename)?;
    }
    Ok(())
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args = parse_args(argv)?.validate()?;
    let db_dir = args.database.as_ref().unwrap();
    let out_dir = args.output.as_ref().unwrap();
    let loci = load_loci(db_dir, out_dir, &args.subset_loci, args.force)?;
    recruit_reads(&loci, &args)?;
    Ok(())
}
