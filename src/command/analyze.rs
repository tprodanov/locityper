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
        kmers::KmerCount,
        recruit::Targets,
        fastx,
    },
    ext::{fmt as fmt_ext, sys as sys_ext},
};
use super::paths;

struct Args {
    input: Vec<PathBuf>,
    database: Option<PathBuf>,
    output: Option<PathBuf>,
    subset_loci: FnvHashSet<String>,

    /// Filter out k-mers with frequency over `max_kmer_freq` across the reference genome.
    max_kmer_freq: KmerCount,
    /// How many k-mers need to be found per read for it to be recruited.
    min_kmer_matches: u16,

    interleaved: bool,
    threads: u16,
    force: bool,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::with_capacity(2),
            database: None,
            output: None,
            subset_loci: FnvHashSet::default(),

            max_kmer_freq: 10,
            min_kmer_matches: 10,

            interleaved: false,
            threads: 4,
            force: false,
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
    println!("    {:KEY$} {:VAL$}  Discard k-mers appearing over {} times in the reference genome [{}].",
        "    --max-freq".green(), "INT".yellow(), "INT".yellow(), defaults.max_kmer_freq);
    println!("    {:KEY$} {:VAL$}  Recruit single-/paired-reads with at least {} k-mers matches [{}].",
        "    --min-matches".green(), "INT".yellow(), "INT".yellow(), defaults.min_kmer_matches);

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

            Long("max-freq") | Long("max-frequency") => args.max_kmer_freq = parser.value()?.parse()?,
            Long("min-matches") => args.min_kmer_matches = parser.value()?.parse()?,

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

fn load_loci(db_path: &Path, subset_loci: &FnvHashSet<String>) -> io::Result<Vec<ContigSet>> {
    log::info!("Loading database.");
    let loci_dir = db_path.join(paths::LOCI_DIR);
    let mut contig_sets = Vec::new();
    let mut total_entries = 0;
    for entry in fs::read_dir(&loci_dir)? {
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
                Ok(set) => contig_sets.push(set),
                Err(e) => log::error!("Could not load locus information from {}: {:?}", fmt_ext::path(&path), e),
            }
        }
    }
    let n = contig_sets.len();
    if n < total_entries {
        log::info!("Loaded {} loci, skipped {} directories", n, total_entries - n);
    } else {
        log::info!("Loaded {} loci", n);
    }
    Ok(contig_sets)
}

/// Returns a vector of tuples (reads_filename, alns_filename) for each locus.
/// Additionally, this function creates missing directories on the path to the files.
fn get_read_aln_filenames(out_dir: &Path, contig_sets: &[ContigSet]) -> io::Result<Vec<(PathBuf, PathBuf)>> {
    let loci_dir = out_dir.join(paths::LOCI_DIR);
    sys_ext::mkdir(&loci_dir)?;
    let mut filenames = Vec::with_capacity(contig_sets.len());
    for set in contig_sets.iter() {
        let locus_dir = loci_dir.join(set.contigs().tag());
        sys_ext::mkdir(&locus_dir)?;
        filenames.push((locus_dir.join("reads.fq.gz"), locus_dir.join("alns.bam")));
    }
    Ok(filenames)
}

fn recruit_reads(contig_sets: &[ContigSet], filenames: &[(PathBuf, PathBuf)], args: &Args) -> io::Result<()> {
    assert_eq!(contig_sets.len(), filenames.len());
    let n_loci = contig_sets.len();
    let (contig_sets, filenames): (Vec<_>, Vec<_>) = contig_sets.iter().zip(filenames)
        .filter_map(|(set, (f1, f2))| if f1.exists() || f2.exists() { None } else { Some((set, f1)) })
        .unzip();
    if contig_sets.is_empty() {
        log::info!("Skipping read recruitment");
        return Ok(());
    }
    if contig_sets.len() < n_loci {
        log::info!("Skipping read recruitment to {} loci", n_loci - contig_sets.len());
    }
    let targets = Targets::new(contig_sets.into_iter(), args.max_kmer_freq, args.min_kmer_matches)?;

    let out_filename = args.output.as_ref().unwrap().join("tmp.fq");
    // Cannot put reader into a box, because `FastxRead` has a type parameter.
    if args.input.len() == 1 && !args.interleaved {
        let reader = fastx::Reader::new(sys_ext::open(&args.input[0])?)?;
        targets.recruit(reader, out_filename, args.threads)
    } else if args.interleaved {
        let reader = fastx::PairedEndInterleaved::new(fastx::Reader::new(sys_ext::open(&args.input[0])?)?);
        targets.recruit(reader, out_filename, args.threads)
    } else {
        let reader = fastx::PairedEndReaders::new(
            fastx::Reader::new(sys_ext::open(&args.input[0])?)?,
            fastx::Reader::new(sys_ext::open(&args.input[1])?)?);
        targets.recruit(reader, out_filename, args.threads)
    }
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args = parse_args(argv)?.validate()?;
    let contig_sets = load_loci(args.database.as_ref().unwrap(), &args.subset_loci)?;
    let read_aln_filenames = get_read_aln_filenames(args.output.as_ref().unwrap(), &contig_sets)?;
    recruit_reads(&contig_sets, &read_aln_filenames, &args)?;
    Ok(())
}
