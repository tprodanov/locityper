use std::{
    fs,
    io::{self, BufRead},
    sync::Arc,
    time::Instant,
    path::{Path, PathBuf},
    cmp::max,
};
use regex::Regex;
use crate::{
    algo::{HashMap, HashSet},
    err::{Error, validate_param, add_path},
    ext::{
        self,
        fmt::{PrettyU32, PrettyU64},
    },
    seq::{
        recruit,
        ContigNames, ContigSet,
        fastx::{self, FastxRecord, SingleRecord},
        kmers::{Kmer, KmerCount, KmerCounts, JfKmerGetter},
    },
    bg::Technology,
};
use colored::Colorize;

struct Args {
    in_files: super::preproc::InputFiles,
    jf_counts: Option<PathBuf>,

    seqs: Vec<String>,
    seqs_all: Option<PathBuf>,
    seqs_list: Option<PathBuf>,
    output: Option<String>,
    regions: Option<PathBuf>,

    minimizer_kw: (u8, u8),
    match_frac: f64,
    match_len: u32,
    preset: Option<String>,
    chunk_length: u64,
    thresh_kmer_count: KmerCount,

    threads: u16,
    jellyfish: PathBuf,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            in_files: Default::default(),
            jf_counts: None,

            seqs: Vec::new(),
            seqs_all: None,
            seqs_list: None,
            output: None,
            regions: None,

            minimizer_kw: (15, 5),
            match_frac: 0.6,
            match_len: 2000,
            preset: None,
            chunk_length: 3_000_000,
            thresh_kmer_count: 10,

            threads: 8,
            jellyfish: PathBuf::from("jellyfish"),
        }
    }
}

impl Args {
    fn update_from_preset(&mut self) -> Result<(), Error> {
        let Some(preset) = self.preset.as_ref() else { return Ok(()) };
        let (tech, paired) = match &preset.to_lowercase() as &str {
            "illumina" | "illumina-pe" | "sr" | "sr-pe" => (Technology::Illumina, true),
            "illumina-se" | "sr-se" => (Technology::Illumina, false),
            "hifi" => (Technology::HiFi, false),
            "pacbio" | "pb" => (Technology::PacBio, false),
            "ont" | "nanopore" => (Technology::Nanopore, false),
            _ => return Err(Error::InvalidInput(format!("Unknown preset `{}`", preset))),
        };
        self.minimizer_kw = tech.default_minim_size();
        self.match_frac = tech.default_match_frac(paired);
        Ok(())
    }

    fn validate(mut self) -> Result<Self, Error> {
        self.in_files.validate(false)?;
        self.threads = max(self.threads, 1);
        if self.in_files.alns.iter().any(|filename| filename.extension()
                .map(|ext| ext == "cram" || ext == "CRAM").unwrap_or(false)) {
            validate_param!(self.in_files.reference.is_some(),
                "Input CRAM file requires a reference file (see -a/--alignment and -r/--reference)");
        }

        validate_param!(self.thresh_kmer_count > 0, "k-mer threshold must not be zero");
        if self.jf_counts.is_some() {
            self.jellyfish = ext::sys::find_exe(self.jellyfish)?;
        }
        Ok(self)
    }
}

/// Locus name in input/output arguments must be replaced with `{}`.
const BRACKETS: &'static str = "{}";

fn print_help() {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Recruit reads to one/multiple loci.".yellow());

    println!("\n{}", "Usage:".bold());
    println!("    {} recruit \\", super::PROGRAM);
    println!("        (-i reads1.fq [reads2.fq] | -a reads.bam [--no-index] | -I in-list) \\");
    println!("        (-s seqs.fa | -S all_seqs.fa) -o out.fastq [args]");

    println!("\n{}  (please see {} for more information on {}/{}/{} arguments)",
        "Input arguments:".bold(), "README".italic(), "-i".green(), "-a".green(), "-I".green());
    println!("    {:KEY$} {:VAL$}  Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Reads in BAM/CRAM format, mutually exclusive with {}.\n\
        {EMPTY}  By default, mapped, sorted and indexed BAM/CRAM file is expected,\n\
        {EMPTY}  please specify {} otherwise.",
        "-a, --alignment".green(), "FILE".yellow(), "-i/--input".green(), "--no-index".green());
    println!("    {:KEY$} {:VAL$}  File with input filenames (see {}).",
        "-I, --in-list".green(), "FILE".yellow(), "README".italic());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file. Required with input CRAM file ({} alns.cram).",
        "-r, --reference".green(), "FILE".yellow(), "-a".green());
    println!("    {:KEY$} {:VAL$}  Canonical k-mer counts across the reference genome, calculated\n\
        {EMPTY}  using Jellyfish. Not required, but highly recommended. k does not\n\
        {EMPTY}  have to match minimizer size. See {} for recommended options.",
        "-j, --jf-counts".green(), "FILE".yellow(), "README".italic());

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Interleaved paired-end reads in single input file.",
        "-^, --interleaved".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Use input BAM/CRAM file ({}) without index: goes over all reads.\n\
        {EMPTY}  Single-end and paired-end interleaved ({}) data is allowed.",
        "    --no-index".green(), super::flag(), "-a".green(), "-^".green());

    println!("\n{}", "Target sequences and output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  FASTA file with target sequences for one locus.\n\
        {EMPTY}  To recruit reads to multiple loci, replace locus name with `{}`.\n\
        {EMPTY}  Then, all matching files will be selected.\n\
        {EMPTY}  Multiple entries allowed, but locus names must not repeat.",
        "-s, --seqs".green(), "FILE".yellow(), BRACKETS);
    println!("    {:KEY$} {:VAL$}  Single FASTA file with target sequences for all loci.\n\
        {EMPTY}  Record names should follow the format `LOCUS*SEQ_NAME`.\n\
        {EMPTY}  Mutually exclusive with {}.",
        "-S, --seqs-all".green(), "FILE".yellow(), "-s".green());
    println!("    {:KEY$} {:VAL$}  Path to the (interleaved) output FASTQ files. If more than\n\
        {EMPTY}  one target locus exists, please replace locus name with `{}`.\n\
        {EMPTY}  {} Will create parent directories, if needed.",
        "-o, --output".green(), "FILE".yellow(), BRACKETS, "Note:".bright_yellow());
    println!("    {:KEY$} {:VAL$}  Two column file with input FASTA and output FASTQ filenames.\n\
        {EMPTY}  Mutually exclusive with {}, {} and {}.",
        "-l, --seqs-list".green(), "FILE".yellow(), "-s".green(), "-S".green(), "-o".green());
    println!("    {:KEY$} {:VAL$}  Recruit unmapped reads and reads with primary alignments to these\n\
        {EMPTY}  regions (BED). Only relevant for mapped and indexed BAM/CRAM files.",
        "-R, --regions".green(), "FILE".yellow());

    println!("\n{}", "Read recruitment:".bold());
    println!("    {}  {}  Use k-mers of size {} (<= {}) with smallest hash\n\
        {EMPTY}  across {} consecutive k-mers [{} {}].",
        "-m, --minimizer".green(), "INT INT".yellow(),
        "INT_1".yellow(), recruit::Minimizer::MAX_KMER_SIZE, "INT_2".yellow(),
        super::fmt_def(defaults.minimizer_kw.0), super::fmt_def(defaults.minimizer_kw.1));
    println!("    {:KEY$} {:VAL$}  Minimal fraction of minimizers that need to match reference [{}].",
        "-M, --match-frac".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.match_frac));
    println!("    {:KEY$} {:VAL$}  Recruit long reads with a matching subregion of this length [{}].",
        "-L, --match-len".green(), "INT".yellow(),
        super::fmt_def(defaults.match_len));
    println!("    {:KEY$} {:VAL$}  Parameter preset (illumina|illumina-SE|hifi|pacbio|ont).\n\
        {EMPTY}  Modifies {} and {}, see {} for default values.",
        "-x, --preset".green(), "STR".yellow(), "-m".green(), "-M".green(), "locityper genotype".underline());
    println!("    {:KEY$} {:VAL$}  Only use k-mers that appear less than {} times in the\n\
        {EMPTY}  reference genome [{}]. Requires {} argument.",
        "-t, --kmer-thresh".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(defaults.thresh_kmer_count),
        "-j".green());
    println!("    {:KEY$} {:VAL$}  Recruit reads in chunks of this sum length [{}].\n\
        {EMPTY}  Impacts runtime in multi-threaded read recruitment.",
        "-c, --chunk-len".green(), "INT".yellow(),
        super::fmt_def(PrettyU64(defaults.chunk_length)));

    println!("\n{}", "Execution arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), super::fmt_def(defaults.threads));
    println!("    {:KEY$} {:VAL$}  Jellyfish executable [{}].",
        "    --jellyfish".green(), "EXE".yellow(), super::fmt_def(defaults.jellyfish.display()));

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
            Short('i') | Long("input") => {
                let mut values = parser.values()?.take(2);
                args.in_files.reads1.push(values.next().expect("First argument is always present").parse()?);
                if let Some(val) = values.next() {
                    args.in_files.reads2.push(val.parse()?);
                }
            }
            Short('a') | Long("aln") | Long("alns") | Long("alignment") | Long("alignments") =>
                args.in_files.alns.push(parser.value()?.parse()?),
            Short('I') | Long("in-list") | Long("input-list") => args.in_files.in_list = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.in_files.reference = Some(parser.value()?.parse()?),
            Short('j') | Long("jf-counts") => args.jf_counts = Some(parser.value()?.parse()?),

            Short('s') | Long("seqs") => {
                for entry in parser.values()? {
                    args.seqs.push(entry.parse()?);
                }
            }
            Short('S') | Long("seqs-all") => args.seqs_all = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),
            Short('l') | Long("seqs-list") => args.seqs_list = Some(parser.value()?.parse()?),
            Short('R') | Long("regions") => args.regions = Some(parser.value()?.parse()?),

            Short('m') | Long("minimizer") | Long("minimizers") =>
                args.minimizer_kw = (parser.value()?.parse()?, parser.value()?.parse()?),
            Short('M') | Long("match-frac") | Long("match-fraction") =>
                args.match_frac = parser.value()?.parse()?,
            Short('x') | Long("preset") =>
                args.preset = Some(parser.value()?.parse()?),
            Short('t') | Long("kmer-thresh") | Long("kmer-threshold") =>
                args.thresh_kmer_count = parser.value()?.parse()?,
            Short('L') | Long("match-len") | Long("match-length") =>
                args.match_len = parser.value()?.parse::<PrettyU32>()?.get(),
            Short('c') | Long("chunk") | Long("chunk-len") =>
                args.chunk_length = parser.value()?.parse::<PrettyU64>()?.get(),

            Short('^') | Long("interleaved") => args.in_files.interleaved = true,
            Long("no-index") => args.in_files.no_index = true,

            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Long("jellyfish") => args.jellyfish = parser.value()?.parse()?,

            Short('V') | Long("version") => {
                super::print_version();
                std::process::exit(0);
            }
            Short('h') | Long("help") | Short('H') | Long("full-help") | Long("hidden-help") => {
                print_help();
                std::process::exit(0);
            }
            _ => Err(arg.unexpected())?,
        }
    }
    Ok(args)
}

type BufFile = io::BufWriter<fs::File>;
type SeqsFiles = (Vec<Vec<Vec<u8>>>, Vec<BufFile>);

/// Creates plain-text output file with large buffer,
fn create_output(filename: impl AsRef<Path>) -> Result<BufFile, Error> {
    let filename = filename.as_ref();
    if let Some(dir) = filename.parent() {
        fs::create_dir_all(dir).map_err(add_path!(dir))?;
    }
    // Output files with a large buffer (4 Mb).
    const BUFFER: usize = 4_194_304;
    fs::File::create(filename).map_err(add_path!(filename))
        .map(|w| io::BufWriter::with_capacity(BUFFER, w))
}

/// Load sequences and output files according to -S and -o arguments.
fn load_seqs_all(seqs_all: &Path, output: &str) -> Result<SeqsFiles, Error> {
    validate_param!(output.contains(BRACKETS),
        "Output path ({}) must contain {}, when using -S", output, BRACKETS);
    log::info!("Loading all sequences from {}", ext::fmt::path(seqs_all));

    // locus -> (seqs, output).
    let mut map = HashMap::default();
    let mut fasta_reader = fastx::Reader::from_path(seqs_all)?;
    let mut record = FastxRecord::default();
    while fasta_reader.read_next_standardized(&mut record)? {
        let name = std::str::from_utf8(record.name())
            .map_err(|_| Error::Utf8("read name", record.name().to_vec()))?;
        let (locus, _) = name.split_once('*').ok_or_else(||
                Error::InvalidData(format!("Invalid record name `{}` for -S {}. Symbol `*` required",
                name, ext::fmt::path(seqs_all))))?;
        let seq = record.seq().to_vec();
        match map.get_mut(locus) {
            None => {
                let file = create_output(output.replace(BRACKETS, locus))?;
                map.insert(locus.to_owned(), (vec![seq], file));
            }
            Some(entry) => entry.0.push(seq),
        }
    }
    Ok(map.into_values().unzip())
}

/// Load sequences and output files according to a single -s and -o arguments.
fn load_seqs_single(
    seqs_filename: &str,
    output_filename: &str,
    seqs: &mut Vec<Vec<Vec<u8>>>,
    files: &mut Vec<BufFile>,
) -> Result<(), Error>
{
    let file = create_output(output_filename)?;
    let mut fasta_reader = fastx::Reader::from_path(seqs_filename)?;
    let mut record = FastxRecord::default();
    let mut curr_seqs = Vec::new();
    while fasta_reader.read_next_standardized(&mut record)? {
        curr_seqs.push(record.seq().to_vec());
    }
    if curr_seqs.is_empty() {
        log::error!("Loaded zero sequences from {}, skipping it", ext::fmt::path(seqs_filename));
    } else {
        seqs.push(curr_seqs);
        files.push(file);
    }
    Ok(())
}

/// Load sequences and output files according to -s and -o arguments with brackets.
fn load_seqs_glob(patterns: &[String], output: &str) -> Result<SeqsFiles, Error> {
    let mut seqs = Vec::new();
    let mut files = Vec::new();
    if patterns.len() == 1 && !patterns[0].contains(BRACKETS) {
        let seq_path = &patterns[0];
        validate_param!(!output.contains(BRACKETS),
            "Single input -s {} does not allow output file with `{}` ({})",
            ext::fmt::path(seq_path), BRACKETS, output);
        log::info!("Loading sequences from {}", ext::fmt::path(seq_path));
        load_seqs_single(seq_path, output, &mut seqs, &mut files)?;
        return Ok((seqs, files));
    }

    let mut loci = HashSet::default();
    validate_param!(output.contains(BRACKETS),
        "Output path ({}) must contain {}, when using -S", output, BRACKETS);
    for pattern in patterns {
        validate_param!(pattern.contains(BRACKETS), "Paths to sequences ({}) must contain brackets `{}`",
            ext::fmt::path(pattern), BRACKETS);
        validate_param!(!pattern.contains("*"), "Paths to sequences ({}) must not contain `*`",
            ext::fmt::path(pattern));
        let can_pattern = fs::canonicalize(pattern).ok()
            .and_then(|p| p.into_os_string().into_string().ok())
            .unwrap_or_else(|| pattern.to_owned());
        let glob_pattern = can_pattern.replace(BRACKETS, "*");
        let regex_pattern = format!("^{}$", can_pattern.replacen(BRACKETS, "([^/\\])+", 1).replace(BRACKETS, "\\1"));
        let regex = Regex::new(&regex_pattern).map_err(|e| Error::RuntimeError(
            format!("Cannot create regex for `{}`: {}", regex_pattern, e)))?;

        log::info!("Searching for files `{}`", glob_pattern);
        let glob_iter = glob::glob(&glob_pattern).map_err(|e| Error::RuntimeError(
            format!("Cannot create glob pattern for `{}`: {}", glob_pattern, e)))?;
        for entry in glob_iter {
            let path = entry.map_err(|e| Error::RuntimeError(
                format!("Glob search failed on pattern `{}`: {}", glob_pattern, e)))?;
            let path_str = path.as_os_str().to_str().ok_or_else(||
                Error::InvalidData(format!("Filename `{}` cannot be represented with Utf-8", path.display())))?;
            // cap.extract() returns tuple (full match, fixed-size array of all captures).
            let Some((_full, [locus])) = regex.captures(path_str).map(|cap| cap.extract()) else { continue };
            if !loci.insert(locus.to_owned()) {
                return Err(Error::InvalidInput(format!("Two sequence files with locus `{}` matched", locus)));
            }
            load_seqs_single(path_str, &output.replace(BRACKETS, locus), &mut seqs, &mut files)?;
        }
    }
    Ok((seqs, files))
}

/// Loads files from
fn load_seqs_list(list_filename: &Path) -> Result<SeqsFiles, Error> {
    let mut out_filenames = HashSet::default();
    let mut seqs = Vec::new();
    let mut files = Vec::new();
    log::info!("Loading sequences based on the list {}", ext::fmt::path(list_filename));
    for line in ext::sys::open(list_filename)?.lines() {
        let line = line.map_err(add_path!(list_filename))?;
        let line_split: smallvec::SmallVec<[&str; 3]> = line.trim().split_whitespace().collect();
        if line_split.len() != 2 {
            return Err(Error::InvalidData(format!("Invalid line `{}` in {}: exactly two filenames required",
                line, ext::fmt::path(list_filename))));
        }
        let seqs_filename = line_split[0];
        let out_filename = line_split[1];
        if !out_filenames.insert(out_filename.to_owned()) {
            return Err(Error::InvalidData(format!("Output filename {} appears twice in the list {}",
                out_filename, ext::fmt::path(list_filename))));
        }
        load_seqs_single(seqs_filename, out_filename, &mut seqs, &mut files)?;
    }
    Ok((seqs, files))
}

/// Based on the input arguments -s/-S, -o and -l, (list of list of sequences, output files).
fn load_seqs_and_outputs(args: &Args) -> Result<SeqsFiles, Error> {
    let (seqs, files) = match &args.seqs_list {
        Some(filename) => {
            validate_param!(args.seqs.is_empty() && args.seqs_all.is_none() && args.output.is_none(),
                "Filename list (-l) is mututally exclusive with -s/-S/-o");
            load_seqs_list(filename)?
        }
        None => {
            let out_filename = args.output.as_ref()
                .ok_or_else(|| Error::InvalidInput("Output path (-o) must be provided".to_string()))?;
            match &args.seqs_all {
                Some(filename) => {
                    validate_param!(args.seqs.is_empty(), "Arguments -s and -S are mutually exclusive");
                    load_seqs_all(filename, out_filename)?
                }
                None => {
                    validate_param!(!args.seqs.is_empty(), "Input sequences (-s or -S) must be provided");
                    load_seqs_glob(&args.seqs, out_filename)?
                }
            }
        }
    };
    assert_eq!(seqs.len(), files.len());
    if seqs.is_empty() {
        Err(Error::InvalidInput("No input sequences loaded".to_string()))
    } else {
        Ok((seqs, files))
    }
}

/// Loads k-mer counts, if needed, and creates recruitment targets.
fn build_targets(
    seqs: Vec<Vec<Vec<u8>>>,
    recr_params: recruit::Params,
    args: &Args,
) -> Result<recruit::Targets, Error>
{
    let mut target_builder = recruit::TargetBuilder::new(recr_params, args.thresh_kmer_count);
    let opt_kmer_getter = args.jf_counts.as_ref()
        .map(|filename| JfKmerGetter::new(args.jellyfish.clone(), filename.to_path_buf())).transpose()?;
    for locus_seqs in seqs.into_iter() {
        let kmer_counts = match &opt_kmer_getter {
            Some(kmer_getter) => kmer_getter.fetch(locus_seqs.clone())?,
            None => KmerCounts::new_zeroed(u32::from(args.minimizer_kw.0), locus_seqs.iter().map(Vec::len)),
        };
        let contig_set = ContigSet::new(Arc::new(ContigNames::empty()), locus_seqs, kmer_counts);
        target_builder.add(&contig_set, 0.0);
    }
    Ok(target_builder.finalize())
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let mut args = parse_args(argv)?;
    args.in_files.fill_from_inlist()?;
    args.update_from_preset()?;
    let args = args.validate()?;
    let recr_params = recruit::Params::new(args.minimizer_kw, args.match_frac, args.match_len)?;

    super::greet();
    let (seqs, files) = load_seqs_and_outputs(&args)?;
    let targets = build_targets(seqs, recr_params, &args)?;

    let timer = Instant::now();
    log::info!("Success! Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
