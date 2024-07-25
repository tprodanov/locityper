use std::{
    fs,
    io::{self, BufRead},
    sync::Arc,
    time::Instant,
    path::{Path, PathBuf},
    cmp::max,
};
use fancy_regex::Regex;
use crate::{
    algo::{HashMap, HashSet},
    err::{Error, error, validate_param, add_path},
    ext::{
        self,
        fmt::{PrettyU32, PrettyU64},
    },
    seq::{
        recruit,
        Interval, ContigNames, ContigSet,
        fastx::{self, FastxRecord, SingleRecord},
        kmers::Kmer,
        counts::{KmerCount, KmerCounts, JfKmerGetter},
    },
    bg::Technology,
};
use colored::Colorize;

const SR_CHUNK_SIZE: usize = 10_000;
const LR_CHUNK_SIZE: usize = 300;

pub const DEFAULT_KMER_THRESH: KmerCount = 10;

struct Args {
    in_files: super::preproc::InputFiles,
    jf_counts: Option<PathBuf>,

    seqs: Vec<String>,
    seqs_all: Option<PathBuf>,
    seqs_list: Option<PathBuf>,
    output: Option<String>,
    regions: Option<PathBuf>,
    alt_contig_len: u32,

    minimizer_kw: (u8, u8),
    match_frac: f64,
    match_len: u32,
    chunk_size: usize,
    thresh_kmer_count: KmerCount,

    subsampling_rate: f64,
    seed: Option<u64>,

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
            alt_contig_len: 10_000_000,

            minimizer_kw: recruit::DEFAULT_MINIM_KW,
            match_frac: 0.5,
            match_len: recruit::DEFAULT_MATCH_LEN,
            chunk_size: SR_CHUNK_SIZE,
            thresh_kmer_count: DEFAULT_KMER_THRESH,

            subsampling_rate: 1.0,
            seed: None,

            threads: 8,
            jellyfish: PathBuf::from("jellyfish"),
        }
    }
}

impl Args {
    fn set_preset(&mut self, preset: String) -> Result<(), lexopt::Error> {
        let (tech, paired) = match &preset.to_lowercase() as &str {
            "illumina" | "illumina-pe" | "sr" | "sr-pe" => (Technology::Illumina, true),
            "illumina-se" | "sr-se" => (Technology::Illumina, false),
            "hifi" => (Technology::HiFi, false),
            "pacbio" | "pb" => (Technology::PacBio, false),
            "ont" | "nanopore" => (Technology::Nanopore, false),
            _ => return Err(lexopt::Error::UnexpectedValue {
                option: "-x/--preset".to_string(),
                value: preset.into(),
            }),
        };
        self.match_frac = tech.default_match_frac(paired);
        if tech != Technology::Illumina {
            self.chunk_size = LR_CHUNK_SIZE;
        }
        Ok(())
    }

    fn validate(mut self) -> crate::Result<Self> {
        self.in_files.validate(false)?;
        self.threads = max(self.threads, 1);

        validate_param!(0.0 < self.subsampling_rate && self.subsampling_rate <= 1.0,
            "Subsampling rate ({}) must be in (0, 1].", self.subsampling_rate);

        validate_param!(self.chunk_size > 0, "Chunk size must be positive");
        if self.jf_counts.is_some() {
            self.jellyfish = ext::sys::find_exe(self.jellyfish)?;
        }

        if self.in_files.has_indexed_alignment() && self.regions.is_none() {
            log::warn!("It is much more efficient to use indexed BAM/CRAM file (-a) together target regions (-R).");
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
    println!("    {:KEY$} {:VAL$}  Recruit reads mapped to contigs shorter than this [{}].\n\
        {EMPTY}  Only used together with {} for mapped and indexed BAM/CRAM files.",
        "    --alt-contig".green(), "INT".yellow(), super::fmt_def(PrettyU32(defaults.alt_contig_len)), "-R".green());

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
    println!("    {:KEY$} {:VAL$}  Only use k-mers that appear less than {} times in the\n\
        {EMPTY}  reference genome [{}]. Requires {} argument.",
        "-t, --kmer-thresh".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(defaults.thresh_kmer_count),
        "-j".green());
    println!("    {:KEY$} {:VAL$}  Recruit reads in chunks of this size [{}].\n\
        {EMPTY}  Impacts runtime in multi-threaded read recruitment.",
        "-c, --chunk-size".green(), "INT".yellow(),
        super::fmt_def(PrettyU64(defaults.chunk_size as u64)));
    println!("    {:KEY$} {:VAL$}  Parameter preset (illumina|illumina-SE|hifi|pacbio|ont).\n\
        {EMPTY}  Modifies {}, {}, see {} for default values.\n\
        {EMPTY}  Additionally, changes {} to {} for long reads.",
        "-x, --preset".green(), "STR".yellow(),
        "-m".green(), "-M".green(), "locityper genotype".underline(), "-c".green(), super::fmt_def(LR_CHUNK_SIZE));

    println!("\n{}", "Subsampling:".bold());
    println!("    {:KEY$} {:VAL$}  Before recruitment, subsample reads at this rate [{}].",
        "    --subsample".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.subsampling_rate));
    println!("    {:KEY$} {:VAL$}  Subsampling seed (optional). Ensures reproducibility\n\
        {EMPTY}  for the same input and program version.",
        "    --seed".green(), "INT".yellow());

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
            Long("alt-contig") => args.alt_contig_len = parser.value()?.parse::<PrettyU32>()?.get(),

            Short('m') | Long("minimizer") | Long("minimizers") =>
                args.minimizer_kw = (parser.value()?.parse()?, parser.value()?.parse()?),
            Short('M') | Long("match-frac") | Long("match-fraction") => args.match_frac = parser.value()?.parse()?,
            Short('x') | Long("preset") => args.set_preset(parser.value()?.parse()?)?,
            Short('t') | Long("kmer-thresh") | Long("kmer-threshold") =>
                args.thresh_kmer_count = parser.value()?.parse()?,
            Short('L') | Long("match-len") | Long("match-length") =>
                args.match_len = parser.value()?.parse::<PrettyU32>()?.get(),
            Short('c') | Long("chunk") | Long("chunk-size") =>
                args.chunk_size = parser.value()?.parse::<PrettyU64>()?.get() as usize,

            Long("subsample") => args.subsampling_rate = parser.value()?.parse()?,
            Long("seed") => args.seed = Some(parser.value()?.parse()?),

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
fn create_output(filename: impl AsRef<Path>) -> crate::Result<BufFile> {
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
fn load_seqs_all(seqs_all: &Path, output: &str) -> crate::Result<SeqsFiles> {
    validate_param!(output.contains(BRACKETS), "When using -S, output path ({}) must contain {}", output, BRACKETS);
    log::info!("Loading all sequences from {}", ext::fmt::path(seqs_all));

    // locus -> (seqs, output).
    let mut map = HashMap::default();
    let mut fasta_reader = fastx::Reader::from_path(seqs_all)?;
    let mut record = FastxRecord::default();
    while fasta_reader.read_next_standardized(&mut record)? {
        let name = std::str::from_utf8(record.name())
            .map_err(|_| Error::Utf8("read name", record.name().to_vec()))?;
        let (locus, _) = name.split_once('*').ok_or_else(||
                error!(InvalidData, "Invalid record name `{}` for -S {}. Symbol `*` required",
                name, ext::fmt::path(seqs_all)))?;
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
    seqs_filename: impl AsRef<Path>,
    output_filename: impl AsRef<Path>,
    seqs: &mut Vec<Vec<Vec<u8>>>,
    files: &mut Vec<BufFile>,
) -> crate::Result<()>
{
    let seqs_filename = seqs_filename.as_ref();
    let output_filename = output_filename.as_ref();
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
fn load_seqs_glob(patterns: &[String], output: &str) -> crate::Result<SeqsFiles> {
    let mut seqs = Vec::new();
    let mut files = Vec::new();
    if patterns.len() == 1 && !patterns[0].contains(BRACKETS) {
        let seq_path = &patterns[0];
        validate_param!(!output.contains(BRACKETS),
            "Single input -s {} does not allow output files with `{}` ({})",
            ext::fmt::path(seq_path), BRACKETS, output);
        log::info!("Loading sequences from {}", ext::fmt::path(seq_path));
        load_seqs_single(seq_path, output, &mut seqs, &mut files)?;
        return Ok((seqs, files));
    }

    let mut loci = HashSet::default();
    validate_param!(output.contains(BRACKETS), "When using -s, output path ({}) must contain {}", output, BRACKETS);
    for pattern in patterns {
        validate_param!(pattern.contains(BRACKETS), "Paths to sequences ({}) must contain brackets `{}`",
            ext::fmt::path(pattern), BRACKETS);
        validate_param!(!pattern.contains("*"), "Paths to sequences ({}) must not contain `*`",
            ext::fmt::path(pattern));
        let can_pattern = fs::canonicalize(pattern).ok()
            .and_then(|p| p.into_os_string().into_string().ok())
            .unwrap_or_else(|| pattern.to_owned());
        let glob_pattern = can_pattern.replace(BRACKETS, "*");
        let regex_pattern = format!("^{}$", can_pattern.replacen(BRACKETS, r"([^/\\]+)", 1).replace(BRACKETS, "\\1"));
        let regex = Regex::new(&regex_pattern).map_err(|e| error!(RuntimeError,
            "Cannot create regex for `{}`: {}", regex_pattern, e))?;

        log::info!("Searching for files `{}`", glob_pattern);
        let glob_iter = glob::glob(&glob_pattern).map_err(|e| error!(RuntimeError,
            "Cannot create glob pattern for `{}`: {}", glob_pattern, e))?;
        for entry in glob_iter {
            let path = entry.map_err(|e| error!(RuntimeError,
                "Glob search failed on pattern `{}`: {}", glob_pattern, e))?;
            let path_str = path.as_os_str().to_str().ok_or_else(||
                error!(InvalidData, "Filename `{}` cannot be represented with Utf-8", path.display()))?;

            let Some(cap1) = regex.captures(path_str)
                .map_err(|e| error!(RuntimeError,
                    "Failed to run regex `{}` on `{}`: {}", regex_pattern, path_str, e))?
                .and_then(|cap| cap.get(1))
                else { continue };
            let locus = cap1.as_str();
            if !loci.insert(locus.to_owned()) {
                return Err(error!(InvalidInput, "Two sequence files with locus `{}` matched", locus));
            }
            load_seqs_single(path_str, &output.replace(BRACKETS, locus), &mut seqs, &mut files)?;
        }
    }
    Ok((seqs, files))
}

/// Loads files from
fn load_seqs_list(list_filename: &Path) -> crate::Result<SeqsFiles> {
    let dirname = list_filename.parent();
    let mut out_filenames = HashSet::default();
    let mut seqs = Vec::new();
    let mut files = Vec::new();
    log::info!("Loading sequences based on the list {}", ext::fmt::path(list_filename));
    for line in ext::sys::open(list_filename)?.lines() {
        let line = line.map_err(add_path!(list_filename))?;
        let line_split: smallvec::SmallVec<[&str; 3]> = line.trim().split_whitespace().collect();
        if line_split.len() != 2 {
            return Err(error!(InvalidData, "Invalid line `{}` in {}: exactly two filenames required",
                line, ext::fmt::path(list_filename)));
        }
        let seqs_filename = ext::sys::add_dir(dirname, line_split[0]);
        let out_filename = ext::sys::add_dir(dirname, line_split[1]);
        if !out_filenames.insert(out_filename.to_owned()) {
            return Err(error!(InvalidData, "Output filename {} appears twice in the list {}",
                ext::fmt::path(out_filename), ext::fmt::path(list_filename)));
        }
        load_seqs_single(seqs_filename, out_filename, &mut seqs, &mut files)?;
    }
    Ok((seqs, files))
}

/// Based on the input arguments -s/-S, -o and -l, (list of list of sequences, output files).
fn load_seqs_and_outputs(args: &Args) -> crate::Result<SeqsFiles> {
    let (seqs, files) = match &args.seqs_list {
        Some(filename) => {
            validate_param!(args.seqs.is_empty() && args.seqs_all.is_none() && args.output.is_none(),
                "Filename list (-l) is mututally exclusive with -s/-S/-o");
            load_seqs_list(filename)?
        }
        None => {
            let out_filename = args.output.as_ref()
                .ok_or_else(|| error!(InvalidInput, "Output path (-o) must be provided"))?;
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
        Err(error!(InvalidInput, "No input sequences loaded"))
    } else {
        Ok((seqs, files))
    }
}

/// Creates recruitment targets based on the input sequences.
///
/// seqs: collection of loci, each with collection of sequences (each Vec<u8>).
///
/// If Jellyfish counts are provided, masks minimizers based on the k-mer content.
fn build_targets(
    seqs: impl IntoIterator<Item = Vec<Vec<u8>>>,
    recr_params: recruit::Params,
    jf_counts: Option<&PathBuf>,
    jellyfish: &PathBuf,
) -> crate::Result<recruit::Targets>
{
    let minimizer_k = recr_params.minimizer_k();
    let mut target_builder = recruit::TargetBuilder::new(recr_params);
    let opt_kmer_getter = match jf_counts {
        Some(filename) => Some(JfKmerGetter::new(jellyfish.clone(), filename.to_path_buf())?),
        None => {
            log::warn!("Consider providing reference k-mer counts (-j)");
            None
        }
    };
    for locus_seqs in seqs.into_iter() {
        let kmer_counts = match &opt_kmer_getter {
            Some(kmer_getter) => kmer_getter.fetch(locus_seqs.clone())?,
            None => KmerCounts::new_zeroed(u32::from(minimizer_k), locus_seqs.iter().map(Vec::len)),
        };
        let contig_set = ContigSet::new(Arc::new(ContigNames::empty()), locus_seqs, kmer_counts);
        target_builder.add(&contig_set, 0.0);
    }
    Ok(target_builder.finalize())
}

fn load_target_regions(contigs: &Arc<ContigNames>, args: &Args) -> crate::Result<Vec<Interval>> {
    match &args.regions {
        Some(filename) => {
            let mut intervals = Vec::new();
            let file = ext::sys::open(filename)?;
            for line in file.lines() {
                let line = line.map_err(add_path!(filename))?;

                match Interval::parse_bed(&mut line.trim().split('\t'), &contigs) {
                    Ok(interv) => intervals.push(interv),
                    Err(Error::ParsingError(e)) => log::error!(
                        "Cannot parse BED line `{}`, the region may be absent from the BAM/CRAM file: {}",
                        line.trim(), e),
                    Err(e) => log::error!("Cannot parse BED line `{}`: {}", line.trim(), e.display()),
                }
            }
            const PADDING: u32 = 1000;
            const MERGE_DISTANCE: u32 = 1000;
            Ok(super::genotype::postprocess_targets(intervals, contigs, PADDING, args.alt_contig_len, MERGE_DISTANCE))
        }
        None => Ok(contigs.ids().map(|id| Interval::full_contig(Arc::clone(contigs), id)).collect()),
    }
}

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let mut args = parse_args(argv)?;
    let timer = Instant::now();
    args.in_files.fill_from_inlist()?;
    let args = args.validate()?;
    super::greet();

    let sampling = if args.subsampling_rate < 1.0 {
        Some((args.subsampling_rate, ext::rand::init_rng(args.seed)))
    } else { None };

    let recr_params = recruit::Params::new(args.minimizer_kw, args.match_frac, args.match_len, args.thresh_kmer_count)?;
    let (seqs, files) = load_seqs_and_outputs(&args)?;
    let targets = build_targets(seqs, recr_params, args.jf_counts.as_ref(), &args.jellyfish)?;
    super::genotype::recruit_to_targets(
        &targets, &args.in_files, files, None, args.threads, args.chunk_size, sampling,
        |contigs| load_target_regions(contigs, &args))?;

    log::info!("Success! Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
