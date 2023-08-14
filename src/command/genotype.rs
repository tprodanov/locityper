use std::{
    fs,
    io::{self, BufRead},
    process::{Command, Stdio},
    cmp::{min, max},
    path::{Path, PathBuf},
    time::Instant,
    sync::Arc,
};
use colored::Colorize;
use const_format::{str_repeat, concatcp};
use fnv::{FnvHashSet, FnvHashMap};
use htslib::bam;
use crate::{
    err::{Error, validate_param, add_path},
    math::Ln,
    seq::{
        recruit, fastx, NamedSeq,
        contigs::{ContigId, ContigNames, ContigSet, Genotype},
        kmers::Kmer,
    },
    bg::{BgDistr, JsonSer, Technology, SequencingInfo},
    ext,
    model::{
        Params as AssgnParams,
        locs::AllAlignments,
        windows::ContigWindows,
        dp_cache::CachedDepthDistrs,
    },
    solvers::scheme::{self, Scheme, SchemeParams},
};
use super::paths;

struct Args {
    input: Vec<PathBuf>,
    database: Option<PathBuf>,
    output: Option<PathBuf>,
    subset_loci: FnvHashSet<String>,
    ploidy: u8,
    priors: Option<PathBuf>,

    interleaved: bool,
    threads: u16,
    rerun: super::Rerun,
    strobealign: PathBuf,
    minimap: PathBuf,
    samtools: PathBuf,
    seed: Option<u64>,
    debug: bool,

    recr_params: recruit::Params,
    assgn_params: AssgnParams,
    scheme_params: SchemeParams,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::new(),
            database: None,
            output: None,
            subset_loci: FnvHashSet::default(),
            ploidy: 2,
            priors: None,

            interleaved: false,
            threads: 8,
            rerun: super::Rerun::None,
            strobealign: PathBuf::from("strobealign"),
            minimap: PathBuf::from("minimap2"),
            samtools: PathBuf::from("samtools"),
            seed: None,
            debug: false,

            recr_params: Default::default(),
            assgn_params: Default::default(),
            scheme_params: Default::default(),
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
        validate_param!(self.ploidy > 0 && self.ploidy <= 11, "Ploidy ({}) must be within [1, 10]", self.ploidy);

        validate_param!(self.database.is_some(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.output.is_some(), "Output directory is not provided (see -o/--output)");
        self.samtools = ext::sys::find_exe(self.samtools)?;

        self.recr_params.validate()?;
        self.assgn_params.validate()?;
        Ok(self)
    }

    fn is_paired_end(&self) -> bool {
        self.input.len() == 2 || self.interleaved
    }
}

fn print_help() {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Genotype complex loci.".yellow());

    println!("\n{} {} genotype -i reads1.fq [reads2.fq] -d db -o out [arguments]",
        "Usage:".bold(), super::PROGRAM);

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Database directory (initialized with {} & {}).",
        "-d, --db".green(), "DIR".yellow(), concatcp!(super::PROGRAM, " create").underline(), "add".underline());
    println!("    {:KEY$} {:VAL$}  Output directory   (initialized with {}).",
        "-o, --output".green(), "DIR".yellow(), concatcp!(super::PROGRAM, " preproc").underline());
    println!("    {:KEY$} {:VAL$}  Input reads are interleaved.",
        "-^, --interleaved".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Optional: only analyze loci with names from this list.",
        "    --subset-loci".green(), "STR+".yellow());
    println!("    {:KEY$} {:VAL$}  Optional: genotype priors. Contains three columns:\n\
        {EMPTY}  <locus>  <genotype (through comma)>  <log10(prior)>.\n\
        {EMPTY}  Missing genotypes are removed from the analysis.",
        "    --priors".green(), "FILE".yellow());

    println!("\n{}", "Read recruitment:".bold());
    println!("    {:KEY$} {:VAL$}  Minimizer k-mer size (no larger than {}) [{}].",
        "-k, --recr-kmer".green(), "INT".yellow(), recruit::Minimizer::MAX_KMER_SIZE,
        super::fmt_def(defaults.recr_params.minimizer_k));
    println!("    {:KEY$} {:VAL$}  Take k-mers with smallest hash across {} consecutive k-mers [{}].",
        "-w, --recr-window".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(defaults.recr_params.minimizer_w));
    println!("    {:KEY$} {:VAL$}  Recruit single-end reads or read pairs with at least this fraction\n\
        {EMPTY}  of minimizers matching one of the targets [{}].",
        "-m, --matches-frac".green(), "FLOAT".yellow(),
        super::fmt_def_f64(f64::from(defaults.recr_params.matches_frac)));
    println!("    {:KEY$} {:VAL$}  Recruit reads in chunks of this size [{}].\n\
        {EMPTY}  May impact runtime in multi-threaded read recruitment.",
        "-c, --chunk-size".green(), "INT".yellow(), super::fmt_def(defaults.recr_params.chunk_size));

    println!("\n{}", "Model parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Solution ploidy [{}]. May be very slow for ploidy over 2.",
        "-p, --ploidy".green(), "INT".yellow(), super::fmt_def(defaults.ploidy));
    println!("    {:KEY$} {:VAL$}\n\
        {EMPTY}  Two p-value thresholds on edit distance:\n\
        {EMPTY}  for good alignments [{}], and for passable alignments [{}].",
        "    --edit-pval".green(), "FLOAT FLOAT".yellow(),
        super::fmt_def_f64(defaults.assgn_params.edit_pvals.0),
        super::fmt_def_f64(defaults.assgn_params.edit_pvals.1));
    println!("    {:KEY$} {:VAL$}  Ignore read alignments that are 10^{} times worse than\n\
        {EMPTY}  the best alignment [{}].",
        "-D, --prob-diff".green(), "FLOAT".yellow(), "FLOAT".yellow(),
        super::fmt_def_f64(Ln::to_log10(defaults.assgn_params.prob_diff)));
    println!("    {:KEY$} {:VAL$}  Unmapped read mate receives 10^{} penalty [{}].",
        "-U, --unmapped".green(), "FLOAT".yellow(), "FLOAT".yellow(),
        super::fmt_def_f64(Ln::to_log10(defaults.assgn_params.unmapped_penalty)));
    println!("    {:KEY$} {}\n\
        {EMPTY}  Calculate window weight based on the fraction of unique k-mers [{} {}].\n\
        {EMPTY}  * {} = breakpoint (0, 1). Weight at breakpoint is 1/2,\n\
        {EMPTY}  * {} = power [0.5, 50]. Regulates sigmoid slope (bigger - steeper).",
        "    --weight".green(), "FLOAT [FLOAT]".yellow(),
        super::fmt_def_f64(defaults.assgn_params.weight_breakpoint),
        super::fmt_def_f64(defaults.assgn_params.weight_power), "FLOAT_1".yellow(), "FLOAT_2".yellow());
    println!("    {:KEY$} {:VAL$}  Ignore reads and windows with weight under this value [{}].",
        "    --min-weight".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.assgn_params.min_weight));
    println!("    {:KEY$} {:VAL$}  Read depth likelihood contribution relative to\n\
        {EMPTY}  read alignment likelihoods [{}].",
        "-C, --dp-contrib".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.assgn_params.depth_contrib));
    println!("    {} {}  Compare window probability to have copy number 1 against two\n\
        {EMPTY}  alternative CN values [{} {}]. First in (0, 1), second > 1.",
        "-A, --alt-cn".green(), "FLOAT FLOAT".yellow(),
        super::fmt_def_f64(defaults.assgn_params.alt_cn.0), super::fmt_def_f64(defaults.assgn_params.alt_cn.1));
    println!("    {:KEY$} {:VAL$}  Use unpaired reads.",
        "    --use-unpaired".green(), super::flag());

    println!("\n{}", "Locus genotyping:".bold());
    println!("    {:KEY$} {:VAL$}  Use at most {} alignments [{}].",
        "    --max-alns".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(defaults.assgn_params.max_alns));
    println!("    {:KEY$} {:VAL$}  Solving stages through comma (see README) [{}].\n\
        {EMPTY}  Possible solvers: {}, {}, {}, {} and {}.",
        "-S, --stages".green(), "STR".yellow(), super::fmt_def(defaults.scheme_params.stages),
        "filter".yellow(), "greedy".yellow(), "anneal".yellow(), "highs".yellow(), "gurobi".yellow());
    println!("    {:KEY$} {:VAL$}  Score threshold for genotype pre-filtering [{}].\n\
        {EMPTY}  Values range from 0 (use all) to 1 (use best-score genotypes).",
        "    --score-thresh".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.assgn_params.score_thresh));
    println!("    {:KEY$} {:VAL$}  After each step, discard genotypes that have\n\
        {EMPTY}  smaller probability than 10^{} to be best [{}].",
        "    --prob-thresh".green(), "FLOAT".yellow(), "FLOAT".yellow(),
        super::fmt_def_f64(Ln::to_log10(defaults.assgn_params.prob_thresh)));
    println!("    {:KEY$} {:VAL$}  Minimum number of genotypes after each step [{}].",
        "    --min-gts".green(), "INT".yellow(), super::fmt_def(defaults.assgn_params.min_gts));
    println!("    {:KEY$} {:VAL$}  Number of attempts per step [{}].",
        "-a, --attempts".green(), "INT".yellow(), super::fmt_def(defaults.assgn_params.attempts));
    println!("    {:KEY$} {:VAL$}  Randomly move read coordinates by at most {} bp [{}].",
        "    --tweak".green(), "INT".yellow(), "INT".yellow(), "auto".cyan());
    println!("        {} {}, {} {}, {} {}, {} {}\n\
        {EMPTY}  Solver parameters (see README).",
        "--greedy".green(), "STR".yellow(),
        "--anneal".green(), "STR".yellow(),
        "--highs".green(), "STR".yellow(),
        "--gurobi".green(), "STR".yellow());

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads.to_string().cyan());
    println!("    {:KEY$} {:VAL$}  Rerun mode [{}]. Rerun all loci ({}); do not rerun\n\
        {EMPTY}  read recruitment ({}); do not rerun completed loci ({}).",
        "    --rerun".green(), "STR".yellow(), defaults.rerun.to_str().cyan(),
        "all".yellow(), "part".yellow(), "none".yellow());
    println!("    {:KEY$} {:VAL$}  Random seed. Ensures reproducibility for the same\n\
        {EMPTY}  input and product version.",
        "-s, --seed".green(), "INT".yellow());
    println!("    {:KEY$} {:VAL$}  Create more files with debug information.",
        "    --debug".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Strobealign executable [{}].",
        "    --strobealign".green(), "EXE".yellow(), super::fmt_def(defaults.strobealign.display()));
    println!("    {:KEY$} {:VAL$}  Minimap2 executable    [{}].",
        "    --minimap".green(), "EXE".yellow(), super::fmt_def(defaults.minimap.display()));
    println!("    {:KEY$} {:VAL$}  Samtools executable    [{}].",
        "    --samtools".green(), "EXE".yellow(), super::fmt_def(defaults.samtools.display()));

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
            Short('i') | Long("input") =>
                args.input = parser.values()?.take(2).map(|s| s.parse()).collect::<Result<_, _>>()?,
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),
            Long("subset-loci") => {
                for val in parser.values()? {
                    args.subset_loci.insert(val.parse()?);
                }
            }
            Long("priors") => args.priors = Some(parser.value()?.parse()?),

            Short('k') | Long("recr-kmer") => args.recr_params.minimizer_k = parser.value()?.parse()?,
            Short('w') | Long("recr-window") => args.recr_params.minimizer_w = parser.value()?.parse()?,
            Short('m') | Long("matches-frac") | Long("matches-fraction") =>
                args.recr_params.matches_frac = parser.value()?.parse()?,
            Short('c') | Long("chunk") | Long("chunk-size") => args.recr_params.chunk_size = parser.value()?.parse()?,

            Short('C') | Long("dp-contrib") | Long("depth-contrib") | Long("dp-contribution") =>
                args.assgn_params.depth_contrib = parser.value()?.parse()?,
            Short('A') | Long("alt-cn") =>
                args.assgn_params.alt_cn = (parser.value()?.parse()?, parser.value()?.parse()?),
            Short('D') | Long("prob-diff") => args.assgn_params.prob_diff = Ln::from_log10(parser.value()?.parse()?),
            Short('U') | Long("unmapped") =>
                args.assgn_params.unmapped_penalty = Ln::from_log10(parser.value()?.parse()?),
            Long("weight") => {
                let mut values = parser.values()?;
                args.assgn_params.weight_breakpoint = values.next()
                    .expect("At least one value must be present").parse()?;
                args.assgn_params.weight_power = values.next()
                    .map(|v| v.parse())
                    .transpose()?
                    .unwrap_or(2.0);
            }
            Long("min-weight") => args.assgn_params.min_weight = parser.value()?.parse()?,
            Long("edit-pval") | Long("edit-pvalue") =>
                args.assgn_params.edit_pvals = (
                    parser.value()?.parse()?,
                    parser.value()?.parse()?,
                ),
            Long("use-unpaired") => args.assgn_params.use_unpaired = true,

            Long("max-alns") | Long("max-alignments") => args.assgn_params.max_alns = parser.value()?.parse()?,
            Short('S') | Long("stages") => args.scheme_params.stages = parser.value()?.parse()?,
            Long("score-thresh") | Long("score-threshold") =>
                args.assgn_params.score_thresh = parser.value()?.parse()?,
            Long("prob-thresh") | Long("prob-threshold") =>
                args.assgn_params.prob_thresh = parser.value()?.parse()?,
            Long("min-gts") | Long("min-genotypes") =>
                args.assgn_params.min_gts = parser.value()?.parse()?,
            Short('a') | Long("attempts") => args.assgn_params.attempts = parser.value()?.parse()?,
            Long("tweak") => {
                let val = parser.value()?;
                args.assgn_params.tweak = if val == "auto" {
                    None
                } else {
                    Some(val.parse()?)
                };
            }
            Long("greedy") => args.scheme_params.greedy_params.push(parser.value()?.parse()?),
            Long("anneal") => args.scheme_params.anneal_params.push(parser.value()?.parse()?),
            Long("highs") => args.scheme_params.highs_params.push(parser.value()?.parse()?),
            Long("gurobi") => args.scheme_params.gurobi_params.push(parser.value()?.parse()?),

            Short('^') | Long("interleaved") => args.interleaved = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Long("rerun") => args.rerun = parser.value()?.parse()?,
            Short('s') | Long("seed") => args.seed = Some(parser.value()?.parse()?),
            Long("debug") => args.debug = true,
            Long("strobealign") => args.strobealign = parser.value()?.parse()?,
            Long("minimap") | Long("minimap2") => args.minimap = parser.value()?.parse()?,
            Long("samtools") => args.samtools = parser.value()?.parse()?,

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

/// Loads priors from a file with three columns: <locus> <genotype> <log10-prior>.
/// Repeated lines <locus> <genotype> are not allowed.
///
/// Output dictionary: `locus -> { genotype -> ln-prior }`.
fn load_priors(path: &Path) -> Result<FnvHashMap<String, FnvHashMap<String, f64>>, Error> {
    let mut res: FnvHashMap<String, FnvHashMap<String, f64>> = FnvHashMap::default();
    for line in ext::sys::open(path)?.lines() {
        let line = line.map_err(add_path!(path))?;
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<_> = line.trim_end().split_whitespace().collect();
        if cols.len() < 3 {
            return Err(Error::InvalidData(format!(
                "Cannot parse genotype priors: invalid line {:?} (expected at least 3 columns)", line)));
        }
        let locus = cols[0];
        let gt = cols[1];
        let prior: f64 = cols[2].parse()
            .map(Ln::from_log10)
            .map_err(|_| Error::InvalidData(format!("Cannot parse genotype priors: offending line {:?}", line)))?;
        if prior > 0.0 {
            return Err(Error::InvalidData(format!(
                "Cannot parse genotype priors: offending line {:?} (priors must be in log-10 space)", line)));
        }

        if let Some(old_prior) = res.entry(locus.to_owned()).or_default().insert(gt.to_owned(), prior) {
            if old_prior != prior {
                return Err(Error::InvalidData(format!(
                    "Cannot parse genotype priors: locus {} genotype {} contains two priors ({} and {})",
                    locus, gt, Ln::from_log10(old_prior), Ln::from_log10(prior))));
            }
        }
    }
    Ok(res)
}

struct LocusData {
    set: ContigSet,
    db_locus_dir: PathBuf,
    /// Output directory with locus data.
    out_dir: PathBuf,
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
    fn new(set: ContigSet, db_locus_dir: &Path, out_loci_dir: &Path) -> Self {
    let out_dir = out_loci_dir.join(set.tag());
        Self {
            db_locus_dir: db_locus_dir.to_owned(),
            tmp_reads_filename: out_dir.join("reads.tmp.fq"),
            reads_filename: out_dir.join("reads.fq"),
            tmp_aln_filename: out_dir.join("aln.tmp.bam"),
            aln_filename: out_dir.join("aln.bam"),
            set, out_dir,
        }
    }
}

/// Loads all loci from the database. If `subset_loci` is not empty, only loads loci that are contained in it.
fn load_loci(
    db_path: &Path,
    out_path: &Path,
    subset_loci: &FnvHashSet<String>,
    rerun: super::Rerun,
) -> Result<Vec<LocusData>, Error>
{
    log::info!("Loading database");
    let db_loci_dir = db_path.join(paths::LOCI_DIR);
    let out_loci_dir = out_path.join(paths::LOCI_DIR);
    ext::sys::mkdir(&out_loci_dir)?;

    let mut loci = Vec::new();
    let mut total_entries = 0;
    for entry in fs::read_dir(&db_loci_dir).map_err(add_path!(db_loci_dir))? {
        let entry = entry.map_err(add_path!(!))?;
        if !entry.file_type().map_err(add_path!(entry.path()))?.is_dir() {
            continue;
        }

        total_entries += 1;
        let path = entry.path();
        if let Some(name) = locus_name_matches(&path, subset_loci) {
            match ContigSet::load(name, &path.join(paths::LOCUS_FASTA), &path.join(paths::KMERS), ()) {
                Ok(set) => {
                    let locus_data = LocusData::new(set, &path, &out_loci_dir);
                    if rerun.need_analysis(&locus_data.out_dir)? {
                        loci.push(locus_data);
                    }
                },
                Err(e) => log::error!("Could not load locus information from {}: {:?}", ext::fmt::path(&path), e),
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
fn recruit_reads(loci: &[LocusData], args: &Args) -> Result<(), Error> {
    let filt_loci: Vec<&LocusData> = loci.iter()
        .filter(|locus| !locus.reads_filename.exists() && !locus.aln_filename.exists())
        .collect();
    if filt_loci.is_empty() {
        log::info!("Skipping read recruitment");
        return Ok(());
    }
    let n_filt_loci = filt_loci.len();
    if n_filt_loci < loci.len() {
        log::info!("Skipping read recruitment to {} loci", loci.len() - n_filt_loci);
    }

    log::info!("Generating recruitment targets");
    let mut targets = recruit::Targets::new(&args.recr_params);
    let mut writers = Vec::with_capacity(n_filt_loci);
    let mut total_seqs = 0;
    for locus in filt_loci.iter() {
        let fasta_path = locus.db_locus_dir.join(paths::LOCUS_FASTA_ALL);
        let mut fasta_reader = fastx::Reader::from_path(&fasta_path)?;
        let locus_all_seqs = fasta_reader.read_all().map_err(add_path!(fasta_path))?;
        total_seqs += locus_all_seqs.len();
        targets.add(locus_all_seqs.iter().map(NamedSeq::seq));
        // Output files with a large buffer (4 Mb).
        const BUFFER: usize = 4_194_304;
        writers.push(fs::File::create(&locus.tmp_reads_filename).map_err(add_path!(&locus.tmp_reads_filename))
            .map(|w| io::BufWriter::with_capacity(BUFFER, w))?);
    }
    log::info!("Collected {} minimizers across {} loci and {} sequences", targets.total_minimizers(),
        n_filt_loci, total_seqs);

    // Cannot put reader into a box, because `FastxRead` has a type parameter.
    if args.input.len() == 1 && !args.interleaved {
        let reader = fastx::Reader::from_path(&args.input[0])?;
        targets.recruit(reader, writers, args.threads).map_err(add_path!(args.input[0]))?;
    } else if args.interleaved {
        let reader = fastx::PairedEndInterleaved::new(fastx::Reader::from_path(&args.input[0])?);
        targets.recruit(reader, writers, args.threads).map_err(add_path!(args.input[0]))?;
    } else {
        let reader = fastx::PairedEndReaders::new(
            fastx::Reader::from_path(&args.input[0])?,
            fastx::Reader::from_path(&args.input[1])?);
        targets.recruit(reader, writers, args.threads).map_err(add_path!(args.input[0], args.input[1]))?;
    }
    for locus in filt_loci.iter() {
        fs::rename(&locus.tmp_reads_filename, &locus.reads_filename)
            .map_err(add_path!(locus.tmp_reads_filename, &locus.reads_filename))?;
    }
    Ok(())
}

/// Map reads with either Strobealign or Minimap2.
/// `n_locs`: number of possible alignment locations to examine for each read.
fn create_mapping_command(
    ref_path: &Path,
    reads_path: &Path,
    seq_info: &SequencingInfo,
    n_locs: &str,
    args: &Args,
) -> Command {
    let mut cmd: Command;
    if seq_info.technology() == Technology::Illumina {
        cmd = Command::new(&args.strobealign);
        // NOTE: Provide single input FASTQ, and do not use --interleaved option.
        // This is because Strobealign starts to connect read pairs in multiple ways,
        // which produces too large output, and takes a lot of time.
        cmd.args(&[
            "--no-progress",
            "-M", n_locs,   // Try as many secondary locations as possible.
            "-N", n_locs,   // Output as many secondary alignments as possible.
            "-S", "0.5",     // Try candidate sites with score >= 0.5 * best score.
            "-f", "0",       // Do not discard repetitive minimizers.
            "-k", "15",      // Use smaller minimizers to get more matches.
            "--eqx",         // Output X/= instead of M operations.
            "-r", &format!("{:.0}", seq_info.mean_read_len()),
            "-t", &args.threads.to_string()]);
    } else {
        cmd = Command::new(&args.minimap);
        cmd.args(&[
            "-a", // Output SAM format,
            "-x", seq_info.technology().minimap_preset(), // Set mapping preset.
            "-N", n_locs,   // Output as many secondary alignments as possible.
            "-f", "0",       // Do not discard repetitive minimizers.
            "--eqx",         // Output X/= instead of M operations.
            "-Y",            // Use soft clipping.
            "-t", &args.threads.to_string()]);
    }
    cmd.arg(&ref_path).arg(&reads_path)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    cmd
}

fn map_reads(locus: &LocusData, seq_info: &SequencingInfo, args: &Args) -> Result<(), Error> {
    if locus.aln_filename.exists() {
        log::info!("    Skipping read mapping");
        return Ok(());
    }

    let in_fasta = locus.db_locus_dir.join(paths::LOCUS_FASTA);
    log::info!("    Mapping reads to {}", ext::fmt::path(&in_fasta));
    // Output at most this number of alignments per read.
    let n_locs = min(30000, locus.set.len() * 30).to_string();
    let start = Instant::now();
    let mut mapping_cmd = create_mapping_command(&in_fasta, &locus.reads_filename, seq_info, &n_locs, args);
    let mut child = mapping_cmd.spawn().map_err(add_path!(!))?;
    let child_stdout = child.stdout.take().unwrap();

    let mut samtools_cmd = Command::new(&args.samtools);
    samtools_cmd.args(&["view", "-b"]) // Output BAM.
        .arg("-o").arg(&locus.tmp_aln_filename)
        .stdin(Stdio::from(child_stdout));

    log::debug!("    {} | {}", ext::fmt::command(&mapping_cmd), ext::fmt::command(&samtools_cmd));
    let samtools_output = samtools_cmd.output().map_err(add_path!(!))?;
    let bwa_output = child.wait_with_output().map_err(add_path!(!))?;
    log::debug!("    Finished in {}", ext::fmt::Duration(start.elapsed()));
    if !bwa_output.status.success() {
        return Err(Error::SubprocessFail(bwa_output));
    } else if !samtools_output.status.success() {
        return Err(Error::SubprocessFail(samtools_output));
    }
    fs::rename(&locus.tmp_aln_filename, &locus.aln_filename)
        .map_err(add_path!(locus.tmp_aln_filename, locus.aln_filename))?;
    fs::remove_file(&locus.reads_filename).map_err(add_path!(locus.reads_filename))?;
    Ok(())
}

/// Generates the list of genotypes and their priors.
fn generate_genotypes(
    contig_ids: &[ContigId],
    contigs: &ContigNames,
    opt_priors: Option<&FnvHashMap<String, f64>>,
    ploidy: usize,
) -> Result<(Vec<Genotype>, Vec<f64>), Error>
{
    if let Some(priors_map) = opt_priors {
        assert_eq!(contig_ids.len(), contigs.len(),
            "All contig IDs must be present when priors are provided");
        let mut genotypes = Vec::with_capacity(priors_map.len());
        let mut priors = Vec::with_capacity(priors_map.len());
        for (s, &prior) in priors_map.iter() {
            let gt = Genotype::parse(s, contigs)?;
            if prior > 0.0 || prior.is_nan() {
                return Err(Error::InvalidInput(format!("Invalid prior {} for genotype {}", prior, s)));
            } else if gt.ploidy() != ploidy {
                return Err(Error::InvalidInput(format!(
                    "Cannot load prior for genotype {} (expected ploidy {})", s, ploidy)));
            } else if prior.is_finite() {
                genotypes.push(gt);
                priors.push(prior);
            }
        }
        Ok((genotypes, priors))
    } else {
        let count = ext::vec::count_combinations_with_repl(contig_ids.len(), ploidy);
        let mut genotypes = Vec::with_capacity(count);
        ext::vec::gen_combinations_with_repl(&contig_ids, ploidy, |ids| genotypes.push(Genotype::new(ids, contigs)));
        Ok((genotypes, vec![0.0; count]))
    }
}

fn analyze_locus(
    locus: &LocusData,
    bg_distr: &BgDistr,
    scheme: &Arc<Scheme>,
    cached_distrs: &CachedDepthDistrs<'_>,
    opt_priors: Option<&FnvHashMap<String, f64>>,
    mut rng: ext::rand::XoshiroRng,
    args: &Args,
) -> Result<(), Error>
{
    log::info!("Analyzing {}", locus.set.tag());
    let timer = Instant::now();
    map_reads(locus, bg_distr.seq_info(), &args)?;

    log::info!("    Calculating read alignment probabilities");
    let bam_reader = bam::Reader::from_path(&locus.aln_filename)?;
    let contigs = locus.set.contigs();

    let mut contig_windows = if args.debug {
        let windows_filename = locus.out_dir.join("windows.bed.gz");
        let windows_writer = ext::sys::create_gzip(&windows_filename)?;
        ContigWindows::new_all(&locus.set, bg_distr.depth(), &args.assgn_params, windows_writer)
            .map_err(add_path!(windows_filename))?
    } else {
        ContigWindows::new_all(&locus.set, bg_distr.depth(), &args.assgn_params, io::sink()).map_err(add_path!(!))?
    };

    let mut all_alns = if args.debug {
        let reads_filename = locus.out_dir.join("reads.csv.gz");
        let reads_writer = ext::sys::create_gzip(&reads_filename)?;
        AllAlignments::load(bam_reader, contigs, bg_distr, &contig_windows, &args.assgn_params, reads_writer)?
    } else {
        AllAlignments::load(bam_reader, contigs, bg_distr, &contig_windows, &args.assgn_params, io::sink())?
    };
    if all_alns.len() > args.assgn_params.max_alns {
        let rate = args.assgn_params.max_alns as f64 / all_alns.len() as f64;
        all_alns.subsample(args.assgn_params.max_alns, &mut rng);
        ContigWindows::define_all_distributions(&mut contig_windows, &cached_distrs.subsample(rate));
    } else {
        ContigWindows::define_all_distributions(&mut contig_windows, cached_distrs);
    }

    let contig_ids: Vec<ContigId> = contigs.ids().collect();
    let (genotypes, priors) = generate_genotypes(&contig_ids, contigs, opt_priors, usize::from(args.ploidy))?;
    if genotypes.is_empty() {
        return Err(Error::RuntimeError(format!("No available genotypes for locus {}", locus.set.tag())));
    }

    let data = scheme::Data {
        scheme: Arc::clone(scheme),
        contigs: Arc::clone(contigs),
        assgn_params: args.assgn_params.clone(),
        debug: args.debug,
        threads: usize::from(args.threads),
        all_alns, genotypes, priors, contig_windows,
    };
    scheme::solve(data, &locus.out_dir, &mut rng)?;
    super::write_success_file(locus.out_dir.join(paths::SUCCESS))?;
    log::info!("    [{}] Successfully finished in {}", locus.set.tag(), ext::fmt::Duration(timer.elapsed()));
    Ok(())
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let mut args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();
    let db_dir = args.database.as_ref().unwrap();
    let out_dir = args.output.as_ref().unwrap();
    let priors = args.priors.as_ref().map(|path| load_priors(path)).transpose()?;

    let bg_path = out_dir.join(paths::BG_DIR).join(paths::BG_DISTR);
    let mut bg_stream = ext::sys::open(&bg_path)?;
    let mut bg_str = String::new();
    bg_stream.read_to_string(&mut bg_str).map_err(add_path!(bg_path))?;
    let mut bg_distr = BgDistr::load(&json::parse(&bg_str)?)?;
    args.assgn_params.set_tweak_size(bg_distr.depth().window_size())?;
    bg_distr.set_edit_pvals(args.assgn_params.edit_pvals);

    validate_param!(bg_distr.insert_distr().is_paired_end() == args.is_paired_end(),
        "Paired-end/Single-end status does not match background data");
    if bg_distr.seq_info().technology() == Technology::Illumina {
        args.strobealign = ext::sys::find_exe(args.strobealign)?;
    } else {
        args.minimap = ext::sys::find_exe(args.minimap)?;
    }

    let loci = load_loci(db_dir, out_dir, &args.subset_loci, args.rerun)?;
    recruit_reads(&loci, &args)?;

    let scheme = Arc::new(Scheme::create(&args.scheme_params)?);
    let cached_distrs = CachedDepthDistrs::new(&bg_distr, &args.assgn_params);
    let mut rng = ext::rand::init_rng(args.seed);
    for locus in loci.iter() {
        // Remove to get ownership of the locus priors.
        let locus_priors = priors.as_ref().and_then(|priors| priors.get(locus.set.tag()));
        if priors.is_some() && locus_priors.is_none() {
            log::warn!("Priors for locus {} are not found. Assuming equal priors", locus.set.tag());
        }

        let rng_clone = rng.clone();
        // Jump over 2^192 random numbers. This way, all loci have independent random numbers.
        rng.long_jump();
        analyze_locus(locus, &bg_distr, &scheme, &cached_distrs, locus_priors, rng_clone, &args)?;
    }
    log::info!("Success. Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
