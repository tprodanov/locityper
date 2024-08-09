use std::{
    fs,
    io::{self, BufRead},
    process::{Command, Stdio},
    cmp::{min, max, Ordering},
    path::{Path, PathBuf},
    time::Instant,
    sync::Arc,
};
use colored::Colorize;
use const_format::{str_repeat, concatcp};
use htslib::bam::{self, Read as BamRead};
use crate::{
    err::{Error, error, validate_param, add_path},
    math::{
        Ln,
        distr::WithQuantile,
    },
    seq::{
        recruit, fastx, div, Interval,
        contigs::{ContigId, ContigNames, ContigSet, Genotype},
        kmers::Kmer,
        counts::KmerCount,
    },
    bg::{
        TECHNOLOGIES, BgDistr, Technology, SequencingInfo,
        err_prof::EditDistCache,
    },
    ext::{
        self,
        rand::XoshiroRng,
        fmt::{PrettyU32, PrettyU64, PrettyUsize},
    },
    model::{
        Params as AssgnParams,
        locs::AllAlignments,
        windows::{ContigInfo, WeightCalculator},
        distr_cache::DistrCache,
    },
    solvers::scheme::{self, Scheme, SchemeParams},
    algo::{HashSet, HashMap},
};
use super::{paths, DebugLvl, preproc::InputFiles};

pub const DEFAULT_CHUNK_LENGTH: u64 = 3_000_000;

struct Args {
    in_files: InputFiles,
    preproc: Option<PathBuf>,
    databases: Vec<PathBuf>,
    output: Option<PathBuf>,

    subset_loci: HashSet<String>,
    ploidy: u8,
    priors: Option<PathBuf>,

    threads: u16,
    rerun: super::Rerun,
    strobealign: PathBuf,
    minimap: PathBuf,
    samtools: PathBuf,
    seed: Option<u64>,
    debug: DebugLvl,

    minimizer_kw: (u8, u8),
    match_frac: Option<f64>,
    match_len: u32,
    thresh_kmer_count: KmerCount,
    chunk_length: u64,

    assgn_params: AssgnParams,
    scheme_params: SchemeParams,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            in_files: Default::default(),
            preproc: None,
            databases: Vec::new(),
            output: None,

            subset_loci: HashSet::default(),
            ploidy: 2,
            priors: None,

            threads: 8,
            rerun: super::Rerun::None,
            strobealign: PathBuf::from("strobealign"),
            minimap: PathBuf::from("minimap2"),
            samtools: PathBuf::from("samtools"),
            seed: None,
            debug: DebugLvl::None,

            minimizer_kw: recruit::DEFAULT_MINIM_KW,
            match_frac: None,
            match_len: recruit::DEFAULT_MATCH_LEN,
            thresh_kmer_count: 5,
            chunk_length: DEFAULT_CHUNK_LENGTH,

            assgn_params: Default::default(),
            scheme_params: Default::default(),
        }
    }
}

impl Args {
    /// Validate arguments, modifying some, if needed.
    fn validate(mut self) -> crate::Result<Self> {
        self.in_files.validate(false)?;
        self.threads = max(self.threads, 1);

        validate_param!(self.preproc.is_some(), "Preprocessing directory is not provided (see -p/--preproc)");
        validate_param!(!self.databases.is_empty(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.output.is_some(), "Output directory is not provided (see -o/--output)");
        self.samtools = ext::sys::find_exe(self.samtools)?;
        self.assgn_params.validate()?;
        Ok(self)
    }
}

fn fmt_def_match_frac() -> String {
    super::describe_defaults(
        TECHNOLOGIES.iter()
            .flat_map(|&tech| [(tech, false), (tech, true)].into_iter())
            .filter(|(tech, paired_end)| !paired_end || tech.paired_end_allowed()),
        |(tech, paired_end)|
            format!("{}{}", tech, if *paired_end { "-PE" } else if tech.paired_end_allowed() { "-SE" } else { "" }),
        |(tech, paired_end)| tech.default_match_frac(*paired_end)
    )
}

fn print_help(extended: bool) {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Genotype complex loci.".yellow());

    println!("\n{}", "Usage:".bold());
    println!("    {} genotype (-i reads1.fq [reads2.fq] | -a reads.bam [--no-index] | -I in-list) \\", super::PROGRAM);
    println!("        -p preproc -d db -o out [args]");
    if !extended {
        println!("\nThis is a {} help message. Please use {} to see the full help.",
            "short".red(), "-H/--full-help".green());
    }

    println!("\n{}  (please see {} for more information on {}/{}/{} arguments)",
        "Input/output arguments:".bold(), "README".italic(), "-i".green(), "-a".green(), "-I".green());
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
    println!("    {:KEY$} {:VAL$}  Preprocessed dataset information (see {}).",
        "-p, --preproc".green(), "DIR".yellow(), concatcp!(super::PROGRAM, " preproc").underline());
    println!("    {:KEY$} {:VAL$}  Database directory (see {}).\n\
        {EMPTY}  Multiple databases allowed, but must contain unique loci names.",
        "-d, --database[s]".green(), "DIR+".yellow(), concatcp!(super::PROGRAM, " add").underline());
    println!("    {:KEY$} {:VAL$}  Output directory.",
        "-o, --output".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Output BAM files for {} best genotypes [{}].",
        "-O, --out-bams".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(defaults.assgn_params.out_bams));

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Interleaved paired-end reads in single input file.",
        "-^, --interleaved".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Use input BAM/CRAM file ({}) without index: goes over all reads.\n\
        {EMPTY}  Single-end and paired-end interleaved ({}) data is allowed.",
        "    --no-index".green(), super::flag(), "-a".green(), "-^".green());
    if extended {
        println!("    {:KEY$} {:VAL$}  Optional: only analyze loci with names from this list.",
            "    --subset-loci".green(), "STR+".yellow());
        println!("    {:KEY$} {:VAL$}  Optional: genotype priors. Contains three columns:\n\
            {EMPTY}  <locus>  <genotype (through comma)>  <log10(prior)>.\n\
            {EMPTY}  Missing genotypes are removed from the analysis.",
            "    --priors".green(), "FILE".yellow());

        println!("\n{}", "Read recruitment:".bold());
        println!("    {}  {}  Use k-mers of size {} (<= {}) with smallest hash\n\
            {EMPTY}  across {} consecutive k-mers [{} {}].",
            "-m, --minimizer".green(), "INT INT".yellow(),
            "INT_1".yellow(), recruit::Minimizer::MAX_KMER_SIZE, "INT_2".yellow(),
            super::fmt_def(defaults.minimizer_kw.0), super::fmt_def(defaults.minimizer_kw.1));
        println!("    {:KEY$} {:VAL$}  Minimal fraction of minimizers that need to match reference.\n\
            {EMPTY}  Default: {}.",
            "-M, --match-frac".green(), "FLOAT".yellow(), fmt_def_match_frac());
        println!("    {:KEY$} {:VAL$}  Recruit long reads with a matching subregion of this length [{}].",
            "-L, --match-len".green(), "INT".yellow(),
            super::fmt_def(defaults.match_len));
        println!("    {:KEY$} {:VAL$}  Only use k-mers that appear less than {} times off target [{}].",
            "    --kmer-thresh".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(defaults.thresh_kmer_count));
        println!("    {:KEY$} {:VAL$}  Recruit reads in chunks of this sum length [{}].\n\
            {EMPTY}  Impacts runtime in multi-threaded read recruitment.",
            "-c, --chunk-len".green(), "INT".yellow(),
            super::fmt_def(PrettyU64(defaults.chunk_length)));

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
            {EMPTY}  * {} = power [0.5, 50]. Regulates sigmoid slope (bigger - steeper).\n\
            {EMPTY}  Use {} to disable weights.",
            "    --kmers-weight".green(), "FLOAT FLOAT | off".yellow(),
            super::fmt_def_f64(defaults.assgn_params.kmers_weight_calc.as_ref().unwrap().breakpoint()),
            super::fmt_def_f64(defaults.assgn_params.kmers_weight_calc.as_ref().unwrap().power()),
            "FLOAT_1".yellow(), "FLOAT_2".yellow(), "off".yellow());
        println!("    {:KEY$} {}\n\
            {EMPTY}  Calculate window weight based on the linguistic complexity\n\
            {EMPTY}  of the window sequence [{} {}]. Same format as {}.",
            "    --compl-weight".green(), "FLOAT FLOAT | off".yellow(),
            super::fmt_def_f64(defaults.assgn_params.compl_weight_calc.as_ref().unwrap().breakpoint()),
            super::fmt_def_f64(defaults.assgn_params.compl_weight_calc.as_ref().unwrap().power()),
            "--kmers-weight".green());
        println!("    {:KEY$} {:VAL$}  Ignore windows with weight under this value [{}].",
            "    --min-weight".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.assgn_params.min_weight));
        println!("    {:KEY$} {:VAL$}  Use reads/read pairs that have at least this number\n\
            {EMPTY}  of non-overlapping k-mers, unique to the target locus [{}].",
            "    --min-kmers".green(), "INT".yellow(), super::fmt_def(defaults.assgn_params.min_unique_kmers));
        println!("    {:KEY$} {:VAL$}  Likelihood skew (-1, 1) [{}]. Negative: alignment probabilities\n\
            {EMPTY}  matter more; Positive: read depth matters more.",
            "    --skew".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.assgn_params.lik_skew));
        println!("    {} {}  Compare window probability to have copy number 1 against two\n\
            {EMPTY}  alternative CN values [{} {}]. First in (0, 1), second > 1.",
            "-A, --alt-cn".green(), "FLOAT FLOAT".yellow(),
            super::fmt_def_f64(defaults.assgn_params.alt_cn.0), super::fmt_def_f64(defaults.assgn_params.alt_cn.1));

        println!("\n{}", "Locus genotyping:".bold());
        println!("    {:KEY$} {:VAL$}  Solving stages through comma (see {}) [{}].\n\
            {EMPTY}  Possible solvers: {}, {}, {}, {} and {}.",
            "-S, --stages".green(), "STR".yellow(), "README".italic(), super::fmt_def(defaults.scheme_params.stages),
            "filter".yellow(), "greedy".yellow(), "anneal".yellow(), "highs".yellow(), "gurobi".yellow());
        println!("    {:KEY$} {:VAL$}  During pre-filtering, discard genotypes that have 10^{}\n\
            {EMPTY}  worse alignment probability than the best genotype [{}].",
            "    --filt-diff".green(), "FLOAT".yellow(), "FLOAT".yellow(),
            super::fmt_def_f64(Ln::to_log10(defaults.assgn_params.filt_diff)));
        println!("    {:KEY$} {:VAL$}  After each step, discard genotypes that have\n\
            {EMPTY}  smaller probability than 10^{} to be best [{}].",
            "    --prob-thresh".green(), "FLOAT".yellow(), "FLOAT".yellow(),
            super::fmt_def_f64(Ln::to_log10(defaults.assgn_params.prob_thresh)));
        println!("    {:KEY$} {:VAL$}  Minimum number of genotypes after each step [{}].",
            "    --min-gts".green(), "INT".yellow(),
            super::fmt_def(PrettyUsize(defaults.assgn_params.min_gts)));
        println!("    {:KEY$} {:VAL$}  Number of attempts per step [{}].",
            "    --attempts".green(), "INT".yellow(), super::fmt_def(defaults.assgn_params.attempts));
        println!("    {:KEY$} {:VAL$}  Randomly move read coordinates by at most {} bp [{}].",
            "-t, --tweak".green(), "INT".yellow(), "INT".yellow(), super::fmt_def("auto"));
        println!("    {:KEY$} {:VAL$}  Normalize depth likelihoods based on sum window weight across\n\
            {EMPTY}  genotype, raised to this power (0 - no normalization) [{}].",
            "    --depth-norm".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.assgn_params.depth_norm_power));
        println!("        {} {}, {} {}, {} {}, {} {}\n\
            {EMPTY}  Solver parameters (see {}).",
            "--greedy".green(), "STR".yellow(), "--anneal".green(), "STR".yellow(),
            "--highs".green(), "STR".yellow(), "--gurobi".green(), "STR".yellow(), "README".italic());
    }

    println!("\n{}", "Execution arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), super::fmt_def(defaults.threads));
    println!("    {:KEY$} {:VAL$}  Rerun mode [{}]. Rerun all loci ({}); do not rerun\n\
        {EMPTY}  read recruitment ({}); do not rerun completed loci ({}).",
        "    --rerun".green(), "STR".yellow(), super::fmt_def(defaults.rerun),
        "all".yellow(), "part".yellow(), "none".yellow());
    println!("    {:KEY$} {:VAL$}  Random seed. Ensures reproducibility for the same\n\
        {EMPTY}  input and program version.",
        "-s, --seed".green(), "INT".yellow());
    println!("    {:KEY$} {:VAL$}  Save debug CSV files: 0 = none, 1 = some, 2 = all.",
        "    --debug".green(), "[INT]".yellow());
    if extended {
        println!("    {:KEY$} {:VAL$}  Strobealign executable [{}].",
            "    --strobealign".green(), "EXE".yellow(), super::fmt_def(defaults.strobealign.display()));
        println!("    {:KEY$} {:VAL$}  Minimap2 executable    [{}].",
            "    --minimap".green(), "EXE".yellow(), super::fmt_def(defaults.minimap.display()));
        println!("    {:KEY$} {:VAL$}  Samtools executable    [{}].",
            "    --samtools".green(), "EXE".yellow(), super::fmt_def(defaults.samtools.display()));
    }

    println!("\n{}", "Other arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Show this help message.", "-h, --help".green(), "");
    println!("    {:KEY$} {:VAL$}  Show {} help message.", "-H, --full-help".green(), "", "extended".red());
    println!("    {:KEY$} {:VAL$}  Show version.", "-V, --version".green(), "");
}

fn parse_args(argv: &[String]) -> crate::Result<Args> {
    if argv.is_empty() {
        print_help(false);
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
            Short('p') | Long("preproc") | Long("preprocessing") => args.preproc = Some(parser.value()?.parse()?),
            Short('d') | Long("db") | Long("database") | Long("databases") => {
                for val in parser.values()? {
                    args.databases.push(val.parse()?);
                }
            }
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),
            Long("subset-loci") | Long("loci-subset") => {
                for val in parser.values()? {
                    args.subset_loci.insert(val.parse()?);
                }
            }
            Long("priors") => args.priors = Some(parser.value()?.parse()?),

            Short('m') | Long("minimizer") | Long("minimizers") =>
                args.minimizer_kw = (parser.value()?.parse()?, parser.value()?.parse()?),
            Short('M') | Long("match-frac") | Long("match-fraction") =>
                args.match_frac = Some(parser.value()?.parse()?),
            Short('L') | Long("match-len") | Long("match-length") =>
                args.match_len = parser.value()?.parse::<PrettyU32>()?.get(),
            Long("kmer-thresh") => args.thresh_kmer_count = parser.value()?.parse()?,
            Short('c') | Long("chunk") | Long("chunk-len") =>
                args.chunk_length = parser.value()?.parse::<PrettyU64>()?.get(),

            Long("skew") => args.assgn_params.lik_skew = parser.value()?.parse()?,
            Short('A') | Long("alt-cn") =>
                args.assgn_params.alt_cn = (parser.value()?.parse()?, parser.value()?.parse()?),
            Short('D') | Long("prob-diff") => args.assgn_params.prob_diff = Ln::from_log10(parser.value()?.parse()?),
            Short('U') | Long("unmapped") =>
                args.assgn_params.unmapped_penalty = Ln::from_log10(parser.value()?.parse()?),
            Long("kmers-weight") => {
                let first = parser.value()?;
                if first == "off" {
                    args.assgn_params.kmers_weight_calc = None;
                } else {
                    args.assgn_params.kmers_weight_calc = Some(
                        WeightCalculator::new(first.parse()?, parser.value()?.parse()?)?);
                }
            }
            Long("compl-weight") => {
                let first = parser.value()?;
                if first == "off" {
                    args.assgn_params.compl_weight_calc = None;
                } else {
                    args.assgn_params.compl_weight_calc = Some(
                        WeightCalculator::new(first.parse()?, parser.value()?.parse()?)?);
                }
            }
            Long("min-weight") => args.assgn_params.min_weight = parser.value()?.parse()?,
            Long("min-kmers") => args.assgn_params.min_unique_kmers = parser.value()?.parse()?,
            Long("edit-pval") | Long("edit-pvalue") =>
                args.assgn_params.edit_pvals = (
                    parser.value()?.parse()?,
                    parser.value()?.parse()?,
                ),

            Short('S') | Long("stages") => args.scheme_params.stages = parser.value()?.parse()?,
            Long("filt-diff") | Long("filt-difference") | Long("filter-diff") =>
                args.assgn_params.filt_diff = Ln::from_log10(parser.value()?.parse()?),
            Long("prob-thresh") | Long("prob-threshold") =>
                args.assgn_params.prob_thresh = Ln::from_log10(parser.value()?.parse()?),
            Long("min-gts") | Long("min-genotypes") =>
                args.assgn_params.min_gts = parser.value()?.parse::<PrettyUsize>()?.get(),
            Long("attempts") => args.assgn_params.attempts = parser.value()?.parse()?,
            Short('t') | Long("tweak") => {
                let val = parser.value()?;
                args.assgn_params.tweak = if val == "auto" {
                    None
                } else {
                    Some(val.parse()?)
                };
            }
            Long("depth-norm") => args.assgn_params.depth_norm_power = parser.value()?.parse()?,
            Long("greedy") => args.scheme_params.greedy_params.push(parser.value()?.parse()?),
            Long("anneal") => args.scheme_params.anneal_params.push(parser.value()?.parse()?),
            Long("highs") => args.scheme_params.highs_params.push(parser.value()?.parse()?),
            Long("gurobi") => args.scheme_params.gurobi_params.push(parser.value()?.parse()?),
            Short('O') | Long("out-bams") => args.assgn_params.out_bams = parser.value()?.parse::<PrettyUsize>()?.get(),

            Short('^') | Long("interleaved") => args.in_files.interleaved = true,
            Long("no-index") => args.in_files.no_index = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Long("rerun") => args.rerun = parser.value()?.parse()?,
            Short('s') | Long("seed") => args.seed = Some(parser.value()?.parse()?),
            Long("debug") => args.debug = DebugLvl::from(parser.value()?.parse::<u8>()?),
            Long("strobealign") => args.strobealign = parser.value()?.parse()?,
            Long("minimap") | Long("minimap2") => args.minimap = parser.value()?.parse()?,
            Long("samtools") => args.samtools = parser.value()?.parse()?,

            Short('V') | Long("version") => {
                super::print_version();
                std::process::exit(0);
            }
            Short('h') | Long("help") => {
                print_help(false);
                std::process::exit(0);
            }
            Short('H') | Long("full-help") | Long("hidden-help") => {
                print_help(true);
                std::process::exit(0);
            }
            _ => Err(arg.unexpected())?,
        }
    }
    Ok(args)
}

fn locus_name_matches<'a>(path: &'a Path, subset_loci: &HashSet<String>) -> Option<&'a str> {
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
fn load_priors(path: &Path) -> crate::Result<HashMap<String, HashMap<String, f64>>> {
    let mut res: HashMap<String, HashMap<String, f64>> = HashMap::default();
    for line in ext::sys::open(path)?.lines() {
        let line = line.map_err(add_path!(path))?;
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<_> = line.trim_end().split_whitespace().collect();
        if cols.len() < 3 {
            return Err(error!(InvalidData,
                "Cannot parse genotype priors: invalid line {:?} (expected at least 3 columns)", line));
        }
        let locus = cols[0];
        let gt = cols[1];
        let prior: f64 = cols[2].parse()
            .map(Ln::from_log10)
            .map_err(|_| error!(InvalidData, "Cannot parse genotype priors: offending line {:?}", line))?;
        if prior > 0.0 {
            return Err(error!(InvalidData,
                "Cannot parse genotype priors: offending line {:?} (priors must be in log-10 space)", line));
        }

        if let Some(old_prior) = res.entry(locus.to_owned()).or_default().insert(gt.to_owned(), prior) {
            if old_prior != prior {
                return Err(error!(InvalidData,
                    "Cannot parse genotype priors: locus {} genotype {} contains two priors ({} and {})",
                    locus, gt, Ln::from_log10(old_prior), Ln::from_log10(prior)));
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

/// Removes all files with `.gz` extension, as well as `alns` directory, if it exists.
fn clean_dir(dir: &Path, n_warnings: &mut usize) -> crate::Result<()> {
    let gz_files = ext::sys::filenames_with_ext(dir, "gz")?;
    let alns_dir = dir.join(paths::ALNS_DIR);
    let alns_exist = alns_dir.exists();
    if !gz_files.is_empty() || alns_exist {
        match (*n_warnings).cmp(&3) {
            Ordering::Less  => log::warn!("    Partially cleaning {}", ext::fmt::path(dir)),
            Ordering::Equal => log::warn!("    More directories cleaned, warnings suppressed"),
            Ordering::Greater => {}
        }
        *n_warnings += 1;
    }
    for filename in gz_files.into_iter() {
        if let Err(e) = fs::remove_file(&filename) {
            log::error!("Cannot remove {}: {}", filename.display(), e);
        }
    }
    if alns_exist {
        if let Err(e) = fs::remove_dir_all(&alns_dir) {
            log::error!("Cannot remove directory {}: {}", alns_dir.display(), e);
        }
    }
    Ok(())
}

/// Loads all loci from the database. If `subset_loci` is not empty, only loads loci that are contained in it.
fn load_loci(
    databases: &[PathBuf],
    out_path: &Path,
    subset_loci: &HashSet<String>,
    rerun: super::Rerun,
) -> crate::Result<Vec<LocusData>>
{
    log::info!("Loading database");
    let out_loci_dir = out_path.join(paths::LOCI_DIR);
    ext::sys::mkdir(&out_loci_dir)?;
    let mut loci = Vec::new();
    let mut total_entries = 0;
    let mut n_warnings = 0;
    let mut loci_names = HashSet::default();

    for db_path in databases {
        let db_path = db_path.join(paths::LOCI_DIR);
        if !db_path.exists() {
            log::error!("Database directory {} does not exist", ext::fmt::path(&db_path));
            continue;
        }

        for entry in fs::read_dir(&db_path).map_err(add_path!(db_path))? {
            let entry = entry.map_err(add_path!(!))?;
            let file_type = entry.file_type().map_err(add_path!(entry.path()))?;
            if !(file_type.is_dir() || file_type.is_symlink()) {
                continue;
            }

            total_entries += 1;
            let path = entry.path();
            if let Some(name) = locus_name_matches(&path, subset_loci) {
                if !path.join(paths::SUCCESS).exists() {
                    log::error!("Skipping directory {} (success file missing)", ext::fmt::path(&path));
                    continue;
                }
                if !loci_names.insert(name.to_owned()) {
                    log::error!("Duplicate locus {} in the database, ignoring second instance", name);
                    continue;
                }
                match ContigSet::load(name, &path.join(paths::LOCUS_FASTA), &path.join(paths::KMERS)) {
                    Ok(set) => {
                        let locus_data = LocusData::new(set, &path, &out_loci_dir);
                        if rerun.prepare_and_clean_dir(&locus_data.out_dir, |path| clean_dir(path, &mut n_warnings))? {
                            loci.push(locus_data);
                        }
                    },
                    Err(e) => log::error!("Could not load locus information from {}: {}",
                        ext::fmt::path(&path), e.display()),
                }
            }
        }
    }
    let n = loci.len();
    if n == 0 {
        log::warn!("Zero loci to analyze");
    } else if n < total_entries {
        log::info!("Loaded {} loci, skipped {} directories", n, total_entries - n);
    } else {
        log::info!("Loaded {} loci", n);
    }
    Ok(loci)
}

/// Add padding to regions, add full short contigs and merge everything within `merge_distance`.
pub(super) fn postprocess_targets(
    mut intervals: Vec<Interval>,
    contigs: &Arc<ContigNames>,
    padding: u32,
    alt_contig_len: u32,
    merge_distance: u32,
) -> Vec<Interval>
{
    for interv in intervals.iter_mut() {
        *interv = interv.add_padding(padding);
    }
    for (id, &len) in contigs.ids().zip(contigs.lengths()) {
        if len < alt_contig_len {
            intervals.push(Interval::full_contig(Arc::clone(&contigs), id));
        }
    }
    intervals.sort();
    let merged = Interval::merge(&intervals, merge_distance);
    log::debug!("    Fetch reads from {} regions (sum length {:.1} Mb) + unmapped reads", merged.len(),
        1e-6 * f64::from(merged.iter().map(Interval::len).sum::<u32>()));
    merged
}

/// Prepare fetch regions for a list of regions.
/// Additionally, add all contigs from
fn create_fetch_targets(
    contigs: &Arc<ContigNames>,
    bg_distr: &BgDistr,
    filt_loci: &[&LocusData],
) -> crate::Result<Vec<Interval>>
{
    let mut intervals = Vec::new();
    for locus in filt_loci {
        let bed_filename = locus.db_locus_dir.join(paths::LOCUS_BED);
        let bed_str = fs::read_to_string(&bed_filename).map_err(add_path!(bed_filename))?;
        match Interval::parse_bed(&mut bed_str.split('\t'), &contigs) {
            Ok(interv) => intervals.push(interv),
            Err(Error::ParsingError(e)) => log::error!(
                "[{}] Cannot parse locus coordinates: {}, the region may be absent from the BAM/CRAM file",
                locus.set.tag(), e),
            Err(e) => log::error!("[{}] Cannot parse locus coordinates: {}", locus.set.tag(), e.display()),
        }
    }

    let padding = if let Some(distr) = bg_distr.insert_distr().distr() {
        distr.quantile(0.999).ceil().min(10_000.0) as u32 + 100
    } else {
        1000
    };
    // Recruit reads from all contigs under 10 Mb.
    const ALT_CONTIG_LEN: u32 = 10_000_000;
    const MERGE_DISTANCE: u32 = 1000;
    Ok(postprocess_targets(intervals, contigs, padding, ALT_CONTIG_LEN, MERGE_DISTANCE))
}

/// Check first read in the BAM/CRAM file, and check if it is paired.
fn identify_pairedness(bam_reader: &mut bam::IndexedReader) -> crate::Result<bool> {
    bam_reader.fetch(bam::FetchDefinition::All)?;
    let mut record = bam::Record::new();
    if bam_reader.read(&mut record).transpose()?.is_none() {
        Err(error!(InvalidData, "Input BAM/CRAM file is empty"))
    } else {
        log::debug!("Input BAM/CRAM file is identified as {}-end",
            if record.is_paired() { "paired" } else { "single"});
        Ok(record.is_paired())
    }
}

/// Recruit reads to targets, depending on the reader type.
///
/// is_paired_end - needed for BAM/CRAM files. If unknown, `identify_pairedness` is called.
pub(super) fn recruit_to_targets(
    targets: &recruit::Targets,
    in_files: &InputFiles,
    writers: Vec<impl io::Write>,
    is_paired_end: Option<bool>,
    threads: u16,
    chunk_size: usize,
    subsampling: Option<(f64, XoshiroRng)>,
    get_targets: impl FnOnce(&Arc<ContigNames>) -> crate::Result<Vec<Interval>>,
) -> crate::Result<recruit::Progress>
{
    if in_files.has_indexed_alignment() {
        assert!(in_files.alns.len() == 1);
        let bam_filename = in_files.alns[0].to_path_buf();
        let mut bam_reader = bam::IndexedReader::from_path(&bam_filename)?;
        fastx::set_reference(&bam_filename, &mut bam_reader, &in_files.reference, None)?;
        let contigs = Arc::new(ContigNames::from_bam_header("bam", bam_reader.header())?);
        let fetch_regions = get_targets(&contigs)?;
        let is_paired_end = is_paired_end.map(Ok).unwrap_or_else(|| identify_pairedness(&mut bam_reader))?;

        let reader = fastx::IndexedBamReader::new(bam_filename, bam_reader, fetch_regions)?;
        // Need a lot of if-elses to have compile-time optimizations based on input data and the number of targets :]
        if is_paired_end {
            targets.recruit(fastx::PairedBamReader::new(reader), writers, threads, chunk_size, subsampling)
        } else {
            targets.recruit(reader, writers, threads, chunk_size, subsampling)
        }
    } else {
        fastx::process_readers!(in_files, None; let {} reader;
            { targets.recruit(reader, writers, threads, chunk_size, subsampling) })
    }
}

/// Calculate chunk size based on the chunk length (sum bp, processed in one chunk) and mean read length.
pub fn calculate_chunk_size(chunk_length: u64, mean_read_len: f64, is_paired_end: bool) -> usize {
    max(1, (chunk_length as f64 / mean_read_len * (if is_paired_end { 0.5 } else { 1.0 })) as usize)
}

/// Recruits reads to all loci, where neither reads nor alignments are available.
fn recruit_reads(
    loci: &[LocusData],
    bg_distr: &BgDistr,
    args: &Args,
    recr_params: recruit::Params,
) -> crate::Result<()>
{
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
    let mut target_builder = recruit::TargetBuilder::new(recr_params);
    let mut writers = Vec::with_capacity(n_filt_loci);
    let mean_read_len = bg_distr.seq_info().mean_read_len();

    for locus in filt_loci.iter() {
        target_builder.add(&locus.set, mean_read_len);
        // Output files with a large buffer (4 Mb).
        const BUFFER: usize = 4_194_304;
        writers.push(fs::File::create(&locus.tmp_reads_filename).map_err(add_path!(&locus.tmp_reads_filename))
            .map(|w| io::BufWriter::with_capacity(BUFFER, w))?);
    }
    let targets = target_builder.finalize();
    let is_paired_end = bg_distr.insert_distr().is_paired_end();
    let chunk_size = calculate_chunk_size(args.chunk_length, mean_read_len, is_paired_end);
    recruit_to_targets(&targets, &args.in_files, writers, Some(is_paired_end), args.threads, chunk_size, None,
        |contigs| create_fetch_targets(contigs, bg_distr, &filt_loci))?;

    for locus in filt_loci.iter() {
        fs::rename(&locus.tmp_reads_filename, &locus.reads_filename)
            .map_err(add_path!(locus.tmp_reads_filename, &locus.reads_filename))?;
    }
    Ok(())
}

/// Map reads with either Strobealign or Minimap2.
fn create_mapping_command(
    ref_path: &Path,
    reads_path: &Path,
    seq_info: &SequencingInfo,
    n_seqs: usize,
    args: &Args,
) -> Command
{
    let n_locs = min(25000, n_seqs * 4).to_string();

    let mut cmd: Command;
    if seq_info.technology() == Technology::Illumina {
        // Force strobealign v0.13 to use more minimizers.
        let rescue_param = (n_seqs * 4 / 100).clamp(3, 25000);

        cmd = Command::new(&args.strobealign);
        // NOTE: Provide single input FASTQ, and do not use --interleaved option.
        // This is because Strobealign starts to connect read pairs in multiple ways,
        // which produces too large output, and takes a lot of time.
        cmd.args(&[
            "--no-progress",
            "-M", &n_locs,   // Try as many secondary locations as possible.
            "-N", &n_locs,   // Output as many secondary alignments as possible.
            "-R", &rescue_param.to_string(),
            "-S", "0.5",     // Try candidate sites with minimizer hits >= FLOAT * best hits.
            "-f", "0.001",
            "-k", "15",      // Use smaller minimizers to get more matches.
            "--eqx",         // Output X/= instead of M operations.
            "-r", &format!("{:.0}", seq_info.mean_read_len()),
            "-t", &args.threads.to_string()]);
    } else {
        cmd = Command::new(&args.minimap);
        cmd.args(&[
            "-a", // Output SAM format,
            "-x", seq_info.technology().minimap_preset(), // Set mapping preset.
            "-N", &n_locs,   // Output as many secondary alignments as possible.
            "-f", "0.05",   // Filter out 5% minimizers. Should be still enough minimizers for long reads.
            "--eqx",        // Output X/= instead of M operations.
            "-t", &args.threads.to_string()]);
    }
    cmd.arg(&ref_path).arg(&reads_path)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    cmd
}

fn map_reads(locus: &LocusData, bg_distr: &BgDistr, args: &Args) -> crate::Result<()> {
    if locus.aln_filename.exists() {
        log::info!("    Skipping read mapping");
        return Ok(());
    }

    let in_fasta = locus.db_locus_dir.join(paths::LOCUS_FASTA);
    log::info!("    Mapping reads to {}", ext::fmt::path(&in_fasta));
    let start = Instant::now();
    let mut mapping_cmd = create_mapping_command(&in_fasta, &locus.reads_filename, bg_distr.seq_info(),
        locus.set.len(), args);
    let mapping_exe = PathBuf::from(mapping_cmd.get_program().to_owned());
    let mut mapping_child = mapping_cmd.spawn().map_err(add_path!(mapping_exe))?;
    let child_stdout = mapping_child.stdout.take().unwrap();
    let mut pipe_guard = ext::sys::PipeGuard::new(mapping_exe, mapping_child);

    let mut samtools_cmd = Command::new(&args.samtools);
    samtools_cmd.args(&["view", "-b"]); // Output BAM.
    // if !bg_distr.insert_distr().is_paired_end() {
    //     samtools_cmd.arg("-F4"); // ignore unmapped reads in case of single-end reads.
    // }
    samtools_cmd.arg("-o").arg(&locus.tmp_aln_filename);
    samtools_cmd.stdin(Stdio::from(child_stdout)).stdout(Stdio::piped()).stderr(Stdio::piped());

    log::debug!("    {} | {}", ext::fmt::command(&mapping_cmd), ext::fmt::command(&samtools_cmd));
    let samtools_child = samtools_cmd.spawn().map_err(add_path!(args.samtools))?;
    pipe_guard.push(args.samtools.clone(), samtools_child);
    pipe_guard.wait()?;

    log::debug!("    Finished in {}", ext::fmt::Duration(start.elapsed()));
    fs::rename(&locus.tmp_aln_filename, &locus.aln_filename)
        .map_err(add_path!(locus.tmp_aln_filename, locus.aln_filename))?;
    Ok(())
}

/// Generates the list of genotypes and their priors.
fn generate_genotypes(
    contig_ids: &[ContigId],
    contigs: &ContigNames,
    opt_priors: Option<&HashMap<String, f64>>,
    ploidy: usize,
) -> crate::Result<(Vec<Genotype>, Vec<f64>)>
{
    if let Some(priors_map) = opt_priors {
        assert_eq!(contig_ids.len(), contigs.len(),
            "All contig IDs must be present when priors are provided");
        let mut genotypes = Vec::with_capacity(priors_map.len());
        let mut priors = Vec::with_capacity(priors_map.len());
        for (s, &prior) in priors_map.iter() {
            let gt = Genotype::parse(s, contigs)?;
            if prior > 0.0 || prior.is_nan() {
                return Err(error!(InvalidInput, "Invalid prior {} for genotype {}", prior, s));
            } else if gt.ploidy() != ploidy {
                return Err(error!(InvalidInput,
                    "Cannot load prior for genotype {} (expected ploidy {})", s, ploidy));
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
    edit_dist_cache: &EditDistCache,
    scheme: &Arc<Scheme>,
    distr_cache: &Arc<DistrCache>,
    opt_priors: Option<&HashMap<String, f64>>,
    mut rng: ext::rand::XoshiroRng,
    args: &Args,
) -> crate::Result<()>
{
    log::info!("{} {}", "Analyzing".bold(), locus.set.tag().bold());
    let timer = Instant::now();
    map_reads(locus, bg_distr, &args)?;
    let is_paired_end = bg_distr.insert_distr().is_paired_end();

    log::info!("    Calculating read alignment probabilities");
    let bam_reader = bam::Reader::from_path(&locus.aln_filename)?;
    let contigs = locus.set.contigs();

    let all_contig_infos = if args.debug >= DebugLvl::Some {
        let windows_filename = locus.out_dir.join("windows.bed.gz");
        let windows_writer = ext::sys::create_gzip(&windows_filename)?;
        ContigInfo::new_all(&locus.set, bg_distr.depth(), &args.assgn_params, windows_writer)?
    } else {
        ContigInfo::new_all(&locus.set, bg_distr.depth(), &args.assgn_params, io::sink())?
    };

    let all_alns = if args.debug >= DebugLvl::Full {
        let reads_writer = ext::sys::create_gzip(&locus.out_dir.join("reads.csv.gz"))?;
        let read_kmer_writer = ext::sys::create_gzip(&locus.out_dir.join("read_kmers.csv.gz"))?;
        AllAlignments::load(bam_reader, &locus.set, bg_distr, edit_dist_cache,
            &args.assgn_params, reads_writer, read_kmer_writer)?
    } else {
        AllAlignments::load(bam_reader, &locus.set, bg_distr, edit_dist_cache,
            &args.assgn_params, io::sink(), io::sink())?
    };

    if locus.reads_filename.exists() {
        // Alignments are succesfully loaded, now we can remove file with reads.
        fs::remove_file(&locus.reads_filename).map_err(add_path!(locus.reads_filename))?;
    }

    let genotyping = if all_alns.reads().is_empty() {
        log::error!("[{}] No available reads", locus.set.tag());
        scheme::Genotyping::empty_result(locus.set.tag().to_string(), vec![scheme::GenotypingWarning::NoReads])
    } else {
        if is_paired_end && args.debug >= DebugLvl::Full {
            let read_pairs_filename = locus.out_dir.join("read_pairs.csv.gz");
            let pairs_writer = ext::sys::create_gzip(&read_pairs_filename)?;
            all_alns.write_read_pair_info::<false>(pairs_writer, contigs).map_err(add_path!(read_pairs_filename))?;
        }

        let contig_ids: Vec<ContigId> = contigs.ids().collect();
        let (genotypes, priors) = generate_genotypes(&contig_ids, contigs, opt_priors, usize::from(args.ploidy))?;
        if genotypes.is_empty() {
            return Err(error!(RuntimeError, "No available genotypes for locus {}", locus.set.tag()));
        }

        let dist_filename = locus.db_locus_dir.join(paths::DISTANCES);
        let contig_distances = if dist_filename.exists() {
            let dist_file = io::BufReader::new(fs::File::open(&dist_filename)
                .map_err(add_path!(dist_filename))?);
            let (_k, _w, dists) = div::load_divergences(dist_file, &dist_filename, contigs.len())?;
            Some(dists)
        } else {
            None
        };

        let data = scheme::Data {
            scheme: Arc::clone(scheme),
            contigs: Arc::clone(contigs),
            contig_distances,
            distr_cache: Arc::clone(distr_cache),
            assgn_params: args.assgn_params.clone(),
            debug: args.debug,
            threads: usize::from(args.threads),
            all_alns, genotypes, priors, all_contig_infos, is_paired_end,
        };
        scheme::solve(data, bg_distr, &locus.out_dir, &mut rng)?
    };

    let res_filename = locus.out_dir.join(paths::RES_JSON);
    let mut res_writer = ext::sys::create_gzip(&res_filename)?;
    genotyping.to_json().write_pretty(&mut res_writer, 4).map_err(add_path!(res_filename))?;
    super::write_success_file(locus.out_dir.join(paths::SUCCESS))?;
    log::info!("    [{}] Successfully finished in {}", locus.set.tag(), ext::fmt::Duration(timer.elapsed()));
    Ok(())
}

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let mut args = parse_args(argv)?;
    args.in_files.fill_from_inlist()?;
    let mut args = args.validate()?;

    super::greet();
    let timer = Instant::now();
    let mut rng = ext::rand::init_rng(args.seed);
    let out_dir = args.output.as_ref().unwrap();
    ext::sys::mkdir(out_dir)?;
    let preproc_dir = args.preproc.as_ref().unwrap();

    let bg_distr = BgDistr::load_from(&preproc_dir.join(paths::BG_DISTR), Some(&preproc_dir.join(paths::SUCCESS)))?;
    args.assgn_params.set_tweak_size(bg_distr.depth().window_size())?;
    let tech = bg_distr.seq_info().technology();
    let recr_params = recruit::Params::new(
        args.minimizer_kw,
        args.match_frac.unwrap_or_else(|| tech.default_match_frac(bg_distr.insert_distr().is_paired_end())),
        args.match_len, args.thresh_kmer_count)?;

    // Add 1 to good edit distance.
    const GOOD_DISTANCE_ADD: u32 = 1;
    let edit_dist_cache = EditDistCache::new(bg_distr.error_profile(), GOOD_DISTANCE_ADD, args.assgn_params.edit_pvals);
    edit_dist_cache.describe_with(bg_distr.seq_info().mean_read_len());

    validate_param!(args.in_files.has_indexed_alignment()
        || bg_distr.insert_distr().is_paired_end() == args.in_files.is_paired_end(),
        "Paired-end/Single-end status does not match background data");
    if tech == Technology::Illumina {
        args.strobealign = ext::sys::find_exe(args.strobealign)?;
    } else {
        args.minimap = ext::sys::find_exe(args.minimap)?;
    }

    let priors = args.priors.as_ref().map(|path| load_priors(path)).transpose()?;
    let loci = load_loci(&args.databases, out_dir, &args.subset_loci, args.rerun)?;
    if loci.is_empty() {
        return Ok(());
    }
    recruit_reads(&loci, &bg_distr, &args, recr_params)?;

    let scheme = Arc::new(Scheme::create(&args.scheme_params)?);
    let distr_cache = Arc::new(DistrCache::new(&bg_distr, args.assgn_params.alt_cn));
    let mut successes = 0;
    for locus in loci.iter() {
        // Remove to get ownership of the locus priors.
        let locus_priors = priors.as_ref().and_then(|priors| priors.get(locus.set.tag()));
        if priors.is_some() && locus_priors.is_none() {
            log::warn!("Priors for locus {} are not found. Assuming equal priors", locus.set.tag());
        }

        let rng_clone = rng.clone();
        // Jump over 2^192 random numbers. This way, all loci have independent random numbers.
        rng.long_jump();
        match analyze_locus(locus, &bg_distr, &edit_dist_cache, &scheme, &distr_cache, locus_priors, rng_clone, &args) {
            Err(e) => log::error!("Error in locus {}: {}", locus.set.tag(), e.display()),
            Ok(()) => successes += 1,
        }
    }
    let nloci = loci.len();
    if successes == 0 {
        log::error!("Failed at {} loci", nloci);
    } else if successes < nloci {
        log::warn!("Successfully analysed {} loci, failed at {} loci", successes, nloci - successes);
    } else {
        log::info!("Successfully analysed {} loci", successes);
    }
    log::info!("Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
