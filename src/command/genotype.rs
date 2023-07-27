use std::{
    fs, fmt,
    io::{self, Write, BufRead},
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
use lexopt::ValueExt;
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
        locs::{self, AllAlignments},
        windows::ContigWindows,
        dp_cache::CachedDepthDistrs,
    },
    solvers::scheme::{self, Scheme},
};
use super::paths;


#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct UsizeOrInf(pub usize);

impl UsizeOrInf {
    const INF: Self = Self(usize::MAX);
}

impl std::str::FromStr for UsizeOrInf {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "inf" || s == "Inf" || s == "INF" {
            Ok(Self::INF)
        } else {
            usize::from_str(s)
                .map_err(|_| format!("Cannot parse {:?}, possible values: integer or 'inf'", s))
                .map(Self)
        }
    }
}

impl fmt::Display for UsizeOrInf {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.0 == usize::MAX {
            f.write_str("inf")
        } else {
            write!(f, "{}", self.0)
        }
    }
}

pub struct FilterSubparams {
    /// Boundaries on the smallest and largest number of haplotypes/genotypes.
    pub min_size: UsizeOrInf,
    pub max_size: UsizeOrInf,
    /// Relative score threshold (between the minimum and maximum score).
    pub score_thresh: f64,
}

impl FilterSubparams {
    pub fn validate(&self, prefix: &'static str) -> Result<(), Error> {
        validate_param!(0.0 <= self.score_thresh && self.score_thresh <= 1.0,
            "{}: score threshold ({}) must be within [0, 1]", prefix, self.score_thresh);
        validate_param!(self.min_size.0 > 1,
            "{}: minimal size ({}) must be at least 2", prefix, self.min_size);
        validate_param!(self.min_size <= self.max_size,
            "{}: min size ({}) must be <= the max size ({})",
            prefix, self.min_size, self.max_size);
        Ok(())
    }

    pub fn set_sizes(&mut self, mut values: lexopt::ValuesIter) -> Result<(), lexopt::Error> {
        self.min_size = values.next().expect("First value must be present").parse()?;
        self.max_size = values.next()
            .map(|v| v.parse())
            .transpose()?
            .unwrap_or(UsizeOrInf::INF);
        Ok(())
    }

    pub fn disable(&mut self) {
        self.min_size = UsizeOrInf::INF;
        self.max_size = UsizeOrInf::INF;
    }

    /// Finds partition point and score threshold based on the decreasing list of scores.
    pub fn get_partition<T>(&self, scores: &[(T, f64)]) -> (usize, f64) {
        let n = scores.len();
        let best = scores[0].1;
        let worst = scores[n - 1].1;
        let range = best - worst;

        let mut thresh = worst + range * self.score_thresh;
        let mut m = scores.partition_point(|(_, score)| *score >= thresh);
        if self.min_size.0 >= n {
            thresh = worst;
            m = n;
        } else if m < self.min_size.0 {
            thresh = scores[self.min_size.0 - 1].1;
            m = scores.partition_point(|(_, score)| *score >= thresh);
        } else if m > self.max_size.0 {
            thresh = scores[self.max_size.0 - 1].1;
            m = scores.partition_point(|(_, score)| *score >= thresh);
        }
        // m can still be over max_size, but we do not want to discard entries with equal scores.
        (m, thresh)
    }
}

/// Filter haplotypes and genotypes before assigning reads.
pub struct FilterParams {
    /// Rank haplotypes based on the top X read alignments, where X = haplotype_aln_coef / ploidy`.
    pub nreads_mult: f64,
    pub haplotype: FilterSubparams,
    pub genotype: FilterSubparams,
}

impl Default for FilterParams {
    fn default() -> Self {
        Self {
            nreads_mult: 0.8,
            haplotype: FilterSubparams {
                min_size: UsizeOrInf(5),
                max_size: UsizeOrInf(100),
                score_thresh: 0.9,
            },
            genotype: FilterSubparams {
                min_size: UsizeOrInf(10),
                max_size: UsizeOrInf(1000),
                score_thresh: 0.9,
            },
        }
    }
}

impl FilterParams {
    pub fn validate(&self) -> Result<(), Error> {
        validate_param!(0.0 < self.nreads_mult && self.nreads_mult < 2.0,
            "Number of reads multiplier ({}) must be reasonably close to 1",
            self.nreads_mult);
        self.haplotype.validate("Haplotype filtering")?;
        self.genotype.validate("Genotype filtering")?;
        Ok(())
    }
}

struct Args {
    input: Vec<PathBuf>,
    database: Option<PathBuf>,
    output: Option<PathBuf>,
    subset_loci: FnvHashSet<String>,
    ploidy: u8,
    solvers: Option<PathBuf>,
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
    filt_params: FilterParams,
    assgn_params: AssgnParams,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::new(),
            database: None,
            output: None,
            subset_loci: FnvHashSet::default(),
            ploidy: 2,
            solvers: None,
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
            filt_params: Default::default(),
            assgn_params: Default::default(),
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
        self.filt_params.validate()?;
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
    println!("    {:KEY$} {:VAL$}  Ignore read alignments that are 10^{} times worse than\n\
        {EMPTY}  the best alignment [{}].",
        "-D, --prob-diff".green(), "FLOAT".yellow(), "FLOAT".yellow(),
        super::fmt_def_f64(Ln::to_log10(defaults.assgn_params.prob_diff)));
    println!("    {:KEY$} {:VAL$}  Unmapped read mate receives 10^{} penalty [{}].",
        "-U, --unmapped".green(), "FLOAT".yellow(), "FLOAT".yellow(),
        super::fmt_def_f64(Ln::to_log10(defaults.assgn_params.unmapped_penalty)));
    println!("    {:KEY$} {} \n\
        {EMPTY}  Contig windows receive different weight depending on average k-mer\n\
        {EMPTY}  frequency. Windows with values under {} [{}] receive full weight.\n\
        {EMPTY}  Windows with values equal to {} [{}] receive half weight.",
        "    --rare-kmer".green(), "FLOAT FLOAT".yellow(),
        "FLOAT_1".yellow(), super::fmt_def(defaults.assgn_params.rare_kmer),
        "FLOAT_2".yellow(), super::fmt_def(defaults.assgn_params.semicommon_kmer));
    println!("    {:KEY$} {:VAL$}  Read depth likelihood contribution relative to\n\
        {EMPTY}  read alignment likelihoods [{}].",
        "-C, --dp-contrib".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.assgn_params.depth_contrib));
    println!("    {} {}  Compare window probability to have copy number 1 against two\n\
        {EMPTY}  alternative CN values [{} {}]. First in (0, 1), second > 1.",
        "-A, --alt-cn".green(), "FLOAT FLOAT".yellow(),
        super::fmt_def_f64(defaults.assgn_params.alt_cn.0), super::fmt_def_f64(defaults.assgn_params.alt_cn.1));

    let filt_defaults = &defaults.filt_params;
    println!("\n{}", "Pre-filtering:".bold());
    println!("    {:KEY$} {:VAL$}  Skip pre-filtering (same as {} or {}).\n\
        {EMPTY}  Note: haplotype filtering is always skipped if priors are set.",
        "    --no-filtering".green(), super::flag(), "--n-haps inf".underline(), "--n-gts inf".underline());
    println!("    {:KEY$} {:VAL$}  Rank haplotypes based on the top {}/ploidy read alignments [{}].",
        "    --nreads-mult".green(), "FLOAT".yellow(), "FLOAT".yellow(),
        super::fmt_def_f64(filt_defaults.nreads_mult));
    println!("    {:KEY$} {:VAL$}\n\
        {EMPTY}  Score threshold for haplotype [{}] and genotype [{}] filtering.\n\
        {EMPTY}  Values range from 0 (use all) to 1 (only entries with the best score).",
        "    --score-thresh".green(), "FLOAT [FLOAT]".yellow(),
        super::fmt_def_f64(filt_defaults.haplotype.score_thresh),
        super::fmt_def_f64(filt_defaults.genotype.score_thresh));
    println!("    {}    {}  Minimum and maximum number of haplotypes after filtering [{} {}].",
        "    -n-haps".green(), "INT [INT]".yellow(),
        super::fmt_def(filt_defaults.haplotype.min_size), super::fmt_def(filt_defaults.haplotype.max_size));
    println!("    {}     {}  Minimum and maximum number of genotypes after filtering [{} {}].",
        "    -n-gts".green(), "INT [INT]".yellow(),
        super::fmt_def(filt_defaults.genotype.min_size), super::fmt_def(filt_defaults.genotype.max_size));

    println!("\n{}", "Locus genotyping:".bold());
    println!("    {:KEY$} {:VAL$}  Optional: describe sequence of solvers in a JSON file.\n\
        {EMPTY}  Please see README for information on the file content.",
        "-S, --solvers".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Randomly move read coordinates by at most {} bp [{}].",
        "    --tweak".green(), "INT".yellow(), "INT".yellow(), "auto".cyan());

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

            Short('S') | Long("solvers") => args.solvers = Some(parser.value()?.parse()?),
            Short('C') | Long("dp-contrib") | Long("depth-contrib") | Long("dp-contribution") =>
                args.assgn_params.depth_contrib = parser.value()?.parse()?,
            Short('A') | Long("alt-cn") =>
                args.assgn_params.alt_cn = (parser.value()?.parse()?, parser.value()?.parse()?),
            Short('D') | Long("prob-diff") => args.assgn_params.prob_diff = Ln::from_log10(parser.value()?.parse()?),
            Short('U') | Long("unmapped") =>
                args.assgn_params.unmapped_penalty = Ln::from_log10(parser.value()?.parse()?),
            Long("rare-kmer") => {
                args.assgn_params.rare_kmer = parser.value()?.parse()?;
                args.assgn_params.semicommon_kmer = parser.value()?.parse()?;
            }
            Long("tweak") => {
                let val = parser.value()?;
                args.assgn_params.tweak = if val == "auto" {
                    None
                } else {
                    Some(val.parse()?)
                };
            }

            Long("no-filt") | Long("no-filter") | Long("no-filtering") => {
                args.filt_params.haplotype.disable();
                args.filt_params.genotype.disable();
            }
            Long("nreads-mult") | Long("nreads-multiplier") =>
                args.filt_params.nreads_mult = parser.value()?.parse()?,
            Long("score-thresh") | Long("score-threshold") => {
                let mut values = parser.values()?;
                args.filt_params.haplotype.score_thresh = values.next()
                    .expect("At least one value must be present").parse()?;
                args.filt_params.genotype.score_thresh = values.next()
                    .map(|v| v.parse())
                    .transpose()?
                    .unwrap_or(args.filt_params.haplotype.score_thresh);
            }
            Long("n-haps") | Long("n-haplotypes") =>
                args.filt_params.haplotype.set_sizes(parser.values()?)?,
            Long("n-gts") | Long("n-genotypes") =>
                args.filt_params.genotype.set_sizes(parser.values()?)?,

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
    /// File with haplotype-pair likelihoods.
    lik_filename: PathBuf,
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
            lik_filename: out_dir.join("lik.csv.gz"),
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
) -> Result<Vec<(Genotype, f64)>, Error>
{
    if let Some(priors_map) = opt_priors {
        assert_eq!(contig_ids.len(), contigs.len(),
            "All contig IDs must be present when priors are provided");
        let mut gt_priors = Vec::with_capacity(priors_map.len());
        for (s, &prior) in priors_map.iter() {
            let gt = Genotype::parse(s, contigs)?;
            if prior > 0.0 || prior.is_nan() {
                return Err(Error::InvalidInput(format!("Invalid prior {} for genotype {}", prior, s)));
            } else if gt.ploidy() != ploidy {
                return Err(Error::InvalidInput(format!(
                    "Cannot load prior for genotype {} (expected ploidy {})", s, ploidy)));
            } else if prior.is_finite() {
                gt_priors.push((gt, prior));
            }
        }
        Ok(gt_priors)
    } else {
        let mut gt_priors = Vec::with_capacity(ext::vec::count_combinations_with_repl(contig_ids.len(), ploidy));
        ext::vec::gen_combinations_with_repl(&contig_ids, ploidy,
            |ids| gt_priors.push((Genotype::new(ids, contigs), 0.0)));
        Ok(gt_priors)
    }
}

// ln(10^-50).
const MIN_SCORE_RANGE: f64 = -115.12925464970229;

// fn select_haplotypes(
//     contigs: &ContigNames,
//     all_alns: &AllAlignments,
//     mut lik_writer: impl Write,
//     params: &SelectParams,
//     ploidy: u8,
//     debug: bool,
// ) -> io::Result<Vec<ContigId>>
// {
//     let n = contigs.len();
//     if params.haplotype.min_size >= n && !debug {
//         return Ok(contigs.ids().collect());
//     }
//     let n_reads = all_alns.n_reads();
//     assert!(n_reads != 0);
//     let n_best_reads = (((n_reads as f64) * params.nreads_mult / f64::from(ploidy)).ceil() as usize)
//         .clamp(1, n_reads);
//     log::info!("    Filtering haplotypes based on the best {} reads out of {}", n_best_reads, n_reads);

//     // Vector of vectors: contig id -> read pair -> best aln prob.
//     let mut best_aln_probs = vec![Vec::<f64>::with_capacity(n_reads); n];
//     for read_alns in all_alns.iter() {
//         for (contig_probs, best_prob) in best_aln_probs.iter_mut().zip(read_alns.best_for_each_contig(contigs)) {
//             contig_probs.push(best_prob);
//         }
//     }

//     let mut scores = Vec::with_capacity(n);
//     for (id, contig_probs) in contigs.ids().zip(best_aln_probs.iter_mut()) {
//         contig_probs.sort_unstable_by(|a, b| b.total_cmp(&a));
//         let score: f64 = contig_probs[..n_best_reads].iter().sum();
//         writeln!(lik_writer, "0a\t{}\t{:.3}", contigs.get_name(id), score)?;
//         scores.push((id, score));
//     }

//     // Decreasing sort by score.
//     scores.sort_unstable_by(|a, b| b.1.total_cmp(&a.1));
//     let best = scores[0];
//     let worst = scores[n - 1];
//     let range = best.1 - worst.1;
//     if range < MIN_SCORE_RANGE {
//         log::warn!("        Difference between haplotypes is too small ({:.1}), keeping all", Ln::to_log10(range));
//         return Ok(contigs.ids().collect());
//     }
//     let (m, thresh) = params.haplotype.get_partition(&scores);
//     log::debug!("        Selected {} haplotypes out of {}, relative score thresh {:.3}", m, n,
//         (thresh - worst.1) / range);
//     log::debug!("        Scores: worst {} ({:.1}), best {} ({:.1}), threshold: {:.1}",
//         contigs.get_name(worst.0), Ln::to_log10(worst.1),
//         contigs.get_name(best.0), Ln::to_log10(best.1), Ln::to_log10(thresh));
//     Ok(scores[..m].iter().map(|&(id, _)| id).collect())
// }

/// Filter genotypes based on read alignments alone (without accounting for read depth).
fn filter_genotypes(
    contig_ids: &[ContigId],
    gt_priors: Vec<(Genotype, f64)>,
    all_alns: &AllAlignments,
    mut lik_writer: impl Write,
    params: &FilterParams,
    debug: bool,
) -> io::Result<Vec<(Genotype, f64)>>
{
    let n = gt_priors.len();
    // Need to run anyway in case of `debug`, so that we output all values.
    if params.genotype.min_size.0 >= n && !debug {
        return Ok(gt_priors);
    }
    log::info!("    Filtering genotypes based on read alignment likelihood");

    let best_aln_matrix = all_alns.best_aln_matrix(contig_ids);
    // Vector (genotype index, score).
    let mut scores = Vec::with_capacity(n);
    let mut gt_best_probs = vec![0.0_f64; all_alns.n_reads()];
    for (i, (gt, prior)) in gt_priors.iter().enumerate() {
        let gt_ids = gt.ids();
        gt_best_probs.copy_from_slice(&best_aln_matrix[gt_ids[0].ix()]);
        for id in gt_ids[1..].iter() {
            gt_best_probs.iter_mut().zip(&best_aln_matrix[id.ix()])
                .for_each(|(best_prob, &curr_prob)| *best_prob = best_prob.max(curr_prob));
        }
        let score = *prior + gt_best_probs.iter().sum::<f64>();
        writeln!(lik_writer, "0\t{}\t{:.3}", gt, Ln::to_log10(score))?;
        scores.push((i, score));
    }

    // Decreasing sort by score.
    scores.sort_unstable_by(|a, b| b.1.total_cmp(&a.1));
    let best = scores[0];
    let worst = scores[n - 1];
    let range = best.1 - worst.1;
    if range < MIN_SCORE_RANGE {
        log::warn!("        Difference between genotypes is too small ({:.1}), keeping all", Ln::to_log10(range));
        return Ok(gt_priors);
    }
    let (m, thresh) = params.genotype.get_partition(&scores);
    log::debug!("        Selected {} genotypes out of {}, relative score thresh {:.3}", m, n,
        (thresh - worst.1) / range);
    log::debug!("        Worst {} ({:.1}), best {} ({:.1}), threshold: {:.1}",
        gt_priors[worst.0].0, Ln::to_log10(worst.1), gt_priors[best.0].0, Ln::to_log10(best.1), Ln::to_log10(thresh));
    Ok(scores[..m].iter().map(|&(i, _)| gt_priors[i].clone()).collect())
}

fn analyze_locus(
    locus: &LocusData,
    bg_distr: &BgDistr,
    scheme: &Scheme,
    cached_distrs: &CachedDepthDistrs,
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

    let contig_windows = ContigWindows::new_all(&locus.set, cached_distrs, &args.assgn_params);
    if args.debug || scheme.has_dbg_output() {
        let windows_filename = locus.out_dir.join("windows.bed.gz");
        let mut windows_writer = ext::sys::create_gzip(&windows_filename)?;
        writeln!(windows_writer, "#{}", ContigWindows::BED_HEADER).map_err(add_path!(windows_filename))?;
        for curr_windows in contig_windows.iter() {
            curr_windows.write_to(&mut windows_writer, &contigs).map_err(add_path!(windows_filename))?;
        }
    }

    let all_alns = if args.debug {
        let reads_filename = locus.out_dir.join("reads.csv.gz");
        let mut reads_writer = ext::sys::create_gzip(&reads_filename)?;
        writeln!(reads_writer, "{}", locs::CSV_HEADER).map_err(add_path!(reads_filename))?;
        AllAlignments::load(bam_reader, contigs, bg_distr, &args.assgn_params, reads_writer)?
    } else {
        AllAlignments::load(bam_reader, contigs, bg_distr, &args.assgn_params, io::sink())?
    };

    let mut lik_writer = ext::sys::create_gzip(&locus.lik_filename)?;
    writeln!(lik_writer, "stage\tgenotype\tlik").map_err(add_path!(locus.lik_filename))?;

    let contig_ids: Vec<ContigId> = contigs.ids().collect();
    let gt_priors = generate_genotypes(&contig_ids, contigs, opt_priors, usize::from(args.ploidy))?;
    if gt_priors.is_empty() {
        return Err(Error::RuntimeError(format!("No available genotypes for locus {}", locus.set.tag())));
    }

    let gt_priors = filter_genotypes(&contig_ids, gt_priors, &all_alns, &mut lik_writer,
        &args.filt_params, args.debug).map_err(add_path!(locus.lik_filename))?;
    let data = scheme::Data {
        scheme: scheme.clone(),
        contigs: Arc::clone(&contigs),
        contig_windows: Arc::new(contig_windows),
        params: args.assgn_params.clone(),
        debug: args.debug,
        all_alns, gt_priors,
    };
    scheme::solve(data, lik_writer, &locus.out_dir, &mut rng, args.threads)?;
    super::write_success_file(locus.out_dir.join(paths::SUCCESS))?;
    log::info!("    [{}] Successfully finished in {}", locus.set.tag(), ext::fmt::Duration(timer.elapsed()));
    Ok(())
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let timer = Instant::now();
    let mut args = parse_args(argv)?.validate()?;
    let db_dir = args.database.as_ref().unwrap();
    let out_dir = args.output.as_ref().unwrap();
    let priors = args.priors.as_ref().map(|path| load_priors(path)).transpose()?;

    let bg_path = out_dir.join(paths::BG_DIR).join(paths::BG_DISTR);
    let mut bg_stream = ext::sys::open(&bg_path)?;
    let mut bg_str = String::new();
    bg_stream.read_to_string(&mut bg_str).map_err(add_path!(bg_path))?;
    let bg_distr = BgDistr::load(&json::parse(&bg_str)?)?;
    args.assgn_params.set_tweak_size(bg_distr.depth().window_size())?;
    validate_param!(bg_distr.insert_distr().is_paired_end() == args.is_paired_end(),
        "Paired-end/Single-end status does not match background data");
    if bg_distr.seq_info().technology() == Technology::Illumina {
        args.strobealign = ext::sys::find_exe(args.strobealign)?;
    } else {
        args.minimap = ext::sys::find_exe(args.minimap)?;
    }

    let loci = load_loci(db_dir, out_dir, &args.subset_loci, args.rerun)?;
    recruit_reads(&loci, &args)?;

    let scheme = match args.solvers.as_ref() {
        Some(filename) => Scheme::from_json(&ext::sys::load_json(filename)?)?,
        None => Scheme::default(),
    };
    let cached_distrs = CachedDepthDistrs::new(&bg_distr, args.assgn_params.alt_cn);

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
