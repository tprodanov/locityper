//! Preprocess WGS dataset.

use std::{
    fs, thread,
    io::{self, Write},
    fmt::Write as FmtWrite,
    cmp::{min, max},
    path::{Path, PathBuf},
    process::{Stdio, Command, ChildStdin},
    time::Instant,
    ops::Deref,
    sync::Arc,
    str::FromStr,
};
use colored::Colorize;
use const_format::str_repeat;
use htslib::bam;
use regex::{Regex, RegexBuilder};
use crate::{
    ext::{
        self,
        rand::{init_rng, XoshiroRng},
        fmt::PrettyU32,
    },
    err::{Error, validate_param, add_path},
    seq::{
        Interval,
        kmers::{JfKmerGetter, KmerCounts},
        fastx::{self, FastxRead},
        cigar::{self, Cigar},
        aln::{NamedAlignment, Alignment, ReadEnd},
        contigs::{ContigNames, GenomeVersion},
    },
    bg::{
        self, BgDistr, Technology, SequencingInfo,
        insertsz::{self, InsertDistr},
        err_prof::{ErrorProfile, SingleEditDistCache},
        depth::ReadDepth,
        ser::{JsonSer, json_get},
    },
};
use super::paths;

fn check_filename(
    filename: &Path,
    should_match: &Regex,
    exp_format: &'static str,
    shouldnt_match: &Regex,
    wrong_format: &'static str,
) -> Result<(), Error> {
    if !filename.exists() {
        if filename == std::ffi::OsStr::new("!") || filename == std::ffi::OsStr::new("-")
                || filename.starts_with("/dev") {
            return Ok(());
        } else {
            log::error!("Input file {} does not exist. Continuing for now", ext::fmt::path(filename));
        }
    }

    if let Some(s) = filename.to_str() {
        if shouldnt_match.is_match(s) {
            return Err(Error::InvalidInput(format!(
                "Incorrect file format for {} (expected {}, found {}). Please check -i and -a arguments",
                ext::fmt::path(filename), exp_format, wrong_format)));
        } else if !should_match.is_match(s) {
            log::warn!("Could not guess file format for {} (expected {}). Can lead to problems later",
                ext::fmt::path(filename), exp_format);
        }
    } else {
        log::warn!("Could not guess file format for {} (expected {}). Can lead to problems later",
            ext::fmt::path(filename), exp_format);
    }
    Ok(())
}

/// Checks if input reads have fastq/fasta extensions, or input alignments have bam/cram extensions.
pub(super) fn check_input_filenames(input: &[PathBuf], alns: &Option<PathBuf>) -> Result<(), Error> {
    // Should only be run once, so there is no need for lazy static.
    let re_fastx = RegexBuilder::new(r"\.f(ast)?[aq](\.[^.]{1,3})?$").case_insensitive(true).build().unwrap();
    let re_bam = RegexBuilder::new(r"\.(bam|cram)$").case_insensitive(true).build().unwrap();
    let fastx_descr = "fasta/fastq[.gz]";
    let bam_descr = "bam/cram";

    for filename in input {
        check_filename(filename, &re_fastx, fastx_descr, &re_bam, bam_descr)?;
    }
    if let Some(filename) = alns {
        check_filename(filename, &re_bam, bam_descr, &re_fastx, fastx_descr)?;
    }
    Ok(())
}

struct Args {
    input: Vec<PathBuf>,
    alns: Option<PathBuf>,
    reference: Option<PathBuf>,
    jf_counts: Option<PathBuf>,
    output: Option<PathBuf>,
    similar_dataset: Option<PathBuf>,
    bg_region: Option<String>,

    technology: Technology,
    min_mapq: u8,
    subsampling_rate: f64,
    seed: Option<u64>,

    interleaved: bool,
    no_index: bool,
    threads: u16,
    rerun: super::Rerun,
    strobealign: PathBuf,
    minimap: PathBuf,
    samtools: PathBuf,
    jellyfish: PathBuf,
    debug: bool,
    /// Was technology explicitely provided?
    explicit_technology: bool,
    /// When calculating insert size distributions and read error profiles,
    /// ignore reads with `clipping > max_clipping * read_len`.
    max_clipping: f64,

    bg_params: bg::Params,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::new(),
            alns: None,
            reference: None,
            jf_counts: None,
            output: None,
            similar_dataset: None,
            bg_region: None,

            technology: Technology::Illumina,
            min_mapq: 20,
            subsampling_rate: 0.1,
            seed: None,

            interleaved: false,
            no_index: false,
            threads: 8,
            rerun: super::Rerun::None,
            debug: false,
            strobealign: PathBuf::from("strobealign"),
            minimap: PathBuf::from("minimap2"),
            samtools: PathBuf::from("samtools"),
            jellyfish: PathBuf::from("jellyfish"),
            explicit_technology: false,
            max_clipping: 0.02,
            bg_params: bg::Params::default(),
        }
    }
}

impl Args {
    /// Validate arguments, modifying some, if needed.
    fn validate(mut self) -> Result<Self, Error> {
        self.threads = max(self.threads, 1);
        let n_input = self.input.len();
        validate_param!(n_input > 0 || self.alns.is_some(),
            "Neither read files, nor alignment files are not provided (see -i and -a)");
        validate_param!(n_input == 0 || self.alns.is_none(),
            "Read files (-i) and an alignment file (-a) cannot be provided together");
        validate_param!(n_input != 2 || !self.interleaved,
            "Two read files (-i/--input) are provided, however, --interleaved is specified");

        if !self.has_indexed_alignment() {
            let paired_end_allowed = self.technology.paired_end_allowed();
            validate_param!(!self.is_paired_end() || paired_end_allowed,
                "Paired end reads are not supported by {}", self.technology.long_name());
            if !self.is_paired_end() && paired_end_allowed {
                log::warn!("Running in single-end mode.");
            }
        } else {
            validate_param!(self.similar_dataset.is_none(),
                "Similar dataset (-~) can only be used together with input reads (-i) \
                or unindexed alignments (-a ... --no-index)");
        }
        check_input_filenames(&self.input, &self.alns)?;

        validate_param!(self.jf_counts.is_some(), "Jellyfish counts are not provided (see -j/--jf-counts)");
        validate_param!(self.reference.is_some(), "Reference fasta file is not provided (see -r/--reference)");
        validate_param!(self.output.is_some(), "Output directory is not provided (see -o/--output)");

        validate_param!(0.0 < self.subsampling_rate && self.subsampling_rate <= 1.0,
            "Subsample rate ({}) must be within (0, 1]", self.subsampling_rate);
        if self.subsampling_rate > 0.99 {
            self.subsampling_rate = 1.0;
        }

        if self.technology == Technology::Illumina {
            self.strobealign = ext::sys::find_exe(self.strobealign)?;
        } else {
            self.minimap = ext::sys::find_exe(self.minimap)?;
        }
        self.samtools = ext::sys::find_exe(self.samtools)?;
        self.jellyfish = ext::sys::find_exe(self.jellyfish)?;

        validate_param!(0.0 <= self.max_clipping && self.max_clipping <= 1.0,
            "Max clipping ({:.5}) must be within [0, 1]", self.max_clipping);
        if self.has_indexed_alignment() {
            self.subsampling_rate = 1.0;
            self.seed = None;
        }
        self.bg_params.validate()?;
        Ok(self)
    }

    /// Returns true if need to process indexed alignment file.
    fn has_indexed_alignment(&self) -> bool {
        self.alns.is_some() && !self.no_index
    }

    /// Returns true if input reads are paired-end. Panics for indexed alignment file.
    fn is_paired_end(&self) -> bool {
        assert!(!self.has_indexed_alignment());
        self.input.len() == 2 || self.interleaved
    }
}

fn print_help(extended: bool) {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Preprocess WGS dataset.".yellow());

    println!("\n{}", "Usage:".bold());
    println!("    {} preproc -i reads1.fq [reads2.fq]  -j counts.jf -r reference.fa -o out [args]", super::PROGRAM);
    println!("    {} preproc -a reads.bam [--no-index] -j counts.jf -r reference.fa -o out [args]", super::PROGRAM);
    println!("    {} preproc -i/-a <input> -~ similar  -j counts.jf -r reference.fa -o out [args]", super::PROGRAM);
    if !extended {
        println!("\nThis is a {} help message. Please use {} to see the full help.",
            "short".red(), "-H/--full-help".green());
    }

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Reads in indexed BAM/CRAM format, already mapped to the whole genome.\n\
        {EMPTY}  Mutually exclusive with {}.",
        "-a, --alignment".green(), "FILE".yellow(), "-i/--input".green());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file. Must contain FAI index.",
        "-r, --reference".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Jellyfish k-mer counts (see README).",
        "-j, --jf-counts".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Output directory.",
        "-o, --output".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  This dataset is similar to already preprocessed dataset.\n\
        {EMPTY}  {}. Only utilizes difference in the number of reads.",
        "-~, --like".green(), "DIR".yellow(), "Use with care".red());

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Interleaved paired-end reads in single input file.",
        "-^, --interleaved".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Use input full BAM/CRAM file ({}) without index.\n\
        {EMPTY}  Single-end and paired-end interleaved ({}) data is allowed.",
        "    --no-index".green(), super::flag(), "-a".green(), "-^".green());
    println!("    {:KEY$} {:VAL$}  Sequencing technology [{}]:\n\
        {EMPTY}  sr  | illumina : short-read sequencing,\n\
        {EMPTY}    hifi         : PacBio HiFi,\n\
        {EMPTY}  pb  | pacbio   : PacBio CLR,\n\
        {EMPTY}  ont | nanopore : Oxford Nanopore.",
        "-t, --technology".green(), "STR".yellow(), super::fmt_def(defaults.technology));
    println!("    {:KEY$} {:VAL$}  Preprocess WGS data based on this background region,\n\
        {EMPTY}  preferably >3 Mb and without many duplications.\n\
        {EMPTY}  Default regions are defined for CHM13, GRCh38 and GRCh37.",
        "-b, --bg-region".green(), "STR".yellow());

    if extended {
        println!("\n{}", "Insert size and error profile estimation:".bold());
        println!("    {:KEY$} {:VAL$}  Ignore reads with mapping quality less than {} [{}].",
            "-q, --min-mapq".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(defaults.min_mapq));
        println!("    {:KEY$} {:VAL$}  Ignore reads with soft/hard clipping over {} * read length [{}].",
            "-c, --max-clipping".green(), "FLOAT".yellow(), "FLOAT".yellow(),
            super::fmt_def_f64(defaults.max_clipping));
        println!("    {:KEY$} {:VAL$}\n\
            {EMPTY}  Two p-value thresholds for filtering recruited reads:\n\
            {EMPTY}  on insert size [{}], and on edit distance [{}].",
            "    --pval-thresh".green(), "FLOAT FLOAT".yellow(),
            super::fmt_def_f64(defaults.bg_params.insert_pval),
            super::fmt_def_f64(defaults.bg_params.edit_pval));

        println!("\n{}", "Background read depth estimation:".bold());
        println!("    {:KEY$} {:VAL$}  Specie ploidy [{}].",
            "-p, --ploidy".green(), "INT".yellow(), super::fmt_def(defaults.bg_params.depth.ploidy));
        println!("    {:KEY$} {:VAL$}  Subsample input reads by this factor [{}].",
            "-S, --subsample".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.subsampling_rate));
        println!("    {:KEY$} {:VAL$}  Count read depth in windows of this size [{}].\n\
            {EMPTY}  Default: half of the mean read length.",
            "-w, --window".green(), "INT".yellow(), "auto".cyan());
        println!("    {:KEY$} {:VAL$}  Skip {} bp near the edge of the background region [{}].",
            "    --boundary".green(), "INT".yellow(), "INT".yellow(),
            super::fmt_def(PrettyU32(defaults.bg_params.depth.boundary_size)));
        println!("    {:KEY$} {:VAL$}  Ignore windows, where less than {}% k-mers are unique [{}].",
            "    --kmer-perc".green(), "FLOAT".yellow(), "FLOAT".yellow(),
            super::fmt_def_f64(defaults.bg_params.depth.uniq_kmer_perc));
        println!("    {:KEY$} {:VAL$}  This fraction of all windows is used in LOESS during read depth\n\
            {EMPTY}  estimation [{}]. Smaller values lead to less robust estimates,\n\
            {EMPTY}  larger values - to similar estimates across different GC-contents.",
            "    --frac-windows".green(), "FLOAT".yellow(),
            super::fmt_def_f64(defaults.bg_params.depth.frac_windows));
        println!("    {:KEY$} {}\n\
            {EMPTY}  Read depth estimates are blured for windows with extreme GC-content\n\
            {EMPTY}  (less than {} windows with smaller/larger GC). There, read depth\n\
            {EMPTY}  is set to the last non-extreme depth, while variance is increased\n\
            {EMPTY}  by a {} factor for each addition GC value [{} {}].",
            "    --blur-extreme".green(), "INT FLOAT".yellow(), "INT".yellow(), "FLOAT".yellow(),
            super::fmt_def(defaults.bg_params.depth.min_tail_obs),
            super::fmt_def_f64(defaults.bg_params.depth.tail_var_mult));
    }

    println!("\n{}", "Execution arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), super::fmt_def(defaults.threads));
    println!("    {:KEY$} {:VAL$}  Rerun mode [{}]. Rerun everything ({}); do not rerun\n\
        {EMPTY}  read mapping ({}); do not rerun ({}).",
        "    --rerun".green(), "STR".yellow(), super::fmt_def(defaults.rerun),
        "all".yellow(), "part".yellow(), "none".yellow());
    println!("    {:KEY$} {:VAL$}  Random seed. Ensures reproducibility for the same\n\
        {EMPTY}  input and program version.",
        "-s, --seed".green(), "INT".yellow());
    println!("    {:KEY$} {:VAL$}  Save debug CSV files.",
        "    --debug".green(), "INT".yellow());
    if extended {
        println!("    {:KEY$} {:VAL$}  Strobealign executable [{}].",
            "    --strobealign".green(), "EXE".yellow(), super::fmt_def(defaults.strobealign.display()));
        println!("    {:KEY$} {:VAL$}  Minimap2 executable    [{}].",
            "    --minimap".green(), "EXE".yellow(), super::fmt_def(defaults.minimap.display()));
        println!("    {:KEY$} {:VAL$}  Samtools executable    [{}].",
            "    --samtools".green(), "EXE".yellow(), super::fmt_def(defaults.samtools.display()));
        println!("    {:KEY$} {:VAL$}  Jellyfish executable   [{}].",
            "    --jellyfish".green(), "EXE".yellow(), super::fmt_def(defaults.jellyfish.display()));
    }

    println!("\n{}", "Other arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Show short help message.", "-h, --help".green(), "");
    println!("    {:KEY$} {:VAL$}  Show {} help message.", "-H, --full-help".green(), "", "extended".red());
    println!("    {:KEY$} {:VAL$}  Show version.", "-V, --version".green(), "");
}

fn parse_args(argv: &[String]) -> Result<Args, lexopt::Error> {
    if argv.is_empty() {
        print_help(false);
        std::process::exit(1);
    }
    use lexopt::prelude::*;
    let mut args = Args::default();
    let mut parser = lexopt::Parser::from_args(argv);

    while let Some(arg) = parser.next()? {
        match arg {
            Short('i') | Long("input") =>
                args.input = parser.values()?.take(2).map(|s| s.parse()).collect::<Result<_, _>>()?,
            Short('a') | Long("aln") | Long("alns") | Long("alignment") | Long("alignments") =>
                args.alns = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),
            Short('t') | Long("tech") | Long("technology") => {
                args.explicit_technology = true;
                args.technology = parser.value()?.parse()?;
            }
            Short('j') | Long("jf-counts") => args.jf_counts = Some(parser.value()?.parse()?),
            Short('~') | Long("like") => args.similar_dataset = Some(parser.value()?.parse()?),
            Short('b') | Long("bg") | Long("bg-region") => args.bg_region = Some(parser.value()?.parse()?),

            Short('q') | Long("min-mapq") | Long("min-mq") => args.min_mapq = parser.value()?.parse()?,
            Short('c') | Long("max-clip") | Long("max-clipping") => args.max_clipping = parser.value()?.parse()?,
            Long("pval-thresh") | Long("pval-threshold") | Long("pvalue-threshold") => {
                args.bg_params.insert_pval = parser.value()?.parse()?;
                args.bg_params.edit_pval = parser.value()?.parse()?;
            }

            Short('p') | Long("ploidy") => args.bg_params.depth.ploidy = parser.value()?.parse()?,
            Short('S') | Long("subsample") => args.subsampling_rate = parser.value()?.parse()?,
            Short('w') | Long("window") => {
                let val = parser.value()?;
                args.bg_params.depth.window_size = if val == "auto" {
                    None
                } else {
                    Some(val.parse()?)
                };
            }
            Long("boundary") => args.bg_params.depth.boundary_size = parser.value()?.parse::<PrettyU32>()?.get(),
            Long("kmer-perc") | Long("kmer-percentage") | Long("kmer-percentile") =>
                args.bg_params.depth.uniq_kmer_perc = parser.value()?.parse()?,
            Long("frac-windows") | Long("fraction-windows") =>
                args.bg_params.depth.frac_windows = parser.value()?.parse()?,
            Long("blur-extreme") => {
                args.bg_params.depth.min_tail_obs = parser.value()?.parse()?;
                args.bg_params.depth.tail_var_mult = parser.value()?.parse()?;
            }

            Short('^') | Long("interleaved") => args.interleaved = true,
            Long("no-index") => args.no_index = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Long("rerun") => args.rerun = parser.value()?.parse()?,
            Long("debug") => args.debug = true,
            Long("strobealign") => args.strobealign = parser.value()?.parse()?,
            Long("minimap") | Long("minimap2") => args.minimap = parser.value()?.parse()?,
            Long("samtools") => args.samtools = parser.value()?.parse()?,
            Long("jellyfish") => args.jellyfish = parser.value()?.parse()?,
            Short('s') | Long("seed") => args.seed = Some(parser.value()?.parse()?),

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
            format!("Reference genome {} does not contain region {}", ext::fmt::path(ref_filename), s)))?;
        return if contigs.in_bounds(&region) {
            Ok(region)
        } else {
            Err(Error::InvalidInput(format!("Region {} is out of bounds", s)))
        };
    }

    let genome = GenomeVersion::guess(contigs)
        .ok_or_else(|| Error::RuntimeError(
            "Could not recognize reference genome. Please provide background region (-b) explicitely, \
            preferably >3 Mb long and without significant duplications".to_owned()))?;
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
    Err(Error::InvalidInput(format!("Reference genome {} does not contain any of the default background regions. \
        Please provide background region (-b) explicitely, preferably >3 Mb long and without significant duplications",
        ext::fmt::path(ref_filename))))
}

/// Returns handle, which returns the total number of consumed reads/read pairs after finishing.
fn set_mapping_stdin(
    args: &Args,
    contigs: &ContigNames,
    child_stdin: ChildStdin,
    rng: &mut XoshiroRng,
) -> Result<thread::JoinHandle<Result<u64, Error>>, Error>
{
    let rate = args.subsampling_rate;
    let mut local_rng = rng.clone();
    rng.long_jump();
    let mut writer = io::BufWriter::new(child_stdin);
    let handle = fastx::process_readers!(args, Some(contigs); let {mut} reader; {
        thread::spawn(move || {
            Ok(if rate < 1.0 {
                reader.subsample(&mut writer, rate, &mut local_rng)?
            } else {
                reader.copy(&mut writer)?
            })
        })
    });
    Ok(handle)
}

fn create_mapping_command(args: &Args, seq_info: &SequencingInfo, ref_filename: &Path) -> Command {
    let mut command;
    if seq_info.technology() == Technology::Illumina {
        // Strobealign for short reads.
        command = Command::new(&args.strobealign);
        command.args(&[
            "-N", "0",     // Retain 0 additional alignments,
            "-R", "0",     // Do not rescue reads,
            "-U",          // Do not output unmapped reads,
            "-f", "0.001", // Discard more minimizers to speed up alignment,
            "--eqx",       // Output X/= instead of M operations,
            "-t", &args.threads.to_string(), // Specify the number of threads,
            "-r", &format!("{:.0}", seq_info.mean_read_len()), // Provide mean read length.
            "--no-progress",
            ]);
        if args.is_paired_end() {
            command.arg("--interleaved");
        }
    } else {
        // Minimap2 for long reads.
        command = Command::new(&args.minimap);
        command.args(&[
            "-a", // Output SAM format,
            "-x", seq_info.technology().minimap_preset(), // Set mapping preset.
            "-N", "2",     // Examine at most 2 additional alignments,
            "--secondary=no", // Do not output secondary alignments,
            "-f", "0.001", // Discard more minimizers to speed up alignment,
            "--eqx",       // Output X/= instead of M operations,
            "-t", &args.threads.to_string(), // Specify the number of threads,
            ]);
    }
    // Provide paths to the reference and pipe reads.
    command.arg(&ref_filename).arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    command
}

fn first_step_str(args: &Args) -> String {
    let mut s = String::new();
    if args.subsampling_rate == 1.0 {
        s.push_str("_interleave_");
    } else {
        write!(s, "_subsample_ --rate {}", args.subsampling_rate).unwrap();
        if let Some(seed) = args.seed {
            write!(s, " --seed {}", seed).unwrap();
        }
    }
    s.push_str(" -i ");
    if let Some(bam_filename) = &args.alns {
        s.push_str(&ext::fmt::path(bam_filename));
    } else {
        s.push_str(&ext::fmt::path(&args.input[0]));
    }
    if args.input.len() > 1 {
        write!(s, " {}", ext::fmt::path(&args.input[1])).unwrap();
    }
    s.push_str(" | ");
    s
}

/// Mapping parameters will be stored at `analysis/bg/MAPPING_PARAMS`.
const MAPPING_PARAMS: &'static str = "mapping.json";

struct MappingParams {
    technology: Technology,
    min_mapq: u8,
    subsampling: f64,
    bg_region: String,
}

impl MappingParams {
    fn new(args: &Args, bg_region: &Interval) -> Self {
        MappingParams {
            technology: args.technology,
            min_mapq: args.min_mapq,
            subsampling: args.subsampling_rate,
            bg_region: bg_region.to_string(),
        }
    }

    /// Compares mapping parameters against parameters from the previous run.
    /// Returns `Ok` if two parameters match, and there is no reason for remapping.
    /// Otherwise, returns `Err` with explanation.
    fn compare(&self, path: &Path) -> Result<(), String> {
        if !path.exists() {
            log::error!("Cannot find previous mapping parameters at {}. Continuing anyway", ext::fmt::path(path));
            return Ok(())
        }
        let old = match ext::sys::load_json(&path).map_err(Error::from).and_then(|json| MappingParams::load(&json)) {
            Err(e) => {
                return Err(format!("Cannot load previous mapping parameters: {}", e.display()));
            }
            Ok(val) => val,
        };
        if self.technology != old.technology {
            Err(format!("Sequencing technologies do not match ({} -> {})", old.technology, self.technology))
        } else if self.min_mapq < old.min_mapq {
            Err(format!("Minimal mapping quality decreased ({} -> {})", old.min_mapq, self.min_mapq))
        } else if (self.subsampling - old.subsampling).abs() > 1e-8 {
            Err(format!("Subsampling rate changed ({} -> {})", old.subsampling, self.subsampling))
        } else if self.bg_region != old.bg_region {
            Err(format!("Background region has changed ({} -> {})", old.bg_region, self.bg_region))
        } else {
            Ok(())
        }
    }
}

impl JsonSer for MappingParams {
    fn save(&self) -> json::JsonValue {
        json::object!{
            technology: self.technology.to_str(),
            min_mapq: self.min_mapq,
            subsampling: self.subsampling,
            bg_region: &self.bg_region as &str,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj => technology (as_str), min_mapq (as_u8), subsampling (as_f64), bg_region (as_str));
        let technology = Technology::from_str(technology).map_err(|e| Error::ParsingError(e))?;
        Ok(Self {
            technology, min_mapq, subsampling,
            bg_region: bg_region.to_owned(),
        })
    }
}

/// Returns true if read mapping can be skipped.
/// Panics if output BAM file exists, but cannot be used.
fn need_mapping(args: &Args, out_dir: &Path, out_bam: &Path, bg_region: &Interval) -> Result<bool, Error> {
    let params_path = out_dir.join(MAPPING_PARAMS);
    let curr_params = MappingParams::new(args, bg_region);
    if !out_bam.exists() {
        let mut params_file = ext::sys::create_file(&params_path)?;
        curr_params.save().write_pretty(&mut params_file, 4).map_err(add_path!(params_path))?;
        return Ok(true);
    }
    match curr_params.compare(&params_path) {
        Ok(()) => {
            log::warn!("BAM file {} exists, skipping read mapping", ext::fmt::path(out_bam));
            Ok(false)
        }
        Err(s) => {
            log::error!("Problem comparing read mapping parameters against {}", ext::fmt::path(&params_path));
            log::error!("    {}", s);
            log::error!("Please use `--rerun full`, or remove file {} to force partial analysis",
                ext::fmt::path(&params_path));
            std::process::exit(1)
        }
    }
}

/// Map reads to the whole reference genome, and then take only reads mapped to the corresponding BED file.
/// Subsample reads if the corresponding rate is less than 1.
/// Return the total number of reads/read pairs.
fn run_mapping(
    args: &Args,
    seq_info: &mut SequencingInfo,
    ref_filename: &Path,
    out_dir: &Path,
    out_bam: &Path,
    bg_region: &BgRegion,
    rng: &mut XoshiroRng,
) -> Result<(), Error>
{
    if !need_mapping(args, out_dir, out_bam, &bg_region.interval)? {
        return Ok(());
    }

    let start = Instant::now();
    let mut mapping = create_mapping_command(args, seq_info, ref_filename);
    let mapping_exe = PathBuf::from(mapping.get_program().to_owned());
    let mut mapping_child = mapping.spawn().map_err(add_path!(mapping_exe))?;
    let mapping_stdin = mapping_child.stdin.take();
    let mapping_stdout = mapping_child.stdout.take();
    let mut pipe_guard = ext::sys::PipeGuard::new(mapping_exe, mapping_child);
    let handle = set_mapping_stdin(args, &bg_region.ref_contigs, mapping_stdin.unwrap(), rng)?;

    let tmp_bed = out_dir.join("tmp.bed");
    let mut bed_file = ext::sys::create_file(&tmp_bed)?;
    writeln!(bed_file, "{}", bg_region.interval .bed_fmt()).map_err(add_path!(tmp_bed))?;
    std::mem::drop(bed_file);

    let tmp_bam = out_dir.join("tmp.bam");
    let mut samtools = Command::new(&args.samtools);
    // See SAM flags here: https://broadinstitute.github.io/picard/explain-flags.html.
    samtools.args(&["view",
            "-b", // Output BAM.
            // Ignore reads where any of the mates is unmapped,
            // + ignore secondary & supplementary alignments + ignore failed checks.
            "-F", "3852",
            "-q", &args.min_mapq.to_string(),
            ])
        .arg("-L").arg(&tmp_bed)
        .arg("-o").arg(&tmp_bam)
        .stdin(Stdio::from(mapping_stdout.unwrap()))
        .stdout(Stdio::piped()).stderr(Stdio::piped());
    log::debug!("    {}{} | {}", first_step_str(&args), ext::fmt::command(&mapping), ext::fmt::command(&samtools));
    let samtools_child = samtools.spawn().map_err(add_path!(args.samtools))?;
    pipe_guard.push(args.samtools.clone(), samtools_child);
    pipe_guard.wait()?;
    log::debug!("    Finished in {}", ext::fmt::Duration(start.elapsed()));
    let total_reads = handle.join()
        .map_err(|e| Error::RuntimeError(format!("Read mapping failed: {:?}", e)))??;
    seq_info.set_total_reads(total_reads);
    fs::rename(&tmp_bam, out_bam).map_err(add_path!(tmp_bam, out_bam))?;
    fs::remove_file(&tmp_bed).map_err(add_path!(tmp_bed))?;
    Ok(())
}

/// Loads records and converts them to Alignments.
/// All returned alignments are primary, have MAPQ over the `args.min_mapq`,
/// and clipping under the threshold `bg_params.max_clipping`.
///
/// Ignore reads, for which `cigar_getter` returns None.
///
/// Second returned argument is true if records are paired-end.
///
/// Alignment probabilities are not set.
fn load_alns(
    reader: &mut impl bam::Read,
    cigar_getter: impl Fn(&bam::Record) -> Option<Cigar>,
    contigs: &Arc<ContigNames>,
    args: &Args
) -> Result<(Vec<NamedAlignment>, bool), Error>
{
    let min_mapq = args.min_mapq;
    let max_clipping = args.max_clipping;
    let mut paired_counts = [0_u64, 0];

    let mut alns = Vec::new();
    let mut record = bam::Record::new();
    let mut ignored = 0;
    let mut wo_cigar = 0;
    while let Some(()) = reader.read(&mut record).transpose()? {
        if record.flags() & 3844 == 0 && record.mapq() >= min_mapq && cigar::clipping_rate(&record) <= max_clipping {
            if let Some(cigar) = cigar_getter(&record) {
                alns.push(NamedAlignment::new(&record, cigar, ReadEnd::from_record(&record),
                    Arc::clone(contigs), f64::NAN));
                paired_counts[usize::from(record.is_paired())] += 1;
            } else {
                wo_cigar += 1;
            }
        } else {
            ignored += 1;
        }
    }
    if paired_counts[0] > 0 && paired_counts[1] > 0 {
        return Err(Error::InvalidData(format!("BAM file contains both paired and unpaired reads")));
    }
    if alns.is_empty() {
        return Err(Error::InvalidData(format!("BAM file contains no reads in the target region")));
    }
    log::debug!("    Loaded {} alignments, discarded {}", alns.len(), ignored);
    if wo_cigar > 0 {
        log::warn!("    Could not load CIGAR for {} records", wo_cigar);
    }
    Ok((alns, paired_counts[1] > 0))
}

const MEAN_LEN_RECORDS: usize = 5000;

/// Calculate mean read length from existing alignments.
fn read_len_from_alns(alns: &[NamedAlignment]) -> f64 {
    let n = min(alns.len(), MEAN_LEN_RECORDS);
    alns[..n].iter()
        .map(|aln| f64::from(aln.cigar().query_len()))
        .sum::<f64>() / n as f64
}

/// Calculate mean read length from input reads.
fn read_len_from_reads(args: &Args, contigs: Option<&ContigNames>) -> Result<f64, Error> {
    log::info!("Calculating mean read length");
    fastx::process_readers!(args, contigs; let {mut} reader; { fastx::mean_read_len(&mut reader, MEAN_LEN_RECORDS) })
}

fn estimate_bg_from_paired(
    alns: Vec<NamedAlignment>,
    seq_info: SequencingInfo,
    args: &Args,
    opt_out_dir: Option<&Path>,
    bg_region: &BgRegion,
) -> Result<BgDistr, Error>
{
    // Group reads into pairs, and estimate insert size from them.
    let pair_ixs = insertsz::group_mates(&alns)?;
    let insert_distr = InsertDistr::estimate(&alns, &pair_ixs, opt_out_dir)?;
    let conf_lvl = 1.0 - args.bg_params.insert_pval;
    let (min_insert, max_insert) = insert_distr.confidence_interval(conf_lvl);
    log::info!("    Allowed insert size: [{}, {}]  ({}%-confidence interval)",
        min_insert, max_insert, crate::math::fmt_signif(100.0 * conf_lvl, 5));

    let interval = &bg_region.interval;
    let windows = bg::Windows::create(interval, &bg_region.sequence, &bg_region.kmer_counts, &seq_info,
        &args.bg_params.depth, opt_out_dir)?;

    // Estimate error profile from read pairs with appropriate insert size.
    let mut errprof_alns = Vec::with_capacity(pair_ixs.len() * 2);
    for &(i, j) in pair_ixs.iter() {
        let first = &alns[i];
        let second = &alns[j];
        let insert = first.insert_size(second.deref());
        if min_insert <= insert && insert <= max_insert {
            errprof_alns.push(first);
            errprof_alns.push(second);
        }
    }
    let err_prof = ErrorProfile::estimate(&errprof_alns, interval, &windows, opt_out_dir)?;
    let edit_dist_cache = SingleEditDistCache::new(&err_prof, args.bg_params.edit_pval);
    edit_dist_cache.print_log(seq_info.mean_read_len());

    // Estimate backgorund read depth from read pairs with both good probabilities and good insert size prob.
    let mut depth_alns: Vec<&Alignment> = Vec::with_capacity(errprof_alns.len());
    for chunk in errprof_alns.chunks_exact(2) {
        let aln1 = chunk[0];
        let aln2 = chunk[1];
        let dist1 = aln1.count_region_operations(interval).edit_distance();
        let dist2 = aln2.count_region_operations(interval).edit_distance();
        if dist1.edit() <= edit_dist_cache.get(dist1.read_len())
                && dist2.edit() <= edit_dist_cache.get(dist2.read_len()) {
            depth_alns.push(aln1.deref());
            depth_alns.push(aln2.deref());
        }
    }
    let depth_distr = ReadDepth::estimate(&depth_alns, &windows,
        &args.bg_params.depth, args.subsampling_rate, true, &seq_info, opt_out_dir)?;
    Ok(BgDistr::new(seq_info, insert_distr, err_prof, depth_distr))
}

fn estimate_bg_from_unpaired(
    alns: Vec<NamedAlignment>,
    seq_info: SequencingInfo,
    args: &Args,
    opt_out_dir: Option<&Path>,
    bg_region: &BgRegion,
) -> Result<BgDistr, Error>
{
    let insert_distr = InsertDistr::undefined();
    let interval = &bg_region.interval;
    let windows = bg::Windows::create(interval, &bg_region.sequence, &bg_region.kmer_counts, &seq_info,
        &args.bg_params.depth, opt_out_dir)?;
    let err_prof = ErrorProfile::estimate(&alns, interval, &windows, opt_out_dir)?;
    let edit_dist_cache = SingleEditDistCache::new(&err_prof, args.bg_params.edit_pval);
    edit_dist_cache.print_log(seq_info.mean_read_len());

    let filt_alns: Vec<&Alignment> = alns.iter()
        .filter(|aln| {
            let dist = aln.count_region_operations(interval).edit_distance();
            dist.edit() <= edit_dist_cache.get(dist.read_len())
        })
        .map(Deref::deref)
        .collect();
    let depth_distr = ReadDepth::estimate(&filt_alns, &windows,
        &args.bg_params.depth, args.subsampling_rate, true, &seq_info, opt_out_dir)?;
    Ok(BgDistr::new(seq_info, insert_distr, err_prof, depth_distr))
}

/// Estimate background distributions from input reads or existing alignments.
fn estimate_bg_distrs(
    args: &Args,
    out_dir: &Path,
    bg_region: &BgRegion,
) -> Result<BgDistr, Error>
{
    let opt_out_dir = if args.debug { Some(out_dir) } else { None };
    let ref_filename = args.reference.as_ref().unwrap();
    let aln_is_paired_end: bool;
    let alns: Vec<NamedAlignment>;
    let mut seq_info: SequencingInfo;

    if args.has_indexed_alignment() {
        let alns_filename = args.alns.as_ref().unwrap();
        log::debug!("Loading mapped reads into memory ({})", ext::fmt::path(alns_filename));
        let mut bam_reader = bam::IndexedReader::from_path(alns_filename)?;
        fastx::set_reference(alns_filename, &mut bam_reader, &args.reference, Some(&bg_region.ref_contigs))?;

        let interval = &bg_region.interval;
        let interval_start = interval.start();
        let ref_seq = &bg_region.sequence;
        bam_reader.fetch((interval.contig_name(), i64::from(interval_start), i64::from(interval.end())))?;
        (alns, aln_is_paired_end) = load_alns(&mut bam_reader,
            |record| Cigar::infer_ext_cigar(record, ref_seq, interval_start),
            &bg_region.ref_contigs, args)?;
        if aln_is_paired_end && !args.technology.paired_end_allowed() {
            return Err(Error::InvalidInput(format!("Paired end reads are not supported by {}",
                args.technology.long_name())));
        }
        seq_info = SequencingInfo::new(read_len_from_alns(&alns), args.technology, args.explicit_technology)?;
    } else {
        seq_info = SequencingInfo::new(read_len_from_reads(args, Some(&bg_region.ref_contigs))?,
            args.technology, args.explicit_technology)?;
        log::info!("Mean read length = {:.1}", seq_info.mean_read_len());

        let bam_filename = out_dir.join("aln.bam");
        log::info!("Mapping reads to the reference");
        let mut rng = init_rng(args.seed);
        run_mapping(args, &mut seq_info, &ref_filename, &out_dir, &bam_filename, bg_region, &mut rng)?;

        log::debug!("Loading mapped reads into memory ({})", ext::fmt::path(&bam_filename));
        let mut bam_reader = bam::Reader::from_path(&bam_filename)?;
        let get_cigar = |record: &bam::Record| -> Option<Cigar> {
            let cigar = Cigar::from_raw(record.raw_cigar());
            assert!(!cigar.has_hard_clipping(), "Cannot process primary alignments with hard clipping");
            Some(cigar)
        };
        (alns, aln_is_paired_end) = load_alns(&mut bam_reader, get_cigar, &bg_region.ref_contigs, args)?;
        if aln_is_paired_end != args.is_paired_end() {
            return Err(Error::RuntimeError(format!(
                "Input data is {}-end, while alignment file is {}-end. \
                Perhaps preprocessing was run multiple times with different inputs?",
                if args.is_paired_end() { "paired" } else { "single" },
                if aln_is_paired_end { "paired" } else { "single" },
            )));
        }
    }

    if aln_is_paired_end {
        estimate_bg_from_paired(alns, seq_info, args, opt_out_dir, bg_region)
    } else {
        estimate_bg_from_unpaired(alns, seq_info, args, opt_out_dir, bg_region)
    }
}

/// Assume that this WGS dataset is similar to another dataset, and use its parameters.
fn estimate_like(args: &Args, like_dir: &Path) -> Result<BgDistr, Error> {
    let similar_path = like_dir.join(paths::BG_DISTR);
    log::info!("Loading distribution parameters from {}", ext::fmt::path(&similar_path));
    let similar_distr = BgDistr::load_from(&similar_path, &like_dir.join(paths::SUCCESS))?;
    let similar_seq_info = similar_distr.seq_info();
    let similar_n_reads = match similar_seq_info.total_reads() {
        Some(count) => count,
        None => return Err(Error::InvalidInput(format!(
            "Cannot use similar dataset {}: total number of reads was not estimated", ext::fmt::path(like_dir)))),
    };

    let seq_info = SequencingInfo::new(read_len_from_reads(args, None)?, args.technology, args.explicit_technology)?;
    if similar_seq_info.technology() != seq_info.technology() {
        return Err(Error::InvalidInput(format!(
            "Cannot use similar dataset {}: different sequencing technology ({} and {})",
            ext::fmt::path(like_dir), similar_seq_info.technology(), seq_info.technology())));
    } else if similar_distr.insert_distr().is_paired_end() != args.is_paired_end() {
        return Err(Error::InvalidInput(format!(
            "Cannot use similar dataset {}: paired-end status does not match", ext::fmt::path(like_dir))));
    } else if !seq_info.technology().is_read_len_similar(seq_info.mean_read_len(), similar_seq_info.mean_read_len()) {
        return Err(Error::InvalidInput(format!(
            "Cannot use similar dataset {}: read lengths are different ({:.0} and {:.0})",
            ext::fmt::path(like_dir), seq_info.mean_read_len(), similar_seq_info.mean_read_len())));
    }
    let mut total_reads = if let Some(bam_filename) = args.alns.as_ref() {
        log::info!("Counting reads in {}", ext::fmt::path(bam_filename));
        fastx::count_reads_bam(bam_filename, &args.samtools, &args.reference, args.threads)?
    } else {
        log::info!("Counting reads in {}", ext::fmt::path(&args.input[0]));
        fastx::count_reads_fastx(&args.input[0])?
    };
    if args.interleaved {
        total_reads /= 2;
    }
    log::debug!("    Input contains {} reads", total_reads);

    let mut new_distr = similar_distr.clone();
    new_distr.set_seq_info(seq_info);
    // NOTE: total reads are deliberately not provided to `seq_info`, so as not to propagate errors.
    let rate = total_reads as f64 / similar_n_reads as f64;
    if rate < 0.1 {
        return Err(Error::InvalidInput(format!(
            "Read depth changed too much (by a factor of {:.4}), please estimate parameters anew", rate)))
    } else if (rate - 1.0).abs() < 0.01 {
        log::debug!("    Almost identical read depth");
    } else {
        log::debug!("    Adapt read depth to new dataset: average changed by {:.4}", rate);
        new_distr.depth_mut().mul_depth(rate);
    }
    Ok(new_distr)
}

/// Information about the background region.
struct BgRegion {
    ref_contigs: Arc<ContigNames>,
    /// Background region.
    interval: Interval,
    /// Sequence of the background region.
    sequence: Vec<u8>,
    /// k-mer counts on the background region.
    kmer_counts: KmerCounts,
}

impl BgRegion {
    fn new(args: &Args) -> Result<Self, Error> {
        let ref_filename = args.reference.as_ref().unwrap();
        let (ref_contigs, mut ref_fasta) = ContigNames::load_indexed_fasta("ref", &ref_filename)?;
        let ref_contigs = Arc::new(ref_contigs);
        let interval = select_bg_interval(&ref_filename, &ref_contigs, &args.bg_region)?;

        let kmer_getter = JfKmerGetter::new(args.jellyfish.clone(), args.jf_counts.clone().unwrap())?;
        let sequence = interval.fetch_seq(&mut ref_fasta)?;
        log::info!("Calculating k-mer counts on the background region");
        let kmer_counts = kmer_getter.fetch([sequence.clone()])?;
        Ok(Self { ref_contigs, interval, sequence, kmer_counts })
    }
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();
    let out_dir = args.output.as_ref().unwrap();
    if !args.rerun.prepare_dir(&out_dir)? {
        std::process::exit(0);
    }
    ext::sys::mkdir(out_dir)?;

    let bg_distr = if let Some(similar_dataset) = args.similar_dataset.as_ref() {
        estimate_like(&args, similar_dataset)?
    } else {
        let bg_region = BgRegion::new(&args)?;
        estimate_bg_distrs(&args, &out_dir, &bg_region)?
    };

    let distr_filename = out_dir.join(paths::BG_DISTR);
    let mut distr_file = ext::sys::create_gzip(&distr_filename)?;
    bg_distr.save().write_pretty(&mut distr_file, 4).map_err(add_path!(distr_filename))?;
    log::info!("Success. Total time: {}", ext::fmt::Duration(timer.elapsed()));
    super::write_success_file(out_dir.join(paths::SUCCESS))?;
    Ok(())
}
