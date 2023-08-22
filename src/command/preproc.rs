//! Preprocess WGS dataset.

use std::{
    fs, io, thread,
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
use const_format::{str_repeat, concatcp};
use htslib::bam;
use crate::{
    ext::{
        self,
        rand::{init_rng, XoshiroRng},
        fmt::PrettyU32,
    },
    err::{Error, validate_param, add_path},
    seq::{
        ContigNames, Interval,
        kmers::KmerCounts,
        fastx::{self, FastxRead},
        cigar::{self, Cigar},
        aln::{Alignment, LightAlignment, ReadEnd},
    },
    bg::{
        self, BgDistr, Technology, SequencingInfo,
        insertsz::{self, InsertDistr},
        err_prof::ErrorProfile,
        depth::ReadDepth,
        ser::{JsonSer, json_get},
    },
};
use super::paths;

/// Parameters, that need to be saved between preprocessing executions.
/// On a new run, rerun is recommended if the parameters have changed.
struct Params {
    technology: Technology,
    min_mapq: u8,
    subsampling_rate: f64,
    seed: Option<u64>,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            technology: Technology::Illumina,
            min_mapq: 20,
            subsampling_rate: 0.1,
            seed: None,
        }
    }
}

impl Params {
    fn validate(&mut self) -> Result<(), Error> {
        validate_param!(0.0 < self.subsampling_rate && self.subsampling_rate <= 1.0,
            "Subsample rate ({}) must be within (0, 1]", self.subsampling_rate);
        if self.subsampling_rate > 0.99 {
            self.subsampling_rate = 1.0;
        }
        Ok(())
    }

    /// Loads parameters from the previous run. If the parameters do not match, or cannot be loaded, returns true.
    fn need_rerun(&self, from_reads: bool, path: &Path) -> bool {
        if !path.exists() {
            log::error!("Cannot find old parameters at {}", ext::fmt::path(path));
            return true;
        }
        let old_params = match ext::sys::load_json(&path)
                .map_err(|e| Error::from(e))
                .and_then(|json| Params::load(&json)) {
            Err(_) => {
                log::error!("Cannot load old parameters from {}", ext::fmt::path(path));
                return true;
            }
            Ok(val) => val,
        };
        if self.technology != old_params.technology {
            log::error!("Sequencing technology has changed ({} -> {})", old_params.technology, self.technology);
            true
        } else if from_reads && self.min_mapq != old_params.min_mapq {
            log::error!("Minimal MAPQ has changed ({} -> {})", old_params.min_mapq, self.min_mapq);
            true
        } else if from_reads && self.subsampling_rate != old_params.subsampling_rate {
            log::error!("Subsampling rate has changed ({} -> {})", old_params.subsampling_rate, self.subsampling_rate);
            true
        } else {
            // Not critical.
            if from_reads && self.seed != old_params.seed {
                log::warn!("Subsampling seed has changed ({:?} -> {:?})",
                    old_params.seed, self.seed);
            }
            false
        }
    }
}

impl JsonSer for Params {
    fn save(&self) -> json::JsonValue {
        json::object!{
            technology: self.technology.to_str(),
            min_mapq: self.min_mapq,
            subsampling_rate: self.subsampling_rate,
            subsampling_seed: self.seed,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> technology (as_str), min_mapq (as_u8), subsampling_rate (as_f64), subsampling_seed? (as_u64));
        let technology = Technology::from_str(technology).map_err(|e| Error::ParsingError(e))?;
        Ok(Self {
            technology, min_mapq, subsampling_rate,
            seed: subsampling_seed,
        })
    }
}

struct Args {
    input: Vec<PathBuf>,
    alns: Option<PathBuf>,
    database: Option<PathBuf>,
    reference: Option<PathBuf>,
    output: Option<PathBuf>,

    interleaved: bool,
    threads: u16,
    rerun: super::Rerun,
    strobealign: PathBuf,
    minimap: PathBuf,
    samtools: PathBuf,
    debug: bool,
    /// Was technology explicitely provided?
    explicit_technology: bool,

    /// When calculating insert size distributions and read error profiles,
    /// ignore reads with `clipping > max_clipping * read_len`.
    max_clipping: f64,

    params: Params,
    bg_params: bg::Params,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::new(),
            alns: None,
            database: None,
            reference: None,
            output: None,

            interleaved: false,
            threads: 8,
            rerun: super::Rerun::None,
            debug: false,
            strobealign: PathBuf::from("strobealign"),
            minimap: PathBuf::from("minimap2"),
            samtools: PathBuf::from("samtools"),
            explicit_technology: false,

            max_clipping: 0.02,
            params: Params::default(),
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
        if n_input > 0 {
            let paired_end_allowed = self.params.technology.paired_end_allowed();
            validate_param!(self.is_single_end() || paired_end_allowed,
                "Paired end reads are not supported by {}", self.params.technology.long_name());
            if self.is_single_end() && paired_end_allowed {
                log::warn!("Running in single-end mode.");
            }
        }

        validate_param!(self.database.is_some(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.reference.is_some(), "Reference fasta file is not provided (see -r/--reference)");
        validate_param!(self.output.is_some(), "Output directory is not provided (see -o/--output)");

        if self.params.technology == Technology::Illumina {
            self.strobealign = ext::sys::find_exe(self.strobealign)?;
        } else {
            self.minimap = ext::sys::find_exe(self.minimap)?;
        }
        self.samtools = ext::sys::find_exe(self.samtools)?;

        validate_param!(0.0 <= self.max_clipping && self.max_clipping <= 1.0,
            "Max clipping ({:.5}) must be within [0, 1]", self.max_clipping);
        if self.alns.is_some() {
            self.params.subsampling_rate = 1.0;
            self.params.seed = None;
        }

        self.params.validate()?;
        self.bg_params.validate()?;
        Ok(self)
    }

    fn is_paired_end(&self) -> bool {
        self.input.len() == 2 || self.interleaved
    }

    fn is_single_end(&self) -> bool {
        !self.is_paired_end()
    }
}

fn print_help(extended: bool) {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Preprocess WGS dataset.".yellow());

    println!("\n{} {} preproc -i reads1.fq [reads2.fq] -d db -r reference.fa -o out [arguments]",
        "Usage:".bold(), super::PROGRAM);
    if !extended {
        println!("\nThis is a {} help message. Please use {} to see the full help.",
            "short".red(), "-H/--full-help".green());
    }

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Reads in indexed BAM/CRAM format, already mapped to the whole genome.\n\
        {EMPTY}  Mutually exclusive with {}.\n\
        {EMPTY}  This method is much faster, but {}.",
        "-a, --alignment".green(), "FILE".yellow(), "-i/--input".green(), "less accurate".red());
    println!("    {:KEY$} {:VAL$}  Database directory (initialized with {} & {}).",
        "-d, --db".green(), "DIR".yellow(), concatcp!(super::PROGRAM, " create").underline(), "add".underline());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file. Must contain FAI index.",
        "-r, --reference".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Output directory.",
        "-o, --output".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Input reads are interleaved.",
        "-^, --interleaved".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Sequencing technology [{}]:\n\
        {EMPTY}  sr  | illumina : short-read sequencing,\n\
        {EMPTY}    hifi         : PacBio HiFi,\n\
        {EMPTY}  pb  | pacbio   : PacBio CLR,\n\
        {EMPTY}  ont | nanopore : Oxford Nanopore.",
        "-t, --technology".green(), "STR".yellow(), super::fmt_def(defaults.params.technology));

    if extended {
        println!("\n{}", "Insert size and error profile estimation:".bold());
        println!("    {:KEY$} {:VAL$}  Ignore reads with mapping quality less than {} [{}].",
            "-q, --min-mapq".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(defaults.params.min_mapq));
        println!("    {:KEY$} {:VAL$}  Ignore reads with soft/hard clipping over {} * read length [{}].",
            "-c, --max-clipping".green(), "FLOAT".yellow(), "FLOAT".yellow(),
            super::fmt_def_f64(defaults.max_clipping));
        println!("    {:KEY$} {:VAL$}\n\
            {EMPTY}  Two p-value thresholds for filtering recruited reads:\n\
            {EMPTY}  on insert size [{}], and on edit distance [{}].",
            "    --pval-thresh".green(), "FLOAT FLOAT".yellow(),
            super::fmt_def_f64(defaults.bg_params.insert_pval),
            super::fmt_def_f64(defaults.bg_params.edit_pval));
        println!("    {:KEY$} {:VAL$}  Multiply error rates by this factor, in order to correct for\n\
            {EMPTY}  read mappings missed due to higher error rate [{}].",
            "-m, --err-mult".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.bg_params.err_rate_mult));

        println!("\n{}", "Background read depth estimation:".bold());
        println!("    {:KEY$} {:VAL$}  Specie ploidy [{}].",
            "-p, --ploidy".green(), "INT".yellow(), super::fmt_def(defaults.bg_params.depth.ploidy));
        println!("    {:KEY$} {:VAL$}  Subsample input reads by this factor [{}].",
            "    --subsample".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.params.subsampling_rate));
        println!("    {:KEY$} {:VAL$}  Count read depth in windows of this size [{}].\n\
            {EMPTY}  Default: half of the mean read length.",
            "-w, --window".green(), "INT".yellow(), "auto".cyan());
        println!("    {:KEY$} {:VAL$}  Skip {} bp near the edge of the background region [{}].",
            "    --boundary".green(), "INT".yellow(), "INT".yellow(),
            super::fmt_def(PrettyU32(defaults.bg_params.depth.boundary_size)));
        println!("    {:KEY$} {:VAL$}  Ignore windows, where less than {}% k-mers are unique [{}].",
            "    --kmer-perc".green(), "FLOAT".yellow(), "FLOAT".yellow(),
            super::fmt_def_f64(defaults.bg_params.depth.kmer_perc));
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

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), super::fmt_def(defaults.threads));
    println!("    {:KEY$} {:VAL$}  Rerun mode [{}]. Rerun everything ({}); do not rerun\n\
        {EMPTY}  read mapping ({}); do not rerun ({}).",
        "    --rerun".green(), "STR".yellow(), super::fmt_def(defaults.rerun),
        "all".yellow(), "part".yellow(), "none".yellow());
    println!("    {:KEY$} {:VAL$}  Random seed. Ensures reproducibility for the same\n\
        {EMPTY}  input and program version.",
        "-s, --seed".green(), "INT".yellow());
    println!("    {:KEY$} {:VAL$}  Create more files with debug information.",
        "    --debug".green(), super::flag());
    if extended {
        println!("    {:KEY$} {:VAL$}  Strobealign executable [{}].",
            "    --strobealign".green(), "EXE".yellow(), super::fmt_def(defaults.strobealign.display()));
        println!("    {:KEY$} {:VAL$}  Minimap2 executable    [{}].",
            "    --minimap".green(), "EXE".yellow(), super::fmt_def(defaults.minimap.display()));
        println!("    {:KEY$} {:VAL$}  Samtools executable    [{}].",
            "    --samtools".green(), "EXE".yellow(), super::fmt_def(defaults.samtools.display()));
    }

    println!("\n{}", "Other parameters:".bold());
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
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),
            Short('t') | Long("technology") => {
                args.explicit_technology = true;
                args.params.technology = parser.value()?.parse()?;
            }

            Short('q') | Long("min-mapq") | Long("min-mq") => args.params.min_mapq = parser.value()?.parse()?,
            Short('c') | Long("max-clip") | Long("max-clipping") => args.max_clipping = parser.value()?.parse()?,
            Long("pval-thresh") | Long("pval-threshold") | Long("pvalue-threshold") => {
                args.bg_params.insert_pval = parser.value()?.parse()?;
                args.bg_params.edit_pval = parser.value()?.parse()?;
            }
            Short('m') | Long("err-mult") | Long("err-multiplier") =>
                args.bg_params.err_rate_mult = parser.value()?.parse()?,

            Short('p') | Long("ploidy") => args.bg_params.depth.ploidy = parser.value()?.parse()?,
            Long("subsample") => args.params.subsampling_rate = parser.value()?.parse()?,
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
                args.bg_params.depth.kmer_perc = parser.value()?.parse()?,
            Long("frac-windows") | Long("fraction-windows") =>
                args.bg_params.depth.frac_windows = parser.value()?.parse()?,
            Long("blur-extreme") => {
                args.bg_params.depth.min_tail_obs = parser.value()?.parse()?;
                args.bg_params.depth.tail_var_mult = parser.value()?.parse()?;
            }

            Short('^') | Long("interleaved") => args.interleaved = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Long("rerun") => args.rerun = parser.value()?.parse()?,
            Long("debug") => args.debug = true,
            Long("strobealign") => args.strobealign = parser.value()?.parse()?,
            Long("minimap") | Long("minimap2") => args.minimap = parser.value()?.parse()?,
            Long("samtools") => args.samtools = parser.value()?.parse()?,
            Short('s') | Long("seed") => args.params.seed = Some(parser.value()?.parse()?),

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

fn create_out_dir(args: &Args) -> Result<PathBuf, Error> {
    let out_dir = args.output.as_ref().unwrap();
    ext::sys::mkdir(out_dir)?;

    let bg_dir = out_dir.join(paths::BG_DIR);
    let params_path = bg_dir.join(paths::PREPROC_PARAMS);
    if !args.rerun.need_analysis(&bg_dir)? {
        if args.params.need_rerun(!args.input.is_empty(), &params_path) {
            log::error!("Please rerun with {} or {}", "--rerun part".red(), "--rerun all".red());
            std::process::exit(1);
        }
        std::process::exit(0);
    }

    let mut params_file = ext::sys::create_file(&params_path)?;
    args.params.save().write_pretty(&mut params_file, 4).map_err(add_path!(params_path))?;
    Ok(bg_dir)
}

fn set_mapping_stdin(
    args: &Args,
    child_stdin: ChildStdin,
    rng: &mut XoshiroRng,
) -> Result<thread::JoinHandle<Result<(), io::Error>>, Error>
{
    let subsampling_rate = args.params.subsampling_rate;
    let mut local_rng = rng.clone();
    rng.long_jump();
    let mut writer = io::BufWriter::new(child_stdin);
    // Cannot put reader into a box, because `FastxRead` has a type parameter.
    if args.input.len() == 1 && !args.interleaved {
        let mut reader = fastx::Reader::from_path(&args.input[0])?;
        Ok(thread::spawn(move || reader.subsample(&mut writer, subsampling_rate, &mut local_rng)))
    } else if args.interleaved {
        let mut reader = fastx::PairedEndInterleaved::new(fastx::Reader::from_path(&args.input[0])?);
        Ok(thread::spawn(move || reader.subsample(&mut writer, subsampling_rate, &mut local_rng)))
    } else {
        let mut reader = fastx::PairedEndReaders::new(
            fastx::Reader::from_path(&args.input[0])?,
            fastx::Reader::from_path(&args.input[1])?);
        Ok(thread::spawn(move || reader.subsample(&mut writer, subsampling_rate, &mut local_rng)))
    }
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
            ]);
        if (args.params.subsampling_rate == 1.0 && args.interleaved) || args.is_paired_end() {
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
            "-t", &args.threads.to_string() // Specify the number of threads,
            ]);
    }

    // Providing paths to the reference and reads.
    command.arg(&ref_filename);
    if args.params.subsampling_rate == 1.0 {
        command.args(&args.input);
    } else {
        command.arg("-").stdin(Stdio::piped());
    }
    command.stdout(Stdio::piped());
    command
}

fn first_step_str(args: &Args) -> String {
    let mut s = String::new();
    if args.params.subsampling_rate == 1.0 {
        return s;
    }

    write!(s, "_subsample_ --rate {}", args.params.subsampling_rate).unwrap();
    if let Some(seed) = args.params.seed {
        write!(s, " --seed {}", seed).unwrap();
    }
    write!(s, " -i {}", ext::fmt::path(&args.input[0])).unwrap();
    if args.input.len() > 1 {
        write!(s, " {}", ext::fmt::path(&args.input[1])).unwrap();
    }
    write!(s, " | ").unwrap();
    s
}

/// Map reads to the whole reference genome, and then take only reads mapped to the corresponding BED file.
/// Subsample reads if the corresponding rate is less than 1.
fn run_mapping(
    args: &Args,
    seq_info: &SequencingInfo,
    ref_filename: &Path,
    bed_target: &Path,
    out_bam: &Path,
    rng: &mut XoshiroRng,
) -> Result<(), Error>
{
    let tmp_bam = out_bam.with_extension("tmp.bam");
    if out_bam.exists() {
        log::warn!("    BAM file {} exists, skipping read mapping", ext::fmt::path(out_bam));
        return Ok(());
    }

    let start = Instant::now();
    let mut mapping = create_mapping_command(args, seq_info, ref_filename);
    let mut mapping_child = mapping.spawn().map_err(add_path!(!))?;
    let mapping_stdin = mapping_child.stdin.take();
    let mapping_stdout = mapping_child.stdout.take();
    let mut guard = ext::sys::ChildGuard::new(mapping_child);
    let handle = mapping_stdin.map(|stdin| set_mapping_stdin(args, stdin, rng)).transpose()?;

    let mut samtools = Command::new(&args.samtools);
    // See SAM flags here: https://broadinstitute.github.io/picard/explain-flags.html.
    samtools.args(&["view",
            "-b", // Output BAM.
            // Ignore reads where any of the mates is unmapped,
            // + ignore secondary & supplementary alignments + ignore failed checks.
            "-F", "3852",
            "-q", &args.params.min_mapq.to_string(),
            ])
        .arg("-L").arg(&bed_target)
        .arg("-o").arg(&tmp_bam)
        .stdin(Stdio::from(mapping_stdout.unwrap()));
    log::debug!("    {}{} | {}", first_step_str(&args), ext::fmt::command(&mapping), ext::fmt::command(&samtools));
    let output = samtools.output().map_err(add_path!(!))?;
    log::debug!("");
    log::debug!("    Finished in {}", ext::fmt::Duration(start.elapsed()));
    if !output.status.success() {
        return Err(Error::Subprocess(output));
    }
    if let Some(handle) = handle {
        // handle.join() returns Result<Result<(), crate::Error>, Any>.
        // expect unwraps the outer Err, then ? returns the inner Err, if any.
        handle.join().expect("Process failed for unknown reason").map_err(add_path!(!))?;
    }
    fs::rename(&tmp_bam, out_bam).map_err(add_path!(tmp_bam, out_bam))?;
    guard.disarm();
    Ok(())
}

/// Loads records and converts them to Alignments.
/// All returned alignments are primary, have MAPQ over the `args.min_mapq`,
/// and clipping under the threshold `argsbg_params.max_clipping`.
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
) -> Result<(Vec<Alignment>, bool), Error>
{
    let min_mapq = args.params.min_mapq;
    let max_clipping = args.max_clipping;
    let mut paired_counts = [0_u64, 0];

    let mut alns = Vec::new();
    let mut record = bam::Record::new();
    let mut discarded = 0;
    while let Some(()) = reader.read(&mut record).transpose()? {
        if record.flags() & 3844 == 0 && record.mapq() >= min_mapq && cigar::clipping_rate(&record) <= max_clipping {
            if let Some(cigar) = cigar_getter(&record) {
                alns.push(Alignment::new(&record, cigar, ReadEnd::from_record(&record), Arc::clone(contigs), f64::NAN));
                paired_counts[usize::from(record.is_paired())] += 1;
            }
        } else {
            discarded += 1;
        }
    }
    if paired_counts[0] > 0 && paired_counts[1] > 0 {
        return Err(Error::InvalidData(format!("BAM file contains both paired and unpaired reads")));
    }
    if alns.is_empty() {
        return Err(Error::InvalidData(format!("BAM file contains no reads in the target region")));
    }
    log::debug!("    Loaded {} alignments, discarded {}", alns.len(), discarded);
    Ok((alns, paired_counts[1] > 0))
}

const MEAN_LEN_RECORDS: usize = 10000;

/// Calculate mean read length from existing alignments.
fn read_len_from_alns(alns: &[Alignment]) -> f64 {
    let n = min(alns.len(), MEAN_LEN_RECORDS);
    alns[..n].iter()
        .map(|aln| f64::from(aln.cigar().query_len()))
        .sum::<f64>() / n as f64
}

/// Calculate mean read length from input reads.
fn read_len_from_reads(input: &[PathBuf]) -> Result<f64, Error> {
    log::info!("Calculating mean read length");
    let lengths: Vec<f64> = input.iter()
        .map(|path| {
            let mut reader = fastx::Reader::from_path(path)?;
            reader.mean_read_len(MEAN_LEN_RECORDS / input.len()).map_err(add_path!(path))
        })
        .collect::<Result<_, Error>>()?;
    Ok(ext::vec::F64Ext::mean(&lengths))
}

fn estimate_bg_from_paired(
    alns: Vec<Alignment>,
    seq_info: SequencingInfo,
    args: &Args,
    opt_out_dir: Option<&Path>,
    data: &RefData,
) -> Result<BgDistr, Error>
{
    // Group reads into pairs, and estimate insert size from them.
    let pair_ixs = insertsz::group_mates(&alns)?;
    let insert_distr = InsertDistr::estimate(&alns, &pair_ixs, opt_out_dir)?;
    let conf_lvl = 1.0 - args.bg_params.insert_pval;
    let (min_insert, max_insert) = insert_distr.confidence_interval(conf_lvl);
    log::info!("    Allowed insert size: [{}, {}]  ({}%-confidence interval)",
        min_insert, max_insert, crate::math::fmt_signif(100.0 * conf_lvl, 5));

    let interval = &data.interval;
    let windows = bg::Windows::create(interval, &data.sequence, &data.kmer_counts, &seq_info,
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
    let err_prof = ErrorProfile::estimate(&errprof_alns, interval, &windows, seq_info.mean_read_len(),
        &args.bg_params, opt_out_dir)?;

    // Estimate backgorund read depth from read pairs with both good probabilities and good insert size prob.
    let mut depth_alns: Vec<&LightAlignment> = Vec::with_capacity(errprof_alns.len());
    for chunk in errprof_alns.chunks_exact(2) {
        let aln1 = chunk[0];
        let aln2 = chunk[1];
        let (edit1, len1) = aln1.count_region_operations(interval).edit_and_read_len();
        let (edit2, len2) = aln2.count_region_operations(interval).edit_and_read_len();
        if edit1 <= err_prof.allowed_edit_dist(len1).0 && edit2 <= err_prof.allowed_edit_dist(len2).0 {
            depth_alns.push(aln1.deref());
            depth_alns.push(aln2.deref());
        }
    }
    let depth_distr = ReadDepth::estimate(&depth_alns, &windows,
        &args.bg_params.depth, args.params.subsampling_rate, true, &seq_info, opt_out_dir)?;
    Ok(BgDistr::new(seq_info, insert_distr, err_prof, depth_distr))
}

fn estimate_bg_from_unpaired(
    alns: Vec<Alignment>,
    seq_info: SequencingInfo,
    args: &Args,
    opt_out_dir: Option<&Path>,
    data: &RefData,
) -> Result<BgDistr, Error>
{
    let insert_distr = InsertDistr::undefined();
    let interval = &data.interval;
    let windows = bg::Windows::create(interval, &data.sequence, &data.kmer_counts, &seq_info,
        &args.bg_params.depth, opt_out_dir)?;
    let err_prof = ErrorProfile::estimate(&alns, interval, &windows, seq_info.mean_read_len(),
        &args.bg_params, opt_out_dir)?;
    let filt_alns: Vec<&LightAlignment> = alns.iter()
        .filter(|aln| {
            let (edit, len) = aln.count_region_operations(interval).edit_and_read_len();
            edit <= err_prof.allowed_edit_dist(len).0
        })
        .map(Deref::deref)
        .collect();
    let depth_distr = ReadDepth::estimate(&filt_alns, &windows,
        &args.bg_params.depth, args.params.subsampling_rate, true, &seq_info, opt_out_dir)?;
    Ok(BgDistr::new(seq_info, insert_distr, err_prof, depth_distr))
}

/// Estimate background distributions from input reads or existing alignments.
fn estimate_bg_distrs(
    args: &Args,
    out_dir: &Path,
    data: &RefData,
    rng: &mut XoshiroRng,
) -> Result<BgDistr, Error>
{
    let opt_out_dir = if args.debug { Some(out_dir) } else { None };
    let is_paired_end: bool;
    let alns: Vec<Alignment>;
    let seq_info: SequencingInfo;

    if let Some(alns_filename) = args.alns.as_ref() {
        log::debug!("Loading mapped reads into memory ({})", ext::fmt::path(alns_filename));
        let mut bam_reader = bam::IndexedReader::from_path(alns_filename)?;
        bam_reader.set_reference(&data.ref_filename)?;
        let interval = &data.interval;
        let interval_start = interval.start();
        let ref_seq = &data.sequence;
        bam_reader.fetch((interval.contig_name(), i64::from(interval_start), i64::from(interval.end())))?;
        (alns, is_paired_end) = load_alns(&mut bam_reader,
            |record| Cigar::infer_ext_cigar(record, ref_seq, interval_start),
            &data.ref_contigs, args)?;
        if is_paired_end && !args.params.technology.paired_end_allowed() {
            return Err(Error::InvalidInput(format!("Paired end reads are not supported by {}",
                args.params.technology.long_name())));
        }
        seq_info = SequencingInfo::new(read_len_from_alns(&alns), args.params.technology, args.explicit_technology)?;
    } else {
        seq_info = SequencingInfo::new(read_len_from_reads(&args.input)?, args.params.technology,
            args.explicit_technology)?;
        log::info!("Mean read length = {:.1}", seq_info.mean_read_len());

        let bam_filename = out_dir.join("aln.bam");
        log::info!("Mapping reads to the reference");
        run_mapping(args, &seq_info, &data.ref_filename, &data.bed_filename, &bam_filename, rng)?;

        log::debug!("Loading mapped reads into memory ({})", ext::fmt::path(&bam_filename));
        let mut bam_reader = bam::Reader::from_path(&bam_filename)?;
        (alns, is_paired_end) = load_alns(&mut bam_reader, |record| Some(Cigar::from_raw(record)),
            &data.ref_contigs, args)?;
        assert!(is_paired_end == args.is_paired_end());
    }

    if is_paired_end {
        estimate_bg_from_paired(alns, seq_info, args, opt_out_dir, data)
    } else {
        estimate_bg_from_unpaired(alns, seq_info, args, opt_out_dir, data)
    }
}

struct RefData {
    ref_filename: PathBuf,
    ref_contigs: Arc<ContigNames>,
    bed_filename: PathBuf,
    interval: Interval,
    sequence: Vec<u8>,
    kmer_counts: KmerCounts,
}

impl RefData {
    fn new(ref_filename: PathBuf, db_path: &Path) -> Result<Self, Error> {
        let (ref_contigs, mut ref_fasta) = ContigNames::load_indexed_fasta("ref", &ref_filename)?;

        // Load and parse background region coordinates.
        let bed_filename = db_path.join(paths::BG_BED);
        let bed_contents = fs::read_to_string(&bed_filename).map_err(add_path!(bed_filename))?;
        let bg_intervals: Vec<_> = bed_contents.split('\n')
            .filter(|s| !s.starts_with('#') && s.contains('\t')).collect();
        if bg_intervals.len() != 1 {
            return Err(Error::InvalidData(format!("Incorrect number of regions in {}", ext::fmt::path(&bed_filename))));
        }
        let interval = Interval::parse_bed(&mut bg_intervals[0].split('\t'), &ref_contigs)?;
        let sequence = interval.fetch_seq(&mut ref_fasta).map_err(add_path!(ref_filename))?;

        // Load k-mer counts for the background region.
        let kmers_filename = db_path.join(paths::KMERS);
        let kmer_counts = KmerCounts::load(ext::sys::open(&kmers_filename)?, &[interval.len()])?;
        Ok(Self { ref_filename, ref_contigs, bed_filename, interval, sequence, kmer_counts })
    }
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();
    let mut rng = init_rng(args.params.seed);
    // `args.output`/bg
    let out_bg_dir = create_out_dir(&args)?;

    log::info!("Loading background non-duplicated region into memory");
    let db_path = args.database.as_ref().unwrap().join(paths::BG_DIR);
    let ref_data = RefData::new(args.reference.clone().unwrap(), &db_path)?;

    let bg_distr = estimate_bg_distrs(&args, &out_bg_dir, &ref_data, &mut rng)?;
    let distr_filename = out_bg_dir.join(paths::BG_DISTR);
    let mut distr_file = ext::sys::create_gzip(&distr_filename)?;
    bg_distr.save().write_pretty(&mut distr_file, 4).map_err(add_path!(distr_filename))?;
    log::info!("Success. Total time: {}", ext::fmt::Duration(timer.elapsed()));
    super::write_success_file(out_bg_dir.join(paths::SUCCESS))?;
    Ok(())
}
