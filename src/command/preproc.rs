//! Preprocess WGS dataset.

use std::{
    fs, io, thread,
    fmt::Write as FmtWrite,
    cmp::max,
    path::{Path, PathBuf},
    process::{Stdio, Command, ChildStdin},
    time::Instant,
    ops::Deref,
    sync::Arc,
};
use colored::Colorize;
use const_format::{str_repeat, concatcp};
use htslib::bam;
use crate::{
    err::{Error, validate_param},
    math::Ln,
    ext,
    seq::{
        ContigNames, Interval,
        kmers::KmerCounts,
        fastx::{self, FastxRead},
        cigar::{self, Cigar},
        aln::{Alignment, LightAlignment, ReadEnd},
    },
    bg::{
        self, BgDistr,
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
    min_mapq: u8,
    subsampling_rate: f64,
    subsampling_seed: Option<u64>,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            min_mapq: 20,
            subsampling_rate: 0.1,
            subsampling_seed: None,
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
    fn need_rerun(&self, path: &Path) -> bool {
        if !path.exists() {
            log::warn!("Cannot find old parameters at {}", ext::fmt::path(path));
            return true;
        }
        let old_params = match ext::sys::load_json(&path)
                .map_err(|e| Error::from(e))
                .and_then(|json| Params::load(&json)) {
            Err(_) => {
                log::warn!("Cannot load old parameters from {}", ext::fmt::path(path));
                return true;
            }
            Ok(val) => val,
        };
        if self.min_mapq != old_params.min_mapq || self.subsampling_rate != old_params.subsampling_rate {
            log::warn!("Parameters has changed");
            true
        } else {
            if self.subsampling_seed != old_params.subsampling_seed {
                log::warn!("Subsampling seed has changed");
            }
            false
        }
    }
}

impl JsonSer for Params {
    fn save(&self) -> json::JsonValue {
        json::object!{
            min_mapq: self.min_mapq,
            subsampling_rate: self.subsampling_rate,
            subsampling_seed: self.subsampling_seed,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, Error> {
        json_get!(obj -> min_mapq (as_u8), subsampling_rate (as_f64), subsampling_seed? (as_u64));
        Ok(Self { min_mapq, subsampling_rate, subsampling_seed })
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
    force: bool,
    strobealign: PathBuf,
    samtools: PathBuf,
    debug: bool,

    /// When calculating insert size distributions and read error profiles,
    /// ignore reads with `clipping > max_clipping * read_len`.
    pub max_clipping: f64,
    /// Do not use reads with alignment probability < `min_aln_prob` (ln-space) for read depth calculation.
    min_aln_prob: f64,

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
            force: false,
            debug: false,
            strobealign: PathBuf::from("strobealign"),
            samtools: PathBuf::from("samtools"),

            max_clipping: 0.02,
            min_aln_prob: Ln::from_log10(-3.0),

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
        if !self.interleaved && n_input == 1 {
            log::warn!("Running in single-end mode.");
        }

        validate_param!(self.database.is_some(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.reference.is_some(), "Reference fasta file is not provided (see -r/--reference)");
        validate_param!(self.output.is_some(), "Output directory is not provided (see -o/--output)");

        self.strobealign = ext::sys::find_exe(self.strobealign)?;
        self.samtools = ext::sys::find_exe(self.samtools)?;

        validate_param!(0.0 <= self.max_clipping && self.max_clipping <= 1.0,
            "Max clipping ({:.5}) must be within [0, 1]", self.max_clipping);
        validate_param!(self.min_aln_prob < 0.0,
            "Minimum alignment probability ({:.5}) must be under 0.", Ln::to_log10(self.min_aln_prob));
        if self.alns.is_some() {
            self.params.subsampling_rate = 1.0;
            self.params.subsampling_seed = None;
        }

        self.params.validate()?;
        self.bg_params.validate()?;
        Ok(self)
    }

    fn is_paired_end(&self) -> bool {
        self.input.len() == 2 || self.interleaved
    }
}

fn print_help(extended: bool) {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Preprocess WGS dataset.".yellow());

    println!("\n{} {} preproc (-i reads1.fq [reads2.fq] | -a reads.bam) -d db -r reference.fa -o out [arguments]",
        "Usage:".bold(), super::PKG_NAME);
    if !extended {
        println!("\nThis is a short help message. Please use {} to see the full help.",
            "-H/--full-help".green());
    }

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Reads in indexed BAM/CRAM format, already mapped to the whole genome.\n\
        {EMPTY}  Mutually exclusive with {}.",
        "-a, --alignment".green(), "FILE".yellow(), "-i/--input".green());
    println!("    {:KEY$} {:VAL$}  Database directory (initialized with {} & {}).",
        "-d, --db".green(), "DIR".yellow(), concatcp!(super::PKG_NAME, " create").underline(), "add".underline());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file. Must contain FAI index.",
        "-r, --reference".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Output directory.",
        "-o, --output".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Input reads are interleaved.",
        "-^, --interleaved".green(), super::flag());

    if extended {
        println!("\n{}", "Insert size and error profile estimation:".bold());
        println!("    {:KEY$} {:VAL$}  Ignore reads with mapping quality less than {} [{}].",
            "-q, --min-mapq".green(), "INT".yellow(), "INT".yellow(), defaults.params.min_mapq);
        println!("    {:KEY$} {:VAL$}  Ignore reads with soft/hard clipping over {} * read length [{}].",
            "-c, --max-clipping".green(), "FLOAT".yellow(), "FLOAT".yellow(), defaults.max_clipping);
        println!("    {:KEY$} {}\n\
            {EMPTY}  Max allowed insert size is calculated as {} multiplied by\n\
            {EMPTY}  the {}-th insert size quantile [{}, {}].",
            "-I, --ins-quantile".green(), "FLOAT FLOAT".yellow(),
            "FLOAT_1".yellow(), "FLOAT_2".yellow(),
            defaults.bg_params.ins_quantile_mult, defaults.bg_params.ins_quantile);
        println!("    {:KEY$} {:VAL$}  Multiply error rates by this factor, in order to correct for\n\
            {EMPTY}  read mappings missed due to higher error rate [{}].",
            "-m, --err-mult".green(), "FLOAT".yellow(), defaults.bg_params.err_rate_mult);

        println!("\n{}", "Background read depth estimation:".bold());
        println!("    {:KEY$} {:VAL$}  Specie ploidy [{}].",
            "-p, --ploidy".green(), "INT".yellow(), defaults.bg_params.depth.ploidy);
        println!("    {:KEY$} {:VAL$}\n\
            {EMPTY}  Subsample input reads by a factor of {} [{}]\n\
            {EMPTY}  Use all reads for {} or if alignment file ({}) is provided.\n\
            {EMPTY}  Second value sets the subsampling seed (optional).",
            "-s, --subsample".green(), "FLOAT [INT]".yellow(), "FLOAT".yellow(), defaults.params.subsampling_rate,
            "-s 1".green(), "-a".green());
        println!("    {:KEY$} {:VAL$}  Count read depth per {} bp windows [{}].",
            "-w, --window".green(), "INT".yellow(), "INT".yellow(), defaults.bg_params.depth.window_size);
        println!("    {:KEY$} {:VAL$}  Extend window by {} bp to both sides in order to calculate\n\
            {EMPTY}  GC-content and average k-mer frequency [{}].",
            "    --window-padd".green(), "INT".yellow(), "INT".yellow(), defaults.bg_params.depth.window_padding);
        println!("    {:KEY$} {:VAL$}  Skip {} bp near the edge of the background region.\n\
            {EMPTY}  Must not be smaller than {} [{}].",
            "    --boundary".green(), "INT".yellow(), "INT".yellow(), "--window-padd".green(),
            defaults.bg_params.depth.boundary_size);
        println!("    {:KEY$} {:VAL$}  Single-end alignment likelihood threshold in log10-space [{:.1}].",
            "    --min-aln-lik".green(), "FLOAT".yellow(), Ln::to_log10(defaults.min_aln_prob));
        println!("    {:KEY$} {:VAL$}  Ignore windows with average k-mer frequency over {} [{}].",
            "    --kmer-freq".green(), "FLOAT".yellow(), "FLOAT".yellow(), defaults.bg_params.depth.max_kmer_freq);
        println!("    {:KEY$} {:VAL$}  This fraction of all windows is used to estimate read depth for\n\
            {EMPTY}  each GC-content [{}]. Smaller values lead to less robust estimates,\n\
            {EMPTY}  larger values - to similar estimates across different GC-contents.",
            "    --frac-windows".green(), "FLOAT".yellow(), defaults.bg_params.depth.frac_windows);
        println!("    {:KEY$} {}\n\
            {EMPTY}  Read depth estimates are blured for windows with extreme GC-content\n\
            {EMPTY}  (less than {} windows with smaller/larger GC). There, read depth\n\
            {EMPTY}  is set to the last non-extreme depth, while variance is increased\n\
            {EMPTY}  by a {} factor for each addition GC value [{} {}].",
            "    --blur-extreme".green(), "INT FLOAT".yellow(), "INT".yellow(), "FLOAT".yellow(),
            defaults.bg_params.depth.min_tail_obs, defaults.bg_params.depth.tail_var_mult);
    }

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads);
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Create more files with debug information.",
        "    --debug".green(), super::flag());
    if extended {
        println!("    {:KEY$} {:VAL$}  Strobealign executable [{}].",
            "    --strobealign".green(), "EXE".yellow(), defaults.strobealign.display());
        println!("    {:KEY$} {:VAL$}  Samtools executable [{}].",
            "    --samtools".green(), "EXE".yellow(), defaults.samtools.display());
    }

    println!("\n{}", "Other parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Show short help message.", "-h, --help".green(), "");
    println!("    {:KEY$} {:VAL$}  Show extended help message.", "-H, --full-help".green(), "");
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
            Short('a') | Long("aln") | Long("alns") | Long("alignment") | Long("alignments") =>
                args.alns = Some(parser.value()?.parse()?),
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),

            Short('q') | Long("min-mapq") | Long("min-mq") => args.params.min_mapq = parser.value()?.parse()?,
            Short('c') | Long("max-clip") | Long("max-clipping") => args.max_clipping = parser.value()?.parse()?,
            Short('I') | Long("ins-quant") | Long("ins-quantile") => {
                args.bg_params.ins_quantile_mult = parser.value()?.parse()?;
                args.bg_params.ins_quantile = parser.value()?.parse()?;
            }
            Short('m') | Long("err-mult") | Long("err-multiplier") =>
                args.bg_params.err_rate_mult = parser.value()?.parse()?,

            Short('p') | Long("ploidy") => args.bg_params.depth.ploidy = parser.value()?.parse()?,
            Short('s') | Long("subsample") => {
                let mut values = parser.values()?;
                args.params.subsampling_rate = values.next()
                    .ok_or_else(|| lexopt::Error::MissingValue { option: Some("subsample".to_owned()) })?
                    .parse()?;
                args.params.subsampling_seed = values.next().map(|v| v.parse()).transpose()?;
            }
            Short('w') | Long("window") => args.bg_params.depth.window_size = parser.value()?.parse()?,
            Long("window-padd") | Long("window-padding") => args.bg_params.depth.window_padding = parser.value()?.parse()?,
            Long("boundary") => args.bg_params.depth.boundary_size = parser.value()?.parse()?,
            Long("min-aln-lik") => args.min_aln_prob = Ln::from_log10(parser.value()?.parse()?),
            Long("kmer-freq") | Long("kmer-frequency") => args.bg_params.depth.max_kmer_freq = parser.value()?.parse()?,
            Long("frac-windows") | Long("fraction-windows") =>
                args.bg_params.depth.frac_windows = parser.value()?.parse()?,
            Long("blur-extreme") => {
                args.bg_params.depth.min_tail_obs = parser.value()?.parse()?;
                args.bg_params.depth.tail_var_mult = parser.value()?.parse()?;
            }

            Short('^') | Long("interleaved") => args.interleaved = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,
            Long("debug") => args.debug = true,
            Long("strobealign") => args.strobealign = parser.value()?.parse()?,
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

fn create_out_dir(args: &Args) -> Result<PathBuf, Error> {
    let out_dir = args.output.as_ref().unwrap();
    ext::sys::mkdir(out_dir)?;

    let bg_dir = out_dir.join(paths::BG_DIR);
    let params_path = bg_dir.join(paths::PREPROC_PARAMS);
    if bg_dir.exists() {
        if args.force || !bg_dir.join(paths::BG_DISTR).exists() {
            log::warn!("Clearing output directory {}", ext::fmt::path(&bg_dir));
            fs::remove_dir_all(&bg_dir)?;
        } else if args.alns.is_none() && args.params.need_rerun(&params_path) {
            log::error!("Please rerun with -F/--force");
            std::process::exit(1);
        }
    }
    ext::sys::mkdir(&bg_dir)?;

    let mut params_file = io::BufWriter::new(fs::File::create(&params_path)?);
    args.params.save().write_pretty(&mut params_file, 4)?;
    Ok(bg_dir)
}

fn set_strobealign_stdin(
    args: &Args,
    child_stdin: ChildStdin,
) -> Result<thread::JoinHandle<Result<(), io::Error>>, Error>
{
    fn create_job(args: &Args, mut writer: impl io::Write + 'static, mut reader: impl FastxRead + 'static,
    ) -> impl FnOnce() -> Result<(), io::Error> + 'static {
        let subsampling_rate = args.params.subsampling_rate;
        let subsampling_seed = args.params.subsampling_seed.clone();
        move || reader.subsample(&mut writer, subsampling_rate, &mut ext::rand::init_rng(subsampling_seed))
    }

    let writer = io::BufWriter::new(child_stdin);
    // Cannot put reader into a box, because `FastxRead` has a type parameter.
    if args.input.len() == 1 && !args.interleaved {
        let reader = fastx::Reader::from_path(&args.input[0])?;
        Ok(thread::spawn(create_job(args, writer, reader)))
    } else if args.interleaved {
        let reader = fastx::PairedEndInterleaved::new(fastx::Reader::from_path(&args.input[0])?);
        Ok(thread::spawn(create_job(args, writer, reader)))
    } else {
        let reader = fastx::PairedEndReaders::new(
            fastx::Reader::from_path(&args.input[0])?,
            fastx::Reader::from_path(&args.input[1])?);
        Ok(thread::spawn(create_job(args, writer, reader)))
    }
}

fn first_step_str(args: &Args) -> String {
    let mut s = String::new();
    if args.params.subsampling_rate == 1.0 {
        return s;
    }

    write!(s, "_subsample_ --rate {}", args.params.subsampling_rate).unwrap();
    if let Some(seed) = args.params.subsampling_seed {
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
fn run_strobealign(args: &Args, ref_filename: &Path, bed_target: &Path, out_bam: &Path) -> Result<(), Error> {
    let tmp_bam = out_bam.with_extension("tmp.bam");
    if out_bam.exists() {
        log::warn!("    BAM file {} exists, skipping read mapping", ext::fmt::path(out_bam));
        return Ok(());
    }

    let start = Instant::now();
    let mut strobealign = Command::new(&args.strobealign);
    strobealign.args(&[
        "-N0", "-R0", "-U", // Retain 0 additional alignments, do not rescue reads, do not output unmapped reads.
        "-f", "0.001", // Discard more minimizers to speed up alignment
        "-t", &args.threads.to_string()]) // Specify the number of threads.
        .stdout(Stdio::piped());
    if args.params.subsampling_rate == 1.0 {
        // Provide input files as they are.
        if args.interleaved {
            strobealign.arg("--interleaved");
        }
        strobealign.arg(&ref_filename).args(&args.input);
    } else {
        if args.is_paired_end() {
            strobealign.arg("--interleaved");
        }
        // Subsample or take first N records from the input files, and pass them through stdin.
        strobealign.arg(&ref_filename).arg("-").stdin(Stdio::piped());
    }
    let mut strobealign_child = strobealign.spawn()?;
    let strobealign_stdin = strobealign_child.stdin.take();
    let strobealign_stdout = strobealign_child.stdout.take();
    let mut guard = ext::sys::ChildGuard::new(strobealign_child);
    let handle = strobealign_stdin.map(|stdin| set_strobealign_stdin(args, stdin)).transpose()?;

    let mut samtools = Command::new(&args.samtools);
    // See SAM flags here: https://broadinstitute.github.io/picard/explain-flags.html.
    samtools.args(&["view",
            "-b", // Output BAM.
            // Ignore reads where any of the mates is unmapped,
            // + ignore secondary & supplementary alignments + ignore failed checks.
            "-F3852",
            "-q", &args.params.min_mapq.to_string(),
            ])
        .arg("-L").arg(&bed_target)
        .arg("-o").arg(&tmp_bam)
        .stdin(Stdio::from(strobealign_stdout.unwrap()));
    log::debug!("    {}{} | {}", first_step_str(&args), ext::fmt::command(&strobealign), ext::fmt::command(&samtools));
    let output = samtools.output()?;
    log::debug!("");
    log::debug!("    Finished in {}", ext::fmt::Duration(start.elapsed()));
    if !output.status.success() {
        return Err(Error::SubprocessFail(output));
    }
    if let Some(handle) = handle {
        // handle.join() returns Result<Result<(), crate::Error>, Any>.
        // expect unwraps the outer Err, then ? returns the inner Err, if any.
        handle.join().expect("Process failed for unknown reason")?;
    }
    fs::rename(&tmp_bam, out_bam)?;
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
    while let Some(()) = reader.read(&mut record).transpose()? {
        if record.flags() & 3844 == 0 && record.mapq() >= min_mapq && cigar::clipping_rate(&record) <= max_clipping {
            if let Some(cigar) = cigar_getter(&record) {
                alns.push(Alignment::new(&record, cigar, ReadEnd::from_record(&record), Arc::clone(contigs), f64::NAN));
                paired_counts[usize::from(record.is_paired())] += 1;
            }
        }
    }
    if paired_counts[0] > 0 && paired_counts[1] > 0 {
        return Err(Error::InvalidData(format!("BAM file contains both paired and unpaired reads")));
    }
    if alns.is_empty() {
        return Err(Error::InvalidData(format!("BAM file contains no reads in the target region")));
    }
    Ok((alns, paired_counts[1] > 0))
}

fn estimate_bg_from_paired(
    alns: Vec<Alignment>,
    args: &Args,
    opt_out_dir: Option<&Path>,
    data: &RefData,
) -> Result<BgDistr, Error>
{
    // Group reads into pairs, and estimate insert size from them.
    let pair_ixs = insertsz::group_mates(&alns)?;
    let insert_distr = InsertDistr::estimate(&alns, &pair_ixs, &args.bg_params, opt_out_dir)?;

    // Estimate error profile from read pairs with appropriate insert size.
    let max_insert_size = insert_distr.max_size();
    let mut errprof_alns = Vec::with_capacity(pair_ixs.len() * 2);
    for &(i, j) in pair_ixs.iter() {
        let first = &alns[i];
        let second = &alns[j];
        if first.insert_size(second.deref()) <= max_insert_size {
            errprof_alns.push(first);
            errprof_alns.push(second);
        }
    }
    let err_prof = ErrorProfile::estimate(errprof_alns.iter().map(|aln| aln.cigar()), args.bg_params.err_rate_mult);

    // Estimate backgorund read depth from read pairs with both good probabilities and good insert size prob.
    let mut depth_alns: Vec<&LightAlignment> = Vec::with_capacity(errprof_alns.len());
    let thresh_prob = 2.0 * args.min_aln_prob + insert_distr.mode_prob();
    log::info!("    Using paired-end alignment likelihood threshold {:.2}", Ln::to_log10(thresh_prob));
    for chunk in errprof_alns.chunks_exact(2) {
        let first = chunk[0];
        let second = chunk[1];
        let prob = err_prof.ln_prob(first.cigar()) + err_prof.ln_prob(second.cigar())
            + first.insert_size_prob(second.deref(), &insert_distr);
        if prob >= thresh_prob {
            depth_alns.push(first.deref());
            depth_alns.push(second.deref());
        }
    }
    let depth_distr = ReadDepth::estimate(depth_alns.into_iter(), &data.interval, &data.sequence, &data.kmer_counts,
        &args.bg_params.depth, args.params.subsampling_rate, true, opt_out_dir)?;
    Ok(BgDistr::new(insert_distr, err_prof, depth_distr))
}

fn estimate_bg_from_unpaired(
    alns: Vec<Alignment>,
    args: &Args,
    opt_out_dir: Option<&Path>,
    data: &RefData,
) -> Result<BgDistr, Error>
{
    let insert_distr = InsertDistr::undefined();
    let err_prof = ErrorProfile::estimate(alns.iter().map(|aln| aln.cigar()), args.bg_params.err_rate_mult);
    let filt_alns = alns.iter()
        .filter(|aln| err_prof.ln_prob(aln.cigar()) >= args.min_aln_prob)
        .map(Deref::deref);
    let depth_distr = ReadDepth::estimate(filt_alns, &data.interval, &data.sequence, &data.kmer_counts,
        &args.bg_params.depth, args.params.subsampling_rate, false, opt_out_dir)?;
    Ok(BgDistr::new(insert_distr, err_prof, depth_distr))
}

/// Map reads to the genome and to the background genome to estimate background distributions based solely on WGS reads.
fn estimate_bg_distrs(
    args: &Args,
    out_dir: &Path,
    data: &RefData,
) -> Result<BgDistr, Error>
{
    let opt_out_dir = if args.debug { Some(out_dir) } else { None };
    let (alns, is_paired_end) = if let Some(alns_filename) = args.alns.as_ref() {
        log::debug!("Loading mapped reads into memory ({})", ext::fmt::path(alns_filename));
        let mut bam_reader = bam::IndexedReader::from_path(alns_filename)?;
        bam_reader.set_reference(&data.ref_filename)?;
        let interval = &data.interval;
        let interval_start = interval.start();
        let ref_seq = &data.sequence;
        bam_reader.fetch((interval.contig_name(), i64::from(interval_start), i64::from(interval.end())))?;
        load_alns(&mut bam_reader, |record| Cigar::infer_ext_cigar(record, ref_seq, interval_start),
            &data.ref_contigs, args)?
    } else {
        let bam_filename = out_dir.join("aln.bam");
        log::info!("Mapping reads to the reference");
        run_strobealign(args, &data.ref_filename, &data.bed_filename, &bam_filename)?;

        log::debug!("Loading mapped reads into memory ({})", ext::fmt::path(&bam_filename));
        let mut bam_reader = bam::Reader::from_path(&bam_filename)?;
        load_alns(&mut bam_reader, |record| Some(Cigar::from_raw(record)), &data.ref_contigs, args)?
    };
    assert!(args.alns.is_some() || (is_paired_end == args.is_paired_end()));

    if is_paired_end {
        estimate_bg_from_paired(alns, args, opt_out_dir, data)
    } else {
        estimate_bg_from_unpaired(alns, args, opt_out_dir, data)
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
        let bed_contents = fs::read_to_string(&bed_filename)?;
        let bg_intervals: Vec<_> = bed_contents.split('\n')
            .filter(|s| !s.starts_with('#') && s.contains('\t')).collect();
        if bg_intervals.len() != 1 {
            return Err(Error::InvalidData(format!("Incorrect number of regions in {}", ext::fmt::path(&bed_filename))));
        }
        let interval = Interval::parse_bed(&mut bg_intervals[0].split('\t'), &ref_contigs)?;
        let sequence = interval.fetch_seq(&mut ref_fasta)?;

        // Load k-mer counts for the background region.
        let kmers_filename = db_path.join(paths::KMERS);
        let kmer_counts = KmerCounts::load(ext::sys::open(&kmers_filename)?, &[interval.len()])?;
        Ok(Self { ref_filename, ref_contigs, bed_filename, interval, sequence, kmer_counts })
    }
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let timer = Instant::now();
    let args = parse_args(argv)?.validate()?;
    let out_dir = create_out_dir(&args)?;

    log::info!("Loading background non-duplicated region into memory");
    let db_path = args.database.as_ref().unwrap().join(paths::BG_DIR);
    let ref_data = RefData::new(args.reference.clone().unwrap(), &db_path)?;

    let bg_distr = estimate_bg_distrs(&args, &out_dir, &ref_data)?;
    let mut distr_file = ext::sys::create_gzip(&out_dir.join(paths::BG_DISTR))?;
    bg_distr.save().write_pretty(&mut distr_file, 4)?;
    log::info!("Success. Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
