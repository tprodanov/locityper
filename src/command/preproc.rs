//! Preprocess WGS dataset.

use std::{
    fs, io,
    fmt::Write as FmtWrite,
    cmp::max,
    path::{Path, PathBuf},
    process::{Stdio, Command, ChildStdin},
    time::Instant,
    rc::Rc,
    thread,
};
use colored::Colorize;
use const_format::{str_repeat, concatcp};
use htslib::bam::{
    self,
    Read as BamRead,
    record::Record as BamRecord,
};
use crate::{
    err::{Error, validate_param},
    math::Ln,
    algo::parse_int,
    ext,
    seq::{
        cigar, ContigNames, ContigSet, ContigId, Interval,
        fastx::{self, FastxRead},
    },
    bg::{
        BgDistr, JsonSer, Params as BgParams,
        insertsz::{ReadMateGrouping, InsertDistr},
        err_prof::ErrorProfile,
        depth::ReadDepth,
    },
};
use super::paths;

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

    n_reads: usize,
    subsampling_rate: f64,
    subsampling_seed: Option<u64>,
    params: BgParams,
}

const DEF_N_READS: &'static str = "500k";

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::new(),
            alns: None,
            database: None,
            reference: None,
            output: None,

            interleaved: false,
            threads: 4,
            force: false,
            strobealign: PathBuf::from("strobealign"),
            samtools: PathBuf::from("samtools"),

            n_reads: parse_int(DEF_N_READS).unwrap(),
            subsampling_rate: 0.25,
            subsampling_seed: None,
            params: BgParams::default(),
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
        validate_param!(0.0 < self.subsampling_rate && self.subsampling_rate <= 1.0,
            "Subsample rate ({}) must be within (0, 1]", self.subsampling_rate);
        if self.subsampling_rate > 0.99 {
            self.subsampling_rate = 1.0;
        }

        self.params.validate()?;
        self.strobealign = ext::sys::find_exe(self.strobealign)?;
        self.samtools = ext::sys::find_exe(self.samtools)?;
        Ok(self)
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
        println!("    {:KEY$} {:VAL$}  Use first {} reads (or read pairs)\n\
            {EMPTY} to estimate profiles [{}].",
            "-n, --n-reads".green(), "INT".yellow(), "INT".yellow(), DEF_N_READS);
        println!("    {:KEY$} {:VAL$}  Ignore reads with mapping quality less than {} [{}].",
            "-q, --min-mapq".green(), "INT".yellow(), "INT".yellow(), defaults.params.min_mapq);
        println!("    {:width$}  Please rerun with {}, if {} or {} values have changed.",
            "WARNING".red(), "--force".red(), "-n".green(), "-q".green(), width = KEY + VAL + 1);
        println!("    {:KEY$} {:VAL$}  Ignore reads with soft/hard clipping over {} * read length [{}].",
            "-c, --max-clipping".green(), "FLOAT".yellow(), "FLOAT".yellow(), defaults.params.max_clipping);
        println!("    {:KEY$} {}\n\
            {EMPTY}  Max allowed insert size is calculated as {} multiplied by\n\
            {EMPTY}  the {}-th insert size quantile [{}, {}].",
            "-I, --ins-quantile".green(), "FLOAT FLOAT".yellow(),
            "FLOAT_1".yellow(), "FLOAT_2".yellow(), defaults.params.ins_quantile_mult, defaults.params.ins_quantile);
        println!("    {:KEY$} {:VAL$}  Multiply error rates by this factor, in order to correct for\n\
            {EMPTY}  read mappings missed due to higher error rate [{}].",
            "-m, --err-mult".green(), "FLOAT".yellow(), defaults.params.err_rate_mult);

        println!("\n{}", "Background read depth estimation:".bold());
        println!("    {:KEY$} {:VAL$}  Specie ploidy [{}].",
            "-p, --ploidy".green(), "INT".yellow(), defaults.params.depth.ploidy);
        println!("    {:KEY$} {:VAL$}\n\
            {EMPTY}  Subsample input reads by a factor of {} [{}]\n\
            {EMPTY}  Use all reads for {} or if alignment file ({}) is provided.\n\
            {EMPTY}  Second value sets the subsampling seed (optional).",
            "-s, --subsample".green(), "FLOAT [INT]".yellow(), "FLOAT".yellow(), defaults.subsampling_rate,
            "-s 1".green(), "-a".green());
        println!("    {:KEY$} {:VAL$}  Count read depth per {} bp windows [{}].",
            "-w, --window".green(), "INT".yellow(), "INT".yellow(), defaults.params.depth.window_size);
        println!("    {:KEY$} {:VAL$}  Extend window by {} bp to both sides in order to calculate\n\
            {EMPTY}  GC-content and average k-mer frequency [{}].",
            "    --window-padd".green(), "INT".yellow(), "INT".yellow(), defaults.params.depth.window_padding);
        println!("    {:KEY$} {:VAL$}  Skip {} bp near the edge of the background region.\n\
            {EMPTY}  Must not be smaller than {} [{}].",
            "    --boundary".green(), "INT".yellow(), "INT".yellow(), "--window-padd".green(),
            defaults.params.depth.boundary_size);
        println!("    {:KEY$} {:VAL$}  Ignore reads with alignment likelihood under 10^{} [{:.1}]",
            "    --min-aln-lik".green(), "FLOAT".yellow(), "FLOAT".yellow(),
            Ln::to_log10(defaults.params.depth.min_aln_prob));
        println!("    {:KEY$} {:VAL$}  Ignore windows with average k-mer frequency over {} [{}].",
            "    --kmer-freq".green(), "FLOAT".yellow(), "FLOAT".yellow(), defaults.params.depth.max_kmer_freq);
        println!("    {:KEY$} {:VAL$}  This fraction of all windows is used to estimate read depth for\n\
            {EMPTY}  each GC-content [{}]. Smaller values lead to less robust estimates,\n\
            {EMPTY}  larger values - to similar estimates across different GC-contents.",
            "    --frac-windows".green(), "FLOAT".yellow(), defaults.params.depth.frac_windows);
        println!("    {:KEY$} {}\n\
            {EMPTY}  Read depth estimates are blured for windows with extreme GC-content\n\
            {EMPTY}  (less than {} windows with smaller/larger GC). There, read depth\n\
            {EMPTY}  is set to the last non-extreme depth, while variance is increased\n\
            {EMPTY}  by a {} factor for each addition GC value [{} {}].",
            "    --blur-extreme".green(), "INT FLOAT".yellow(), "INT".yellow(), "FLOAT".yellow(),
            defaults.params.depth.min_tail_obs, defaults.params.depth.tail_var_mult);
    }

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads);
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());
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
                args.input = parser.values()?.take(2).map(|s| s.parse()).collect::<Result<Vec<_>, _>>()?,
            Short('a') | Long("aln") | Long("alns") | Long("alignment") | Long("alignments") =>
                args.alns = Some(parser.value()?.parse()?),
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),

            Short('n') | Long("n-reads") => args.n_reads = parser.value()?.parse_with(parse_int)?,
            Short('q') | Long("min-mapq") | Long("min-mq") => args.params.min_mapq = parser.value()?.parse()?,
            Short('c') | Long("max-clip") | Long("max-clipping") => args.params.max_clipping = parser.value()?.parse()?,
            Short('I') | Long("ins-quant") | Long("ins-quantile") => {
                args.params.ins_quantile_mult = parser.value()?.parse()?;
                args.params.ins_quantile = parser.value()?.parse()?;
            }
            Short('m') | Long("err-mult") | Long("err-multiplier") =>
                args.params.err_rate_mult = parser.value()?.parse()?,

            Short('p') | Long("ploidy") => args.params.depth.ploidy = parser.value()?.parse()?,
            Short('s') | Long("subsample") => {
                let mut values = parser.values()?;
                args.subsampling_rate = values.next()
                    .ok_or_else(|| lexopt::Error::MissingValue { option: Some("subsample".to_owned()) })?
                    .parse()?;
                args.subsampling_seed = values.next().map(|v| v.parse()).transpose()?;
            }
            Short('w') | Long("window") => args.params.depth.window_size = parser.value()?.parse()?,
            Long("window-padd") | Long("window-padding") => args.params.depth.window_padding = parser.value()?.parse()?,
            Long("boundary") => args.params.depth.boundary_size = parser.value()?.parse()?,
            Long("min-aln-lik") => args.params.depth.min_aln_prob = Ln::from_log10(parser.value()?.parse()?),
            Long("kmer-freq") | Long("kmer-frequency") => args.params.depth.max_kmer_freq = parser.value()?.parse()?,
            Long("frac-windows") | Long("fraction-windows") =>
                args.params.depth.frac_windows = parser.value()?.parse()?,
            Long("blur-extreme") => {
                args.params.depth.min_tail_obs = parser.value()?.parse()?;
                args.params.depth.tail_var_mult = parser.value()?.parse()?;
            }

            Short('^') | Long("interleaved") => args.interleaved = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,
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
    if bg_dir.exists() && args.force {
        log::warn!("Clearing output directory {}", ext::fmt::path(&bg_dir));
        fs::remove_dir_all(&bg_dir)?;
    }
    ext::sys::mkdir(&bg_dir)?;
    Ok(bg_dir)
}

fn set_strobealign_stdin(
    args: &Args,
    to_bg: bool,
    child_stdin: ChildStdin,
) -> Result<thread::JoinHandle<Result<(), io::Error>>, Error>
{
    fn create_job(args: &Args, to_bg: bool, child_stdin: ChildStdin, mut reader: impl FastxRead + 'static,
    ) -> impl FnOnce() -> Result<(), io::Error> + 'static {
        let n_reads = args.n_reads;
        let subsampling_rate = args.subsampling_rate;
        let subsampling_seed = args.subsampling_seed.clone();
        let mut writer = io::BufWriter::new(child_stdin);
        move || {
            match to_bg {
                true => {
                    let mut rng = ext::rand::init_rng(subsampling_seed);
                    reader.subsample(&mut writer, subsampling_rate, &mut rng)
                }
                false => reader.write_next_n(&mut writer, n_reads).map(|_| ()),
            }
        }
    }

    // Cannot put reader into a box, because `FastxRead` has a type parameter.
    if args.input.len() == 1 && !args.interleaved {
        let reader = fastx::Reader::new(ext::sys::open(&args.input[0])?)?;
        Ok(thread::spawn(create_job(args, to_bg, child_stdin, reader)))
    } else if args.interleaved {
        let reader = fastx::PairedEndInterleaved::new(fastx::Reader::new(ext::sys::open(&args.input[0])?)?);
        Ok(thread::spawn(create_job(args, to_bg, child_stdin, reader)))
    } else {
        let reader = fastx::PairedEndReaders::new(
            fastx::Reader::new(ext::sys::open(&args.input[0])?)?,
            fastx::Reader::new(ext::sys::open(&args.input[1])?)?);
        Ok(thread::spawn(create_job(args, to_bg, child_stdin, reader)))
    }
}

fn create_samtools_command(args: &Args, stdin: Stdio, out_file: &Path) -> Command {
    let mut cmd = Command::new(&args.samtools);
    // See SAM flags here: https://broadinstitute.github.io/picard/explain-flags.html.
    cmd.args(&["view",
            "--bam", // Output BAM.
            // Ignore reads where any of the mates is unmapped,
            // + ignore secondary & supplementary alignments + ignore failed checks.
            "--excl-flags", "3852",
            "--min-MQ", &args.params.min_mapq.to_string()
            ])
        .arg("-o").arg(out_file)
        .stdin(stdin);
    cmd
}

fn first_step_str(args: &Args, to_bg: bool) -> String {
    let mut s = String::new();
    if to_bg && args.subsampling_rate == 1.0 {
        return s;
    }

    if to_bg {
        write!(s, "_subsample_ --rate {}", args.subsampling_rate).unwrap();
        if let Some(seed) = args.subsampling_seed {
            write!(s, " --seed {}", seed).unwrap();
        }
    } else {
        write!(s, "_head_ --n-reads {}", args.n_reads).unwrap();
    }
    write!(s, " -i {}", ext::fmt::path(&args.input[0])).unwrap();
    if args.input.len() > 1 {
        write!(s, " {}", ext::fmt::path(&args.input[1])).unwrap();
    }
    write!(s, " | ").unwrap();
    s
}

/// Run strobealign for the input WGS reads.
/// If `to_bg` is true, reads are mapped to a single background region. If needed, subsampling is performed.
/// If `to_bg` is false, first `args.n_reads` reads are mapped to the whole genome.
// fn run_strobealign(args: &Args, ref_filename: &Path, out_bam: &Path, to_bg: bool) -> Result<(), Error> {
fn run_strobealign(args: &Args, ref_filename: &Path, out_bam: &Path, to_bg: bool) -> Result<(), Error> {
    let tmp_bam = out_bam.with_extension("tmp.bam");
    if tmp_bam.exists() {
        fs::remove_file(&tmp_bam)?;
    }
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
    if to_bg && args.subsampling_rate == 1.0 {
        // Provide input files as they are.
        if args.interleaved {
            strobealign.arg("--interleaved");
        }
        strobealign.arg(&ref_filename).args(&args.input);
    } else {
        // Subsample or take first N records from the input files, and pass them through stdin.
        strobealign.arg("--interleaved").arg(&ref_filename).arg("-").stdin(Stdio::piped());
    }
    let mut strobealign_child = strobealign.spawn()?;
    let strobealign_stdin = strobealign_child.stdin.take();
    let strobealign_stdout = strobealign_child.stdout.take();
    let mut guard = ext::sys::ChildGuard::new(strobealign_child);

    let handle = strobealign_stdin.map(|stdin| set_strobealign_stdin(args, to_bg, stdin)).transpose()?;
    let mut samtools = create_samtools_command(args, Stdio::from(strobealign_stdout.unwrap()), &tmp_bam);
    log::debug!("    {}{} | {}", first_step_str(&args, to_bg), ext::fmt::command(&strobealign),
        ext::fmt::command(&samtools));

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

/// Map reads to the genome and to the background genome to estimate background distributions based solely on WGS reads.
fn estimate_bg_from_reads(
    args: &Args,
    bg_fasta_filename: &Path,
    out_dir: &Path,
    contig_set: &ContigSet,
) -> Result<BgDistr, Error>
{
    let bam_filename1 = out_dir.join("to_ref.bam");
    log::info!("Mapping reads to the whole reference");
    run_strobealign(args, args.reference.as_ref().unwrap(), &bam_filename1, false)?;

    let bam_filename2 = out_dir.join("to_bg.bam");
    let contig0 = ContigId::new(0);
    let interval = Interval::full_contig(Rc::clone(contig_set.contigs()), contig0);
    log::info!("Mapping reads to non-duplicated region {}", interval);
    run_strobealign(args, &bg_fasta_filename, &bam_filename2, true)?;

    log::debug!("Loading mapped reads into memory ({})", ext::fmt::path(&bam_filename1));
    let mut bam_reader1 = bam::Reader::from_path(&bam_filename1)?;
    let params = &args.params;
    let records1: Vec<BamRecord> = bam_reader1.records()
        .filter(|res: &Result<BamRecord, _>| match res {
            Ok(record) => record.mapq() >= params.min_mapq && cigar::clipping_rate(record) <= params.max_clipping,
            Err(_) => true,
        })
        .collect::<Result<_, _>>()?;
    let insert_distr = if args.input.len() > 1 || args.interleaved {
        let pairings = ReadMateGrouping::from_unsorted_bam(records1.iter(), Some(records1.len()));
        InsertDistr::estimate(&pairings, params)?
    } else {
        // Single-end input.
        InsertDistr::undefined()
    };

    let max_insert_size = i64::from(insert_distr.max_size());
    let err_prof = ErrorProfile::estimate(records1.iter(), |record| cigar::Cigar::from_raw(record.raw_cigar()),
        max_insert_size, params);

    log::debug!("Loading mapped reads into memory ({})", ext::fmt::path(&bam_filename2));
    let mut bam_reader2 = bam::Reader::from_path(&bam_filename2)?;
    let records2: Vec<BamRecord> = bam_reader2.records()
        .filter(|res: &Result<BamRecord, _>| match res {
            Ok(record) => record.mapq() >= params.min_mapq,
            Err(_) => true,
        })
        .collect::<Result<_, _>>()?;
    let depth_distr = ReadDepth::estimate(records2.iter(), &interval, contig_set.get_seq(contig0),
        contig_set.kmer_counts(), &err_prof, max_insert_size, &params.depth, args.subsampling_rate);
    Ok(BgDistr::new(insert_distr, err_prof, depth_distr))
}

/// Estimate background distributions based on WGS read mappings to the genome.
fn estimate_bg_from_alns(
    bam_filename: &Path,
    args: &Args,
    bg_fasta_filename: &Path,
    bg_contig_set: &ContigSet,
    bg_interval_name: &str,
) -> Result<BgDistr, Error>
{
    log::debug!("Comparing reference sequence");
    let ref_fasta_filename = args.reference.as_ref().unwrap();
    let (ref_contigs, mut ref_fasta) = ContigNames::load_indexed_fasta("reference", ref_fasta_filename)?;
    let bg_interval = Interval::parse(bg_interval_name, &ref_contigs)?;
    let ref_seq = bg_interval.fetch_seq(&mut ref_fasta)?;
    let bg_seq = bg_contig_set.get_seq(ContigId::new(0));
    if ref_seq != bg_seq {
        return Err(Error::InvalidInput(format!(
            "Fasta file {} contains interval {}, but it has different sequence than the background fasta file {}",
            ext::fmt::path(ref_fasta_filename), bg_interval_name, ext::fmt::path(bg_fasta_filename))));
    }

    log::debug!("Loading mapped reads into memory ({})", ext::fmt::path(&bam_filename));
    let mut bam_reader = bam::IndexedReader::from_path(bam_filename)?;
    bam_reader.set_reference(ref_fasta_filename)?;
    let params = &args.params;
    bam_reader.fetch((bg_interval.contig_name(), i64::from(bg_interval.start()), i64::from(bg_interval.end())))?;
    let records: Vec<BamRecord> = bam_reader.records()
        .filter(|res: &Result<BamRecord, _>| match res {
            Ok(record) => record.mapq() >= params.min_mapq && cigar::clipping_rate(record) <= params.max_clipping,
            Err(_) => true,
        })
        .collect::<Result<_, _>>()?;
    let pairings = ReadMateGrouping::from_mixed_bam(&records)?;

    let insert_distr = InsertDistr::estimate(&pairings, params)?;
    let seq_shift = bg_interval.start();
    let max_insert_size = i64::from(insert_distr.max_size());
    let err_prof = ErrorProfile::estimate(records.iter(),
        |record| cigar::Cigar::infer_ext_cigar(record, &bg_seq, seq_shift), max_insert_size, params);
    let depth_distr = ReadDepth::estimate(records.iter(), &bg_interval, &bg_seq, bg_contig_set.kmer_counts(), &err_prof,
        max_insert_size, &params.depth, 1.0);
    Ok(BgDistr::new(insert_distr, err_prof, depth_distr))
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args = parse_args(argv)?.validate()?;
    let out_dir = create_out_dir(&args)?;

    log::info!("Loading background non-duplicated region into memory");
    let db_path = args.database.as_ref().unwrap();
    let db_bg_dir = db_path.join(paths::BG_DIR);
    let bg_fasta_filename = db_bg_dir.join(paths::BG_FASTA);
    let kmers_filename = db_bg_dir.join(paths::KMERS);
    let mut descriptions = Vec::with_capacity(1);
    let contig_set = ContigSet::load("bg", &bg_fasta_filename, &kmers_filename, &mut descriptions)?;
    if contig_set.len() != 1 {
        return Err(Error::InvalidData(format!("File {:?} contains more than one region", bg_fasta_filename)));
    }
    let bg_interval_name = descriptions[0].as_ref().ok_or_else(|| Error::InvalidData(format!(
        "First contig in {:?} does not have a valid description", bg_fasta_filename)))?;

    let bg_distr = if let Some(alns_filename) = args.alns.as_ref() {
        estimate_bg_from_alns(alns_filename, &args, &bg_fasta_filename, &contig_set, bg_interval_name)?
    } else {
        estimate_bg_from_reads(&args, &bg_fasta_filename, &out_dir, &contig_set)?
    };
    let mut bg_file = ext::sys::create_gzip(&out_dir.join(paths::SAMPLE_PARAMS))?;
    bg_distr.save().write_pretty(&mut bg_file, 4)?;
    log::info!("Success!");
    Ok(())
}
