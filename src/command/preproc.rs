//! Preprocess WGS dataset.

use std::{
    fs::{self, File},
    io::BufReader,
    fmt::Write as FmtWrite,
    cmp::max,
    path::{Path, PathBuf},
    process::{Stdio, Command},
    time::Instant,
    rc::Rc,
};
use colored::Colorize;
use htslib::bam::{
    self,
    Read as BamRead,
    record::Record as BamRecord,
};
use crate::{
    Error,
    algo::parse_int,
    seq::{
        ContigNames, ContigId, Interval,
        kmers::KmerCounts,
        cigar,
    },
    bg::{
        BgDistr, JsonSer, Params as BgParams,
        insertsz::{ReadMateGrouping, InsertDistr},
        err_prof::ErrorProfile,
        depth::ReadDepth,
    },
};
use super::common;

struct Args {
    input: Vec<PathBuf>,
    database: Option<PathBuf>,
    reference: Option<PathBuf>,
    output: Option<PathBuf>,

    interleaved: bool,
    threads: u16,
    force: bool,
    strobealign: PathBuf,
    samtools: PathBuf,

    n_alns: u64,
    min_mapq: u8,
    params: BgParams,
}

const DEF_N_ALNS: &'static str = "1M";

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::with_capacity(2),
            database: None,
            reference: None,
            output: None,

            interleaved: false,
            threads: 4,
            force: false,
            strobealign: PathBuf::from("strobealign"),
            samtools: PathBuf::from("samtools"),

            n_alns: parse_int(DEF_N_ALNS).unwrap(),
            min_mapq: 20,
            params: BgParams::default(),
        }
    }
}

fn print_help() {
    const KEY: usize = 17;
    const VAL: usize = 4;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Preprocess WGS dataset.".yellow());

    println!("\n{} {} preproc -i reads1 [reads2] -d db -r reference.fa -o out [arguments]",
        "Usage:".bold(), env!("CARGO_PKG_NAME"));

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {} Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Database directory.",
        "-d, --db".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file. Must contain FAI index.",
        "-r, --reference".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Output directory for the sample.",
        "-o, --output".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Input reads are interleaved.",
        "    --interleaved".green(), super::flag());

    println!("\n{}", "Optional parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Estimate insert size and error profiles from the first {}\n\
        {EMPTY}  alignments to the reference genome [{}].",
        "-n, --n-alns".green(), "INT".yellow(), "INT".yellow(), DEF_N_ALNS);
    println!("    {:KEY$} {:VAL$}  Ignore reads with mapping quality less than {} [{}].",
        "-q, --min-mapq".green(), "INT".yellow(), "INT".yellow(), defaults.min_mapq);
    println!("    {:width$}  Please rerun with {}, if {} or {} values have changed.",
        "WARNING".red(), "--force".red(), "-n".green(), "-q".green(), width = KEY + VAL + 1);

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads);
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Strobealign executable [{}].",
        "    --strobealign".green(), "EXE".yellow(), defaults.strobealign.display());
    println!("    {:KEY$} {:VAL$}  Samtools executable [{}].",
        "    --samtools".green(), "EXE".yellow(), defaults.samtools.display());

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
                .extend(parser.values()?.map(|s| s.parse()).collect::<Result<Vec<_>, _>>()?),
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),

            Short('n') | Long("n-alns") => args.n_alns = parser.value()?.parse_with(parse_int)?,
            Short('q') | Long("min-mapq") | Long("min-mq") => args.min_mapq = parser.value()?.parse()?,
            Long("inter") | Long("interleaved") => args.interleaved = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,
            Long("strobealign") => args.strobealign = parser.value()?.parse()?,
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

fn process_args(mut args: Args) -> Result<Args, Error> {
    args.threads = max(args.threads, 1);
    let n_input = args.input.len();
    if n_input == 0 {
        Err(lexopt::Error::from("Read files are not provided (see -i/--input)"))?;
    } else if args.interleaved && n_input != 1 {
        Err(lexopt::Error::from("Two read files (-i/--input) are provided, however, --interleaved is specified"))?;
    } else if n_input > 2 {
        Err(lexopt::Error::from("Too many read files (see -i/--input). Expected 1 or 2 files"))?;
    }

    if args.database.is_none() {
        Err(lexopt::Error::from("Database directory is not provided (see -d/--database)"))?;
    }
    if args.reference.is_none() {
        Err(lexopt::Error::from("Reference fasta file is not provided (see -r/--reference)"))?;
    }
    if args.output.is_none() {
        Err(lexopt::Error::from("Output directory is not provided (see -o/--output)"))?;
    }
    args.params.validate()?;
    args.strobealign = super::find_exe(args.strobealign)?;
    args.samtools = super::find_exe(args.samtools)?;
    Ok(args)
}

fn create_out_dir(args: &Args) -> Result<PathBuf, Error> {
    let out_dir = args.output.as_ref().unwrap();
    super::mkdir(out_dir)?;

    let bg_dir = out_dir.join("bg");
    if bg_dir.exists() && args.force {
        log::warn!("Clearing output directory {}", super::fmt_path(&bg_dir));
        fs::remove_dir_all(&bg_dir)?;
    }
    super::mkdir(&bg_dir)?;
    Ok(bg_dir)
}

/// Run strobealign.
/// If `head_lines` is `Some`, run `head -n INT` before running `samtools view`.
/// In such case, pass additional parameters `-f 0.001 -R 0` to `strobalign`.
fn run_strobealign(args: &Args, ref_filename: &Path, out_bam: &Path, head_lines: Option<u64>) -> Result<(), Error> {
    let tmp_bam = out_bam.with_extension("tmp.bam");
    if tmp_bam.exists() {
        fs::remove_file(&tmp_bam)?;
    }
    if out_bam.exists() {
        log::warn!("    BAM file {} exists, skipping read mapping", super::fmt_path(out_bam));
        return Ok(());
    }

    let mut strobealign = Command::new(&args.strobealign);
    strobealign.args(&[
        "-N0", // Retain 0 additional alignments.
        "-R0", // Do not rescue reads,
        "-U", // Do not output unmapped reads.
        // "--no-progress",
        "-t", &args.threads.to_string()]); // Specify the number of threads.
    if head_lines.is_some() {
        strobealign.args(&[
            "-f", "0.001"]); // Discard more minimizers to speed up alignment
    }
    if args.interleaved {
        strobealign.arg("--interleaved");
    }
    strobealign.arg(&ref_filename).args(&args.input);
    strobealign.stdout(Stdio::piped());
    let mut pipeline_str = super::fmt_cmd(&strobealign);
    let strobealign_child = strobealign.spawn()?;

    let latest_stdout = if let Some(n) = head_lines {
        let mut head = Command::new("head");
        head.args(&["-n", &n.to_string()])
            .stdin(Stdio::from(strobealign_child.stdout.unwrap()))
            .stdout(Stdio::piped());
        write!(pipeline_str, " | {}", super::fmt_cmd(&head)).unwrap();
        let head_child = head.spawn()?;
        head_child.stdout.unwrap()
    } else {
        strobealign_child.stdout.unwrap()
    };

    let mut samtools = Command::new(&args.samtools);
    // See SAM flags here: https://broadinstitute.github.io/picard/explain-flags.html.
    samtools.args(&["view",
            "--bam", // Output BAM.
            // Ignore reads where any of the mates is unmapped,
            // + ignore secondary & supplementary alignments + ignore failed checks.
            "--excl-flags", "3852",
            "--min-MQ", &args.min_mapq.to_string()
            ])
        .arg("-o").arg(&tmp_bam)
        .stdin(Stdio::from(latest_stdout));
    write!(pipeline_str, " | {}", super::fmt_cmd(&samtools)).unwrap();
    log::debug!("    {}", pipeline_str);

    let start = Instant::now();
    let output = samtools.output()?;
    log::debug!("");
    log::debug!("    Finished in {}", common::fmt_duration(start.elapsed()));
    if !output.status.success() {
        return Err(Error::SubprocessFail(output));
    }
    fs::rename(&tmp_bam, out_bam)?;
    Ok(())
}

/// Run strobealign to map some reads to the full reference.
/// Keep approximately the first `n` alignments.
/// Returns the path to the output BAM file.
fn run_strobealign_full_ref(args: &Args, out_dir: &Path) -> Result<PathBuf, Error> {
    let out_bam = out_dir.join("to_ref.bam");
    log::info!("Mapping reads to the whole reference");

    let ref_filename = args.reference.as_ref().unwrap();
    // Need to know the number of contigs in order to
    let fai_filename = common::append_path(ref_filename, ".fai");
    let fai_file = File::open(&fai_filename).map(BufReader::new)?;
    let n_contigs = common::count_lines(fai_file)?;

    // Save the header + n_alns + 10 lines just in case.
    let head_lines = n_contigs + args.n_alns + 10;
    run_strobealign(args, &ref_filename, &out_bam, Some(head_lines))?;
    Ok(out_bam)
}

/// Run strobealign to map all reads to a single background region.
/// Returns the path to the output BAM file.
fn run_strobealign_bg_region(args: &Args, fasta_filename: &Path, out_dir: &Path) -> Result<PathBuf, Error> {
    let out_bam = out_dir.join("to_bg.bam");
    log::info!("Mapping reads to a single non-duplicated region");
    run_strobealign(args, fasta_filename, &out_bam, None)?;
    Ok(out_bam)
}

fn process_first_bam(bam_filename: &Path, params: &BgParams) -> Result<(InsertDistr, ErrorProfile), Error> {
    log::debug!("Loading mapped reads into memory ({})", super::fmt_path(bam_filename));
    let mut bam_reader = bam::Reader::from_path(&bam_filename)?;
    let records: Vec<BamRecord> = bam_reader.records().collect::<Result<_, _>>()?;
    let full_mappings: Vec<_> = records.iter()
        .filter(|record| cigar::clipping_rate(record) <= params.max_clipping)
        .collect();
    let pairings = ReadMateGrouping::from_unsorted_bam(full_mappings.iter().copied(), Some(full_mappings.len()));

    let insert_distr = InsertDistr::estimate(&pairings, params);
    let err_prof = ErrorProfile::estimate(full_mappings.into_iter(),
        |record| cigar::Cigar::from_raw(record.raw_cigar()),
        i64::from(insert_distr.max_size()), params);
    Ok((insert_distr, err_prof))
}

fn process_second_bam(
    bam_filename: &Path,
    interval: &Interval,
    ref_seq: &[u8],
    kmer_counts: &KmerCounts,
    err_prof: &ErrorProfile,
    insert_distr: &InsertDistr,
    params: &BgParams,
) -> Result<ReadDepth, Error> {
    log::debug!("Loading mapped reads into memory ({})", super::fmt_path(bam_filename));
    let mut bam_reader = bam::Reader::from_path(&bam_filename)?;
    let records: Vec<BamRecord> = bam_reader.records().collect::<Result<_, _>>()?;
    let max_insert_size = i64::from(insert_distr.max_size());
    Ok(ReadDepth::estimate(records.iter(), interval, ref_seq, kmer_counts, err_prof, max_insert_size, &params.depth))
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args: Args = process_args(parse_args(argv)?)?;
    let out_dir = create_out_dir(&args)?;

    log::info!("Loading background non-duplicated region into memory");
    let db_path = args.database.as_ref().unwrap();
    let fasta_filename = super::create::bg_fasta_filename(db_path);
    let (contigs, mut ref_seqs) = ContigNames::load_fasta(common::open(&fasta_filename)?, "bg".to_string())?;
    if contigs.len() != 1 {
        return Err(Error::InvalidData(format!("File {:?} contains more than one region", fasta_filename)));
    }
    let ref_seq = ref_seqs.pop().unwrap();
    let interval = Interval::full_contig(Rc::clone(&contigs), ContigId::new(0));

    let kmers_filename = super::create::kmers_filename(fasta_filename.parent().unwrap());
    let kmers_file = common::open(&kmers_filename)?;
    let kmer_counts = KmerCounts::load(kmers_file, contigs.lengths())?;

    let bam_filename1 = run_strobealign_full_ref(&args, &out_dir)?;
    let bam_filename2 = run_strobealign_bg_region(&args, &fasta_filename, &out_dir)?;

    let (insert_distr, err_prof) = process_first_bam(&bam_filename1, &args.params)?;
    process_second_bam(
        &bam_filename2,
        &interval,
        &ref_seq,
        &kmer_counts,
        &err_prof,
        &insert_distr,
        &args.params,
    )?;

    // let mut bam_reader = bam::Reader::from_path(&bam_filename)?;
    // log::debug!("    Loading mapped reads into memory");
    // let records: Vec<BamRecord> = bam_reader.records().collect::<Result<_, _>>()?;

    // let interval = Interval::full_contig(Rc::clone(&contigs), ContigId::new(0));
    // let bg = BgDistr::estimate(&records, &interval, &ref_seq, &kmer_counts, &args.params)?;

    // let params_filename = out_dir.join("params.gz");
    // let mut bg_file = common::create_gzip(&params_filename)?;
    // bg.save().write_pretty(&mut bg_file, 4)?;
    Ok(())
}
