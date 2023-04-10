//! Preprocess WGS dataset.

use std::{
    fs,
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
    bg::{Params as BgParams, BgDistr, JsonSer},
    seq::{
        ContigNames, ContigId, Interval,
        kmers::KmerCounts,
    },
};
use super::common;

struct Args {
    input: Vec<PathBuf>,
    database: Option<PathBuf>,
    output: Option<PathBuf>,

    interleaved: bool,
    threads: u16,
    force: bool,
    strobealign: PathBuf,
    samtools: PathBuf,

    bg_params: BgParams,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::with_capacity(2),
            database: None,
            output: None,

            interleaved: false,
            threads: 4,
            force: false,
            strobealign: PathBuf::from("strobealign"),
            samtools: PathBuf::from("samtools"),

            bg_params: BgParams::default(),
        }
    }
}

fn print_help() {
    const KEY: usize = 17;
    const VAL: usize = 4;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Preprocess WGS dataset.".yellow());

    println!("\n{} {} preproc -i reads1 [reads2] -d db -o out [arguments]",
        "Usage:".bold(), env!("CARGO_PKG_NAME"));

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {} Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Database directory.",
        "-d, --db".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Output directory for the sample.",
        "-o, --output".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Input reads are interleaved.",
        "    --interleaved".green(), super::flag());

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
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),

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
    if args.output.is_none() {
        Err(lexopt::Error::from("Output directory is not provided (see -o/--output)"))?;
    }
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

fn run_strobealign(args: &Args, fasta_filename: &Path, out_dir: &Path) -> Result<PathBuf, Error> {
    let out_bam = out_dir.join("aln.bam");
    if out_bam.exists() {
        log::warn!("BAM file {} exists, skipping read mapping", super::fmt_path(&out_bam));
        return Ok(out_bam);
    }
    log::info!("Mapping reads to background region");

    let mut strobealign = Command::new(&args.strobealign);
    strobealign.args(&[
        "-N", "0", // Retain 0 additional alignments.
        "-U", // Do not output unmapped reads.
        "--no-progress",
        "-t", &args.threads.to_string()]); // Specify the number of threads.
    if args.interleaved {
        strobealign.arg("--interleaved");
    }
    strobealign.arg(&fasta_filename).args(&args.input);
    strobealign.stdout(Stdio::piped()).stderr(Stdio::piped());
    let strobealign_str = super::fmt_cmd(&strobealign);
    let strobealign_child = strobealign.spawn()?;

    let mut samtools = Command::new(&args.samtools);
    // See SAM flags here: https://broadinstitute.github.io/picard/explain-flags.html.
    samtools.args(&["view",
        "--bam", // Output BAM.
        // Ignore reads where any of the mates is unmapped,
        // + ignore secondary & supplementary alignments + ignore failed checks.
        "-F", "3852",
        "--min-MQ", stringify!(crate::bg::MIN_MAPQ), // Ignore reads with low MAPQ.
    ]);
    samtools.arg("-o").arg(&out_bam);
    samtools.stdin(Stdio::from(strobealign_child.stdout.unwrap()));
    let samtools_str = super::fmt_cmd(&samtools);
    log::debug!("    {} | \\\n{:25} {}", strobealign_str, "", samtools_str);

    let start = Instant::now();
    let output = samtools.output()?;
    log::debug!("    Finished in {}", common::fmt_duration(start.elapsed()));
    if !output.status.success() {
        return Err(Error::SubprocessFail(output));
    }
    Ok(out_bam)
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args: Args = process_args(parse_args(argv)?)?;
    let out_dir = create_out_dir(&args)?;

    let db_path = args.database.as_ref().unwrap();
    let fasta_filename = super::create::bg_fasta_filename(db_path);
    let (contigs, mut ref_seqs) = ContigNames::load_fasta(common::open(&fasta_filename)?, "bg".to_string())?;
    if contigs.len() != 1 {
        return Err(Error::InvalidData(format!("File {:?} contains more than one region", fasta_filename)));
    }
    let ref_seq = ref_seqs.pop().unwrap();

    let kmers_filename = super::create::kmers_filename(fasta_filename.parent().unwrap());
    let kmers_file = common::open(&kmers_filename)?;
    let kmer_counts = KmerCounts::load(kmers_file, contigs.lengths())?;

    let bam_filename = run_strobealign(&args, &fasta_filename, &out_dir)?;
    let mut bam_reader = bam::Reader::from_path(&bam_filename)?;
    log::debug!("    Loading mapped reads into memory");
    let records: Vec<BamRecord> = bam_reader.records().collect::<Result<_, _>>()?;

    let interval = Interval::full_contig(Rc::clone(&contigs), ContigId::new(0));
    let bg = BgDistr::estimate(&records, &interval, &ref_seq, &kmer_counts, &args.bg_params)?;

    let params_filename = out_dir.join("params.gz");
    let mut bg_file = common::create_gzip(&params_filename)?;
    bg.save().write_pretty(&mut bg_file, 4)?;
    Ok(())
}
