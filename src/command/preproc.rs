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
    alns: Option<PathBuf>,
    database: Option<PathBuf>,
    reference: Option<PathBuf>,
    output: Option<PathBuf>,

    interleaved: bool,
    threads: u16,
    force: bool,
    strobealign: PathBuf,
    samtools: PathBuf,

    n_alns: u64,
    params: BgParams,
}

const DEF_N_ALNS: &'static str = "1M";

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::with_capacity(2),
            alns: None,
            database: None,
            reference: None,
            output: None,

            interleaved: false,
            threads: 4,
            force: false,
            strobealign: PathBuf::from("strobealign"),
            samtools: PathBuf::from("samtools"),

            n_alns: parse_int(DEF_N_ALNS).unwrap(),
            params: BgParams::default(),
        }
    }
}

fn print_help(extended: bool) {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Preprocess WGS dataset.".yellow());

    println!("\n{} {} preproc (-i reads1.fq [reads2.fq] | -a reads.bam) -d db -r reference.fa -o out [arguments]",
        "Usage:".bold(), env!("CARGO_PKG_NAME"));

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Reads in indexed BAM/CRAM format, already mapped to the whole genome.\n\
        {EMPTY}  Mutually exclusive with {}.",
        "-a, --alignment".green(), "FILE".yellow(), "-i/--input".green());
    println!("    {:KEY$} {:VAL$}  Database directory.",
        "-d, --db".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file. Must contain FAI index.",
        "-r, --reference".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Output directory for the sample.",
        "-o, --output".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Input reads are interleaved.",
        "    --interleaved".green(), super::flag());

    if extended {
        println!("\n{}", "Insert size and error profile estimation:".bold());
        println!("    {:KEY$} {:VAL$}  Use first {} alignments to estimate profiles [{}].",
            "-n, --n-alns".green(), "INT".yellow(), "INT".yellow(), DEF_N_ALNS);
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
        println!("    {:KEY$} {}\n\
            {EMPTY}  Min allowed alignment log-likelihood is calculated as {}\n\
            {EMPTY}  multiplied by the {}-th log-likelihood quantile [{}, {}].",
            "-E, --err-quantile".green(), "FLOAT FLOAT".yellow(),
            "FLOAT_1".yellow(), "FLOAT_2".yellow(), defaults.params.err_quantile_mult, defaults.params.err_quantile);
        println!("    {:KEY$} {:VAL$}  Multiply error rates by this factor, in order to correct for\n\
            {EMPTY}  read mappings missed due to higher error rate [{}].",
            "-m, --err-mult".green(), "FLOAT".yellow(), defaults.params.err_rate_mult);

        println!("\n{}", "Background read depth estimation:".bold());
        println!("    {:KEY$} {:VAL$}  Specie ploidy [{}].",
            "-p, --ploidy".green(), "INT".yellow(), defaults.params.depth.ploidy);
        println!("    {:KEY$} {:VAL$}  Count read depth per {} bp windows [{}].",
            "-w, --window".green(), "INT".yellow(), "INT".yellow(), defaults.params.depth.window_size);
        println!("    {:KEY$} {:VAL$}  Extend window by {} bp to both sides in order to calculate\n\
            {EMPTY}  GC-content and average k-mer frequency [{}].",
            "    --window-padd".green(), "INT".yellow(), "INT".yellow(), defaults.params.depth.window_padding);
        println!("    {:KEY$} {:VAL$}  Skip {} bp near the edge of the background region.\n\
            {EMPTY}  Must not be smaller than {} [{}].",
            "    --edge-padd".green(), "INT".yellow(), "INT".yellow(), "--window-padd".green(),
            defaults.params.depth.edge_padding);
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
            Short('i') | Long("input") => args.input
                .extend(parser.values()?.take(2).map(|s| s.parse()).collect::<Result<Vec<_>, _>>()?),
            Short('a') | Long("aln") | Long("alns") | Long("alignment") | Long("alignments") =>
                args.alns = Some(parser.value()?.parse()?),
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),

            Short('n') | Long("n-alns") => args.n_alns = parser.value()?.parse_with(parse_int)?,
            Short('q') | Long("min-mapq") | Long("min-mq") => args.params.min_mapq = parser.value()?.parse()?,
            Short('c') | Long("max-clip") | Long("max-clipping") => args.params.max_clipping = parser.value()?.parse()?,
            Short('I') | Long("ins-quant") | Long("ins-quantile") => {
                args.params.ins_quantile_mult = parser.value()?.parse()?;
                args.params.ins_quantile = parser.value()?.parse()?;
            }
            Short('E') | Long("err-quant") | Long("err-quantile") => {
                args.params.err_quantile_mult = parser.value()?.parse()?;
                args.params.err_quantile = parser.value()?.parse()?;
            }
            Short('m') | Long("err-mult") | Long("err-multiplier") =>
                args.params.err_rate_mult = parser.value()?.parse()?,

            Short('p') | Long("ploidy") => args.params.depth.ploidy = parser.value()?.parse()?,
            Short('w') | Long("window") => args.params.depth.window_size = parser.value()?.parse()?,
            Long("window-padd") | Long("window-padding") => args.params.depth.window_padding = parser.value()?.parse()?,
            Long("edge-padd") | Long("edge-padding") => args.params.depth.edge_padding = parser.value()?.parse()?,
            Long("kmer-freq") | Long("kmer-frequency") => args.params.depth.max_kmer_freq = parser.value()?.parse()?,
            Long("frac-windows") | Long("fraction-windows") =>
                args.params.depth.frac_windows = parser.value()?.parse()?,
            Long("blur-extreme") => {
                args.params.depth.min_tail_obs = parser.value()?.parse()?;
                args.params.depth.tail_var_mult = parser.value()?.parse()?;
            }

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

fn process_args(mut args: Args) -> Result<Args, Error> {
    args.threads = max(args.threads, 1);
    let n_input = args.input.len();
    if n_input == 0 && args.alns.is_none() {
        Err(lexopt::Error::from("Neither read files, nor alignment files are not provided (see -i and -a)"))?;
    } else if n_input > 0 && args.alns.is_some() {
        Err(lexopt::Error::from("Read files (-i) and an alignment file (-a) cannot be provided together"))?;
    }

    // Ignore error if interleaved and n_input == 0.
    if args.interleaved && n_input == 2 {
        Err(lexopt::Error::from("Two read files (-i/--input) are provided, however, --interleaved is specified"))?;
    } else if n_input == 1 {
        log::warn!("Running in single-end mode.");
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
            "--min-MQ", &args.params.min_mapq.to_string()
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

/// Map reads to the genome and to the background genome to estimate background distributions based solely on WGS reads.
fn estimate_bg_from_reads(
    args: &Args,
    bg_fasta_filename: &Path,
    out_dir: &Path,
    interval: &Interval,
    ref_seq: &[u8],
    kmer_counts: &KmerCounts,
) -> Result<BgDistr, Error>
{
    let bam_filename1 = run_strobealign_full_ref(&args, &out_dir)?;
    let bam_filename2 = run_strobealign_bg_region(&args, &bg_fasta_filename, &out_dir)?;

    log::debug!("Loading mapped reads into memory ({})", super::fmt_path(&bam_filename1));
    let mut bam_reader1 = bam::Reader::from_path(&bam_filename1)?;
    let params = &args.params;
    let records1: Vec<BamRecord> = bam_reader1.records()
        .filter(|res: &Result<BamRecord, _>| match res {
            Ok(record) => record.mapq() >= params.min_mapq && cigar::clipping_rate(record) <= params.max_clipping,
            Err(_) => true,
        })
        .collect::<Result<_, _>>()?;
    let pairings = ReadMateGrouping::from_unsorted_bam(records1.iter(), Some(records1.len()));

    let insert_distr = InsertDistr::estimate(&pairings, params)?;
    let err_prof = ErrorProfile::estimate(records1.iter(),
        |record| cigar::Cigar::from_raw(record.raw_cigar()),
        i64::from(insert_distr.max_size()), params);

    log::debug!("Loading mapped reads into memory ({})", super::fmt_path(&bam_filename2));
    let mut bam_reader2 = bam::Reader::from_path(&bam_filename2)?;
    let records2: Vec<BamRecord> = bam_reader2.records()
        .filter(|res: &Result<BamRecord, _>| match res {
            Ok(record) => record.mapq() >= params.min_mapq,
            Err(_) => true,
        })
        .collect::<Result<_, _>>()?;
    let depth_distr = ReadDepth::estimate(records2.iter(), interval, ref_seq, kmer_counts, &err_prof,
        i64::from(insert_distr.max_size()), &params.depth);
    Ok(BgDistr::new(insert_distr, err_prof, depth_distr))
}

/// Estimate background distributions based on WGS read mappings to the genome.
fn estimate_bg_from_alns(
    bam_filename: &Path,
    args: &Args,
    bg_fasta_filename: &Path,
    bg_seq: &[u8],
    bg_interval_name: &str,
) -> Result<BgDistr, Error>
{
    let ref_fasta_filename = args.reference.as_ref().unwrap();
    let (ref_contigs, mut ref_fasta) = ContigNames::load_indexed_fasta(ref_fasta_filename, "full_ref".to_owned())?;
    let bg_interval = Interval::parse(bg_interval_name, &ref_contigs)?;
    let ref_seq = bg_interval.fetch_seq(&mut ref_fasta)?;
    if ref_seq != bg_seq {
        return Err(Error::InvalidInput(format!(
            "Fasta file {} contains interval {}, but it has different sequence than the background fasta file {}",
            super::fmt_path(ref_fasta_filename), bg_interval_name, super::fmt_path(bg_fasta_filename))));
    }

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

    unimplemented!()
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args: Args = process_args(parse_args(argv)?)?;
    let out_dir = create_out_dir(&args)?;

    log::info!("Loading background non-duplicated region into memory");
    let db_path = args.database.as_ref().unwrap();
    let bg_fasta_filename = super::create::bg_fasta_filename(db_path);
    let mut descriptions = Vec::with_capacity(1);
    let (contigs, mut ref_seqs) = ContigNames::load_fasta(common::open(&bg_fasta_filename)?, "bg".to_owned(),
        &mut descriptions)?;
    if contigs.len() != 1 {
        return Err(Error::InvalidData(format!("File {:?} contains more than one region", bg_fasta_filename)));
    }
    let bg_interval_name = descriptions[0].as_ref().ok_or_else(|| Error::InvalidData(format!(
        "First contig in {:?} does not have a valid description", bg_fasta_filename)))?;
    let ref_seq = ref_seqs.pop().unwrap();
    let interval = Interval::full_contig(Rc::clone(&contigs), ContigId::new(0));

    let kmers_filename = super::create::kmers_filename(bg_fasta_filename.parent().unwrap());
    let kmers_file = common::open(&kmers_filename)?;
    let kmer_counts = KmerCounts::load(kmers_file, contigs.lengths())?;

    let bg_distr = if let Some(alns_filename) = args.alns.as_ref() {
        estimate_bg_from_alns(alns_filename, &args, &bg_fasta_filename, &ref_seq, bg_interval_name)?
    } else {
        estimate_bg_from_reads(&args, &bg_fasta_filename, &out_dir, &interval, &ref_seq, &kmer_counts)?
    };
    let mut bg_file = common::create_gzip(&out_dir.join("params.gz"))?;
    bg_distr.save().write_pretty(&mut bg_file, 4)?;
    log::info!("Success!");
    Ok(())
}
