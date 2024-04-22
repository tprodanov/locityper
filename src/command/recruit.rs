use std::{
    time::Instant,
    path::PathBuf,
};
use crate::{
    Error,
    ext::{
        self,
        fmt::{PrettyU32, PrettyU64},
    },
    seq::{
        recruit,
        kmers::Kmer,
    },
    bg::TECHNOLOGIES,
};
use colored::Colorize;

struct Args {
    in_files: super::preproc::InputFiles,

    seqs: Vec<String>,
    seqs_all: Option<PathBuf>,
    seqs_list: Option<PathBuf>,
    output: Option<String>,
    regions: Option<PathBuf>,

    minimizer_kw: Option<(u8, u8)>,
    match_frac: Option<f64>,
    match_len: u32,
    chunk_length: u64,

    threads: u16,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            in_files: Default::default(),

            seqs: Vec::new(),
            seqs_all: None,
            seqs_list: None,
            output: None,
            regions: None,

            minimizer_kw: None,
            match_frac: None,
            match_len: 2000,
            chunk_length: 3_000_000,

            threads: 8,
        }
    }
}

impl Args {
    fn validate(self) -> Result<Self, Error> {
        Ok(self)
    }
}

pub(crate) fn fmt_def_minimizers() -> String {
    super::describe_defaults(
        TECHNOLOGIES.iter().copied(),
        |tech| tech.to_string(),
        |tech| {
            let (k, w) = tech.default_minim_size();
            format!("{} {}", k, w)
        }
    )
}

pub(crate) fn fmt_def_match_frac() -> String {
    super::describe_defaults(
        TECHNOLOGIES.iter()
            .flat_map(|&tech| [(tech, false), (tech, true)].into_iter())
            .filter(|(tech, paired_end)| !paired_end || tech.paired_end_allowed()),
        |(tech, paired_end)|
            format!("{}{}", tech, if *paired_end { "-PE" } else if tech.paired_end_allowed() { "-SE" } else { "" }),
        |(tech, paired_end)| tech.default_match_frac(*paired_end)
    )
}

fn print_help() {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Recruit reads to one/multiple loci.".yellow());

    println!("\n{}", "Usage:".bold());
    println!("    {} recruit \\", super::PROGRAM);
    println!("        (-i reads1.fq [reads2.fq] | -a reads.bam [--no-index] | -I in-list) \\");
    println!("        (-s seqs.fa | -S all_seqs.fa) -o out.fastq [args]");

    println!("\n{}  (please see {} for more information on {}/{}/{} arguments)",
        "Input arguments:".bold(), "README".italic(), "-i".green(), "-a".green(), "-I".green());
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

    println!("\n{}", "Target sequences and output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  FASTA file with target sequences for one locus.\n\
        {EMPTY}  To recruit reads to multiple loci, replace locus name with `{{}}`.\n\
        {EMPTY}  Then, all matching files will be selected.\n\
        {EMPTY}  Multiple entries allowed, but locus names must not repeat.",
        "-s, --seqs".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Single FASTA file with target sequences for all loci.\n\
        {EMPTY}  Record names should follow the format `LOCUS*SEQ_NAME`.\n\
        {EMPTY}  Mutually exclusive with {}.",
        "-S, --seqs-all".green(), "FILE".yellow(), "-s".green());
    println!("    {:KEY$} {:VAL$}  Path to the (interleaved) output FASTQ files. If more than\n\
        {EMPTY}  one target locus exists, please replace locus name with `{{}}`.\n\
        {EMPTY}  {} Will create parent directories, if needed.",
        "-o, --output".green(), "FILE".yellow(), "Note:".bright_yellow());
    println!("    {:KEY$} {:VAL$}  Two column file with input FASTA and output FASTQ filenames.\n\
        {EMPTY}  Mutually exclusive with {}, {} and {}.",
        "-l, --seqs-list".green(), "FILE".yellow(), "-s".green(), "-S".green(), "-o".green());
    println!("    {:KEY$} {:VAL$}  Recruit unmapped reads and reads with primary alignments to these\n\
        {EMPTY}  regions (BED). Only relevant for mapped and indexed BAM/CRAM files.",
        "-R, --regions".green(), "FILE".yellow());

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Interleaved paired-end reads in single input file.",
        "-^, --interleaved".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Use input BAM/CRAM file ({}) without index: goes over all reads.\n\
        {EMPTY}  Single-end and paired-end interleaved ({}) data is allowed.",
        "    --no-index".green(), super::flag(), "-a".green(), "-^".green());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads.to_string().cyan());

    println!("\n{}", "Read recruitment:".bold());
    println!("    {}  {}  Use k-mers of size {} (<= {}) with smallest hash\n\
        {EMPTY}  across {} consecutive k-mers.\n\
        {EMPTY}  Default: {}.",
        "-m, --minimizer".green(), "INT INT".yellow(),
        "INT_1".yellow(), recruit::Minimizer::MAX_KMER_SIZE, "INT_2".yellow(), fmt_def_minimizers());
    println!("    {:KEY$} {:VAL$}  Minimal fraction of minimizers that need to match reference.\n\
        {EMPTY}  Default: {}.",
        "-M, --match-frac".green(), "FLOAT".yellow(), fmt_def_match_frac());
    println!("    {:KEY$} {:VAL$}  Recruit long reads with a matching subregion of this length [{}].",
        "-L, --match-len".green(), "INT".yellow(),
        super::fmt_def(defaults.match_len));
    println!("    {:KEY$} {:VAL$}  Recruit reads in chunks of this sum length [{}].\n\
        {EMPTY}  Impacts runtime in multi-threaded read recruitment.",
        "-c, --chunk-len".green(), "INT".yellow(),
        super::fmt_def(PrettyU64(defaults.chunk_length)));

    println!("\n{}", "Other arguments:".bold());
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

            Short('s') | Long("seqs") => {
                for entry in parser.values()? {
                    args.seqs.push(entry.parse()?);
                }
            }
            Short('S') | Long("seqs-all") => args.seqs_all = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),
            Short('l') | Long("seqs-list") => args.seqs_list = Some(parser.value()?.parse()?),
            Short('R') | Long("regions") => args.regions = Some(parser.value()?.parse()?),

            Short('m') | Long("minimizer") | Long("minimizers") =>
                args.minimizer_kw = Some((parser.value()?.parse()?, parser.value()?.parse()?)),
            Short('M') | Long("match-frac") | Long("match-fraction") =>
                args.match_frac = Some(parser.value()?.parse()?),
            Short('L') | Long("match-len") | Long("match-length") =>
                args.match_len = parser.value()?.parse::<PrettyU32>()?.get(),
            Short('c') | Long("chunk") | Long("chunk-len") =>
                args.chunk_length = parser.value()?.parse::<PrettyU64>()?.get(),

            Short('^') | Long("interleaved") => args.in_files.interleaved = true,
            Long("no-index") => args.in_files.no_index = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,

            Short('V') | Long("version") => {
                super::print_version();
                std::process::exit(0);
            }
            Short('h') | Long("help") | Short('H') | Long("full-help") | Long("hidden-help") => {
                print_help();
                std::process::exit(0);
            }
            _ => Err(arg.unexpected())?,
        }
    }
    Ok(args)
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let mut args = parse_args(argv)?;
    args.in_files.fill_from_inlist()?;
    let mut args = args.validate()?;

    super::greet();
    let timer = Instant::now();
    log::info!("Success! Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
