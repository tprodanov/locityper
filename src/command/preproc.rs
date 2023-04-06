//! Preprocess WGS dataset.

use std::{
    cmp::max,
    path::PathBuf,
};
use colored::Colorize;
use crate::Error;

struct Args {
    input: Vec<PathBuf>,
    database: Option<PathBuf>,
    output: Option<PathBuf>,

    threads: u16,
    force: bool,
    strobealign: PathBuf,
    interleaved: bool,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::with_capacity(2),
            database: None,
            output: None,

            threads: 4,
            force: false,
            strobealign: PathBuf::from("strobealign"),
            interleaved: false,
        }
    }
}

fn print_help() {
    const KEY: usize = 17;
    const VAL: usize = 5;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Preprocess WGS dataset.".yellow());

    println!("\n{} {} preproc -i reads1 [reads2] -d db -o out [arguments]",
        "Usage:".bold(), env!("CARGO_PKG_NAME"));

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Database directory.",
        "-d, --db".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Output directory for the sample.",
        "-o, --output".green(), "DIR".yellow());

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads);
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Strobealign executable [{}].",
        "    --strobealign".green(), "EXE".yellow(), defaults.strobealign.display());
    println!("    {:KEY$} {:VAL$}  Input reads are interleaved.",
        "    --interleaved".green(), super::flag());

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

            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,
            Long("strobealign") => args.strobealign = parser.value()?.parse()?,
            Long("inter") | Long("interleaved") => args.interleaved = true,

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
    if args.input.is_empty() {
        Err(lexopt::Error::from("Read files are not provided (see -i/--input)"))?;
    }
    if args.input.len() > 2 {
        Err(lexopt::Error::from("Too many read files (see -i/--input). Expected 1 or 2 files"))?;
    }
    if args.database.is_none() {
        Err(lexopt::Error::from("Database directory is not provided (see -d/--database)"))?;
    }
    if args.output.is_none() {
        Err(lexopt::Error::from("Output directory is not provided (see -o/--output)"))?;
    }
    args.strobealign = super::find_exe(args.strobealign)?;
    Ok(args)
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args: Args = process_args(parse_args(argv)?)?;
    Ok(())
}
