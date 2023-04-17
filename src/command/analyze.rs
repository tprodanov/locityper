use std::{
    cmp::max,
    path::{Path, PathBuf},
};
use colored::Colorize;
use const_format::{str_repeat, concatcp};
use crate::{
    err::{Error, validate_param},
};

struct Args {
    input: Vec<PathBuf>,
    database: Option<PathBuf>,
    output: Option<PathBuf>,

    interleaved: bool,
    threads: u16,
    force: bool,
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

        validate_param!(self.database.is_some(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.output.is_some(), "Output directory is not provided (see -o/--output)");
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Analyze WGS dataset.".yellow());

    println!("\n{} {} analyze -i reads1.fq [reads2.fq] -d db -o out [arguments]",
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
        "    --interleaved".green(), super::flag());

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads);
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());

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
                .extend(parser.values()?.take(2).map(|s| s.parse()).collect::<Result<Vec<_>, _>>()?),
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),

            Long("interleaved") => args.interleaved = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,

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

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args = parse_args(argv)?.validate()?;
    Ok(())
}
