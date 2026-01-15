use std::{
    path::{Path, PathBuf},
    time::Instant,
    io::{Write, BufRead},
    ffi::OsStr,
};
use colored::Colorize;
use crate::{
    ext,
    err::{self, validate_param, error},
};

struct Args {
    input: Option<PathBuf>,
    output: Option<PathBuf>,
    alignments: PathBuf,
    subset_loci: Option<PathBuf>,

    skip_tree: bool,
    div_field: String,
    div_thresh: f64,

    // all_pairs: bool,
    // pairs: Vec<String>,
    // pairs_file: Option<PathBuf>,
    // threads: u16,

    // params: dist::Params,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: None,
            output: None,
            alignments: PathBuf::from("haplotypes.paf.gz"),
            subset_loci: None,

            skip_tree: false,
            div_field: "dv".to_string(),
            div_thresh: 0.005,

            // all_pairs: false,
            // pairs: Vec::new(),
            // pairs_file: None,
            // threads: 8,

            // params: Default::default(),
        }
    }
}

impl Args {
    fn validate(mut self) -> crate::Result<Self> {
        validate_param!(self.input.is_some(), "Input database is not provided (see -i/--input)");
        validate_param!(self.output.is_some(), "Output database is not provided (see -o/--output)");

        // validate_param!(u8::from(self.all_pairs)
        //     + u8::from(!self.pairs.is_empty())
        //     + u8::from(self.pairs_file.is_some())
        //     == 1, "Exactly one argument -p/-P/-A is required");

        // self.params.validate()?;
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 17;
    const VAL: usize = 4;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Remove similar target haplotypes.".yellow());

    print!("\n{}", "Usage:".bold());
    println!(" {} -i db -o pruned_db [args]", super::PROGRAM);

    println!("\n{}", "Input arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Input database directory.",
        "-i, --input".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Output pruned database directory.",
        "-o, --output".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Path to alignment .paf[.gz] files [{}].\n\
        {EMPTY}  Should either contain {{}}, which are then replaced with locus names,\n\
        {EMPTY}  or direct to files located in {}/loci/<locus>/{}.\n\
        {EMPTY}  Alignments can be constructed using {}.",
        "-a, --alignments".green(), "PATH".yellow(), super::fmt_def(&defaults.alignments.display()),
        "INPUT".yellow(), "PATH".yellow(), const_format::concatcp!(super::PROGRAM, " align").underline());

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  PAF field with divergence values [{}].",
        "    --div-field".green(), "STR".yellow(), super::fmt_def(&defaults.div_field));
    println!("    {:KEY$} {:VAL$}  Divergence threshold for pruning [{}].",
        "    --div-thresh".green(), "NUM".yellow(), super::fmt_def_f64(defaults.div_thresh));
    println!("    {:KEY$} {:VAL$}  Limit the pruning to loci from this file.",
        "    --subset-loci".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Do not write trees in the output directory.",
        "    --skip-tree".green(), super::flag());

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
            Short('i') | Long("input") => args.input = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),

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

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();



    log::info!("Success! Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
