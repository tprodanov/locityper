use std::{
    time::Instant,
};
use colored::Colorize;
use crate::{
    err::{error, add_path, validate_param},
    ext::{self},
};

struct Args {
    threads: u16,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            threads: 8,
        }
    }
}

impl Args {
    fn validate(self) -> crate::Result<Self> {
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 16;
    const VAL: usize = 4;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Augment database.".yellow());

    print!("\n{}", "Usage:".bold());
    println!(" {} [TODO]", super::PROGRAM);

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Database with loci haplotypes.",
        "-d, --database".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Limit augmentation to these loci.",
        "    --subset-loci".green(), "STR+".yellow());
    // println!("    {:KEY$} {:VAL$}")

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), super::fmt_def(defaults.threads));

    println!("\n{}", "Other arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Show this help message.", "-h, --help".green(), "");
    println!("    {:KEY$} {:VAL$}  Show version.", "-V, --version".green(), "");
}

fn parse_args(argv: &[String]) -> crate::Result<Args> {
    if argv.is_empty() {
        print_help();
        std::process::exit(1);
    }
    use lexopt::prelude::*;
    let mut args = Args::default();
    let mut parser = lexopt::Parser::from_args(argv);

    while let Some(arg) = parser.next()? {
        match arg {
            // Short('p') | Long("paf") => args.paf = Some(parser.value()?.parse()?),
            // Short('f') | Long("fasta") => args.fasta = Some(parser.value()?.parse()?),
            // Short('d') | Long("discarded") => args.disc_filename = parser.value()?.parse()?,
            // Short('o') | Long("output") => {
            //     let mut values = parser.values()?.take(2);
            //     args.out_merged = Some(values.next().expect("First argument is always present").parse()?);
            //     if let Some(val) = values.next() {
            //         args.out_separate = Some(val.parse()?);
            //     }
            // }
            // Short('r') | Long("ref-hap") => args.ref_hap = Some(parser.value()?.parse()?),
            // Short('R') | Long("region") => args.region = parser.value()?.parse()?,

            Short('V') | Long("version") => {
                super::print_version();
                std::process::exit(0);
            }
            Short('h') | Long("help") | Long("full-help") | Long("hidden-help") => {
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
