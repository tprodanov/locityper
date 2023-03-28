//! Create a new database.

use colored::Colorize;

#[derive(PartialEq, Eq)]
enum Input {
    Fasta(String),
    Agc(String, String),
    None,
}

struct Args {
    input: Input,
    database: String,
    kmer_size: u8,
    jellyfish: String,
    bg_region: String,
    threads: u8,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Input::None,
            database: "".to_string(),
            kmer_size: 25,
            jellyfish: "jellyfish".to_string(),
            bg_region: "chr17:72062001-76562000".to_string(),
            threads: 4,
        }
    }
}

fn print_help() {
    let defaults = Args::default();
    println!("{}", "Create an empty database of complex loci.\n".bold().bright_yellow());
    println!("{} {} create (-f genome.fa | -a col.agc[@sample]) -d out_dir [arguments]",
        "Usage:".bold().red(), env!("CARGO_PKG_NAME"));

    println!("\n{}", "Input/output arguments:".bold().red());
    println!("    {:<15} {:<10}  Input FASTA file.",
        "-f, --fasta".green(), "FILE".cyan());
    println!("    {:<15} {:<10}  Genome collection in AGC format.\n\
        {:<30}  Optional sample name after '@' (use first sample by default).",
        "-a, --agc".green(), "FILE[@STR]".cyan(), "");
    println!("    {:<15} {:<10}  Output database directory.", "-d, --db".green(), "DIR".cyan());

    println!("\n{}", "Counting k-mers:".bold().red());
    println!("    {:<15} {:<10}  k-mer size [{}].",
        "-k, --kmer".green(), "INT".cyan(), defaults.kmer_size);
    println!("    {:<15} {:<10}  Jellyfish executable [{}].",
        "-j, --jellyfish".green(), "EXE".cyan(), defaults.jellyfish);

    println!("\n{}", "Background regions:".bold().red());
    println!("    {:<15} {:<10}  Calculate background distributions based on reads,\n\
        {:<30}  mapped to this region [{}].",
        "-r, --region".green(), "REGION".cyan(), "", defaults.bg_region);

    println!("\n{}", "Execution parameters:".bold().red());
    println!("    {:<15} {:<10}  Number of threads [{}].",
        "-@, --threads".green(), "INT".cyan(), defaults.threads);

    println!("\n{}", "Other parameters:".bold().red());
    println!("    {:<15} {:<10}  Show this help message.", "-h, --help".green(), "");
    println!("    {:<15} {:<10}  Show version.", "-V, --version".green(), "");
}

fn parse_args(argv: &[String]) -> Result<Args, lexopt::Error> {
    use lexopt::prelude::*;
    let mut args = Args::default();
    let mut parser = lexopt::Parser::from_args(argv);

    while let Some(arg) = parser.next()? {
        match arg {
            Short('f') | Long("fasta") => {
                if args.input != Input::None {
                    Err("Input files provided several times!")?;
                }
                args.input = Input::Fasta(parser.value()?.parse()?);
            }
            Short('a') | Long("agc") => {
                if args.input != Input::None {
                    Err("Input files provided several times!")?;
                }
                let val: String = parser.value()?.parse()?;
                args.input = if let Some((first, second)) = val.split_once('@') {
                    Input::Agc(first.to_string(), second.to_string())
                } else {
                    Input::Agc(val, "".to_string())
                };
            }
            Short('d') | Long("database") => {
                args.database = parser.value()?.parse()?;
            }

            Short('k') | Long("kmer") => {
                args.kmer_size = parser.value()?.parse()?;
            }
            Short('j') | Long("jellyfish") => {
                args.jellyfish = parser.value()?.parse()?;
            }
            Short('r') | Long("region") => {
                args.bg_region = parser.value()?.parse()?;
            }
            Short('@') | Long("threads") => {
                args.threads = parser.value()?.parse()?;
            }

            Short('V') | Long("version") => {
                super::print_version();
                std::process::exit(0);
            }
            Short('h') | Long("help") => {
                print_help();
                std::process::exit(0);
            }
            _ => return Err(arg.unexpected()),
        }
    }
    process_args(args)
}

fn process_args(mut args: Args) -> Result<Args, lexopt::Error> {
    if args.threads == 0 {
        args.threads = 1;
    }
    if args.input == Input::None {
        Err("Please provide input file with --fasta OR --agc.")?;
    }
    if args.database.is_empty() {
        Err("Please provide database directory with --db.")?;
    }
    Ok(args)
}

pub(super) fn run(argv: &[String]) {
    parse_args(argv).unwrap();
}