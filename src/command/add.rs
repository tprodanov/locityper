use std::{
    io::{BufRead, Read, Seek},
    rc::Rc,
    collections::HashSet,
    path::{Path, PathBuf},
};
use bio::io::fasta::IndexedReader;
use colored::Colorize;
use const_format::str_repeat;
use crate::{
    Error,
    seq::{Interval, NamedInterval, ContigNames},
};

struct Args {
    database: Option<PathBuf>,
    reference: Option<PathBuf>,
    variants: Option<PathBuf>,

    loci: Vec<String>,
    bed_files: Vec<PathBuf>,

    jellyfish: PathBuf,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            database: None,
            reference: None,
            variants: None,

            loci: Vec::new(),
            bed_files: Vec::new(),

            jellyfish: PathBuf::from("jellyfish"),
        }
    }
}

fn print_help() {
    const KEY: usize = 15;
    const VAL: usize = 4;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Adds complex locus/loci to the database.".yellow());

    println!("\n{} {} add -d db -r reference.fa -v vars.vcf.gz -l/-L loci [arguments]",
        "Usage:".bold(), env!("CARGO_PKG_NAME"));

    println!("\n{}", "Input arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Input database directory.",
        "-d, --db".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file.",
        "-r, --reference".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  PanGenie input VCF file. Encodes variation across pangenome samples.\n\
        {EMPTY}  Must be compressed and indexed with `tabix`.",
        "-v, --vcf".green(), "FILE".yellow());

    println!("\n{}", "Complex loci coordinates:".bold());
    println!("    {:KEY$} {:VAL$}  Complex locus coordinates. Multiple loci are allowed.\n\
        {EMPTY}  Format: 'chrom:start-end' or 'chrom:start-end@name',\n\
        {EMPTY}  where 'name' is the locus name (must be unique).",
        "-l, --locus".green(), "STR+".yellow());
    println!("    {:KEY$} {:VAL$}  BED file with complex loci coordinates. May be repeated multiple times.\n\
        {EMPTY}  If fourth column is present, it is used for the locus name (must be unique).",
        "-L, --loci-bed".green(), "FILE".yellow());

    // println!("\n{}", "Optional parameters:".bold());
    // println!("    {:KEY$} {:VAL$}  Extend loci boundaries by at most {} bp,\n\
    //     {EMPTY}  in order to select the best ")

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Jellyfish executable [{}].",
        "    --jellyfish".green(), "EXE".yellow(), defaults.jellyfish.display());

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
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),
            Short('v') | Long("vcf") => args.variants = Some(parser.value()?.parse()?),

            Short('l') | Long("locus") =>
                args.loci.extend(parser.values()?.map(ValueExt::string).collect::<Result<Vec<_>, _>>()?),
            Short('L') | Long("loci") | Long("loci-bed") => args.bed_files.push(parser.value()?.parse()?),

            Short('j') | Long("jellyfish") => args.jellyfish = parser.value()?.parse()?,

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
    if args.database.is_none() {
        Err(lexopt::Error::from("Database directory is not provided (see -d/--database)"))?;
    }
    if args.reference.is_none() {
        Err(lexopt::Error::from("Reference fasta file is not provided (see -r/--reference)"))?;
    }
    if args.variants.is_none() {
        Err(lexopt::Error::from("Variants VCF file is not provided (see -v/--vcf)"))?;
    }
    if args.loci.is_empty() && args.bed_files.is_empty() {
        Err(lexopt::Error::from("Complex loci are not provided (see -l/--locus and -L/--loci-bed)"))?;
    }
    args.jellyfish = super::find_exe(args.jellyfish)?;
    Ok(args)
}

/// Loads named interval from a list of loci and a list of BED files. Names must not repeat.
fn load_loci(contigs: &Rc<ContigNames>, loci: &[String], bed_files: &[PathBuf]) -> Result<Vec<NamedInterval>, Error> {
    let mut intervals = Vec::new();
    let mut names = HashSet::new();
    for locus in loci.iter() {
        let interval = NamedInterval::parse(locus, contigs)?;
        if !names.insert(interval.name().to_string()) {
            return Err(Error::InvalidInput(format!("Locus name '{}' appears at least twice", interval.name())));
        }
        intervals.push(interval);
    }

    for filename in bed_files.iter() {
        for line in super::common::open(filename)?.lines() {
            let interval = NamedInterval::parse_bed(&mut line?.split('\t'), contigs)?;
            if !names.insert(interval.name().to_string()) {
                return Err(Error::InvalidInput(format!("Locus name '{}' appears at least twice", interval.name())));
            }
            intervals.push(interval);
        }
    }
    Ok(intervals)
}

fn extract_locus<R>(contigs: &ContigNames, fasta: &mut IndexedReader<R>, locus: &NamedInterval) -> Result<(), Error>
where R: Read + Seek,
{
    log::info!("Analyzing locus {}", locus);
    Ok(())
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args = process_args(parse_args(argv)?)?;

    // unwrap as args.reference was previously checked to be Some.
    let ref_filename = args.reference.as_ref().unwrap();
    let (contigs, mut fasta) = ContigNames::load_indexed_fasta(&ref_filename, "reference".to_string())?;
    let loci = load_loci(&contigs, &args.loci, &args.bed_files)?;
    for locus in loci.iter() {
        extract_locus(&contigs, &mut fasta, locus)?;
    }

    Ok(())
}