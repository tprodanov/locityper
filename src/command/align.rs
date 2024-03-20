use std::{
    path::PathBuf,
    time::Instant,
    io::{self, Write, BufRead},
    fmt::Write as FmtWrite,
};
use crate::{
    err::{Error, add_path, validate_param},
    ext::{
        self,
        TriangleMatrix,
        fmt::PrettyU32,
    },
    seq::{
        fastx, div, dist,
        NamedSeq,
        cigar::{Cigar, Operation},
        wfa::{self, Penalties},
        kmers::{self, Kmer},
    },
    algo::{HashMap, Hasher},
};
use colored::Colorize;
use smallvec::SmallVec;

struct Args {
    input: Option<PathBuf>,
    output: Option<PathBuf>,

    all_pairs: bool,
    pairs: Vec<String>,
    pairs_file: Option<PathBuf>,
    threads: u16,

    skip_div: bool,
    div_k: u8,
    div_w: u8,
    max_div: f64,

    penalties: Penalties,
    backbone_ks: String,
    accuracy: u8,
    max_gap: u32,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: None,
            output: None,

            all_pairs: false,
            pairs: Vec::new(),
            pairs_file: None,
            threads: 8,

            skip_div: false,
            div_k: 25,
            div_w: 15,
            max_div: 0.5,

            penalties: Default::default(),
            backbone_ks: "25,51,101".to_string(),
            accuracy: 9,
            max_gap: 500,
        }
    }
}

impl Args {
    fn validate(self) -> Result<Self, Error> {
        validate_param!(self.input.is_some(), "Input FASTA file is not provided (see -i/--input)");
        validate_param!(self.output.is_some(), "Output PAF path is not provided (see -o/--output)");

        validate_param!(u8::from(self.all_pairs)
            + u8::from(!self.pairs.is_empty())
            + u8::from(self.pairs_file.is_some())
            == 1, "Exactly one argument -p/-P/-A is required");

        validate_param!(0 < self.div_k && self.div_k <= u64::MAX_KMER_SIZE,
            "k-mer size ({}) must be between 1 and {}", self.div_k, u64::MAX_KMER_SIZE);
        validate_param!(0 < self.div_w && self.div_w <= kmers::MAX_MINIMIZER_W,
            "Minimizer window ({}) must be between 1 and {}", self.div_w, kmers::MAX_MINIMIZER_W);
        validate_param!(0.0 <= self.max_div && self.max_div <= 1.0,
            "Maximum divergence ({}) must be within [0, 1]", self.max_div);
        validate_param!(1 <= self.accuracy && self.accuracy <= wfa::MAX_ACCURACY,
            "Alignment accuracy level ({}) must be between 0 and {}.", self.accuracy, wfa::MAX_ACCURACY);

        self.penalties.validate()?;
        Ok(self)
    }

    fn get_backbones(&self) -> Result<Vec<u32>, Error> {
        let backbones = self.backbone_ks.split(',').map(str::parse)
            .collect::<Result<Vec<u32>, _>>()
            .map_err(|_| Error::InvalidInput(
                format!("Cannot parse `-k {}`: must be list of integers separated by comma", self.backbone_ks)))?;
        validate_param!(!backbones.is_empty(), "Expect at least one backbone k-mer");
        validate_param!(backbones.iter().all(|&k| k >= 5),
            "Backbone k-mer sizes must be at least 5 ({:?})", backbones);
        Ok(backbones)
    }
}

fn print_help() {
    const KEY: usize = 16;
    const VAL: usize = 5;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Align medium-size sequence to each other.".yellow());

    print!("\n{}", "Usage:".bold());
    println!(" {} -i input.fa -o out.paf (-p name,name | -P pairs.txt | -A) [args]", super::PROGRAM);

    println!("\n{}", "Input arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Input FASTA file.",
        "-i, --input".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Output PAF file.",
        "-o, --output".green(), "FILE".yellow());

    println!("\n{} (mutually exclusive):", "Alignment pairs".bold());
    println!("    {:KEY$} {:VAL$}  Find alignments for these pairs:\n\
        {EMPTY}  comma-separated sequence names, multiple allowed.",
        "-p, --pairs".green(), "STR+".yellow());
    println!("    {:KEY$} {:VAL$}  Find alignments for these pairs (two column file).",
        "-P, --pairs-file".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Find alignments for all pairs.",
        "-A, --all".green(), super::flag());

    println!("\n{}", "Alignment arguments:".bold());
    println!("    {} {} (k,w)-minimizers for sequence divergence calculation [{} {}].",
        "-m, --minimizer".green(), "INT INT".yellow(),
        super::fmt_def(defaults.div_k), super::fmt_def(defaults.div_w));
    println!("    {:KEY$} {:VAL$}  Skip divergence calculation.",
        "    --skip-div".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Do not align sequences with bigger divergence than this [{}].",
        "-D, --max-div".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.max_div));
    println!("    {:KEY$} {:VAL$}  One or more k-mer size for backbone alignment,\n\
        {EMPTY}  separated by comma [{}].",
        "-k, --backbone".green(), "INT".yellow(), super::fmt_def(&defaults.backbone_ks));
    println!("    {:KEY$} {:VAL$}  Do not complete gaps over this size [{}].",
        "-g, --max-gap".green(), "INT".yellow(), super::fmt_def(PrettyU32(defaults.max_gap)));
    println!("    {:KEY$} {:VAL$}  Alignment accuracy level (1-{}) [{}].",
        "-a, --accuracy".green(), "INT".yellow(), wfa::MAX_ACCURACY, super::fmt_def(defaults.accuracy));
    println!("    {:KEY$} {:VAL$}  Penalty for mismatch [{}].",
        "-M, --mismatch".green(), "INT".yellow(), super::fmt_def(defaults.penalties.mismatch));
    println!("    {:KEY$} {:VAL$}  Gap open penalty [{}].",
        "-O, --gap-open".green(), "INT".yellow(), super::fmt_def(defaults.penalties.gap_open));
    println!("    {:KEY$} {:VAL$}  Gap extend penalty [{}].",
        "-E, --gap-extend".green(), "INT".yellow(), super::fmt_def(defaults.penalties.gap_extend));

    println!("\n{}", "Execution arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), super::fmt_def(defaults.threads));

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

            Short('p') | Long("pairs") => {
                for val in parser.values()? {
                    args.pairs.push(val.parse()?);
                }
            }
            Short('P') | Long("pairs-file") => args.pairs_file = Some(parser.value()?.parse()?),
            Short('A') | Long("all") | Long("all-pairs") => args.all_pairs = true,

            Short('m') | Long("minimizer") | Long("minimizers") =>
            {
                args.div_k = parser.value()?.parse()?;
                args.div_w = parser.value()?.parse()?;
            }
            Long("skip-div") => args.skip_div = true,
            Short('D') | Long("max-div") => args.max_div = parser.value()?.parse()?,
            Short('k') | Long("backbone") | Long("backbone-ks") => args.backbone_ks = parser.value()?.parse()?,
            Short('g') | Long("max-gap") => args.max_gap = parser.value()?.parse::<PrettyU32>()?.get(),
            Short('a') | Long("accuracy") => args.accuracy = parser.value()?.parse()?,
            Short('M') | Long("mismatch") => args.penalties.mismatch = parser.value()?.parse()?,
            Short('O') | Long("gap-open") | Long("gap-opening") =>
                args.penalties.gap_open = parser.value()?.parse()?,
            Short('E') | Long("gap-extend") | Long("gap-extension") =>
                args.penalties.gap_extend = parser.value()?.parse()?,

            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
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

fn parse_pair(split_pair: &[&str], name2id: &HashMap<&str, u32>) -> Result<(u32, u32), Error> {
    if split_pair.len() != 2 {
        return Err(Error::InvalidInput(format!("Cannot parse pair `{:?}`: exactly two names required",
            split_pair)));
    }
    let id1 = *name2id.get(split_pair[0]).ok_or_else(|| Error::InvalidInput(
        format!("Cannot find sequence `{}`", split_pair[0])))?;
    let id2 = *name2id.get(split_pair[1]).ok_or_else(|| Error::InvalidInput(
        format!("Cannot find sequence `{}`", split_pair[1])))?;
    if id1 == id2 {
        return Err(Error::InvalidInput(
            format!("Cannot align sequence to itself ({})", split_pair[0])));
    }
    Ok((id2, id1))
}

fn load_pairs(args: &Args, seqs: &[NamedSeq]) -> Result<Vec<(u32, u32)>, Error> {
    if args.all_pairs {
        return Ok(TriangleMatrix::indices(seqs.len()).map(|(i, j)| (i as u32, j as u32)).collect())
    }

    let mut name2id = HashMap::<&str, u32>::with_capacity_and_hasher(seqs.len(), Hasher::default());
    for (i, entry) in seqs.iter().enumerate() {
        if name2id.insert(entry.name(), i as u32).is_some() {
            return Err(Error::InvalidInput(
                format!("Duplicate sequence {} in the input FASTA file", entry.name())));
        }
    }
    let mut pairs = Vec::new();
    for pair in args.pairs.iter() {
        let split: SmallVec<[&str; 2]> = pair.split(',').collect();
        pairs.push(parse_pair(&split, &name2id)?);
    }

    if let Some(path) = &args.pairs_file {
        let f = ext::sys::open(&path)?;
        for line in f.lines() {
            let line = line.map_err(add_path!(path))?;
            if line.starts_with('#') {
                continue;
            }
            let split: SmallVec<[&str; 2]> = line.split_whitespace().collect();
            pairs.push(parse_pair(&split, &name2id)?);
        }
    }
    Ok(pairs)
}

fn write_paf(
    seqs: &[NamedSeq],
    pairs: &[(u32, u32)],
    calc_aln: &[bool],
    divergences: &Option<Vec<(u32, f64)>>,
    alignments: &[(Cigar, i32)],
    args: &Args,
    mut f: impl Write,
) -> io::Result<()>
{
    writeln!(f, "# minimizers={},{}; max_divergence={:.5}; backbone-ks={}; accuracy={}; max-gap={}",
        args.div_k, args.div_w, args.max_div, args.backbone_ks, args.accuracy, args.max_gap)?;

    let mut alns_iter = alignments.iter();
    let mut cigar_str = String::new();
    for (k, (&(i, j), &has_aln)) in pairs.iter().zip(calc_aln).enumerate() {
        let qentry = &seqs[j as usize];
        let rentry = &seqs[i as usize];
        write!(f, "{}\t{len}\t0\t{len}\t+\t", qentry.name(), len = qentry.seq().len())?;
        write!(f, "{}\t{len}\t0\t{len}\t+\t", rentry.name(), len = rentry.seq().len())?;

        cigar_str.clear();
        if has_aln {
            let (cigar, score) = alns_iter.next().expect("Too few alignments");
            let mut nmatches = 0;
            let mut nerrs = 0;
            for item in cigar.iter() {
                match item.operation() {
                    Operation::Equal => nmatches += item.len(),
                    _ => nerrs += item.len(),
                }
            }
            let aln_len = nmatches + nerrs;
            write!(f, "{}\t{}\t255", nmatches, aln_len)?;
            let dv = f64::from(nerrs) / f64::from(aln_len);
            write!(f, "\tNM:i:{}\tAS:i:{}\tdv:f:{:.9}", nerrs, score, dv)?;
            write!(cigar_str, "{}", cigar).unwrap();
        } else {
            write!(f, "0\t0\t255")?;
        }
        if let Some(divs) = divergences {
            let (uniq_minims, minim_dv) = divs[k];
            write!(f, "\tum:i:{}\tmd:f:{:.9}", uniq_minims, minim_dv)?;
        }
        if !cigar_str.is_empty() {
            write!(f, "\tcg:Z:{}", cigar_str)?;
        }
        writeln!(f)?;
    }
    Ok(())
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let args = parse_args(argv)?.validate()?;
    let backbones = args.get_backbones()?;
    super::greet();
    let timer = Instant::now();

    let mut fasta_reader = fastx::Reader::from_path(args.input.as_ref().unwrap())?;
    let seqs = fasta_reader.read_all()?;
    let pairs = load_pairs(&args, &seqs)?;
    if pairs.is_empty() {
        return Err(Error::InvalidInput("No alignments to compute".to_string()));
    }
    log::info!("Align {} pairs across {} sequences", pairs.len(), seqs.len());

    // Create output file now, this way we are sure it can be opened.
    let out_filename = args.output.as_ref().unwrap();
    let out = ext::sys::create(out_filename)?;

    let divergences = if args.skip_div {
        log::debug!("    Skipping minimizer divergences");
        None
    } else {
        log::debug!("    Calculating minimizer divergences");
        Some(div::minimizer_divergences(&seqs, &pairs, args.div_k, args.div_w, args.threads))
    };
    let n_pairs = pairs.len();
    let (sub_pairs, calc_aln) = if args.skip_div || args.max_div == 1.0 {
        (pairs.clone(), vec![true; n_pairs])
    } else {
        let mut sub_pairs = Vec::with_capacity(n_pairs);
        let mut calc_aln = Vec::with_capacity(n_pairs);
        for (&pair, div) in pairs.iter().zip(divergences.as_ref().unwrap()) {
            if div.1 <= args.max_div {
                sub_pairs.push(pair);
                calc_aln.push(true);
            } else {
                calc_aln.push(false);
            }
        }
        (sub_pairs, calc_aln)
    };
    let alns = if sub_pairs.is_empty() {
        log::warn!("    Skipping all alignments (divergence too high)");
        Vec::new()
    } else {
        log::debug!("    Find alignments for {} pairs", sub_pairs.len());
        dist::align_sequences(&seqs, sub_pairs, &args.penalties, &backbones, args.accuracy,
            args.max_gap, args.threads)?
    };

    write_paf(&seqs, &pairs, &calc_aln, &divergences, &alns, &args, out).map_err(add_path!(out_filename))?;
    log::info!("Success! Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
