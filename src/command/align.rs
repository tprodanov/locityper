use std::{
    path::{Path, PathBuf},
    time::Instant,
    io::{Write, BufRead},
    ffi::OsStr,
};
use colored::Colorize;
use smallvec::SmallVec;
use rand::Rng;
use crate::{
    err::{error, add_path, validate_param},
    ext::{
        self,
        TriangleMatrix,
        fmt::PrettyU32,
    },
    seq::{
        fastx, dist, wfa,
        NamedSeq,
        kmers::Kmer,
    },
    algo::{HashMap, Hasher},
};

struct Args {
    input: Option<PathBuf>,
    output: Option<PathBuf>,
    prefix: Option<PathBuf>,

    all_pairs: bool,
    pairs: Vec<String>,
    pairs_file: Option<PathBuf>,
    threads: u16,

    params: dist::Params,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: None,
            output: None,
            prefix: None,

            all_pairs: false,
            pairs: Vec::new(),
            pairs_file: None,
            threads: 8,

            params: Default::default(),
        }
    }
}

impl Args {
    fn validate(mut self) -> crate::Result<Self> {
        validate_param!(self.input.is_some(), "Input FASTA file is not provided (see -i/--input)");
        validate_param!(self.output.is_some(), "Output PAF path is not provided (see -o/--output)");

        validate_param!(u8::from(self.all_pairs)
            + u8::from(!self.pairs.is_empty())
            + u8::from(self.pairs_file.is_some())
            == 1, "Exactly one argument -p/-P/-A is required");

        self.params.validate()?;
        Ok(self)
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
    println!("    {:KEY$} {:VAL$}  Output PAF file. In multithreaded calls, output will be unsorted\n\
        {EMPTY}  for performance purposes.",
        "-o, --output".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Prefix for temporary files. Only necessary if multiple threads\n\
        {EMPTY}  are used and the main output goes to stdout.",
        "    --prefix".green(), "PATH".yellow());

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
        super::fmt_def(defaults.params.div_k), super::fmt_def(defaults.params.div_w));
    println!("    {:KEY$} {:VAL$}  Skip divergence calculation.",
        "-s, --skip-div".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Do not align sequences with minimizer divergence >= {} [{}].\n\
        {EMPTY}  Use {} to align everything.",
        "-D, --thresh-div".green(), "NUM".yellow(), "NUM".yellow(), super::fmt_def_f64(defaults.params.thresh_div),
        "-D 1".green());
    println!("    {:KEY$} {:VAL$}  One or more k-mer sizes (5 <= k <= {}) for backbone alignment,\n\
        {EMPTY}  separated by comma [{}].",
        "-k, --backbone".green(), "INT".yellow(), ruint::aliases::U256::MAX_KMER_SIZE,
        super::fmt_def(defaults.params.backbone_str()));
    println!("    {:KEY$} {:VAL$}  Do not complete gaps over this size [{}].",
        "-g, --max-gap".green(), "INT".yellow(), super::fmt_def(PrettyU32(defaults.params.max_gap)));
    println!("    {:KEY$} {:VAL$}  Alignment accuracy level (1-{}) [{}].",
        "-a, --accuracy".green(), "INT".yellow(), wfa::MAX_ACCURACY, super::fmt_def(defaults.params.accuracy));
    println!("    {:KEY$} {:VAL$}  Penalty for mismatch [{}].",
        "-M, --mismatch".green(), "INT".yellow(), super::fmt_def(defaults.params.penalties.mismatch));
    println!("    {:KEY$} {:VAL$}  Gap open penalty [{}].",
        "-O, --gap-open".green(), "INT".yellow(), super::fmt_def(defaults.params.penalties.gap_open));
    println!("    {:KEY$} {:VAL$}  Gap extend penalty [{}].",
        "-E, --gap-extend".green(), "INT".yellow(), super::fmt_def(defaults.params.penalties.gap_extend));

    println!("\n{}", "Execution arguments:".bold());
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
            Short('i') | Long("input") => args.input = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),
            Long("prefix") => args.prefix = Some(parser.value()?.parse()?),

            Short('p') | Long("pairs") => {
                for val in parser.values()? {
                    args.pairs.push(val.parse()?);
                }
            }
            Short('P') | Long("pairs-file") => args.pairs_file = Some(parser.value()?.parse()?),
            Short('A') | Long("all") | Long("all-pairs") => args.all_pairs = true,

            Short('m') | Long("minimizer") | Long("minimizers") =>
            {
                args.params.div_k = parser.value()?.parse()?;
                args.params.div_w = parser.value()?.parse()?;
            }
            Short('s') | Long("skip-div") => args.params.skip_div = true,
            Short('D') | Long("thresh-div") => args.params.thresh_div = parser.value()?.parse()?,
            Short('k') | Long("backbone") | Long("backbone-ks") => {
                let backbone_str: String = parser.value()?.parse()?;
                args.params.backbone_ks = backbone_str.split(',').map(str::parse)
                    .collect::<Result<Vec<u8>, _>>()
                    .map_err(|_| error!(InvalidInput,
                    "Cannot parse `-k {}`: must be list of integers separated by comma", backbone_str))?;
            }
            Short('g') | Long("max-gap") => args.params.max_gap = parser.value()?.parse::<PrettyU32>()?.get(),
            Short('a') | Long("accuracy") => args.params.accuracy = parser.value()?.parse()?,
            Short('M') | Long("mismatch") => args.params.penalties.mismatch = parser.value()?.parse()?,
            Short('O') | Long("gap-open") | Long("gap-opening") =>
                args.params.penalties.gap_open = parser.value()?.parse()?,
            Short('E') | Long("gap-extend") | Long("gap-extension") =>
                args.params.penalties.gap_extend = parser.value()?.parse()?,

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

fn parse_pair(split_pair: &[&str], name2id: &HashMap<&str, u32>) -> crate::Result<(u32, u32)> {
    if split_pair.len() != 2 {
        return Err(error!(InvalidInput, "Cannot parse pair `{:?}`: exactly two names required",
            split_pair));
    }
    let id1 = *name2id.get(split_pair[0])
        .ok_or_else(|| error!(InvalidInput, "Cannot find sequence `{}`", split_pair[0]))?;
    let id2 = *name2id.get(split_pair[1])
        .ok_or_else(|| error!(InvalidInput, "Cannot find sequence `{}`", split_pair[1]))?;
    if id1 == id2 {
        return Err(error!(InvalidInput, "Cannot align sequence to itself ({})", split_pair[0]));
    }
    Ok((id2, id1))
}

fn load_pairs(args: &Args, seqs: &[NamedSeq]) -> crate::Result<Vec<(u32, u32)>> {
    if args.all_pairs {
        return Ok(TriangleMatrix::indices(seqs.len()).map(|(i, j)| (i as u32, j as u32)).collect())
    }

    let mut name2id = HashMap::<&str, u32>::with_capacity_and_hasher(seqs.len(), Hasher::default());
    for (i, entry) in seqs.iter().enumerate() {
        if name2id.insert(entry.name(), i as u32).is_some() {
            return Err(error!(InvalidInput, "Duplicate sequence {} in the input FASTA file", entry.name()));
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

struct TempFilenames(Vec<PathBuf>);

impl TempFilenames {
    fn new(
        output: &Path,
        prefix_arg: &Option<PathBuf>,
        threads: u16,
        rng: &mut impl Rng,
    ) -> crate::Result<Self>
    {
        if threads <= 1 {
            return Ok(Self(Vec::new()));
        }

        let (mut prefix, middle, extension) = if ext::sys::is_tty(output) {
            let Some(prefix) = prefix_arg.as_ref() else {
                return Err(error!(InvalidInput, "Must provide --prefix when using stdout and multiple threads"));
            };
            (prefix.to_path_buf(), format!("{:X}.", rng.random::<u32>()), "paf")
        } else {
            let prefix = match prefix_arg {
                Some(val) => val.to_path_buf(),
                None => output.with_extension(""),
            };
            let extension = match output.extension().and_then(OsStr::to_str) {
                Some("gz") => "paf.gz",
                Some("lz4") => "paf.lz4",
                _ => "paf",
            };
            (prefix, String::new(), extension)
        };

        if prefix.is_dir() {
            prefix.push("locityper-align");
        }
        Ok(Self((1..threads).map(|i| prefix.with_added_extension(format!("{}{}.{}", middle, i, extension))).collect()))
    }

    /// Skip file deletion.
    fn disarm(&mut self) {
        self.0.clear();
    }
}

impl Drop for TempFilenames {
    fn drop(&mut self) {
        for filename in &self.0 {
            if filename.exists() {
                if let Err(_) = std::fs::remove_file(filename) {
                    log::error!("Could not remove {}", filename.display());
                }
            }
        }
    }
}

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();

    let mut fasta_reader = fastx::Reader::from_path(args.input.as_ref().unwrap())?;
    let entries = fasta_reader.read_all()?;
    let pairs = load_pairs(&args, &entries)?;
    if pairs.is_empty() {
        return Err(error!(InvalidInput, "No alignments to compute"));
    }
    log::info!("Align {} pairs across {} sequences", pairs.len(), entries.len());
    let threads = usize::from(args.threads).min(pairs.len()) as u16;

    let mut files = Vec::with_capacity(usize::from(threads));
    let out_filename = args.output.as_ref().unwrap();
    files.push(ext::sys::create(out_filename)?);
    writeln!(files[0], "# minimizers={},{}; max_divergence={:.5}; backbone-ks={}; accuracy={}; max-gap={}",
        args.params.div_k, args.params.div_w, args.params.thresh_div,
        args.params.backbone_str(), args.params.accuracy, args.params.max_gap).map_err(add_path!(out_filename))?;

    let mut rng = ext::rand::init_rng(None, false);
    // Create output file now, this way we are sure it can be opened.
    let mut temp_filenames = TempFilenames::new(out_filename, &args.prefix, threads, &mut rng)?;
    for filename in &temp_filenames.0 {
        files.push(ext::sys::create(filename)?);
    }

    dist::align_sequences(entries, pairs, &args.params, threads, files, &mut rng)?;
    ext::sys::merge_files(out_filename, &temp_filenames.0)?;
    temp_filenames.disarm();

    log::info!("Success! Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
