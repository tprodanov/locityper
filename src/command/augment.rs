use std::{
    time::Instant,
    path::PathBuf,
};
use colored::Colorize;
use crate::{
    err::{error, add_path, validate_param},
    ext::{
        self,
        fmt::PrettyU32,
    },
    seq::{
        dist, wfa,
        kmers::Kmer,
    },
    algo::HashSet,
};

struct Args {
    database: Option<PathBuf>,
    subset_loci: HashSet<String>,
    ref_name: Option<String>,

    skip_basis: bool,
    basis_tag: Option<String>,
    make_default: bool,
    basis_div: f64,
    div_window: u32,
    basis_leaveout: Vec<String>,
    basis_iters: usize,

    threads: u16,
    aln_params: dist::Params,
    rerun: super::Rerun,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            database: None,
            subset_loci: Default::default(),
            ref_name: None,

            skip_basis: false,
            basis_tag: None,
            make_default: true,
            basis_div: 0.01,
            div_window: 1000,
            basis_leaveout: Vec::new(),
            basis_iters: 100,

            threads: 8,
            aln_params: Default::default(),
            rerun: super::Rerun::None,
        }
    }
}

impl Args {
    fn validate(self) -> crate::Result<Self> {
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 18;
    const VAL: usize = 4;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{} For each locus:", "Augment database.".yellow());
    println!("  - constructs pairwise haplotype alignments (locityper align),");
    // println!("  - constructs local VCF file (locityper paf-vcf),");
    println!("  - extracts basis haplotypes for faster read-to-haplotype alignment.");
    println!("Multiple instances can be run over the same output directory at the same time.");

    print!("\n{}", "Usage:".bold());
    println!(" {} augment -d db -n ref_name", super::PROGRAM);

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Database with loci haplotypes.",
        "-d, --database".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Limit augmentation to these loci.",
        "    --subset-loci".green(), "STR+".yellow());

    println!("\n{}", "Alignment parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Reference haplotype name.\n\
        {EMPTY}  Alignments to this haplotype will be constructed in any case.",
        "-n, --ref-name".green(), "STR".yellow());
    println!("    {}  {} (k,w)-minimizers for sequence divergence calculation [{} {}].",
        "-m, --minimizer".green(), "INT INT".yellow(),
        super::fmt_def(defaults.aln_params.div_k), super::fmt_def(defaults.aln_params.div_w));
    println!("    {:KEY$} {:VAL$}  Do not align sequences with minimizer divergence >= {} [{}].\n\
        {EMPTY}  Use {} to align everything.",
        "-D, --thresh-div".green(), "NUM".yellow(), "NUM".yellow(), super::fmt_def_f64(defaults.aln_params.thresh_div),
        "-D 1".green());
    println!("    {:KEY$} {:VAL$}  One or more k-mer sizes (5 <= k <= {}) for backbone alignment,\n\
        {EMPTY}  separated by comma [{}].",
        "-k, --backbone".green(), "INT".yellow(), ruint::aliases::U256::MAX_KMER_SIZE,
        super::fmt_def(defaults.aln_params.backbone_str()));
    println!("    {:KEY$} {:VAL$}  Do not complete gaps over this size [{}].",
        "-g, --max-gap".green(), "INT".yellow(), super::fmt_def(PrettyU32(defaults.aln_params.max_gap)));
    println!("    {:KEY$} {:VAL$}  Alignment accuracy level (1-{}) [{}].",
        "-a, --accuracy".green(), "INT".yellow(), wfa::MAX_ACCURACY, super::fmt_def(defaults.aln_params.accuracy));
    println!("    {:KEY$} {:VAL$}  Penalty for mismatch [{}].",
        "-M, --mismatch".green(), "INT".yellow(), super::fmt_def(defaults.aln_params.penalties.mismatch));
    println!("    {:KEY$} {:VAL$}  Gap open penalty [{}].",
        "-O, --gap-open".green(), "INT".yellow(), super::fmt_def(defaults.aln_params.penalties.gap_open));
    println!("    {:KEY$} {:VAL$}  Gap extend penalty [{}].",
        "-E, --gap-extend".green(), "INT".yellow(), super::fmt_def(defaults.aln_params.penalties.gap_extend));

    println!("\n{} (used for faster read mapping)", "Basis haplotypes selection:".bold());
    println!("    {:KEY$} {:VAL$}  Skip basis haplotypes selection.",
        "    --skip-basis".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Custom tag for the basis haplotypes (haplotypes-basis.{}.fa.gz).",
        "-t, --tag".green(), "STR".yellow(), "STR".yellow());
    println!("    {:KEY$} {:VAL$}  Do not make these basis haplotypes default.",
        "    --not-default".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Maximum seq. divergence from the basis haplotypes [{}].",
        "    --basis-div".green(), "NUM".yellow(), super::fmt_def_f64(defaults.basis_div));
    println!("    {:KEY$} {:VAL$}  Calculate divergence across {} bp moving windows [{}].\n\
        {EMPTY}  Use \"inf\" for global divergence over the whole alignment.",
        "    --div-window".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(PrettyU32(defaults.div_window)));
    println!("    {:KEY$} {:VAL$}  Remove sample(s) from the basis haplotypes.\n\
        {EMPTY}  Removed haplotypes can be replaced by identical haplotypes with another name.\n\
        {EMPTY}  Needed for subsequent {}.",
        "    --basis-lo".green(), "STR+".yellow(), "locityper genotype --lo".underline());
    println!("    {:KEY$} {:VAL$}  Find the smallest basis set in {} iterations [{}].",
        "    --basis-iters".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(defaults.basis_iters));

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Rerun everything ({}); do not rerun haplotype alignment ({});\n\
        {EMPTY}  or do not rerun completed loci ({}, default).",
        "    --rerun".green(), "STR".yellow(), "all".yellow(), "part".yellow(), "none".yellow());
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
            Short('d') | Long("db") | Long("database") => args.database = Some(parser.value()?.parse()?),
            Long("subset-loci") | Long("loci-subset") => {
                for val in parser.values()? {
                    args.subset_loci.insert(val.parse()?);
                }
            }

            Short('n') | Long("ref-name") => args.ref_name = Some(parser.value()?.parse()?),
            Short('m') | Long("minimizer") | Long("minimizers") =>
            {
                args.aln_params.div_k = parser.value()?.parse()?;
                args.aln_params.div_w = parser.value()?.parse()?;
            }
            Short('D') | Long("thresh-div") => args.aln_params.thresh_div = parser.value()?.parse()?,
            Short('k') | Long("backbone") | Long("backbone-ks") => {
                let backbone_str: String = parser.value()?.parse()?;
                args.aln_params.backbone_ks = backbone_str.split(',').map(str::parse)
                    .collect::<Result<Vec<u8>, _>>()
                    .map_err(|_| error!(InvalidInput,
                    "Cannot parse `-k {}`: must be list of integers separated by comma", backbone_str))?;
            }
            Short('g') | Long("max-gap") => args.aln_params.max_gap = parser.value()?.parse::<PrettyU32>()?.get(),
            Short('a') | Long("accuracy") => args.aln_params.accuracy = parser.value()?.parse()?,
            Short('M') | Long("mismatch") => args.aln_params.penalties.mismatch = parser.value()?.parse()?,
            Short('O') | Long("gap-open") | Long("gap-opening") =>
                args.aln_params.penalties.gap_open = parser.value()?.parse()?,
            Short('E') | Long("gap-extend") | Long("gap-extension") =>
                args.aln_params.penalties.gap_extend = parser.value()?.parse()?,

            Long("skip-basis") => args.skip_basis = true,
            Short('t') | Long("tag") => args.basis_tag = Some(parser.value()?.parse()?),
            Long("not-default") => args.make_default = false,
            Long("basis-div") => args.basis_div = parser.value()?.parse()?,
            Long("div-window") => args.div_window = parser.value()?.parse::<PrettyU32>()?.get(),
            Long("basis-lo") => {
                for val in parser.values()? {
                    args.basis_leaveout.push(val.parse()?);
                }
            }
            Long("basis-iters") => args.basis_iters = parser.value()?.parse()?,

            Long("rerun") => args.rerun = parser.value()?.parse()?,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,

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
