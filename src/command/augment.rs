use std::{
    time::Instant,
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
};

enum HaplotypeSelection {
    Dominating,
    Clusters,
}

struct Args {
    tag: Option<String>,
    make_default: bool,
    haplotype_selection: HaplotypeSelection,
    max_div: f64,
    div_window: u32,
    threads: u16,

    align_params: dist::Params,
    rerun: super::Rerun,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            tag: None,
            make_default: true,
            haplotype_selection: HaplotypeSelection::Dominating,
            max_div: 0.01,
            div_window: 1000,
            threads: 8,

            align_params: Default::default(),
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
    println!("  - constructs local VCF file (locityper paf-vcf),");
    println!("  - extracts basis haplotypes for faster read-to-haplotype alignment.");
    println!("Multiple instances can be run over the same output directory at the same time.");

    print!("\n{}", "Usage:".bold());
    println!(" {} [TODO]", super::PROGRAM);

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Database with loci haplotypes.",
        "-d, --database".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Limit augmentation to these loci.",
        "    --subset-loci".green(), "STR+".yellow());

    println!("\n{}", "Alignment parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Reference haplotype name. Needed for the VCF construction.",
        "    --ref-name".green(), "STR".yellow());
    println!("    {}  {} (k,w)-minimizers for sequence divergence calculation [{} {}].",
        "-m, --minimizer".green(), "INT INT".yellow(),
        super::fmt_def(defaults.align_params.div_k), super::fmt_def(defaults.align_params.div_w));
    println!("    {:KEY$} {:VAL$}  Do not align sequences with minimizer divergence >= {} [{}].\n\
        {EMPTY}  Use {} to align everything.",
        "-D, --thresh-div".green(), "NUM".yellow(), "NUM".yellow(), super::fmt_def_f64(defaults.align_params.thresh_div),
        "-D 1".green());
    println!("    {:KEY$} {:VAL$}  One or more k-mer sizes (5 <= k <= {}) for backbone alignment,\n\
        {EMPTY}  separated by comma [{}].",
        "-k, --backbone".green(), "INT".yellow(), ruint::aliases::U256::MAX_KMER_SIZE,
        super::fmt_def(defaults.align_params.backbone_str()));
    println!("    {:KEY$} {:VAL$}  Do not complete gaps over this size [{}].",
        "-g, --max-gap".green(), "INT".yellow(), super::fmt_def(PrettyU32(defaults.align_params.max_gap)));
    println!("    {:KEY$} {:VAL$}  Alignment accuracy level (1-{}) [{}].",
        "-a, --accuracy".green(), "INT".yellow(), wfa::MAX_ACCURACY, super::fmt_def(defaults.align_params.accuracy));
    println!("    {:KEY$} {:VAL$}  Penalty for mismatch [{}].",
        "-M, --mismatch".green(), "INT".yellow(), super::fmt_def(defaults.align_params.penalties.mismatch));
    println!("    {:KEY$} {:VAL$}  Gap open penalty [{}].",
        "-O, --gap-open".green(), "INT".yellow(), super::fmt_def(defaults.align_params.penalties.gap_open));
    println!("    {:KEY$} {:VAL$}  Gap extend penalty [{}].",
        "-E, --gap-extend".green(), "INT".yellow(), super::fmt_def(defaults.align_params.penalties.gap_extend));

    println!("\n{} (used for faster read mapping)", "Basis haplotypes selection:".bold());
    println!("    {:KEY$} {:VAL$}  Skip basis haplotypes selection.",
        "    --skip-basis".green(), super::flag())
    println!("    {:KEY$} {:VAL$}  Custom tag for the basis haplotypes (haplotypes-{}.fa.gz).",
        "-t, --tag".green(), "STR".yellow(), "STR".yellow());
    println!("    {:KEY$} {:VAL$}  Do not make these basis haplotypes default.",
        "    --not-default".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Selection of basis haplotypes:\n\
        {EMPTY}  {} set ({}; default) or hierarchical {} ({}).",
        "    --sel-basis".green(), "STR".yellow(),
        "dominating".yellow(), "d".yellow(), "clusters".yellow(), "c".yellow());
    println!("    {:KEY$} {:VAL$}  Maximum sequence divergence between collapsed haplotypes [{}].",
        "    --collapse-div".green(), "NUM".yellow(), super::fmt_def_f64(defaults.max_div));
    println!("    {:KEY$} {:VAL$}  Calculate divergence across {} bp moving windows [{}].\n\
        {EMPTY}  Use \"inf\" for global divergence over the whole alignment.",
        "    --div-window".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(PrettyU32(defaults.div_window)));
    // println!("    {:KEY$} {:VAL$}  Only for ")

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
