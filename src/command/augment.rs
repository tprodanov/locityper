use std::{
    fs,
    time::Instant,
    path::{Path, PathBuf},
    fmt::Write as FmtWrite,
    cmp::{min, max},
};
use colored::Colorize;
use crate::{
    err::{error, add_path, validate_param},
    ext::{
        self,
        fmt::PrettyU32,
    },
    seq::{
        dist, wfa, fastx,
        contigs::{ContigId, ContigSet},
        kmers::Kmer,
    },
    algo::{self, HashSet, IntSet, TwoU32},
};
use super::{paths, Rerun};

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

    threads: u16,
    aln_params: dist::Params,
    rerun: Rerun,
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

            threads: 8,
            aln_params: Default::default(),
            rerun: Rerun::None,
        }
    }
}

impl Args {
    fn validate(self) -> crate::Result<Self> {
        validate_param!(self.database.is_some(), "Input database is not provided (see -i/--input)");
        validate_param!(self.ref_name.is_some(), "Reference haplotype name must be provided");
        validate_param!(self.basis_div >= 0.0 && self.basis_div <= 1.0,
            "Basis divergence ({}) must be between 0 and 1", self.basis_div);
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
            Short('d') | Long("database") | Short('i') | Long("input") => args.database = Some(parser.value()?.parse()?),
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

/// Construct basis tag as "d{div}-{window}" with suffix "-loNAME,NAME,NAME" if necessary.
fn construct_basis_tag(args: &Args) -> crate::Result<String> {
    let mut tag = format!("d{}", crate::math::fmt_signif(args.basis_div, 5));
    if args.div_window == u32::MAX {
        write!(tag, "-global").unwrap();
    } else {
        write!(tag, "-{}", ext::fmt::PrettyU32(args.div_window)).unwrap();
    }
    if !args.basis_leaveout.is_empty() {
        write!(tag, "-lo{}", args.basis_leaveout.join(",")).unwrap();
    }
    if tag.len() >= 128 {
        Err(error!(
            RuntimeError, "Automatic tag name is too long ({} chars.), please provide tag using --tag", tag.len()))
    } else {
        Ok(tag)
    }
}

fn construct_dominant_set(
    locus: &str,
    dir: &Path,
    args: &Args,
    tag: &str,
) -> crate::Result<bool> {
    let basis_filename = dir.join(format!("haplotypes-basis.{}.fa.gz", tag));
    if basis_filename.exists() {
        if args.rerun == Rerun::Part || args.rerun == Rerun::All {
            log::debug!("    Overwritting basis file `{}`", ext::fmt::path(&basis_filename));
        } else {
            log::debug!("    Basis file `{}` already exists, skipping it", ext::fmt::path(&basis_filename));
            return Ok(false);
        }
    }
    let fasta_fname = dir.join(paths::LOCUS_FASTA);
    let contig_set = ContigSet::load(locus, &fasta_fname)?;
    let contigs = contig_set.contigs();

    // [TODO] Leave-out

    let mut adjacencies = vec![Vec::new(); contig_set.len()];
    let paf_filename = dir.join(paths::LOCUS_PAF);
    let mut paf_file = ext::sys::open(paf_filename).map(dist::PafFile::new)?;
    let max_window_edit = (f64::from(args.div_window) * args.basis_div).floor() as u32;
    let mut pairs = IntSet::default();
    let mut n_edges = 0;
    while let Some(entry) = paf_file.next().transpose()? {
        let Some(i) = contigs.try_get_id(entry.query_name()) else { continue };
        let Some(j) = contigs.try_get_id(entry.target_name()) else { continue };
        if !pairs.insert(TwoU32(u32::from(min(i, j).get()), u32::from(max(i, j).get()))) {
            log::warn!("Pair {} - {} appears twice in the PAF file", entry.query_name(), entry.target_name());
            continue;
        }
        // Can't do anything if there is no CIGAR.
        let Some(cigar) = entry.cigar().transpose()? else { continue };
        let aln_len = entry.aln_len()?;
        let global_edit = aln_len - entry.n_matches()?;

        let add_edge = if aln_len <= args.div_window {
            f64::from(global_edit) <= args.basis_div * f64::from(aln_len)
        } else {
            global_edit <= max_window_edit || cigar.max_local_edit(args.div_window) <= max_window_edit
        };
        if add_edge {
            n_edges += 1;
            adjacencies[i.ix()].push(j.ix());
            adjacencies[j.ix()].push(i.ix());
        }
    }

    log::info!("    In total, {} similar haplotype pairs ({:.1}%)",
        n_edges, 100.0 * f64::from(n_edges) / ext::TriangleMatrix::calc_linear_len(contig_set.len()) as f64);
    let dominant_set = algo::dom_set::find_dominating_set(&adjacencies);
    log::info!("    Identified a basis set of {}/{} haplotypes ({:.1}% reduction)",
        dominant_set.len(), contig_set.len(), 100.0 - dominant_set.len() as f64 / contig_set.len() as f64 * 100.0);
    let tmp_filename = basis_filename.with_extension(".tmp.gz");
    let mut fasta_writer = ext::sys::create_gzip(&tmp_filename)?;
    for i in dominant_set {
        let id = ContigId::new(i);
        fastx::write_fasta(&mut fasta_writer, contigs.get_name(id).as_bytes(), contig_set.get_seq(id))
            .map_err(add_path!(tmp_filename))?;
    }
    fs::rename(&tmp_filename, &basis_filename).map_err(add_path!(tmp_filename, &basis_filename))?;
    Ok(true)
}

fn process_locus(
    locus: &str,
    dir: &Path,
    args: &Args,
    tag: &str,
) -> crate::Result<bool> {
    let Some(lock_file) = ext::sys::LockFile::try_create(dir.join("augment.lock"))? else { return Ok(false) };

    let mut did_anything = false;
    // [TODO] Construct alignments.
    did_anything |= construct_dominant_set(locus, dir, args, tag)?;

    lock_file.release()?;
    Ok(did_anything)
}

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();

    let tag = match args.basis_tag.clone() {
        Some(v) => v,
        None => construct_basis_tag(&args)?,
    };

    let loci_subdirs = super::add::load_loci_subdirs(
        args.database.as_ref().expect("Database directory must be provided"))?;
    for (locus, locus_dir) in loci_subdirs {
        if !args.subset_loci.is_empty() && !args.subset_loci.contains(&locus) {
            log::trace!("Skipping locus {} (not in the subset loci)", locus);
            continue;
        }
        process_locus(&locus, &locus_dir, &args, &tag)?;
    }

    log::info!("Success! Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
