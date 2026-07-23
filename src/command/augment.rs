use std::{
    fs,
    io::BufRead,
    time::Instant,
    path::{Path, PathBuf},
    fmt::Write as FmtWrite,
    cmp::max,
    borrow::Cow,
    sync::Arc,
};
use colored::Colorize;
use crate::{
    err::{error, add_path, validate_param},
    ext::{
        self, TriangleMatrix, LazyResult,
        fmt::{PrettyU32, YesNo},
        vec::BitArray,
    },
    seq::{
        self, wfa, fastx,
        paf::{PafFile, PafParseResult},
        cigar::Cigar,
        contigs::{ContigId, ContigNames, ContigSet, DiscardedHaplotypes},
        kmers::Kmer,
    },
    algo::{self, HashSet},
    math::RoundDiv,
};
use super::{paths, Rerun};

struct Args {
    database: Option<PathBuf>,
    subset_loci: HashSet<String>,
    ref_name: Option<String>,

    divergence: f64,
    window: u32,
    step: Option<u32>,
    skip_basis: bool,
    basis_tag: Option<String>,
    make_default: Option<bool>,
    basis_leaveout: Vec<String>,

    skip_vcf: bool,
    threads: u16,
    aln_params: seq::align::Params,
    rerun: Rerun,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            database: None,
            subset_loci: Default::default(),
            ref_name: None,

            divergence: 0.01,
            window: 250,
            step: None,
            skip_basis: false,
            basis_tag: None,
            make_default: None,
            basis_leaveout: Vec::new(),

            skip_vcf: false,
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
        validate_param!(self.divergence >= 0.0 && self.divergence <= 1.0,
            "Basis divergence ({}) must be between 0 and 1", self.divergence);
        if let Some(step) = self.step {
            validate_param!(0 < step && step <= self.window,
                "Step size ({}) must be over 0 and not greater than the window size ({})", step, self.window);
        }
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
    println!("Multiple instances can be run over the same output directory at the same time,");
    println!("unless --rerun is used.");

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
    println!("    {:KEY$} {:VAL$}  Transitively speed up alignment by using existing alignments with\n\
        {EMPTY}  divergence under this value (use 0 to disable) [{}].",
        "    --tr-div".green(), "NUM".yellow(), super::fmt_def_f64(defaults.aln_params.transitive_div));
    println!("    {:KEY$} {:VAL$}  Anchor size during transitive alignment construction [{}].\n\
        {EMPTY}  Larger values increase accuracy but decrease speed.",
        "    --tr-anchor".green(), "NUM".yellow(), super::fmt_def(defaults.aln_params.transitive_anchor));
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

    println!("\n{}", "Basis haplotypes selection:".bold());
    println!("    Each original haplotype will be similar ({}) to one of the basis haplotypes", "-x".green());
    println!("    at each moving window ({}, {}). Basis haplotypes will be used for faster read mapping.",
        "-w".green(), "-s".green());
    println!("    {:KEY$} {:VAL$}  Maximum seq. divergence from the basis haplotypes [{}].",
        "-x, --divergence".green(), "NUM".yellow(), super::fmt_def_f64(defaults.divergence));
    println!("    {:KEY$} {:VAL$}  Calculate divergence across {} bp moving windows [{}].\n\
        {EMPTY}  Use \"inf\" for global divergence over the whole alignment.",
        "-w, --window".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(PrettyU32(defaults.window)));
    println!("    {:KEY$} {:VAL$}  Moving window step [half of {}].",
        "-s, --step".green(), "INT".yellow(), "-w".green());
    println!("    {:KEY$} {:VAL$}  Custom tag for the basis haplotypes (haplotypes-basis.{}.fa.gz).",
        "-t, --tag".green(), "STR".yellow(), "STR".yellow());
    println!("    {:KEY$} {:VAL$}  Make these basis haplotypes default (yes/no).\n\
        {EMPTY}  Default: {}, unless {} is used.",
        "    --default".green(), "y|n".yellow(), super::fmt_def("yes"), "--basis-lo".green());
    println!("    {:KEY$} {:VAL$}  Remove sample(s) from the basis haplotypes.\n\
        {EMPTY}  Removed haplotypes can be replaced by identical haplotypes with another name.\n\
        {EMPTY}  Needed for subsequent {}.",
        "    --basis-lo".green(), "STR+".yellow(), "locityper genotype --lo".underline());
    println!("    {:KEY$} {:VAL$}  Skip basis haplotypes selection.",
        "    --skip-basis".green(), super::flag());

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Do not construct local VCF files.",
        "    --skip-vcf".green(), super::flag());
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
            Long("tr-div") => args.aln_params.transitive_div = parser.value()?.parse()?,
            Long("tr-anchor") => args.aln_params.transitive_anchor = parser.value()?.parse()?,
            Short('g') | Long("max-gap") => args.aln_params.max_gap = parser.value()?.parse::<PrettyU32>()?.get(),
            Short('a') | Long("accuracy") => args.aln_params.accuracy = parser.value()?.parse()?,
            Short('M') | Long("mismatch") => args.aln_params.penalties.mismatch = parser.value()?.parse()?,
            Short('O') | Long("gap-open") | Long("gap-opening") =>
                args.aln_params.penalties.gap_open = parser.value()?.parse()?,
            Short('E') | Long("gap-extend") | Long("gap-extension") =>
                args.aln_params.penalties.gap_extend = parser.value()?.parse()?,

            Short('x') | Long("divergence") => args.divergence = parser.value()?.parse()?,
            Short('w') | Long("window") => args.window = parser.value()?.parse::<PrettyU32>()?.get(),
            Short('s') | Long("step") => args.step = Some(parser.value()?.parse::<PrettyU32>()?.get()),
            Long("skip-basis") => args.skip_basis = true,
            Short('t') | Long("tag") => args.basis_tag = Some(parser.value()?.parse()?),
            Long("default") => args.make_default = Some(parser.value()?.parse::<YesNo>()?.into()),
            Long("basis-lo") => {
                for val in parser.values()? {
                    args.basis_leaveout.push(val.parse()?);
                }
            }

            Long("skip-vcf") => args.skip_vcf = true,
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
    let mut tag = String::new();
    write!(tag, "x{}", crate::math::fmt_signif(args.divergence, 5)).unwrap();
    if args.window == u32::MAX {
        write!(tag, "-global").unwrap();
    } else {
        write!(tag, "-w{}", ext::fmt::PrettyU32(args.window)).unwrap();
        if let Some(step) = args.step {
            write!(tag, "-s{}", ext::fmt::PrettyU32(step)).unwrap();
        }
    };
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

fn create_symlink(path1: &Path, path2: &Path) -> crate::Result<()> {
    log::debug!("        Soft linking {} -> {}", ext::fmt::path(path1), ext::fmt::path(path2));
    // Second check is needed if the symlink exists but points to a missing file.
    if path2.exists() || path2.is_symlink() {
        std::fs::remove_file(path2).map_err(add_path!(path2))?;
    }
    std::os::unix::fs::symlink(path1.file_name().expect("File has an empty basename"), path2)
        .map_err(add_path!(path1, path2))
}

fn update_bitarray<const IN_QUERY: bool>(
    cigar: &Cigar,
    global_div: f64,
    bitarrays: &mut [BitArray],
    oth_index: usize,
    buf_indices: &mut Vec<usize>,
    step: u32,
    max_window_edit: u32,
    args: &Args,
) {
    buf_indices.clear();
    if (if IN_QUERY { cigar.query_len() } else { cigar.ref_len() }) <= args.window {
        if global_div <= args.divergence {
            buf_indices.push(0);
        }
    } else {
        cigar.locally_similar::<IN_QUERY>(args.window, step, max_window_edit, buf_indices);
    }
    for &w in buf_indices as &[usize] {
        bitarrays[w].set_true(oth_index);
    }
}

fn inner_construct_dominant_set(
    mut paf_file: PafFile<impl BufRead>,
    contigs: &ContigNames,
    args: &Args,
) -> crate::Result<Vec<usize>> {
    let max_window_edit = (f64::from(args.window) * args.divergence).floor() as u32;
    let step = args.step.unwrap_or(max(args.window >> 1, 1));
    // For each contig and each window, store similarities to other contigs.
    let mut bitarrays: Vec<_> = contigs.lengths().iter().enumerate().map(|(i, &l)| {
        let n_windows = (l - args.window).fast_ceil_div(step) as usize + 1;
        let mut bitarr = BitArray::new(contigs.len());
        bitarr.set_true(i);
        vec![bitarr; n_windows as usize]
    }).collect();

    let mut indices = Vec::new();
    while let Some(entry) = paf_file.next(contigs)? {
        let PafParseResult::Entry(entry) = entry else { continue };
        let i = entry.query_id().ix();
        let j = entry.target_id().ix();
        if i == j { continue };
        let Some(cigar) = entry.cigar() else { continue };
        let global_div = entry.divergence();
        update_bitarray::<true >(&cigar, global_div, &mut bitarrays[i], j, &mut indices, step, max_window_edit, args);
        update_bitarray::<false>(&cigar, global_div, &mut bitarrays[j], i, &mut indices, step, max_window_edit, args);
    }

    let mut constraints = HashSet::default();
    for bitarray in bitarrays.into_iter().flatten() {
        constraints.insert(bitarray);
    }
    log::debug!("        Using ILP to find a basis set over {} haplotypes and {} constraints",
        contigs.len(), constraints.len());
    Ok(algo::dom_set::find_dominating_set(contigs.len(), constraints.iter()))
}

fn construct_dominant_set(
    dir: &Path,
    contig_set: &ContigSet,
    disc_haps: &DiscardedHaplotypes,
    args: &Args,
    tag: &str,
) -> crate::Result<bool> {
    let basis_filename = dir.join(format!("haplotypes-basis.{}.fa.gz", tag));
    let default_filename = dir.join(paths::DEFAULT_BASIS_FASTA);
    if basis_filename.exists() {
        if args.rerun == Rerun::None {
            log::debug!("    Basis file `{}` already exists, skipping it", ext::fmt::path(&basis_filename));
            if args.make_default.expect("--default must not be None") {
                create_symlink(&basis_filename, &default_filename)?;
            }
            return Ok(false);
        } else {
            log::debug!("    Overwritting basis file `{}`", ext::fmt::path(&basis_filename));
        }
    }

    log::info!("    Constructing basis haplotypes set");
    let contig_set: Cow<ContigSet> = match args.basis_leaveout.is_empty() {
        true => Cow::Borrowed(contig_set),
        false => Cow::Owned(contig_set.extract_subset(&args.basis_leaveout.iter().cloned().collect(), disc_haps)?.1),
    };
    let contigs = contig_set.contigs();
    let paf_filename = dir.join(paths::LOCUS_PAF);
    let paf_file = ext::sys::open(paf_filename).map(PafFile::new)?;
    let dominant_set = inner_construct_dominant_set(paf_file, contigs, args)?;

    log::info!("        Identified a basis set of {}/{} haplotypes ({:.1}% reduction)",
        dominant_set.len(), contig_set.len(), 100.0 - dominant_set.len() as f64 / contig_set.len() as f64 * 100.0);
    let tmp_filename = basis_filename.with_extension(".tmp.gz");
    let mut fasta_writer = ext::sys::create_gzip(&tmp_filename)?;
    for i in dominant_set {
        let id = ContigId::new(i);
        fastx::write_fasta(&mut fasta_writer, contigs.get_name(id).as_bytes(), contig_set.get_seq(id))
            .map_err(add_path!(tmp_filename))?;
    }
    fs::rename(&tmp_filename, &basis_filename).map_err(add_path!(tmp_filename, &basis_filename))?;

    if args.make_default.expect("--default must not be None") {
        create_symlink(&basis_filename, &default_filename)?;
    }
    Ok(true)
}

fn find_ref_haplotype(
    contigs: &ContigNames,
    disc_haps: &DiscardedHaplotypes,
    ref_name: &str,
) -> crate::Result<Option<ContigId>> {
    if let Some(id) = contigs.try_get_id(ref_name) {
        return Ok(Some(id));
    }
    for (&id, other_haps) in disc_haps.by_contig() {
        if let Some((_, is_identical)) = other_haps.iter().filter(|(name, _)| name == ref_name).next() {
            if !is_identical { break };
            log::debug!("Replacing reference haplotype {} with {}", ref_name, contigs.get_name(id));
            return Ok(Some(id));
        }
    }
    log::warn!("Could not find reference haplotype `{}`", ref_name);
    Ok(None)
}

fn process_locus(
    locus: &str,
    dir: &Path,
    args: &Args,
    tag: &str,
) -> crate::Result<bool> {
    let lock_filename = dir.join("augment.lock");
    let Some(lock_file) = ext::sys::LockFile::try_create(lock_filename.clone())? else {
        log::debug!("Skipping locus {} (lock exists {})", locus, ext::fmt::path(&lock_filename));
        return Ok(false)
    };
    log::info!("Processing locus {}", locus);
    let mut did_anything = false;
    let fasta_filename = dir.join(paths::LOCUS_FASTA);
    // Lazily load contig set and discarded haplotypes.
    let mut lazy_data = LazyResult::new(|| -> crate::Result<_> {
        let contig_set = ContigSet::load(locus, &fasta_filename).map(Arc::new)?;
        let disc_haps = DiscardedHaplotypes::load_if_present(&dir.join(paths::DISCARDED_HAPS), contig_set.contigs())?;
        let ref_id = find_ref_haplotype(contig_set.contigs(), &disc_haps,
            args.ref_name.as_ref().expect("Reference name must be provided"))?;
        Ok((contig_set, disc_haps, ref_id))
    });

    let aln_filename = dir.join(paths::LOCUS_PAF);
    if !aln_filename.exists() || args.rerun == Rerun::All {
        let (contig_set, _, ref_id) = lazy_data.get()?;
        let pairs = TriangleMatrix::indices(contig_set.len())
            .map(|(i, j)| (ContigId::new(i), ContigId::new(j))).collect();
        log::info!("    Running pairwise haplotype alignment");
        let mut against_contig = vec![false; contig_set.len()];
        if let Some(ref_id) = ref_id {
            against_contig[ref_id.ix()] = true;
        }
        super::align::align(Arc::clone(contig_set), pairs, against_contig,
            &aln_filename, &None, args.threads, &args.aln_params)?;
        did_anything = true;
    } else {
        log::debug!("    Skipping alignments (already constructed)");
    }

    let vcf_filename = dir.join("haplotypes.vcf.gz");
    if args.skip_vcf || (vcf_filename.exists() && args.rerun == Rerun::None) {
        log::debug!("    Skipping local VCF file");
    } else {
        log::info!("    Constructing a local VCF file");
        let (contig_set, disc_haps, ref_id) = lazy_data.get()?;
        let contigs = contig_set.contigs();
        if let &Some(ref_id) = ref_id
            && let Some((chrom, shift)) = super::paf_vcf::load_region("auto", &fasta_filename, contigs.get_len(ref_id))?
        {
            super::paf_vcf::convert_to_vcf(
                &aln_filename, contig_set, disc_haps, ref_id, &chrom, shift, &vcf_filename, None)?;
            did_anything = true;
        } else {
            log::error!("Cannot construct local VCF: reference contig not found or \
                could not identify reference coordinates");
        }
    }

    if args.skip_basis {
        log::debug!("    Skipping basis construction");
    } else {
        let (contig_set, disc_haps, _) = lazy_data.get()?;
        did_anything |= construct_dominant_set(dir, contig_set, disc_haps, args, tag)?;
    }

    lock_file.release()?;
    Ok(did_anything)
}

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let mut args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();

    let tag = match args.basis_tag.clone() {
        Some(v) => v,
        None => construct_basis_tag(&args)?,
    };
    // If --default was not explicitely specified, set it to true if --basis-lo is empty.
    args.make_default = Some(args.make_default.unwrap_or(args.basis_leaveout.is_empty()));
    if !args.skip_basis {
        log::info!("Basis haplotypes will have tag \"{}\" and will{} be saved as default",
            tag, if args.make_default.unwrap() { Default::default() } else { " not".bold() });
    }

    let loci_subdirs = super::add::load_loci_subdirs(
        args.database.as_ref().expect("Database directory must be provided"))?;
    let mut completed = 0;
    let mut skipped = 0;
    let mut errors = 0;
    for (locus, locus_dir) in loci_subdirs {
        if !args.subset_loci.is_empty() && !args.subset_loci.contains(&locus) {
            log::trace!("Skipping locus {} (not in the subset loci)", locus);
            skipped += 1;
            continue;
        }
        match process_locus(&locus, &locus_dir, &args, &tag) {
            Ok(true) => completed += 1,
            Ok(false) => skipped += 1,
            Err(e) => {
                log::error!("Error in locus {}: {}", locus, e.display());
                errors += 1;
            }
        }
    }
    log::info!("Completed {} loci, skipped {}, errors in {} loci", completed, skipped, errors);
    log::info!("Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
