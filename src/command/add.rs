use std::{
    cmp::max,
    io::{BufRead, Read, Seek, Write, BufWriter},
    fs::{self, File},
    rc::Rc,
    process::Command,
    collections::HashSet,
    path::{Path, PathBuf},
};
use bio::io::fasta;
use htslib::bcf::{
    self,
    Read as VcfRead,
    record::Record as VcfRecord,
};
use colored::Colorize;
use crate::{
    Error,
    algo::{
        bisect,
        vec_ext::{VecExt, F64Ext, IterExt},
    },
    seq::{
        self, NamedInterval, ContigNames,
        kmers::JfKmerGetter,
    },
};
use super::common;

struct Args {
    database: Option<PathBuf>,
    reference: Option<PathBuf>,
    variants: Option<PathBuf>,

    loci: Vec<String>,
    bed_files: Vec<PathBuf>,

    ref_name: Option<String>,
    max_expansion: u32,
    moving_window: u32,
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

            ref_name: None,
            max_expansion: 2000,
            moving_window: 100,
            jellyfish: PathBuf::from("jellyfish"),
        }
    }
}

fn print_help() {
    const KEY: usize = 15;
    const VAL: usize = 4;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

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
        {EMPTY}  Format: 'chrom:start-end' or 'chrom:start-end=name',\n\
        {EMPTY}  where 'name' is the locus name (must be unique).",
        "-l, --locus".green(), "STR+".yellow());
    println!("    {:KEY$} {:VAL$}  BED file with complex loci coordinates. May be repeated multiple times.\n\
        {EMPTY}  If fourth column is present, it is used for the locus name (must be unique).",
        "-L, --loci-bed".green(), "FILE".yellow());

    println!("\n{}", "Optional parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Reference genome name. If provided, add reference locus sequence\n\
        {EMPTY}  to the list of potential haplotypes.",
        "-g, --genome".green(), "STR".yellow());
    println!("    {:KEY$} {:VAL$}  If needed, expand loci boundaries by at most {} bp outwards [{}].",
        "-e, --expand".green(), "INT".yellow(), "INT".yellow(), defaults.max_expansion);
    println!("    {:KEY$} {:VAL$}  Select best locus boundary based on k-mer frequencies in\n\
        {EMPTY}  moving windows of size {} bp [{}].",
        "-w, --window".green(), "INT".yellow(), "INT".yellow(), defaults.moving_window);

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

            Short('g') | Long("genome") => args.ref_name = Some(parser.value()?.parse()?),
            Short('e') | Long("expand") => args.max_expansion = parser.value()?.parse()?,
            Short('w') | Long("window") => args.moving_window = parser.value()?.parse()?,
            Long("jellyfish") => args.jellyfish = parser.value()?.parse()?,

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
    // Make window size odd.
    args.moving_window += 1 - args.moving_window % 2;
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
        for line in common::open(filename)?.lines() {
            let interval = NamedInterval::parse_bed(&mut line?.split('\t'), contigs)?;
            if !names.insert(interval.name().to_string()) {
                return Err(Error::InvalidInput(format!("Locus name '{}' appears at least twice", interval.name())));
            }
            intervals.push(interval);
        }
    }
    Ok(intervals)
}

/// Best position (across equals) is Minimal/Maximal.
#[derive(Debug, Clone, Copy)]
enum PosLoc { Min, Max }

/// Across the islands `(start, end)`, finds the position with the smallest `aver_kmer_freq` value.
/// If there are several optimal position, returns
fn find_optimal_pos<I>(pos_loc: PosLoc, start: u32, aver_kmer_freqs: &[f64], islands: I) -> u32
where I: Iterator<Item = (u32, u32)>,
{
    let mut best_pos = 0;
    let mut best_freq = f64::INFINITY;
    for (isl_start, isl_end) in islands {
        let subvec = &aver_kmer_freqs[(isl_start - start) as usize..(isl_end - start) as usize];
        let (i, s) = match pos_loc {
            PosLoc::Min => F64Ext::argmin(subvec),
            // Finds last argmin.
            PosLoc::Max => IterExt::arg_optimal(subvec.iter().cloned(), |opt, e| opt <= e),
        };
        if s <= 1.0 {
            // Best frequency is already achieved.
            return isl_start + i as u32;
        } else if s < best_freq {
            best_freq = s;
            best_pos = isl_start + i as u32;
        }
    }
    best_pos
}

/// Find possible islands, where
/// - the variation graph contains no bubbles,
/// - reference sequence contains no Ns,
/// - average k-mer frequency is the smallest,
/// - position is closest to the boundary.
///
/// Returns best position, if found.
fn find_best_boundary(
    pos_loc: PosLoc,
    mov_window: u32,
    seq_shift: u32,
    seq: &[u8],
    vars: &[VcfRecord],
    kmer_getter: &JfKmerGetter,
) -> Result<Option<u32>, Error>
{
    let halfw = mov_window / 2;
    // New position can only be selected between `start` and `end`, where k-mer counts will be defined.
    let start = seq_shift + halfw;
    let end = seq_shift + seq.len() as u32 - halfw;

    /// Skip 20 bp around variants.
    const VAR_MARGIN: u32 = 20;
    let k = kmer_getter.k();
    let nrun_margin = max(VAR_MARGIN, k);
    let nruns = seq::n_runs(seq);

    let mut levels = Vec::with_capacity(1 + vars.len() + nruns.len());
    levels.push((start, end, 1));
    levels.extend(vars.iter().map(|var|
        ((var.pos() as u32).saturating_sub(VAR_MARGIN), var.end() as u32 + VAR_MARGIN, -1)));
    levels.extend(nruns.into_iter().map(|(run_start, run_end)|
        ((seq_shift + run_start).saturating_sub(nrun_margin), seq_shift + run_end + nrun_margin, -1)));
    // Find regions with depth == 1: +1 when inside the region, -1 when there is a bubble, -1 when there is an N run.
    let islands = seq::interv::split_disjoint(&levels, |depth| depth == 1);
    if islands.is_empty() {
        return Ok(None);
    }

    let kmer_counts: Vec<f64> = kmer_getter.fetch(&seq)?.into_iter().map(f64::from).collect();
    let divisor = f64::from(mov_window + 1 - k);
    // Average k-mer counts for each window of size `mov_window` over k-mers of size `k`.
    // Only defined for indices between `start` and `end`.
    let aver_kmer_freqs: Vec<f64> = VecExt::moving_window_sums(&kmer_counts, (mov_window + 1 - k) as usize)
        .into_iter().map(|count| count / divisor).collect();
    debug_assert_eq!(aver_kmer_freqs.len() as u32, end - start);

    // Find position with the smallest average frequency and closest to the boundary.
    let islands_iter = islands.iter().map(|(isl_start, isl_end, _)| (*isl_start, *isl_end));
    Ok(Some(match pos_loc {
        PosLoc::Min => find_optimal_pos(pos_loc, start, &aver_kmer_freqs, islands_iter),
        PosLoc::Max => find_optimal_pos(pos_loc, start, &aver_kmer_freqs, islands_iter.rev()),
    }))
}

fn write_locus(
    locus_dir: &Path,
    locus: &NamedInterval,
    ref_name: &Option<String>,
    seqs: &[(String, Vec<u8>)],
    kmer_getter: &JfKmerGetter,
) -> Result<(), Error>
{
    let mut bed_out = File::create(locus_dir.join("ref.bed"))?;
    bed_out.write_all(format!("{}\n", locus.interval().bed_fmt()).as_bytes())?;
    bed_out.sync_all()?;
    std::mem::drop(bed_out);

    let fasta_filename = locus_dir.join("haplotypes.fa");
    let mut fasta_out = File::create(&fasta_filename).map(BufWriter::new)?;
    for (name, seq) in seqs.iter() {
        let descr = match &ref_name {
            Some(ref_name) if ref_name == name => Some(locus.interval().to_string()),
            _ => None,
        };
        seq::write_fasta(&mut fasta_out, name, descr.as_ref().map(|s| s as &str), seq)?;
    }
    fasta_out.flush()?;
    std::mem::drop(fasta_out);

    log::debug!("    Counting k-mers for {}", locus);
    let contig_lengths: Vec<_> = seqs.iter().map(|(_name, seq)| seq.len() as u32).collect();
    let kmer_counts = kmer_getter.fetch_fasta(&fasta_filename, &contig_lengths)?;
    let mut kmers_out = common::create_gzip(&locus_dir.join("kmers.gz"))?;
    kmer_counts.save(&mut kmers_out)?;

    let fasta_gzip = common::append_path(&fasta_filename, ".gz");
    if fasta_gzip.exists() {
        // Remove file directly, as `gzip --force` is not always available.
        fs::remove_file(fasta_gzip)?;
    }
    let gzip_output = Command::new("gzip").arg(&fasta_filename).output()?;
    if !gzip_output.status.success() {
        return Err(Error::SubprocessFail(gzip_output));
    }
    Ok(())
}

/// Add `locus` to the database.
fn add_locus<R>(
    loci_dir: &Path,
    locus: &NamedInterval,
    fasta_file: &mut fasta::IndexedReader<R>,
    vcf_file: &mut bcf::IndexedReader,
    kmer_getter: &JfKmerGetter,
    args: &Args,
) -> Result<bool, Error>
where R: Read + Seek,
{
    log::info!("Analyzing locus {}", locus);
    let inner_interv = locus.interval();
    // Add extra half-window to each sides.
    let halfw = args.moving_window / 2;
    let outer_interv = inner_interv.expand(args.max_expansion + halfw, args.max_expansion + halfw);
    let outer_seq = outer_interv.fetch_seq(fasta_file)?;

    let vcf_rid = vcf_file.header().name2rid(outer_interv.contig_name().as_bytes())?;
    vcf_file.fetch(vcf_rid, u64::from(outer_interv.start()), Some(u64::from(outer_interv.end())))?;
    let vcf_recs: Vec<_> = vcf_file.records().collect::<Result<_, _>>()?;

    // Best locus coordinates (start, end) would be within
    // outer_start + halfw <= start <= inner_start  <  inner_end <= end <= outer_end - halfw.
    let (inner_start, inner_end) = inner_interv.range();
    let left_var_ix = bisect::right_by(&vcf_recs, |var| var.pos().cmp(&i64::from(inner_start + 1)));
    let right_var_ix = bisect::left_by(&vcf_recs, |var| var.end().cmp(&i64::from(inner_end)));

    // Extend region to the left.
    let outer_start = outer_interv.start();
    let new_start = match find_best_boundary(PosLoc::Max, args.moving_window, outer_start,
            &outer_seq[..(halfw + inner_start + 1 - outer_start) as usize], &vcf_recs[..left_var_ix], kmer_getter)? {
        Some(pos) => pos,
        None => {
            log::error!("Cannot extend locus {} to the left.\n    \
                Try increasing -e/--extend parameter or manually modifying region boundaries.", locus);
            return Ok(false);
        }
    };

    // Extend region to the right.
    let right_shift = inner_end - halfw - 1;
    let new_end = match find_best_boundary(PosLoc::Min, args.moving_window, right_shift,
            &outer_seq[(right_shift - outer_start) as usize..], &vcf_recs[right_var_ix..], kmer_getter)? {
        Some(pos) => pos + 1,
        None => {
            log::error!("Cannot extend locus {} to the right.\n    \
                Try increasing -e/--extend parameter or manually modify region boundaries.", locus);
            return Ok(false);
        }
    };
    let new_locus;
    if new_start != inner_start || new_end != inner_end {
        new_locus = locus.with_new_range(new_start, new_end);
        log::info!("    Extended locus {} by {} bp left and {} bp right. New locus: {}",
            locus, inner_start - new_start, new_end - inner_end, new_locus.interval());
    } else {
        new_locus = locus.clone();
    }

    let dir = loci_dir.join(new_locus.name());
    super::mkdir(&dir)?;
    let seqs = seq::panvcf::reconstruct_sequences(new_start,
        &outer_seq[(new_start - outer_start) as usize..(new_end - outer_start) as usize], &args.ref_name,
        vcf_file.header(), &vcf_recs)?;
    // TODO: Check on Ns within sequences + filter very similar sequences.
    write_locus(&dir, &new_locus, &args.ref_name, &seqs, kmer_getter)?;
    Ok(true)
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let mut args = process_args(parse_args(argv)?)?;
    // unwrap as argsuments were previously checked to be Some.
    let db_path = args.database.as_ref().unwrap();
    let ref_filename = args.reference.as_ref().unwrap();
    let vcf_filename = args.variants.as_ref().unwrap();

    let (contigs, mut fasta_file) = ContigNames::load_indexed_fasta(&ref_filename, "reference".to_string())?;
    let loci = load_loci(&contigs, &args.loci, &args.bed_files)?;

    let mut vcf_file = bcf::IndexedReader::from_path(vcf_filename)?;
    let mut jf_filenames = common::find_filenames(&db_path.join("jf"), "jf".as_ref())?;
    if jf_filenames.len() != 1 {
        return Err(Error::InvalidInput(format!("There are {} files {}/jf/*.jf (expected 1)",
            db_path.display(), jf_filenames.len())));
    }
    // unwrap, as we know that there is exactly one filename.
    let kmer_getter = JfKmerGetter::new(args.jellyfish.clone(), jf_filenames.pop().unwrap())?;
    args.moving_window = max(kmer_getter.k(), args.moving_window);

    let loci_dir = db_path.join("loci");
    for locus in loci.iter() {
        add_locus(&loci_dir, locus, &mut fasta_file, &mut vcf_file, &kmer_getter, &args)?;
    }
    log::info!("Success!");
    Ok(())
}