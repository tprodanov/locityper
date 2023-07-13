use std::{
    cmp::max,
    io::{self, BufRead, Read, Seek, Write},
    fs::{self, File},
    sync::Arc,
    path::{Path, PathBuf},
};
use fnv::FnvHashSet;
use bio::io::fasta;
use htslib::bcf::{
    self,
    Read as VcfRead,
    record::Record as VcfRecord,
};
use colored::Colorize;
use const_format::{str_repeat, concatcp};
use crate::{
    err::{Error, validate_param},
    algo::bisect,
    ext::{
        self,
        vec::{VecExt, F64Ext, IterExt},
    },
    seq::{
        self, NamedInterval, ContigNames, NamedSeq,
        kmers::JfKmerGetter,
        wfa::Penalties,
        dist,
    },
};
use super::paths;

struct Args {
    database: Option<PathBuf>,
    reference: Option<PathBuf>,
    variants: Option<PathBuf>,
    loci: Vec<String>,
    bed_files: Vec<PathBuf>,

    ref_name: String,
    leave_out: FnvHashSet<String>,
    max_expansion: u32,
    moving_window: u32,

    penalties: Penalties,
    max_divergence: f64,

    threads: u16,
    force: bool,
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

            ref_name: "REF".to_owned(),
            leave_out: Default::default(),
            max_expansion: 5000,
            moving_window: 500,

            penalties: Default::default(),
            max_divergence: 0.0001,

            threads: 8,
            force: false,
            jellyfish: PathBuf::from("jellyfish"),
        }
    }
}

impl Args {
    fn validate(mut self) -> Result<Self, Error> {
        validate_param!(self.database.is_some(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.reference.is_some(), "Reference fasta file is not provided (see -r/--reference)");
        validate_param!(self.variants.is_some(), "Variants VCF file is not provided (see -v/--vcf)");
        validate_param!(!self.loci.is_empty() || !self.bed_files.is_empty(),
            "Complex loci are not provided (see -l/--locus and -L/--loci-bed)");

        self.jellyfish = ext::sys::find_exe(self.jellyfish)?;
        // Make window size odd.
        self.moving_window += 1 - self.moving_window % 2;
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 16;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Adds complex locus/loci to the database.".yellow());

    println!("\n{}", "Usage:".bold());
    println!("    {} add -d db -r ref.fa -v vars.vcf.gz -l/-L loci [arguments]", super::PKG_NAME);
    // println!("    {} add -d db -s seqs.fa=name [arguments]", super::PKG_NAME);

    println!("\n{}", "Input arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Input database directory (initialized with {}).",
        "-d, --db".green(), "DIR".yellow(), concatcp!(super::PKG_NAME, " create").underline());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file.",
        "-r, --reference".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  PanGenie input VCF file. Encodes variation across pangenome samples.\n\
        {EMPTY}  Must be compressed and indexed with {}.",
        "-v, --vcf".green(), "FILE".yellow(), "tabix".underline());
    // println!("    {:KEY$} {:VAL$}  Fasta file with haplotype sequences\n
    //     {EMPTY}  Mutually exclusive with {}, {}, {}, {} and {}.",
    //     "-s, --seqs".green(), "FILE".yellow(),
    //     "-r".green(), "-v".green(), "-l".green(), "-L".green(), )

    println!("\n{}", "Complex loci coordinates:".bold());
    println!("    {:KEY$} {:VAL$}  Complex locus coordinates. Multiple loci are allowed.\n\
        {EMPTY}  Format: 'chrom:start-end' or 'chrom:start-end=name',\n\
        {EMPTY}  where 'name' is the locus name (must be unique).",
        "-l, --locus".green(), "STR+".yellow());
    println!("    {:KEY$} {:VAL$}  BED file with complex loci coordinates. May be repeated multiple times.\n\
        {EMPTY}  If fourth column is present, it is used for the locus name (must be unique).",
        "-L, --loci-bed".green(), "FILE".yellow());

    println!("\n{}", "Haplotype extraction parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Reference genome name [{}].",
        "-g, --genome".green(), "STR".yellow(), defaults.ref_name);
    println!("    {:KEY$} {:VAL$}  If needed, expand loci boundaries by at most {} bp outwards [{}].",
        "-e, --expand".green(), "INT".yellow(), "INT".yellow(), defaults.max_expansion);
    println!("    {:KEY$} {:VAL$}  Select best locus boundary based on k-mer frequencies in\n\
        {EMPTY}  moving windows of size {} bp [{}].",
        "-w, --window".green(), "INT".yellow(), "INT".yellow(), defaults.moving_window);
    println!("    {:KEY$} {:VAL$}  Leave out sequences with specified names.",
        "    --leave-out".green(), "STR+".yellow());

    println!("\n{}", "Haplotype clustering parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Penalty for mismatch [{}].",
        "-M, --mismatch".green(), "INT".yellow(), defaults.penalties.mismatch);
    println!("    {:KEY$} {:VAL$}  Gap open penalty [{}].",
        "-O, --gap-open".green(), "INT".yellow(), defaults.penalties.gap_opening);
    println!("    {:KEY$} {:VAL$}  Gap extend penalty [{}].",
        "-E, --gap-extend".green(), "INT".yellow(), defaults.penalties.gap_extension);
    println!("    {:KEY$} {:VAL$}  Sequence divergence threshold, used for haplotypes clustering [{}].",
        "-D, --divergence".green(), "FLOAT".yellow(), defaults.max_divergence);

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads);
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());
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
                args.loci = parser.values()?.map(ValueExt::string).collect::<Result<Vec<_>, _>>()?,
            Short('L') | Long("loci") | Long("loci-bed") => args.bed_files.push(parser.value()?.parse()?),

            Short('g') | Long("genome") => args.ref_name = parser.value()?.parse()?,
            Long("leave-out") | Long("leaveout") => {
                for val in parser.values()? {
                    args.leave_out.insert(val.parse()?);
                }
            }
            Short('e') | Long("expand") => args.max_expansion = parser.value()?.parse()?,
            Short('w') | Long("window") => args.moving_window = parser.value()?.parse()?,

            Short('M') | Long("mismatch") => args.penalties.mismatch = parser.value()?.parse()?,
            Short('O') | Long("gap-open") | Long("gap-opening") =>
                args.penalties.gap_opening = parser.value()?.parse()?,
            Short('E') | Long("gap-extend") | Long("gap-extension") =>
                args.penalties.gap_extension = parser.value()?.parse()?,
            Short('D') | Long("divergence") => args.max_divergence = parser.value()?.parse()?,

            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,
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

/// Loads named interval from a list of loci and a list of BED files. Names must not repeat.
fn load_loci(contigs: &Arc<ContigNames>, loci: &[String], bed_files: &[PathBuf]) -> Result<Vec<NamedInterval>, Error> {
    let mut intervals = Vec::new();
    let mut names = FnvHashSet::default();
    for locus in loci.iter() {
        let interval = NamedInterval::parse(locus, contigs)?;
        if !names.insert(interval.name().to_string()) {
            return Err(Error::InvalidInput(format!("Locus name '{}' appears at least twice", interval.name())));
        }
        intervals.push(interval);
    }

    for filename in bed_files.iter() {
        for line in ext::sys::open(filename)?.lines() {
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
        let (i, opt_freq) = match pos_loc {
            PosLoc::Min => F64Ext::argmin(subvec),
            // Finds last argmin.
            PosLoc::Max => IterExt::arg_optimal(subvec.iter().copied(), |opt, e| opt <= e),
        };
        if opt_freq <= 1.0 {
            // Best frequency is already achieved.
            return isl_start + i as u32;
        } else if opt_freq < best_freq {
            best_freq = opt_freq;
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

    let kmer_counts: Vec<f64> = kmer_getter.fetch_one(seq.to_vec())?
        .take_first().into_iter()
        .map(f64::from).collect();
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

/// Check haplotypes:
/// - warn if there are any Ns in the sequences,
/// - warn if sequence divergence is too high.
fn check_haplotypes(tag: &str, entries: &[NamedSeq], mut divergences: impl Iterator<Item = f64>) {
    for (i, entry1) in entries.iter().enumerate() {
        if entry1.seq().contains(&b'N') {
            log::warn!("    Locus {}, haplotype {} contains Ns!", tag, entry1.name());
        }
        for entry2 in entries[i + 1..].iter() {
            let diverg = divergences.next().unwrap();
            if diverg >= 0.2 {
                log::warn!("    Locus {}, haplotypes {} and {} have high divergence ({:.5})",
                    tag, entry1.name(), entry2.name(), diverg);
            }
        }
    }
}

/// Cluster haplotypes using `kodama` crate
/// (which, in turn, is based on this paper https://arxiv.org/pdf/1109.2378.pdf).
///
/// Clusters with divergence not exceeding the threshold are joined.
/// Returns boolean vector: does the sequence remain after mergin?
fn cluster_haplotypes(
    mut nwk_writer: impl Write,
    entries: &[NamedSeq],
    mut divergences: Vec<f64>,
    thresh: f64,
) -> io::Result<Vec<bool>> {
    let n = entries.len();
    let total_clusters = 2 * n - 1;
    // Use Complete method, meaning that we track maximal distance between two points between two clusters.
    // This is done to cut as little as needed.
    let dendrogram = kodama::linkage(&mut divergences, n, kodama::Method::Complete);
    let mut clusters_nwk = Vec::with_capacity(total_clusters);
    // Cluster representatives.
    let mut cluster_repr = Vec::with_capacity(total_clusters);
    clusters_nwk.extend(entries.iter().map(|entry| entry.name().to_owned()));
    cluster_repr.extend(0..n);

    let steps = dendrogram.steps();
    for step in steps.iter() {
        let i = step.cluster1;
        let j = step.cluster2;
        clusters_nwk.push(format!("({}:{dist},{}:{dist})", &clusters_nwk[i], &clusters_nwk[j],
            dist = 0.5 * step.dissimilarity));
        let size1 = if i < n { 1 } else { steps[i - n].size };
        let size2 = if j < n { 1 } else { steps[j - n].size };
        cluster_repr.push(cluster_repr[if size1 >= size2 { i } else { j }]);
    }
    assert_eq!(clusters_nwk.len(), total_clusters);
    writeln!(nwk_writer, "{};", clusters_nwk.last().unwrap())?;

    let mut queue = vec![total_clusters - 1];
    let mut keep_seqs = vec![false; n];
    while let Some(i) = queue.pop() {
        if i < n {
            keep_seqs[i] = true;
        } else {
            let step = &steps[i - n];
            if step.dissimilarity <= thresh {
                keep_seqs[cluster_repr[i]] = true;
            } else {
                queue.push(step.cluster1);
                queue.push(step.cluster2);
            }
        }
    }
    Ok(keep_seqs)
}

/// Process haplotypes: write FASTA, PAF, kmers files, cluster sequences.
fn process_haplotypes(
    locus_dir: &Path,
    locus: &NamedInterval,
    entries: Vec<NamedSeq>,
    kmer_getter: &JfKmerGetter,
    args: &Args,
) -> Result<(), Error>
{
    let tag = locus.name();
    log::info!("    [{}] Writing {} haplotypes to {}/", tag, entries.len(), ext::fmt::path(locus_dir));
    let mut bed_writer = File::create(locus_dir.join(paths::LOCUS_BED))?;
    writeln!(bed_writer, "{}", locus.bed_fmt())?;
    std::mem::drop(bed_writer);

    log::info!("    [{}] Calculating haplotype divergence", tag);
    let paf_writer = ext::sys::create_gzip(&locus_dir.join(paths::LOCUS_PAF))?;
    let divergences = dist::pairwise_divergences(&entries, paf_writer, &args.penalties, args.threads)?;
    check_haplotypes(tag, &entries, divergences.iter().copied());

    log::info!("    [{}] Clustering haploypes", tag);
    let nwk_writer = io::BufWriter::new(File::create(&locus_dir.join(paths::LOCUS_DENDROGRAM))?);
    let keep_seqs = cluster_haplotypes(nwk_writer, &entries, divergences, args.max_divergence)?;
    let n_filtered = keep_seqs.iter().fold(0, |sum, &keep| sum + usize::from(keep));
    if n_filtered == entries.len() {
        log::info!("        Keep all sequences after clustering");
    } else {
        log::info!("        Discard {} sequences after clustering", entries.len() - n_filtered);
    }

    let filt_fasta_filename = locus_dir.join(paths::LOCUS_FASTA);
    let mut filt_fasta_writer = ext::sys::create_gzip(&filt_fasta_filename)?;
    let mut all_fasta_writer = ext::sys::create_gzip(&locus_dir.join(paths::LOCUS_FASTA_ALL))?;
    let mut filt_seqs = Vec::with_capacity(n_filtered);
    for (&keep, entry) in keep_seqs.iter().zip(entries.into_iter()) {
        seq::write_fasta(&mut all_fasta_writer, entry.name().as_bytes(), entry.seq())?;
        if keep {
            seq::write_fasta(&mut filt_fasta_writer, entry.name().as_bytes(), entry.seq())?;
            filt_seqs.push(entry.take_seq());
        }
    }
    std::mem::drop((filt_fasta_writer, all_fasta_writer));

    log::info!("    [{}] Counting k-mers", tag);
    let kmer_counts = kmer_getter.fetch(filt_seqs)?;
    let mut kmers_writer = ext::sys::create_gzip(&locus_dir.join(paths::KMERS))?;
    kmer_counts.save(&mut kmers_writer)?;
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
        log::info!("    [{}] Extending locus by {} bp left and {} bp right. New locus: {}",
            locus.name(), inner_start - new_start, new_end - inner_end, new_locus.interval());
    } else {
        new_locus = locus.clone();
    }

    let dir = loci_dir.join(new_locus.name());
    if dir.exists() {
        if args.force {
            log::warn!("    Clearing directory {}", ext::fmt::path(&dir));
            fs::remove_dir_all(&dir)?;
        } else {
            log::warn!("    Skipping directory {}", ext::fmt::path(&dir));
            return Ok(true);
        }
    }
    ext::sys::mkdir(&dir)?;
    log::info!("    [{}] Reconstructing haplotypes", new_locus.name());
    let reconstruction = seq::panvcf::reconstruct_sequences(new_start,
        &outer_seq[(new_start - outer_start) as usize..(new_end - outer_start) as usize], &args.ref_name,
        vcf_file.header(), &vcf_recs, &args.leave_out);
    match reconstruction {
        Ok(seqs) => process_haplotypes(&dir, &new_locus, seqs, kmer_getter, &args)?,
        Err(Error::InvalidData(e)) => {
            log::error!("Cannot extract locus {} sequences: {}", locus, e);
            return Ok(false);
        }
        Err(e) => return Err(e),
    }
    Ok(true)
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let mut args = parse_args(argv)?.validate()?;
    // unwrap as argsuments were previously checked to be Some.
    let db_path = args.database.as_ref().unwrap();
    let ref_filename = args.reference.as_ref().unwrap();
    let vcf_filename = args.variants.as_ref().unwrap();

    let (contigs, mut fasta_file) = ContigNames::load_indexed_fasta("reference", &ref_filename)?;
    let loci = load_loci(&contigs, &args.loci, &args.bed_files)?;

    let mut vcf_file = bcf::IndexedReader::from_path(vcf_filename)?;
    let mut jf_filenames = ext::sys::filenames_with_ext(&db_path.join(paths::JF_DIR), "jf")?;
    if jf_filenames.len() != 1 {
        return Err(Error::InvalidInput(format!("There are {} files {}/jf/*.jf (expected 1)",
            db_path.display(), jf_filenames.len())));
    }
    // unwrap, as we know that there is exactly one filename.
    let kmer_getter = JfKmerGetter::new(args.jellyfish.clone(), jf_filenames.pop().unwrap())?;
    args.moving_window = max(kmer_getter.k(), args.moving_window);

    let loci_dir = db_path.join(paths::LOCI_DIR);
    ext::sys::mkdir(&loci_dir)?;
    let mut n_successes = 0;
    for locus in loci.iter() {
        n_successes += u32::from(add_locus(&loci_dir, locus, &mut fasta_file, &mut vcf_file, &kmer_getter, &args)?);
    }
    if n_successes > 0 {
        log::info!("Success!");
    } else {
        log::error!("Could not add any loci.");
    }
    Ok(())
}
