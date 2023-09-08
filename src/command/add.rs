use std::{
    cmp::{min, max},
    io::{self, BufRead, Read, Seek, Write},
    fs::File,
    sync::Arc,
    path::{Path, PathBuf},
    time::Instant,
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
    err::{Error, validate_param, add_path},
    algo::bisect,
    ext::{
        self,
        vec::{VecExt, IterExt},
        fmt::PrettyU32,
    },
    seq::{
        self, NamedInterval, Interval, ContigNames, NamedSeq,
        dist, panvcf, interv, fastx,
        contigs::GenomeVersion,
        kmers::JfKmerGetter,
        wfa::Penalties,
    },
};
use super::paths;

struct Args {
    database: Option<PathBuf>,
    reference: Option<PathBuf>,
    variants: Option<PathBuf>,
    loci: Vec<String>,
    bed_files: Vec<PathBuf>,
    sequences: Vec<String>,

    ref_name: Option<String>,
    leave_out: FnvHashSet<String>,
    max_expansion: u32,
    moving_window: u32,

    penalties: Penalties,
    unknown_frac: f64,
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
            sequences: Vec::new(),

            ref_name: None,
            leave_out: Default::default(),
            max_expansion: 50000,
            moving_window: 500,

            penalties: Default::default(),
            unknown_frac: 0.0001,
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
        if self.sequences.is_empty() {
            validate_param!(self.reference.is_some(), "Reference fasta file is not provided (see -r/--reference)");
            validate_param!(self.variants.is_some(), "Pangenome VCF file is not provided (see -v/--vcf)");
            validate_param!(!self.loci.is_empty() || !self.bed_files.is_empty(),
                "Complex loci are not provided (see -l/--locus and -L/--loci-bed)");
        } else {
            validate_param!(self.reference.is_none(), "Reference (-r) is mutually exclusive with sequences (-s)");
            validate_param!(self.variants.is_none(), "VCF file (-v) is mutually exclusive with sequences (-s)");
            validate_param!(self.loci.is_empty() && self.bed_files.is_empty(),
                "Loci (-l/-L) are mutually exclusive with sequences (-s)");
        }

        validate_param!(0.0 <= self.unknown_frac && self.unknown_frac <= 1.0,
            "Unknown fraction ({}) must be within [0, 1]", self.unknown_frac);
        validate_param!(0.0 <= self.max_divergence && self.max_divergence <= 1.0,
            "Maximum divergence ({}) must be within [0, 1]", self.max_divergence);

        self.jellyfish = ext::sys::find_exe(self.jellyfish)?;
        // Make window size odd.
        self.moving_window += 1 - self.moving_window % 2;
        self.penalties.validate()?;
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
    println!("    {} add -d db -r ref.fa -v vars.vcf.gz -l/-L loci [arguments]", super::PROGRAM);
    println!("    {} add -d db -s seqs.fa=name [seqs2.fa=name2 ...] [arguments]", super::PROGRAM);

    println!("\n{}", "Input arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Input database directory (initialized with {}).",
        "-d, --database".green(), "DIR".yellow(), concatcp!(super::PROGRAM, " create").underline());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file.",
        "-r, --reference".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Input VCF file, encoding variation across pangenome samples.\n\
        {EMPTY}  Must be compressed and indexed with {}.",
        "-v, --vcf".green(), "FILE".yellow(), "tabix".underline());
    println!("    {:KEY$} {:VAL$}\n\
        {EMPTY}  Fasta file(s) with locus alleles. Format: {}={}.\n\
        {EMPTY}  Mutually exclusive with {}, {} and {}/{}.\n\
        {EMPTY}  All sequences {} on the same strand!",
        "-s, --seqs".green(), "FILE=STR".yellow(), "filename".underline(), "locus_name".underline(),
        "-r".green(), "-v".green(), "-l".green(), "-L".green(), "must be".red());

    println!("\n{}", "Complex loci coordinates:".bold());
    println!("    {:KEY$} {:VAL$}  Complex locus coordinates. Multiple loci are allowed.\n\
        {EMPTY}  Format: 'chrom:start-end=name',\n\
        {EMPTY}  where 'name' is the locus name (must be unique).",
        "-l, --locus".green(), "STR+".yellow());
    println!("    {:KEY$} {:VAL$}  BED file with complex loci coordinates. May be repeated multiple times.\n\
        {EMPTY}  Fourth column must contain locus name (all names should be unique).",
        "-L, --loci-bed".green(), "FILE".yellow());

    println!("\n{}", "Haplotype extraction parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Reference genome name, default: tries to guess.",
        "-g, --genome".green(), "STR".yellow());
    println!("    {:KEY$} {:VAL$}  If needed, expand loci boundaries by at most {} bp outwards [{}].",
        "-e, --expand".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(PrettyU32(defaults.max_expansion)));
    println!("    {:KEY$} {:VAL$}  Select best locus boundary based on k-mer frequencies in\n\
        {EMPTY}  moving windows of size {} bp [{}].",
        "-w, --window".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(PrettyU32(defaults.moving_window)));
    println!("    {:KEY$} {:VAL$}  Allow this fraction of unknown nucleotides per haplotype [{}]\n\
        {EMPTY}  (relative to the haplotype length). Variants that have no known\n\
        {EMPTY}  variation in the input VCF pangenome are ignored.",
        "-u, --unknown".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.unknown_frac));
    println!("    {:KEY$} {:VAL$}  Leave out sequences with specified names.",
        "    --leave-out".green(), "STR+".yellow());

    println!("\n{}", "Haplotype clustering:".bold());
    println!("    {:KEY$} {:VAL$}  Penalty for mismatch [{}].",
        "-M, --mismatch".green(), "INT".yellow(), super::fmt_def(defaults.penalties.mismatch));
    println!("    {:KEY$} {:VAL$}  Gap open penalty [{}].",
        "-O, --gap-open".green(), "INT".yellow(), super::fmt_def(defaults.penalties.gap_opening));
    println!("    {:KEY$} {:VAL$}  Gap extend penalty [{}].",
        "-E, --gap-extend".green(), "INT".yellow(), super::fmt_def(defaults.penalties.gap_extension));
    println!("    {:KEY$} {:VAL$}  Sequence divergence threshold,\n\
        {EMPTY}  used to discard almost identical locus alleles [{}].",
        "-D, --divergence".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.max_divergence));

    println!("\n{}", "Execution arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), super::fmt_def(defaults.threads));
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Jellyfish executable [{}].",
        "    --jellyfish".green(), "EXE".yellow(), super::fmt_def(defaults.jellyfish.display()));

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
            Short('d') | Long("db") | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),
            Short('v') | Long("vcf") => args.variants = Some(parser.value()?.parse()?),

            Short('l') | Long("locus") =>
                args.loci = parser.values()?.map(ValueExt::string).collect::<Result<Vec<_>, _>>()?,
            Short('L') | Long("loci") | Long("loci-bed") => args.bed_files.push(parser.value()?.parse()?),
            Short('s') | Long("seqs") | Long("sequences") =>
                args.sequences = parser.values()?.map(|s| s.parse()).collect::<Result<_, _>>()?,

            Short('g') | Long("genome") => args.ref_name = Some(parser.value()?.parse()?),
            Long("leave-out") | Long("leaveout") => {
                for val in parser.values()? {
                    args.leave_out.insert(val.parse()?);
                }
            }
            Short('e') | Long("expand") => args.max_expansion = parser.value()?.parse::<PrettyU32>()?.get(),
            Short('w') | Long("window") => args.moving_window = parser.value()?.parse::<PrettyU32>()?.get(),

            Short('M') | Long("mismatch") => args.penalties.mismatch = parser.value()?.parse()?,
            Short('O') | Long("gap-open") | Long("gap-opening") =>
                args.penalties.gap_opening = parser.value()?.parse()?,
            Short('E') | Long("gap-extend") | Long("gap-extension") =>
                args.penalties.gap_extension = parser.value()?.parse()?,
            Short('D') | Long("divergence") => args.max_divergence = parser.value()?.parse()?,
            Short('u') | Long("unknown") => args.unknown_frac = parser.value()?.parse()?,

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
        let interval = NamedInterval::parse_explicit(locus, contigs)?;
        if !names.insert(interval.name().to_string()) {
            return Err(Error::InvalidInput(format!("Locus name '{}' appears at least twice", interval.name())));
        }
        intervals.push(interval);
    }

    for filename in bed_files.iter() {
        for line in ext::sys::open(filename)?.lines() {
            let line = line.map_err(add_path!(filename))?;
            let interval = NamedInterval::parse_bed(&mut line.split('\t'), contigs)?;
            if !interval.is_name_explicit() {
                return Err(Error::InvalidInput(format!("Interval '{}' has no name", interval)));
            }
            if !names.insert(interval.name().to_string()) {
                return Err(Error::InvalidInput(format!("Locus name '{}' appears at least twice", interval.name())));
            }
            intervals.push(interval);
        }
    }
    Ok(intervals)
}

/// Find possible islands, where
/// - the variation graph contains no bubbles,
/// - average k-mer frequency is the smallest,
/// - position is closest to the boundary.
///
/// Returns best position, if found.
fn find_best_boundary<const LEFT: bool>(
    mov_window: u32,
    seq_shift: u32,
    seq: &[u8],
    vars: &[VcfRecord],
    kmer_getter: &JfKmerGetter,
    expand_size: u32,
) -> Result<Option<u32>, Error>
{
    let halfw = mov_window / 2;
    // New position can only be selected between `start` and `end`, where k-mer counts will be defined.
    let start = seq_shift + halfw;
    let end = seq_shift + seq.len() as u32 - halfw;

    let unique_kmers: Vec<u32> = kmer_getter.fetch_one(seq.to_vec())?
        .take_first().into_iter()
        .map(|count| u32::from(count <= 1))
        .collect();
    let k = kmer_getter.k();
    let divisor = f64::from(mov_window + 1 - k);
    // Initially, set weights to the fraction of unique kmers.
    let mut weights: Vec<f64> = VecExt::moving_window_sums(&unique_kmers, (mov_window + 1 - k) as usize)
        .into_iter().map(|count| f64::from(count) / divisor).collect();
    assert_eq!(weights.len() as u32, end - start);

    // Try to select boundary at least 10 bp from any variant.
    const EFFECT_MARGIN: u32 = 9;
    let effect_divisor = f64::from(EFFECT_MARGIN + 1);
    for var in vars.iter() {
        let var_start = var.pos() as u32;
        let var_end = var.end() as u32;
        // Ignore positions with variants.
        for i in var_start.saturating_sub(start)..min(var_end, end).saturating_sub(start) {
            weights[i as usize] = 0.0;
        }
        // Downgrade positions close to variants.
        // i in 0..EFFECT_MARGIN  &&  start <= var_start - i - 1 < end.
        for i in var_start.saturating_sub(end)..var_start.saturating_sub(start).min(EFFECT_MARGIN) {
            weights[(var_start - start - i - 1) as usize] *= f64::from(EFFECT_MARGIN - i) / effect_divisor;
        }
        // i in 0..EFFECT_MARGIN  &&  start <= var_end + i < end.
        for i in start.saturating_sub(var_end)..end.saturating_sub(var_end).min(EFFECT_MARGIN) {
            weights[(var_end + i - start) as usize] *= f64::from(i + 1) / effect_divisor;
        }
    }

    // `inner_start - expand_size` and `inner_end + expand_size` are penalized by 20%.
    const WEIGHT_DROP: f64 = 0.2;
    let per_bp_drop = WEIGHT_DROP / f64::from(expand_size);
    if LEFT {
        weights.iter_mut().rev().enumerate().for_each(|(i, w)| *w -= *w * per_bp_drop * i as f64);
    } else {
        weights.iter_mut().enumerate().for_each(|(i, w)| *w -= *w * per_bp_drop * i as f64);
    }

    let (i, maxval) = if LEFT {
        // Finds last argmax
        IterExt::arg_optimal(weights.iter().copied(), |opt, e| opt >= e)
    } else {
        IterExt::arg_optimal(weights.iter().copied(), |opt, e| opt > e)
    };
    if maxval == 0.0 {
        Ok(None)
    } else {
        Ok(Some(start + i as u32))
    }
}

/// Check divergencies and warns if they are too high.
fn check_divergencies(tag: &str, entries: &[NamedSeq], mut divergences: impl Iterator<Item = f64>, from_vcf: bool) {
    let n = entries.len();
    let mut count = 0;
    let mut highest = 0.0;
    let mut highest_i = 0;
    let mut highest_j = 0;
    for i in 0..n {
        for j in i + 1..n {
            let diverg = divergences.next().unwrap();
            if diverg >= 0.2 {
                count += 1;
                if diverg > highest {
                    highest = diverg;
                    highest_i = i;
                    highest_j = j;
                }
            }
        }
    }
    if count > 0 {
        log::warn!("    [{}] {} allele pairs with high divergence, highest {:.5} ({} and {})", tag, count,
            highest, entries[highest_i].name(), entries[highest_j].name());
        if highest > 0.5 && !from_vcf {
            log::error!("    Please check that all alleles lie on the same strand");
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
    tag: &str,
    entries: Vec<NamedSeq>,
    kmer_getter: &JfKmerGetter,
    args: &Args,
) -> Result<(), Error>
{
    log::info!("    Writing {} haplotypes to {}/", entries.len(), ext::fmt::path(locus_dir));
    log::info!("    Calculating haplotype divergence");
    let paf_filename = locus_dir.join(paths::LOCUS_PAF);
    let paf_writer = ext::sys::create_gzip(&paf_filename)?;
    let divergences = dist::pairwise_divergences(&entries, paf_writer, &args.penalties, args.threads)
        .map_err(add_path!(paf_filename))?;
    check_divergencies(tag, &entries, divergences.iter().copied(), args.variants.is_some());

    log::info!("    Clustering haploypes");
    let nwk_filename = locus_dir.join(paths::LOCUS_DENDROGRAM);
    let nwk_writer = ext::sys::create_file(&nwk_filename)?;
    let keep_seqs = cluster_haplotypes(nwk_writer, &entries, divergences, args.max_divergence)
        .map_err(add_path!(nwk_filename))?;
    let n_filtered = keep_seqs.iter().fold(0, |sum, &keep| sum + usize::from(keep));
    if n_filtered == entries.len() {
        log::info!("        Keep all sequences after clustering");
    } else {
        log::info!("        Discard {} sequences after clustering", entries.len() - n_filtered);
    }

    let filt_fasta_filename = locus_dir.join(paths::LOCUS_FASTA);
    let mut filt_fasta_writer = ext::sys::create_gzip(&filt_fasta_filename)?;
    let locus_fasta = locus_dir.join(paths::LOCUS_FASTA_ALL);
    let mut all_fasta_writer = ext::sys::create_gzip(&locus_fasta)?;
    let mut filt_seqs = Vec::with_capacity(n_filtered);
    for (&keep, entry) in keep_seqs.iter().zip(entries.into_iter()) {
        seq::write_fasta(&mut all_fasta_writer, entry.name().as_bytes(), entry.seq()).map_err(add_path!(locus_fasta))?;
        if keep {
            seq::write_fasta(&mut filt_fasta_writer, entry.name().as_bytes(), entry.seq())
                .map_err(add_path!(locus_fasta))?;
            filt_seqs.push(entry.take_seq());
        }
    }
    std::mem::drop((filt_fasta_writer, all_fasta_writer));

    log::info!("    Counting k-mers");
    let kmer_counts = kmer_getter.fetch(filt_seqs)?;
    let kmers_filename = locus_dir.join(paths::KMERS);
    let mut kmers_writer = ext::sys::create_lz4_slow(&kmers_filename)?;
    kmer_counts.save(&mut kmers_writer).map_err(add_path!(kmers_filename))?;
    super::write_success_file(locus_dir.join(paths::SUCCESS))?;
    Ok(())
}

/// Checks unknown sequences, and returns max interval that contains `inner_interv` and does not contain Ns.
/// If `inner_interv` contains Ns as well, returns None.
fn process_unknown_seq(
    locus: &str,
    full_seq: &[u8],
    full_interv: &Interval,
    inner_interv: &Interval,
) -> Option<Interval>
{
    let n_runs = seq::n_runs(full_seq);
    if n_runs.is_empty() {
        return Some(full_interv.clone());
    }

    let (full_start, full_end) = full_interv.range();
    let (inner_start, inner_end) = inner_interv.range();
    let mut new_start = full_start;
    let mut new_end = full_end;
    for (mut start, mut end) in n_runs {
        start += full_start;
        end += full_start;
        if start <= inner_start {
            if end <= inner_start {
                new_start = end;
            } else {
                return None;
            }
        } else if end >= inner_end {
            if start >= inner_end {
                new_end = start;
            } else {
                return None;
            }
        }
    }
    if new_start > full_start {
        log::debug!("    [{}] Unknown sequence {} bp to the left ({}:{})", locus, inner_start - new_start,
            inner_interv.contig_name(), new_start + 1);
    }
    if new_end < full_end {
        log::debug!("    [{}] Unknown sequence {} bp to the right ({}:{})", locus, new_end - inner_end + 1,
            inner_interv.contig_name(), new_end + 1);
    }
    Some(Interval::new(Arc::clone(full_interv.contigs()), full_interv.contig_id(), new_start, new_end))
}

/// Add `locus` to the database.
fn add_locus<R>(
    loci_dir: &Path,
    locus: &NamedInterval,
    fasta_file: &mut fasta::IndexedReader<R>,
    vcf_file: &mut bcf::IndexedReader,
    haplotypes: &panvcf::AllHaplotypes,
    kmer_getter: &JfKmerGetter,
    args: &Args,
) -> Result<bool, Error>
where R: Read + Seek,
{
    let locus_dir = loci_dir.join(locus.name());
    if !super::Rerun::from_force(args.force).prepare_dir(&locus_dir)? {
        return Ok(true);
    }
    log::info!("{} {}", "Analyzing locus".bold(), locus);
    let inner_interv = locus.interval();
    // Add extra half-window to each sides.
    let halfw = args.moving_window / 2;
    let full_interv = inner_interv.expand(args.max_expansion + halfw, args.max_expansion + halfw);
    let full_seq = full_interv.fetch_seq(fasta_file).map_err(add_path!(args.reference.as_ref().unwrap()))?;

    let outer_interv = match process_unknown_seq(locus.name(), &full_seq, &full_interv, &inner_interv) {
        Some(interv) => interv,
        None => {
            log::error!("Cannot add locus {}. Interval {} contains unknown sequence", locus, inner_interv);
            return Ok(false);
        }
    };
    let outer_seq = &full_seq[
        (outer_interv.start() - full_interv.start()) as usize..(outer_interv.end() - full_interv.start()) as usize];

    let vcf_rid = vcf_file.header().name2rid(outer_interv.contig_name().as_bytes())?;
    vcf_file.fetch(vcf_rid, u64::from(outer_interv.start()), Some(u64::from(outer_interv.end())))?;
    let vcf_recs = panvcf::filter_variants(vcf_file, haplotypes)?;

    // Best locus coordinates (start, end) would be within
    // outer_start + halfw <= start <= inner_start  <  inner_end <= end <= outer_end - halfw.
    let (inner_start, inner_end) = inner_interv.range();
    let left_var_ix = bisect::right_by(&vcf_recs, |var| var.pos().cmp(&i64::from(inner_start + 1)));
    let right_var_ix = bisect::left_by(&vcf_recs, |var| var.end().cmp(&i64::from(inner_end)));

    // Extend region to the left.
    let outer_start = outer_interv.start();
    let new_start = match find_best_boundary::<true>(args.moving_window, outer_start,
            &outer_seq[..(halfw + inner_start + 1 - outer_start) as usize], &vcf_recs[..left_var_ix],
            kmer_getter, args.max_expansion)? {
        Some(pos) => pos,
        None => {
            log::error!("Cannot expand locus {} to the left due to a variant overlapping boundary.\n    \
                Try increasing -e/--expand parameter or manually modifying region boundaries.", locus.name());
            return Ok(false);
        }
    };

    // Extend region to the right.
    let right_shift = inner_end - halfw - 1;
    let new_end = match find_best_boundary::<false>(args.moving_window, right_shift,
            &outer_seq[(right_shift - outer_start) as usize..], &vcf_recs[right_var_ix..],
            kmer_getter, args.max_expansion)? {
        Some(pos) => pos + 1,
        None => {
            log::error!("Cannot expand locus {} to the right due to a variant overlapping boundary.\n    \
                Try increasing -e/--expand parameter or manually modify region boundaries.", locus.name());
            return Ok(false);
        }
    };
    let new_locus;
    if new_start != inner_start || new_end != inner_end {
        new_locus = locus.with_new_range(new_start, new_end);
        log::info!("    Extending locus by {} bp left and {} bp right -> {}",
            inner_start - new_start, new_end - inner_end, new_locus.interval());
        assert!(new_locus.is_name_explicit());
    } else {
        new_locus = locus.clone();
    }
    let locus_bed = locus_dir.join(paths::LOCUS_BED);
    let mut bed_writer = File::create(&locus_bed).map_err(add_path!(locus_bed))?;
    writeln!(bed_writer, "{}", new_locus.bed_fmt()).map_err(add_path!(locus_bed))?;
    std::mem::drop(bed_writer);

    log::info!("    Reconstructing haplotypes");
    let ref_seq = &outer_seq[(new_start - outer_start) as usize..(new_end - outer_start) as usize];
    let reconstruction = panvcf::reconstruct_sequences(new_start, ref_seq, &vcf_recs, haplotypes,
        vcf_file.header(), args.unknown_frac);
    match reconstruction {
        Ok(seqs) => process_haplotypes(&locus_dir, new_locus.name(), seqs, kmer_getter, &args)?,
        Err(Error::InvalidData(e)) => {
            log::error!("Cannot extract locus {} sequences: {}", new_locus, e);
            return Ok(false);
        }
        Err(e) => return Err(e),
    }
    Ok(true)
}

/// Add loci when pangenome is provided as a VCF file.
/// Returns number of successful and total loci.
fn run_with_vcf(loci_dir: &Path, kmer_getter: &JfKmerGetter, args: &Args) -> Result<(u32, u32), Error> {
    let ref_filename = args.reference.as_ref().unwrap();
    let vcf_filename = args.variants.as_ref().unwrap();

    let (contigs, mut fasta_file) = ContigNames::load_indexed_fasta("reference", &ref_filename)?;
    let ref_name = args.ref_name.as_ref().map(|s| s as &str)
        .or_else(|| GenomeVersion::guess(&contigs).map(GenomeVersion::to_str))
        .ok_or_else(|| Error::RuntimeError("Cannot guess reference genome name, please provide using -g".to_owned()))?;

    let loci = load_loci(&contigs, &args.loci, &args.bed_files)?;
    let mut vcf_file = bcf::IndexedReader::from_path(vcf_filename)?;
    let haplotypes = panvcf::AllHaplotypes::new(&mut vcf_file, &ref_name, &args.leave_out)?;

    let mut succeed = 0;
    let mut total = 0;
    for locus in loci.iter() {
        if add_locus(loci_dir, locus, &mut fasta_file, &mut vcf_file, &haplotypes, kmer_getter, &args)? {
            succeed += 1;
        }
        total += 1;
    }
    Ok((succeed, total))
}

fn check_sequences(seqs: &[NamedSeq], locus: &str) -> Result<(), Error> {
    if seqs.is_empty() {
        return Err(Error::InvalidData(format!("No sequences available for locus {}", locus)));
    }
    let min_size = seqs.iter().map(|s| s.len()).min().unwrap();
    if min_size < 1000 {
        return Err(Error::InvalidData(format!("Locus alleles are too short for locus {} (shortest: {} bp)",
            locus, min_size)));
    } else if min_size < 10000 {
        log::warn!("[{}] Locus alleles may be too short (shortest: {} bp)", locus, min_size);
    }
    if seqs.iter().any(|entry| seq::has_n(entry.seq())) {
        return Err(Error::InvalidData(format!("Locus alleles contain Ns for locus {}", locus)));
    }

    const AFFIX_SIZE: usize = 2;
    let seq0 = seqs[0].seq();
    let prefix = &seq0[..AFFIX_SIZE];
    let suffix = &seq0[seq0.len() - AFFIX_SIZE..];
    if seqs[1..].iter().map(NamedSeq::seq)
            .any(|s| &s[..AFFIX_SIZE] != prefix || &s[s.len() - AFFIX_SIZE..] != suffix) {
        log::warn!("[{}] There are variants on the boundary of the locus", locus);
    }
    Ok(())
}

fn process_locus_from_fasta(
    locus: &str,
    fasta_path: &Path,
    loci_dir: &Path,
    kmer_getter: &JfKmerGetter,
    args: &Args,
) -> Result<(), Error>
{
    let locus_dir = loci_dir.join(locus);
    if !super::Rerun::from_force(args.force).prepare_dir(&locus_dir)? {
        return Ok(())
    }
    let mut fasta_reader = fastx::Reader::from_path(fasta_path)?;
    let seqs = fasta_reader.read_all().map_err(add_path!(fasta_path))?;
    check_sequences(&seqs, locus)?;
    process_haplotypes(&locus_dir, locus, seqs, kmer_getter, &args)?;
    Ok(())
}

/// Add loci based on explicit fasta files.
/// Returns number of successful and failed loci.
fn run_with_fasta(loci_dir: &Path, kmer_getter: &JfKmerGetter, args: &Args) -> Result<(u32, u32), Error> {
    let total = args.sequences.len();
    // First, parse `args.sequences`, so that we do not fail later.
    let mut loci = Vec::with_capacity(total);
    for path_and_name in args.sequences.iter() {
        let (fasta_path, locus) = path_and_name.split_once('=')
            .ok_or_else(|| Error::InvalidInput(format!("Sequence argument {:?} must follow format FASTA=LOCUS_NAME",
                path_and_name)))?;
        interv::check_locus_name(locus)?;
        loci.push((locus, fasta_path));
    }

    let mut succeed = 0;
    for (locus, fasta_path) in loci.into_iter() {
        match process_locus_from_fasta(locus, &Path::new(fasta_path), loci_dir, kmer_getter, args) {
            Ok(()) => succeed += 1,
            Err(e) => log::error!("Error while analyzing locus {}: {}", locus, e.display()),
        }
    }
    Ok((succeed, total as u32))
}

/// Finds exactly one file under `db/jf/*.jf` and creates k-mer getter.
pub(super) fn load_kmer_getter(db_path: &Path, jellyfish_exe: PathBuf) -> Result<JfKmerGetter, Error> {
    let mut jf_filenames = ext::sys::filenames_with_ext(&db_path.join(paths::JF_DIR), "jf")?;
    if jf_filenames.len() != 1 {
        return Err(Error::InvalidInput(format!("There are {} files {}/jf/*.jf (expected 1)",
            db_path.display(), jf_filenames.len())));
    }
    JfKmerGetter::new(jellyfish_exe, jf_filenames.pop().expect("At least one filename must be present"))
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let mut args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();

    let db_path = args.database.as_ref().unwrap();
    let loci_dir = db_path.join(paths::LOCI_DIR);
    ext::sys::mkdir(&loci_dir)?;

    let kmer_getter = load_kmer_getter(db_path, args.jellyfish.clone())?;
    args.moving_window = max(kmer_getter.k(), args.moving_window);

    let (succeed, total) = if args.sequences.is_empty() {
        run_with_vcf(&loci_dir, &kmer_getter, &args)?
    } else {
        run_with_fasta(&loci_dir, &kmer_getter, &args)?
    };

    if succeed == 0 {
        log::error!("Failed at {} loci", total);
    } else if succeed < total {
        log::warn!("Successfully analysed {} loci, failed at {} loci", succeed, total - succeed);
    } else {
        log::info!("Successfully analysed {} loci", succeed);
    }
    log::info!("Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
