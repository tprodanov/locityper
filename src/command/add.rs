use std::{
    cmp::{min, max},
    io::{BufRead, Read, Seek, Write},
    sync::Arc,
    path::{Path, PathBuf},
    time::Instant,
};
use bio::io::fasta;
use htslib::bcf::{
    self,
    Read as VcfRead,
    record::Record as VcfRecord,
};
use colored::Colorize;
use const_format::str_repeat;
use crate::{
    err::{error, validate_param, add_path},
    algo::HashSet,
    ext::{
        self,
        TriangleMatrix,
        vec::IterExt,
        fmt::PrettyU32,
    },
    seq::{
        self, NamedInterval, Interval, ContigNames, NamedSeq,
        panvcf, fastx, div,
        contigs::GenomeVersion,
        kmers::{self, Kmer},
        counts::{JfKmerGetter, KmerCount},
    },
};
use super::paths;

struct Args {
    database: Option<PathBuf>,
    reference: Option<PathBuf>,
    jf_counts: Option<PathBuf>,
    variants: Option<PathBuf>,
    loci: Vec<(String, String, Option<PathBuf>)>,
    bed_files: Vec<PathBuf>,

    ref_name: Option<String>,
    leave_out: HashSet<String>,
    max_expansion: u32,
    moving_window: u32,
    ignore_overlaps: bool,

    div_k: u8,
    div_w: u8,
    unknown_frac: f64,

    threads: u16,
    force: bool,
    jellyfish: PathBuf,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            database: None,
            jf_counts: None,
            reference: None,
            variants: None,
            loci: Vec::new(),
            bed_files: Vec::new(),

            ref_name: None,
            leave_out: Default::default(),
            max_expansion: 50000,
            moving_window: 500,
            ignore_overlaps: false,

            div_k: 15,
            div_w: 15,
            unknown_frac: 0.0001,

            threads: 8,
            force: false,
            jellyfish: PathBuf::from("jellyfish"),
        }
    }
}

impl Args {
    fn validate(mut self) -> crate::Result<Self> {
        validate_param!(self.database.is_some(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.jf_counts.is_some(), "Jellyfish counts are not provided (see -j/--jf-counts)");
        validate_param!(self.reference.is_some(), "Reference fasta file is not provided (see -r/--reference)");
        validate_param!(!self.loci.is_empty() || !self.bed_files.is_empty(),
            "Target loci are not provided (see -l/--locus and -L/--loci-bed)");

        validate_param!(0.0 <= self.unknown_frac && self.unknown_frac <= 1.0,
            "Unknown fraction ({}) must be within [0, 1]", self.unknown_frac);
        validate_param!(0 < self.div_k && self.div_k <= u64::MAX_KMER_SIZE,
            "k-mer size ({}) must be between 1 and {}", self.div_k, u64::MAX_KMER_SIZE);
        validate_param!(0 < self.div_w && self.div_w <= kmers::MAX_MINIMIZER_W,
            "Minimizer window ({}) must be between 1 and {}", self.div_w, kmers::MAX_MINIMIZER_W);
        validate_param!(self.moving_window < u32::from(u16::MAX),
            "Moving window ({}) must fit in two bytes", self.moving_window);

        self.jellyfish = ext::sys::find_exe(self.jellyfish)?;
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Adds target locus/loci to the database.".yellow());

    print!("\n{}", "Usage:".bold());
    println!(" {} add -d db -r ref.fa -j counts.jf [-v vars.vcf.gz] -l/-L loci [args]", super::PROGRAM);

    println!("\n{}", "Input arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Output database directory.",
        "-d, --database".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file.",
        "-r, --reference".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Jellyfish k-mer counts (see documentation).",
        "-j, --jf-counts".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Input VCF file, encoding variation across pangenome samples.",
        "-v, --variants".green(), "FILE".yellow());

    println!("\n{} (may be repeated multiple times)", "Target loci:".bold());
    println!("    {:KEY$} {}\n\
        {EMPTY}  Locus name and coordinates. If VCF ({}) was not provided,\n\
        {EMPTY}  FASTA with locus alleles is required as a third argument.",
        "-l, --locus".green(), "NAME REGION [FILE]".yellow(), "-v".green());
    println!("    {:KEY$} {:VAL$}  BED file with loci coordinates. Fourth column: locus name.\n\
        {EMPTY}  If VCF ({}) was not provided, fifth column is required\n\
        {EMPTY}  with the path to locus alleles.",
        "-L, --loci".green(), "FILE".yellow(), "-v".green());

    println!("\n{}", "Allele extraction parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Name of the reference haplotype [default: tries to guess].",
        "-g, --genome".green(), "STR".yellow());
    println!("    {:KEY$} {:VAL$}  If needed, expand loci boundaries by at most {} bp [{}].",
        "-e, --expand".green(), "INT".yellow(), "INT".yellow(),
        super::fmt_def(PrettyU32(defaults.max_expansion)));
    println!("    {:KEY$} {:VAL$}  Select best locus boundary based on k-mer frequencies in\n\
        {EMPTY}  moving windows of size {} bp [{}].",
        "-w, --window".green(), "INT".yellow(), "INT".yellow(),
        super::fmt_def(PrettyU32(defaults.moving_window)));
    println!("    {:KEY$} {:VAL$}  Allow this fraction of unknown nucleotides per allele [{}]\n\
        {EMPTY}  (relative to the allele length). Variants that have no known\n\
        {EMPTY}  variation in the input VCF pangenome are ignored.",
        "-u, --unknown".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.unknown_frac));
    println!("    {:KEY$} {:VAL$}  Leave out sequences with specified names.",
        "    --leave-out".green(), "STR+".yellow());
    println!("    {}   {} (k,w)-minimizers for sequence divergence calculation [{} {}].",
        "-m, --minimizer".green(), "INT INT".yellow(),
        super::fmt_def(defaults.div_k), super::fmt_def(defaults.div_w));
    println!("    {:KEY$} {:VAL$}  Ignore overlapping variants. Default: fail with error.\n\
        {EMPTY}  Of several overlapping variants only the first one is used.",
        "    --ignore-overl".green(), super::flag());

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
        std::process::exit(1);
    }
    use lexopt::prelude::*;
    let mut args = Args::default();
    let mut parser = lexopt::Parser::from_args(argv);

    while let Some(arg) = parser.next()? {
        match arg {
            Short('d') | Long("db") | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('j') | Long("jf-counts") => args.jf_counts = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.reference = Some(parser.value()?.parse()?),
            Short('v') | Long("vcf") | Long("variants") => args.variants = Some(parser.value()?.parse()?),

            Short('l') | Long("locus") => {
                // There are 2 or 3 values, take first with `value()`.
                let name = parser.value()?.parse()?;
                let mut values = parser.values()?;
                // Collect second value from `values()`, as it requires at least one argument.
                let region = values.next().expect("First value must be present").parse()?;
                let path = values.next().map(|val| val.parse()).transpose()?;
                args.loci.push((name, region, path));
            }
            Short('L') | Long("loci") | Long("loci-bed") => args.bed_files.push(parser.value()?.parse()?),

            Short('g') | Long("genome") => args.ref_name = Some(parser.value()?.parse()?),
            Long("leave-out") | Long("leaveout") => {
                for val in parser.values()? {
                    args.leave_out.insert(val.parse()?);
                }
            }
            Short('e') | Long("expand") => args.max_expansion = parser.value()?.parse::<PrettyU32>()?.get(),
            Short('w') | Long("window") => args.moving_window = parser.value()?.parse::<PrettyU32>()?.get(),
            Short('u') | Long("unknown") => args.unknown_frac = parser.value()?.parse()?,
            Short('m') | Long("minimizer") | Long("minimizers") =>
            {
                args.div_k = parser.value()?.parse()?;
                args.div_w = parser.value()?.parse()?;
            }
            Long("ignore-overl") | Long("ignore-overlaps") => args.ignore_overlaps = true,

            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,
            Long("jellyfish") => args.jellyfish = parser.value()?.parse()?,

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

/// Loads named intervals from a list of loci and a list of BED files. Names must not repeat.
/// Additionally, returns optional FASTA filename for each locus.
fn load_loci(
    contigs: &Arc<ContigNames>,
    loci: &[(String, String, Option<PathBuf>)],
    bed_files: &[PathBuf],
    require_seqs: bool,
) -> crate::Result<Vec<(NamedInterval, Option<PathBuf>)>>
{
    let mut intervals = Vec::new();
    let mut names = HashSet::default();
    for (name, locus, opt_fasta) in loci.iter() {
        if !names.insert(name.clone()) {
            return Err(error!(InvalidInput, "Locus name '{}' appears at least twice", name));
        }
        let interval = Interval::parse(locus, contigs)?;
        intervals.push((NamedInterval::new(name.clone(), interval)?, opt_fasta.clone()));
    }

    for filename in bed_files.iter() {
        let dirname = ext::sys::parent_unless_redirect(filename);
        for line in ext::sys::open(filename)?.lines() {
            let line = line.map_err(add_path!(filename))?;
            let mut split = line.trim_end().split('\t');
            let interval = NamedInterval::parse_bed(&mut split, contigs)?;
            if !names.insert(interval.name().to_string()) {
                return Err(error!(InvalidInput, "Locus name '{}' appears at least twice", interval.name()));
            }
            let opt_fasta = match split.next() {
                None => None,
                Some(s) => {
                    let path = &Path::new(s);
                    // Join with parent directory, if any.
                    Some(dirname.map(|d| d.join(path)).unwrap_or_else(|| path.to_path_buf()))
                },
            };
            intervals.push((interval, opt_fasta));
        }
    }

    // Test FASTA filenames.
    for (locus, opt_fasta) in intervals.iter() {
        if let Some(filename) = opt_fasta {
            if !require_seqs {
                return Err(error!(InvalidInput,
                    "FASTA file with locus alleles cannot be provided if variants (-v) were specified (see locus {})",
                    locus.name()));
            }
            if !filename.exists() || filename.is_dir() {
                return Err(error!(InvalidInput, "FASTA file {} does not exist for locus {}",
                    ext::fmt::path(filename), locus.name()));
            }
        } else if require_seqs {
            return Err(error!(InvalidInput,
                "FASTA file is required if variants (-v) was not provided (see locus {})", locus.name()));
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
    start: u32,
    end: u32,
    vars: &[VcfRecord],
    k: u32,
    kmer_counts: &[KmerCount],
    args: &Args,
) -> crate::Result<Option<u32>>
{
    if start == end {
        if vars.iter().any(|var| var.pos() as u32 <= start && end <= var.end() as u32) {
            return Ok(None);
        } else {
            return Ok(Some(start));
        }
    }

    let cumul_uniq_kmers: Vec<u32> = IterExt::cumul_sums(kmer_counts.iter().map(|&count| u32::from(count <= 1)));
    let kmers_per_window = args.moving_window + 1 - k;
    let divisor = f64::from(kmers_per_window);
    let mut weights: Vec<f64> = cumul_uniq_kmers.iter().zip(&cumul_uniq_kmers[kmers_per_window as usize..])
        .map(|(&lag_sum, &sum)| f64::from(sum - lag_sum) / divisor).collect();
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

    // Furthest point to the region boundary is penalized by 20%.
    const WEIGHT_DROP: f64 = 0.2;
    let per_bp_drop = WEIGHT_DROP / f64::from(args.max_expansion);
    let (i, maxval) = if LEFT {
        weights.iter_mut().rev().enumerate().for_each(|(i, w)| *w -= *w * per_bp_drop * i as f64);
        // Last argmax.
        IterExt::arg_optimal(weights.iter().copied(), |opt, e| opt >= e)
    } else {
        weights.iter_mut().enumerate().for_each(|(i, w)| *w -= *w * per_bp_drop * i as f64);
        // First argmax.
        IterExt::arg_optimal(weights.iter().copied(), |opt, e| opt > e)
    };

    if maxval == 0.0 {
        Ok(None)
    } else {
        Ok(Some(start + i as u32))
    }
}

/// Expands locus boundaries if necessary.
fn expand_locus(
    locus: &NamedInterval,
    fasta_reader: &mut fasta::IndexedReader<impl Read + Seek>,
    vcf_file: &mut bcf::IndexedReader,
    hap_names: &panvcf::HaplotypeNames,
    kmer_getter: &JfKmerGetter,
    args: &Args,
) -> crate::Result<NamedInterval>
{
    assert!(args.max_expansion > 0, "This function should not be called if expansion is forbidden");
    let inner_interval = locus.interval();
    if inner_interval.len() < args.moving_window {
        return Err(error!(RuntimeError, "Locus {} is shorter ({}) than the moving window ({})",
            locus.name(), inner_interval.len(), args.moving_window));
    }
    let contig_name = inner_interval.contig_name();
    let contig_len = inner_interval.contigs().get_len(inner_interval.contig_id());
    let (inner_start, inner_end) = inner_interval.range();
    // Because we need to calculate the average number of unique k-mers per each moving window,
    // we set left end further, so that we have the same number of moving windows as the potential position shift.
    let mut left_start = inner_start.saturating_sub(args.max_expansion);
    let left_end = inner_start + args.moving_window;
    let mut left_seq = seq::fetch_seq(fasta_reader, contig_name, left_start.into(), left_end.into())?;

    // Same with the right boundary, we extend right start further left so that we have enough moving windows.
    let right_start = inner_end.checked_sub(args.moving_window).unwrap();
    let mut right_end = min(inner_end + args.max_expansion, contig_len);
    let mut right_seq = seq::fetch_seq(fasta_reader, contig_name, right_start.into(), right_end.into())?;

    // Crop left region at the last N.
    if let Some(shift) = left_seq.iter().rposition(|&nt| nt == b'N') {
        left_start += shift as u32 + 1;
        if left_start > inner_start {
            return Err(error!(RuntimeError, "Unknown sequence at the locus {}", locus.name()));
        }
        log::warn!("    [{}] Unknown nucleotide {} bp to the left ({}:{})",
            locus.name(), inner_start - left_start + 1, contig_name, left_start + 1);

        let old_len = left_seq.len();
        // Keep only sequence after N.
        left_seq.copy_within(shift + 1.., 0);
        left_seq.truncate(old_len - shift - 1);
        assert_eq!(left_seq.len() as u32, left_end - left_start);
    }
    // Crop right region at the first N.
    if let Some(shift) = right_seq.iter().position(|&nt| nt == b'N') {
        right_end = right_start + shift as u32;
        if right_end < inner_end {
            return Err(error!(RuntimeError, "Unknown sequence at the locus {}", locus.name()));
        }
        log::warn!("    [{}] Unknown nucleotide {} bp to the right ({}:{})",
            locus.name(), right_end - inner_end + 1, contig_name, right_end + 1);
        right_seq.truncate(shift);
        assert_eq!(right_seq.len() as u32, right_end - right_start);
    }

    // Fetch variants to the left and right of the boundary.
    let vcf_rid = vcf_file.header().name2rid(contig_name.as_bytes())?;
    vcf_file.fetch(vcf_rid, u64::from(left_start), Some(u64::from(inner_start + 1)))?;
    let left_vars = panvcf::filter_variants(vcf_file, hap_names)?;
    vcf_file.fetch(vcf_rid, u64::from(inner_end - 1), Some(u64::from(right_end)))?;
    let right_vars = panvcf::filter_variants(vcf_file, hap_names)?;

    let kmer_counts = kmer_getter.fetch([left_seq, right_seq])?;
    // Extend region to the left.
    let Some(new_start) = find_best_boundary::<true>(left_start, inner_start + 1, &left_vars,
            kmer_getter.k(), kmer_counts.get(0), args)? else {
        return Err(error!(RuntimeError,
            "Cannot expand locus {} to the left due to a long variant overlapping boundary.\n    \
            Try increasing -e/--expand parameter or manually modifying region boundaries.", locus.name()))
    };

    // Extend region to the right.
    let Some(mut new_end) = find_best_boundary::<false>(inner_end - 1, right_end, &right_vars,
            kmer_getter.k(), kmer_counts.get(1), args)? else {
        return Err(error!(RuntimeError,
            "Cannot expand locus {} to the right due to a long variant overlapping boundary.\n    \
            Try increasing -e/--expand parameter or manually modifying region boundaries.", locus.name()))
    };
    new_end += 1;
    if new_start != inner_start || new_end != inner_end {
        let new_interval = inner_interval.create_at_same_contig(new_start, new_end);
        log::info!("    Extending locus by {} bp left and {} bp right -> {}",
            inner_start - new_start, new_end - inner_end, new_interval);
        NamedInterval::new(locus.name().to_owned(), new_interval)
    } else {
        Ok(locus.clone())
    }
}

/// Check divergencies and warns if they are too high.
fn check_divergencies(tag: &str, entries: &[NamedSeq], divergences: &TriangleMatrix<(u32, f64)>, from_vcf: bool) {
    let mut count = 0;
    let mut highest = 0.0;
    let mut highest_i = 0;
    let mut highest_j = 0;
    for ((i, j), &(_, diverg)) in TriangleMatrix::indices(divergences.side()).zip(divergences.iter()) {
        if diverg >= 0.2 {
            count += 1;
            if diverg > highest {
                highest = diverg;
                highest_i = i;
                highest_j = j;
            }
        }
    }
    if count > 0 {
        log::warn!("    [{}] {} allele pairs with high divergence, highest {:.5} ({} and {})", tag, count,
            highest, entries[highest_i].name(), entries[highest_j].name());
        if highest > 0.5 && !from_vcf {
            log::warn!("    Extremely high sequence divergence, please check if allele sequences are accurate");
        }
    }
}

// /// Cluster alleles using `kodama` crate
// /// (which, in turn, is based on this paper https://arxiv.org/pdf/1109.2378.pdf).
// ///
// /// Clusters with divergence not exceeding the threshold are joined.
// /// Returns boolean vector: does the sequence remain after mergin?
// fn cluster_alleles(
//     mut nwk_writer: impl Write,
//     entries: &[NamedSeq],
//     divergences: TriangleMatrix<f64>,
//     thresh: f64,
// ) -> io::Result<Vec<bool>> {
//     let n = entries.len();
//     let total_clusters = 2 * n - 1;
//     // Use Complete method, meaning that we track maximal distance between two points between two clusters.
//     // This is done to cut as little as needed.
//     let dendrogram = kodama::linkage(&mut divergences.take_linear(), n, kodama::Method::Complete);
//     let mut clusters_nwk = Vec::with_capacity(total_clusters);
//     // Cluster representatives.
//     let mut cluster_repr = Vec::with_capacity(total_clusters);
//     clusters_nwk.extend(entries.iter().map(|entry| entry.name().to_owned()));
//     cluster_repr.extend(0..n);

//     let steps = dendrogram.steps();
//     for step in steps.iter() {
//         let i = step.cluster1;
//         let j = step.cluster2;
//         clusters_nwk.push(format!("({}:{dist},{}:{dist})", &clusters_nwk[i], &clusters_nwk[j],
//             dist = crate::math::fmt_signif(0.5 * step.dissimilarity, 5)));
//         let size1 = if i < n { 1 } else { steps[i - n].size };
//         let size2 = if j < n { 1 } else { steps[j - n].size };
//         cluster_repr.push(cluster_repr[if size1 >= size2 { i } else { j }]);
//     }
//     assert_eq!(clusters_nwk.len(), total_clusters);
//     writeln!(nwk_writer, "{};", clusters_nwk.last().unwrap())?;

//     let mut queue = vec![total_clusters - 1];
//     let mut keep_seqs = vec![false; n];
//     while let Some(i) = queue.pop() {
//         if i < n {
//             keep_seqs[i] = true;
//         } else {
//             let step = &steps[i - n];
//             if step.dissimilarity == 0.0 || step.dissimilarity <= thresh {
//                 keep_seqs[cluster_repr[i]] = true;
//             } else {
//                 queue.push(step.cluster1);
//                 queue.push(step.cluster2);
//             }
//         }
//     }
//     Ok(keep_seqs)
// }

/// Discards identical haplotypes, and, if any, writes their names into a text files.
fn discard_identical(entries: Vec<NamedSeq>, locus_dir: &Path) -> crate::Result<Vec<NamedSeq>> {
    let mut selected: Vec<NamedSeq> = Vec::with_capacity(entries.len());
    let mut n_discarded = 0;
    let mut disc_names: Vec<Vec<String>> = Vec::with_capacity(entries.len());

    'outer: for entry in entries.into_iter() {
        for (i, entry0) in selected.iter().enumerate() {
            if entry0.seq() == entry.seq() {
                disc_names[i].push(entry.take_name());
                n_discarded += 1;
                continue 'outer;
            }
        }
        selected.push(entry);
        disc_names.push(Vec::new());
    }

    if n_discarded > 0 {
        log::info!("    Discarded {} duplicate haplotypes", n_discarded);
        let filename = locus_dir.join("discarded_haplotypes.txt");
        let mut f = ext::sys::create_file(&filename)?;
        for (entry, discarded) in selected.iter().zip(&disc_names) {
            if discarded.is_empty() {
                continue;
            }
            write!(f, "{} = ", entry.name()).map_err(add_path!(filename))?;
            for (i, name) in discarded.iter().enumerate() {
                if i > 0 {
                    write!(f, ", ").map_err(add_path!(filename))?;
                }
                write!(f, "{}", name).map_err(add_path!(filename))?;
            }
            writeln!(f).map_err(add_path!(filename))?;
        }
    }
    Ok(selected)
}

/// Process alleles: write FASTA, PAF, kmers files, cluster sequences.
fn process_alleles(
    locus: &str,
    locus_dir: &Path,
    mut ref_seq: Vec<u8>,
    entries: Vec<NamedSeq>,
    kmer_getter: &JfKmerGetter,
    args: &Args,
) -> crate::Result<()>
{
    let entries = discard_identical(entries, locus_dir)?;
    let n_entries = entries.len();
    let mut seqs = Vec::with_capacity(n_entries + 1);
    let fasta_filename = locus_dir.join(paths::LOCUS_FASTA);
    let mut fasta_writer = ext::sys::create_gzip(&fasta_filename)?;
    for entry in entries.iter() {
        seq::write_fasta(&mut fasta_writer, entry.name().as_bytes(), entry.seq()).map_err(add_path!(fasta_filename))?;
        seqs.push(entry.seq().to_vec());
    }
    std::mem::drop(fasta_filename);

    log::info!("    Calculating sequence divergence for {} alleles", n_entries);
    let all_pairs: Vec<_> = TriangleMatrix::indices_u32(n_entries).collect();
    let divergences = div::minimizer_divergences(&entries, &all_pairs, args.div_k, args.div_w, args.threads);
    let divergences = TriangleMatrix::from_linear(n_entries, divergences);
    check_divergencies(locus, &entries, &divergences, args.variants.is_some());
    let dist_filename = locus_dir.join(paths::DISTANCES);
    let dist_file = ext::sys::create_file(&dist_filename)?;
    div::write_divergences(args.div_k, args.div_w, &divergences, dist_file).map_err(add_path!(dist_filename))?;

    log::info!("    Counting k-mers");
    // Replace Ns with As in the reference sequence.
    let ref_n_runs = seq::n_runs(&ref_seq);
    for &(start, end) in ref_n_runs.iter() {
        (&mut ref_seq[start as usize..end as usize]).fill(b'A');
    }

    // Calculate k-mer counts for the reference sequence as well.
    seqs.push(ref_seq);
    let mut kmer_counts = kmer_getter.fetch(seqs.clone())?;
    let ref_seq = seqs.pop().unwrap();
    let mut ref_counts = kmer_counts.pop();

    // Replace k-mer counts with 0s where Ns were observed.
    let k = kmer_getter.k();
    let ref_counts_size = ref_counts.len();
    for &(start, end) in ref_n_runs.iter() {
        (&mut ref_counts[(start + 1).saturating_sub(k) as usize..min(end as usize, ref_counts_size)]).fill(0);
    }
    let off_target_counts = kmer_counts.off_target_counts(&seqs, &ref_seq, &ref_counts, ref_n_runs.is_empty());

    let kmers_filename = locus_dir.join(paths::KMERS);
    let mut kmers_writer = ext::sys::create_lz4_slow(&kmers_filename)?;
    // Consecutively save off-target counts and regular counts (if they will sometimes be useful later).
    off_target_counts.save(&mut kmers_writer).map_err(add_path!(kmers_filename))?;
    kmer_counts.save(&mut kmers_writer).map_err(add_path!(kmers_filename))?;
    super::write_success_file(locus_dir.join(paths::SUCCESS))?;
    Ok(())
}

/// Checks sequences for Ns, minimal size and for equal boundaries.
fn check_sequences(seqs: &[NamedSeq], locus: &NamedInterval, ref_seq: Option<&[u8]>) -> crate::Result<()> {
    if seqs.len() < 2 {
        return Err(error!(InvalidData, "Less than two alleles available for locus {}", locus.name()));
    }
    let min_size = seqs.iter().map(|s| s.len()).min().unwrap();
    if min_size < 1000 {
        return Err(error!(InvalidData, "Locus alleles are too short for locus {} (shortest: {} bp)",
            locus, min_size));
    } else if min_size < 10000 {
        log::warn!("[{}] Locus alleles may be too short (shortest: {} bp)", locus.name(), min_size);
    }

    const AFFIX_SIZE: usize = 5;
    let seq0 = seqs[0].seq();
    let prefix = &seq0[..AFFIX_SIZE];
    let suffix = &seq0[seq0.len() - AFFIX_SIZE..];
    if seqs[1..].iter().map(NamedSeq::seq)
            .any(|s| &s[..AFFIX_SIZE] != prefix || &s[s.len() - AFFIX_SIZE..] != suffix) {
        log::warn!("[{}] Allele sequences differ at the boundary", locus.name());
    }

    if let Some(rseq) = ref_seq {
        let prefix = &rseq[..AFFIX_SIZE];
        let suffix = &rseq[rseq.len() - AFFIX_SIZE..];
        let d = seqs.iter().map(NamedSeq::seq)
            .filter(|seq| &seq[..AFFIX_SIZE] != prefix || &seq[seq.len() - AFFIX_SIZE..] != suffix).count();
        let n = seqs.len();
        if d > 0 {
            log::log!(
                if d == n { log::Level::Error } else { log::Level::Warn },
                "[{}] {} of the alleles ({}/{}) do not match the reference sequence at {}. Continuing", locus.name(),
                if d == n { "All" } else if d >= n / 2 { "Most" } else { "Some" },
                d, n, locus.interval());
        }
    }
    Ok(())
}

/// Discard alleles, whose names are in `leave_out`.
/// Additionally, discard sequences `name[._-][0-9]` where `name` is in `leave_out`.
fn discard_leave_out_alleles(alleles: Vec<NamedSeq>, leave_out: &HashSet<String>) -> Vec<NamedSeq> {
    if leave_out.is_empty() {
        return alleles;
    }
    let mut filt_alleles = Vec::with_capacity(alleles.len());
    let mut left_out = Vec::new();
    for allele in alleles.into_iter() {
        let name = allele.name();
        let bytes = name.as_bytes();
        let n = bytes.len();
        if leave_out.contains(name) ||
                (n > 2 && leave_out.contains(&std::str::from_utf8(&bytes[..n - 2]).unwrap() as &str)
                && b'0' <= bytes[n - 1] && bytes[n - 1] <= b'9'
                && (bytes[n - 2] == b'.' || bytes[n - 2] == b'_' || bytes[n - 2] == b'-')) {
            left_out.push(allele.take_name());
            continue;
        }
        filt_alleles.push(allele);
    }
    let discarded = left_out.len();
    if discarded > 0 {
        if left_out.len() > 5 {
            left_out.truncate(5);
            left_out.push("...".to_owned());
        }
        log::warn!("    Leave out {} alleles ({})", discarded, left_out.join(", "));
    } else {
        log::warn!("Zero matches between leave-out and FASTA alleles");
    }
    filt_alleles
}

/// Adds locus to the database.
fn add_locus(
    mut locus: NamedInterval,
    alleles_fasta: &Option<PathBuf>,
    locus_dir: &Path,
    fasta_reader: &mut fasta::IndexedReader<impl Read + Seek>,
    vcf_data: &mut Option<(bcf::IndexedReader, panvcf::HaplotypeNames)>,
    kmer_getter: &JfKmerGetter,
    args: &Args,
) -> crate::Result<()>
{
    log::info!("Analyzing {} ({})", locus.name().bold(), locus.interval());
    // VCF file was provided.
    if let Some((vcf_file, hap_names)) = vcf_data {
        if args.max_expansion > 0 {
            locus = expand_locus(&locus, fasta_reader, vcf_file, hap_names, kmer_getter, args)?;
        }
    }
    let ref_seq = locus.interval().fetch_seq(fasta_reader)?;

    // Write reference coordinates to the BED file.
    let locus_bed = locus_dir.join(paths::LOCUS_BED);
    let mut bed_writer = ext::sys::create_file(&locus_bed)?;
    writeln!(bed_writer, "{}", locus.bed_fmt()).map_err(add_path!(locus_bed))?;

    // Load sequences.
    let allele_seqs = if let Some((vcf_file, hap_names)) = vcf_data {
        panvcf::reconstruct_sequences(locus.interval(), &ref_seq, vcf_file, hap_names,
            args.unknown_frac, args.ignore_overlaps)?
    } else if let Some(fasta_filename) = alleles_fasta {
        let mut fasta_reader = fastx::Reader::from_path(fasta_filename)?;
        let alleles = fasta_reader.read_all()?;
        log::info!("FASTA file contains {} alleles", alleles.len());
        discard_leave_out_alleles(alleles, &args.leave_out)
    } else {
        unreachable!("Either VCF file or alleles FASTA must be specified")
    };
    let obtained_count = allele_seqs.len();
    let allele_seqs: Vec<_> = allele_seqs.into_iter().filter(|entry| !seq::has_n(entry.seq())).collect();
    let without_ns = allele_seqs.len();
    if without_ns < obtained_count {
        log::warn!("    Removed {}/{} alleles with Ns", obtained_count - without_ns, obtained_count);
    }
    check_sequences(&allele_seqs, &locus, if alleles_fasta.is_some() { Some(&ref_seq) } else { None })?;
    process_alleles(locus.name(), locus_dir, ref_seq, allele_seqs, kmer_getter, args)
}

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let mut args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();

    let db_dir = args.database.as_ref().unwrap();
    ext::sys::mkdir(&db_dir)?;
    let loci_dir = db_dir.join(paths::LOCI_DIR);
    ext::sys::mkdir(&loci_dir)?;

    let (contigs, mut fasta_reader) = ContigNames::load_indexed_fasta("REF", args.reference.as_ref().unwrap())?;
    let contigs = Arc::new(contigs);
    let loci = load_loci(&contigs, &args.loci, &args.bed_files, args.variants.is_none())?;
    let mut vcf_data = if let Some(vcf_filename) = &args.variants {
        let mut vcf_file = bcf::IndexedReader::from_path(vcf_filename)?;
        let ref_name = args.ref_name.as_ref().map(String::clone)
            .or_else(|| GenomeVersion::guess(&contigs).map(|ver| ver.to_str().to_owned()))
            .ok_or_else(|| error!(RuntimeError, "Cannot guess reference genome name, please provide using -g"))?;
        let hap_names = panvcf::HaplotypeNames::new(&mut vcf_file, &ref_name, &args.leave_out)?;
        Some((vcf_file, hap_names))
    } else {
        None
    };
    let kmer_getter = JfKmerGetter::new(args.jellyfish.clone(), args.jf_counts.clone().unwrap())?;
    args.moving_window = max(kmer_getter.k(), args.moving_window);

    let total = loci.len();
    let mut failed = 0;
    for (locus, alleles_fasta) in loci.into_iter() {
        let locus_dir = loci_dir.join(locus.name());
        if !super::Rerun::from_force(args.force).prepare_dir(&locus_dir)? {
            continue;
        }
        let res = add_locus(locus.clone(), &alleles_fasta, &locus_dir, &mut fasta_reader,
            &mut vcf_data, &kmer_getter, &args);
        if let Err(e) = res {
            log::error!("Error while analyzing locus {}:\n        {}", locus, e.display());
            failed += 1;
        }
    }

    let succeed = total - failed;
    if succeed == 0 {
        log::error!("Failed to add {} loci", failed);
    } else if failed > 0 {
        log::warn!("Successfully added {} loci, failed to add {} loci", succeed, failed);
    } else {
        log::info!("Successfully added {} loci", succeed);
    }
    log::info!("Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
