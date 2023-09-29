use std::{
    cmp::{min, max},
    io::{self, BufRead, Read, Seek, Write},
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
    err::{Error, validate_param, add_path},
    algo::HashSet,
    ext::{
        self,
        vec::IterExt,
        fmt::PrettyU32,
    },
    seq::{
        self, NamedInterval, Interval, ContigNames, NamedSeq,
        panvcf, fastx,
        contigs::GenomeVersion,
        kmers::{JfKmerGetter, KmerCount},
        wfa::Penalties,
    },
};
#[cfg(feature = "aln")]
use crate::seq::dist;
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

    penalties: Penalties,
    backbone_k: usize,
    /// Alignment accuracy between 0 and 9.
    accuracy: u8,
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
            jf_counts: None,
            reference: None,
            variants: None,
            loci: Vec::new(),
            bed_files: Vec::new(),

            ref_name: None,
            leave_out: Default::default(),
            max_expansion: 50000,
            moving_window: 500,

            penalties: Default::default(),
            backbone_k: 101,
            unknown_frac: 0.0001,

            #[cfg(feature = "aln")]
            accuracy: 2,
            #[cfg(not(feature = "aln"))]
            accuracy: 0,

            #[cfg(feature = "aln")]
            max_divergence: 0.0001,
            #[cfg(not(feature = "aln"))]
            max_divergence: 0.0,

            threads: 8,
            force: false,
            jellyfish: PathBuf::from("jellyfish"),
        }
    }
}

impl Args {
    fn validate(mut self) -> Result<Self, Error> {
        validate_param!(self.database.is_some(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.jf_counts.is_some(), "Jellyfish counts are not provided (see -j/--jf-counts)");
        validate_param!(self.reference.is_some(), "Reference fasta file is not provided (see -r/--reference)");
        validate_param!(self.variants.is_some(), "Pangenome VCF file is not provided (see -v/--variants)");
        validate_param!(!self.loci.is_empty() || !self.bed_files.is_empty(),
            "Target loci are not provided (see -l/--locus and -L/--loci-bed)");

        validate_param!(0.0 <= self.unknown_frac && self.unknown_frac <= 1.0,
            "Unknown fraction ({}) must be within [0, 1]", self.unknown_frac);
        validate_param!(0.0 <= self.max_divergence && self.max_divergence <= 1.0,
            "Maximum divergence ({}) must be within [0, 1]", self.max_divergence);
        validate_param!(self.backbone_k >= 5, "Backbone alignment k-mer size ({}) must be at least 5", self.backbone_k);
        validate_param!(self.accuracy <= crate::seq::wfa::MAX_ACCURACY,
            "Alignment accuracy level ({}) must be between 0 and {}.", self.accuracy, crate::seq::wfa::MAX_ACCURACY);
        #[cfg(not(feature = "aln"))] {
            validate_param!(self.accuracy == 0,
                "Cannot set non-zero alignment accuracy level when `aln` feature is disabled");
            validate_param!(self.max_divergence == 0.0,
                "Cannot set non-zero sequence divergence when `aln` feature is disabled");
        }
        validate_param!(self.moving_window < u32::from(u16::MAX),
            "Moving window ({}) must fit in two bytes", self.moving_window);

        self.jellyfish = ext::sys::find_exe(self.jellyfish)?;
        self.penalties.validate()?;
        Ok(self)
    }
}

fn print_help(extended: bool) {
    const KEY: usize = 16;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Adds target locus/loci to the database.".yellow());

    print!("\n{}", "Usage:".bold());
    println!(" {} add -d db -r ref.fa -j counts.jf [-v vars.vcf.gz] -l/-L loci [args]", super::PROGRAM);
    if !extended {
        println!("\nThis is a {} help message. Please use {} to see the full help.",
            "short".red(), "-H/--full-help".green());
    }

    println!("\n{}", "Input arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Output database directory.",
        "-d, --database".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file.",
        "-r, --reference".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Jellyfish k-mer counts (see README).",
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
        {EMPTY}  with path to locus alleles.",
        "-L, --loci".green(), "FILE".yellow(), "-v".green());

    if extended {
        println!("\n{}", "Allele extraction parameters:".bold());
        println!("    {:KEY$} {:VAL$}  Reference genome name, default: tries to guess.",
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

        println!("\n{}", "Alignment and clustering of alleles:".bold());
        println!("    {:KEY$} {:VAL$}  Penalty for mismatch [{}].",
            "-M, --mismatch".green(), "INT".yellow(), super::fmt_def(defaults.penalties.mismatch));
        println!("    {:KEY$} {:VAL$}  Gap open penalty [{}].",
            "-O, --gap-open".green(), "INT".yellow(), super::fmt_def(defaults.penalties.gap_open));
        println!("    {:KEY$} {:VAL$}  Gap extend penalty [{}].",
            "-E, --gap-extend".green(), "INT".yellow(), super::fmt_def(defaults.penalties.gap_extend));
        println!("    {:KEY$} {:VAL$}  Backbone alignment k-mer size [{}].",
            "-k, --backbone-k".green(), "INT".yellow(), super::fmt_def(defaults.backbone_k));
        println!("    {:KEY$} {:VAL$}  Accuracy level of allele alignments (0-9) [{}].\n\
            {EMPTY}  0: no sequence alignment, 1: fast and inaccurate alignment,\n\
            {EMPTY}  9: slow and accurate alignment.",
            "-a, --accuracy".green(), "INT".yellow(), super::fmt_def(defaults.accuracy));
        println!("    {:KEY$} {:VAL$}  Sequence divergence threshold, used to discard very similar\n\
            {EMPTY}  alleles [{}]. Use 0 to keep all distinct alleles.",
            "-D, --divergence".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.max_divergence));
    }

    println!("\n{}", "Execution arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), super::fmt_def(defaults.threads));
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());
    if extended {
        println!("    {:KEY$} {:VAL$}  Jellyfish executable [{}].",
            "    --jellyfish".green(), "EXE".yellow(), super::fmt_def(defaults.jellyfish.display()));
    }

    println!("\n{}", "Other arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Show this help message.", "-h, --help".green(), "");
    println!("    {:KEY$} {:VAL$}  Show {} help message.", "-H, --full-help".green(), "", "extended".red());
    println!("    {:KEY$} {:VAL$}  Show version.", "-V, --version".green(), "");
}

fn parse_args(argv: &[String]) -> Result<Args, lexopt::Error> {
    if argv.is_empty() {
        print_help(false);
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

            Short('M') | Long("mismatch") => args.penalties.mismatch = parser.value()?.parse()?,
            Short('O') | Long("gap-open") | Long("gap-opening") =>
                args.penalties.gap_open = parser.value()?.parse()?,
            Short('E') | Long("gap-extend") | Long("gap-extension") =>
                args.penalties.gap_extend = parser.value()?.parse()?,
            Short('k') | Long("backbone-k") => args.backbone_k = parser.value()?.parse()?,
            Short('a') | Long("accuracy") => args.accuracy = parser.value()?.parse()?,

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
                print_help(false);
                std::process::exit(0);
            }
            Short('H') | Long("full-help") | Long("hidden-help") => {
                print_help(true);
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
) -> Result<Vec<(NamedInterval, Option<PathBuf>)>, Error>
{
    let mut intervals = Vec::new();
    let mut names = HashSet::default();
    for (name, locus, opt_fasta) in loci.iter() {
        if !names.insert(name.clone()) {
            return Err(Error::InvalidInput(format!("Locus name '{}' appears at least twice", name)));
        }
        let interval = Interval::parse(locus, contigs)?;
        intervals.push((NamedInterval::new(name.clone(), interval)?, opt_fasta.clone()));
    }

    for filename in bed_files.iter() {
        for line in ext::sys::open(filename)?.lines() {
            let line = line.map_err(add_path!(filename))?;
            let mut split = line.trim_end().split('\t');
            let interval = NamedInterval::parse_bed(&mut split, contigs)?;
            if !names.insert(interval.name().to_string()) {
                return Err(Error::InvalidInput(format!("Locus name '{}' appears at least twice", interval.name())));
            }
            let opt_fasta = split.next()
                .map(|s| PathBuf::try_from(s).map_err(|_| Error::InvalidInput(format!("Cannot parse filename {}", s))))
                .transpose()?;
            intervals.push((interval, opt_fasta));
        }
    }

    // Test FASTA filenames.
    for (locus, opt_fasta) in intervals.iter() {
        if let Some(filename) = opt_fasta {
            if !require_seqs {
                return Err(Error::InvalidInput(format!(
                    "FASTA file with locus alleles cannot be provided if variants (-v) were specified (see locus {})",
                    locus.name())));
            }
            if !filename.is_file() {
                return Err(Error::InvalidInput(format!("FASTA file {} does not exist for locus {}",
                    ext::fmt::path(filename), locus.name())));
            }
        } else if require_seqs {
            return Err(Error::InvalidInput(format!(
                "FASTA file is required if variants (-v) was not provided (see locus {})", locus.name())));
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
) -> Result<Option<u32>, Error>
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
        .map(|(&lag_sum, &sum)| f64::from(lag_sum - sum) / divisor).collect();
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
) -> Result<NamedInterval, Error>
{
    assert!(args.max_expansion > 0, "This function should not be called if expansion is forbidden");
    let inner_interval = locus.interval();
    if inner_interval.len() < args.moving_window {
        return Err(Error::RuntimeError(format!("Locus {} is shorter ({}) than the moving window ({})",
            locus.name(), inner_interval.len(), args.moving_window)));
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
            return Err(Error::RuntimeError(format!("Unknown sequence at the locus {}", locus.name())));
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
            return Err(Error::RuntimeError(format!("Unknown sequence at the locus {}", locus.name())));
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
    let new_start = match find_best_boundary::<true>(left_start, inner_start + 1, &left_vars,
            kmer_getter.k(), kmer_counts.get(0), args)? {
        Some(pos) => pos,
        None => return Err(Error::RuntimeError(format!(
            "Cannot expand locus {} to the left due to a long variant overlapping boundary.\n    \
            Try increasing -e/--expand parameter or manually modifying region boundaries.", locus.name()))),
    };

    // Extend region to the right.
    let new_end = match find_best_boundary::<false>(inner_end - 1, right_end, &right_vars,
            kmer_getter.k(), kmer_counts.get(1), args)? {
        Some(pos) => pos + 1,
        None => return Err(Error::RuntimeError(format!(
            "Cannot expand locus {} to the right due to a long variant overlapping boundary.\n    \
            Try increasing -e/--expand parameter or manually modifying region boundaries.", locus.name()))),
    };
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
#[cfg(feature = "aln")]
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
            dist = crate::math::fmt_signif(0.5 * step.dissimilarity, 5)));
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
            if step.dissimilarity == 0.0 || step.dissimilarity <= thresh {
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
    #[allow(unused_variables)] // Not used if feature `aln` is disabled.
    tag: &str,
    entries: Vec<NamedSeq>,
    kmer_getter: &JfKmerGetter,
    args: &Args,
) -> Result<(), Error>
{
    let n_entries = entries.len();
    log::info!("    Writing {} haplotypes to {}/", n_entries, ext::fmt::path(locus_dir));
    log::info!("    Calculating haplotype divergence");
    let divergences: Vec<f64>;
    if args.accuracy > 0 {
        #[cfg(feature = "aln")] {
            let bam_path = locus_dir.join("hap_names.bam");
            divergences = dist::pairwise_divergences(&bam_path, &entries,
                &args.penalties, args.backbone_k, args.accuracy, args.threads)?;
            check_divergencies(tag, &entries, divergences.iter().copied(), args.variants.is_some());
        } #[cfg(not(feature = "aln"))] {
            panic!("Accuracy > 0, but `aln` feature is disabled");
        }
    } else {
        // Need this so that `entries` are not consumed by `move` closure.
        let entries_ref = &entries;
        divergences = (0..n_entries)
            .flat_map(|i| (i + 1..n_entries)
                .map(move |j| if entries_ref[i].seq() == entries_ref[j].seq() { 0.0 } else { 1.0 }))
            .collect();
    }

    log::info!("    Clustering haploypes");
    let nwk_filename = locus_dir.join(paths::LOCUS_DENDROGRAM);
    let nwk_writer = ext::sys::create_file(&nwk_filename)?;
    let keep_seqs = cluster_haplotypes(nwk_writer, &entries, divergences, args.max_divergence)
        .map_err(add_path!(nwk_filename))?;
    let n_filtered = keep_seqs.iter().fold(0, |sum, &keep| sum + usize::from(keep));
    if n_filtered == n_entries {
        log::info!("        Keep all sequences after clustering");
    } else {
        log::info!("        Discard {} sequences after clustering", n_entries - n_filtered);
    }

    let filt_fasta_readername = locus_dir.join(paths::LOCUS_FASTA);
    let mut filt_fasta_writer = ext::sys::create_gzip(&filt_fasta_readername)?;
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

/// Checks sequences for Ns, minimal size and for equal boundaries.
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

/// Adds locus to the database.
fn add_locus(
    mut locus: NamedInterval,
    alleles_fasta: &Option<PathBuf>,
    locus_dir: &Path,
    fasta_reader: &mut fasta::IndexedReader<impl Read + Seek>,
    vcf_data: &mut Option<(bcf::IndexedReader, panvcf::HaplotypeNames)>,
    kmer_getter: &JfKmerGetter,
    args: &Args,
) -> Result<(), Error>
{
    log::info!("Analyzing {} ({})", locus.name().bold(), locus.interval());
    // VCF file was provided.
    if let Some((vcf_file, hap_names)) = vcf_data {
        if args.max_expansion > 0 {
            locus = expand_locus(&locus, fasta_reader, vcf_file, hap_names, kmer_getter, args)?;
        }
    }
    let ref_seq = locus.interval().fetch_seq(fasta_reader)?;
    if seq::has_n(&ref_seq) {
        return Err(Error::RuntimeError(format!("Locus {} sequence contains Ns", locus)));
    }

    // Write reference coordinates to the BED file.
    let locus_bed = locus_dir.join(paths::LOCUS_BED);
    let mut bed_writer = ext::sys::create_file(&locus_bed)?;
    writeln!(bed_writer, "{}", locus.bed_fmt()).map_err(add_path!(locus_bed))?;

    // Load sequences.
    let allele_seqs = if let Some((vcf_file, hap_names)) = vcf_data {
        panvcf::reconstruct_sequences(locus.interval(), &ref_seq, vcf_file, hap_names, args.unknown_frac)?
    } else if let Some(fasta_filename) = alleles_fasta {
        let mut fasta_reader = fastx::Reader::from_path(fasta_filename)?;
        fasta_reader.read_all()?
    } else {
        unreachable!("Either VCF file or alleles FASTA must be specified")
    };
    check_sequences(&allele_seqs, locus.name())?;
    process_haplotypes(locus_dir, locus.name(), allele_seqs, kmer_getter, args)
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let mut args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();

    let db_dir = args.database.as_ref().unwrap();
    ext::sys::mkdir(&db_dir)?;
    let loci_dir = db_dir.join(paths::LOCI_DIR);
    ext::sys::mkdir(&loci_dir)?;

    let (contigs, mut fasta_reader) = ContigNames::load_indexed_fasta("REF", args.reference.as_ref().unwrap())?;
    let loci = load_loci(&contigs, &args.loci, &args.bed_files, args.variants.is_none())?;
    let mut vcf_data = if let Some(vcf_filename) = &args.variants {
        let mut vcf_file = bcf::IndexedReader::from_path(vcf_filename)?;
        let ref_name = args.ref_name.as_ref().map(String::clone)
            .or_else(|| GenomeVersion::guess(&contigs).map(|ver| ver.to_str().to_owned()))
            .ok_or_else(|| Error::RuntimeError(
                "Cannot guess reference genome name, please provide using -g".to_owned()))?;
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
            log::error!("Error while analyzing locus {}: {}", locus, e.display());
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
