use std::{
    path::{Path, PathBuf},
    time::Instant,
    io::{self, Write, BufRead},
    fmt::Write as FmtWrite,
    ffi::OsStr,
};
use colored::Colorize;
use crate::{
    ext::{self, TriangleMatrix},
    algo::HashSet,
    seq::{
        fastx, div, ContigId, ContigNames,
        counts::KmerCounts,
    },
    err::{validate_param, error, add_path},
    math::PowerMean,
};
use super::{
    paths,
    genotype::LocusData,
};

struct Args {
    input: Option<PathBuf>,
    output: Option<PathBuf>,
    alignments: String,
    subset_loci: HashSet<String>,

    skip_tree: bool,
    field: String,
    threshold: f64,
    n_clusters: Option<usize>,
    power: PowerMean,
    force: bool,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: None,
            output: None,
            alignments: "haplotypes.paf.gz".to_string(),
            subset_loci: HashSet::default(),

            skip_tree: false,
            field: "dv".to_string(),
            threshold: 0.0002,
            n_clusters: None,
            power: PowerMean::Pow(2),
            force: false,
        }
    }
}

impl Args {
    fn validate(mut self) -> crate::Result<Self> {
        validate_param!(self.input.is_some(), "Input database is not provided (see -i/--input)");
        validate_param!(self.output.is_some(), "Output database is not provided (see -o/--output)");
        validate_param!(!self.field.contains(':'), "PAF divergence field ({}) must not contain :",
            self.field);
        validate_param!(self.threshold >= 0.0, "Divergence threshold ({}) should be non-negative",
            self.threshold);
        validate_param!(self.n_clusters.unwrap_or(usize::MAX) > 0,
            "Number of clusters (--n-clusters) must be positive");

        if !self.alignments.contains("{}") {
            // Make path to alignments INPUT/loci/{}/ALIGNMENTS
            let mut new_alignments = self.input.clone().unwrap();
            new_alignments.push(paths::LOCI_DIR);
            new_alignments.push("{}");
            new_alignments.push(&self.alignments);
            self.alignments = new_alignments.to_str().ok_or_else(|| error!(InvalidInput,
                "Input database has invalid UTF-8 name `{}`", self.input.as_ref().unwrap().display()))?
                .to_owned();
        }
        Ok(self)
    }
}

fn print_help() {
    const KEY: usize = 17;
    const VAL: usize = 4;
    const EMPTY: &'static str = const_format::str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Remove similar target haplotypes.".yellow());

    print!("\n{}", "Usage:".bold());
    println!(" {} -i db -o pruned_db [args]", super::PROGRAM);

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Input database directory.",
        "-i, --input".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Output pruned database directory.",
        "-o, --output".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  Path to alignment .paf[.gz] files [{}].\n\
        {EMPTY}  Should either contain {{}}, which are then replaced with locus names,\n\
        {EMPTY}  or direct to files located in {}/loci/<locus>/{}.\n\
        {EMPTY}  Alignments can be constructed using {}.",
        "-a, --alignments".green(), "PATH".yellow(), super::fmt_def(&defaults.alignments),
        "INPUT".yellow(), "PATH".yellow(), const_format::concatcp!(super::PROGRAM, " align").underline());

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  PAF field with divergence values [{}].",
        "    --field".green(), "STR".yellow(), super::fmt_def(&defaults.field));
    println!("    {:KEY$} {:VAL$}  Divergence threshold for pruning [{}].",
        "    --threshold".green(), "NUM".yellow(), super::fmt_def_f64(defaults.threshold));
    println!("    {:KEY$} {:VAL$}  Use dynamic threshold to get approximately {} clusters.",
        "    --n-clusters".green(), "INT".yellow(), "INT".yellow());
    println!("    {:KEY$} {:VAL$}  Select cluster representative with the smallest\n\
        {EMPTY}  generalized mean of this power [{}].",
        "    --power".green(), "NUM".yellow(), super::fmt_def(&defaults.power));
    println!("    {:KEY$} {:VAL$}  Limit the pruning to loci from this list.",
        "    --subset-loci".green(), "STR+".yellow());
    println!("    {:KEY$} {:VAL$}  Do not write trees in the output directory.",
        "    --skip-tree".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());

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
            Short('i') | Long("input") => args.input = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),
            Short('a') | Long("aln") | Long("alignments") => args.alignments = parser.value()?.parse()?,

            Long("subset-loci") | Long("loci-subset") => {
                for val in parser.values()? {
                    args.subset_loci.insert(val.parse()?);
                }
            }
            Long("skip-tree") => args.skip_tree = true,
            Long("field") => args.field = parser.value()?.parse()?,
            Long("thresh") | Long("threshold") => args.threshold = parser.value()?.parse()?,
            Long("n-clusters") => args.n_clusters = Some(parser.value()?.parse()?),
            Long("power") => args.power = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,

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

fn load_divergences(
    f: impl BufRead,
    contigs: &ContigNames,
    field: &str,
    repl_missing: f64,
) -> crate::Result<TriangleMatrix<f64>>
{
    let mut n_warnings = 0;
    let prefix = format!("{}:", field);
    let crop_length = prefix.as_bytes().len() + 2;
    let mut divergences = TriangleMatrix::new(contigs.len(), f64::NAN);

    for line in f.lines() {
        let line = line.map_err(add_path!(!))?;
        let split: Vec<_> = line.trim_end().split('\t').collect();
        let Some(id1) = contigs.try_get_id(&split[0]) else { continue };
        let Some(id2) = contigs.try_get_id(&split[5]) else { continue };
        if id1 == id2 { continue }

        let mut opt_val: Option<f64> = None;
        for val_str in &split[12..] {
            if val_str.starts_with(&prefix) {
                opt_val = Some(str::from_utf8(&val_str.as_bytes()[crop_length..]).ok()
                    .and_then(|s| str::parse::<f64>(s).ok())
                    .ok_or_else(|| error!(ParsingError, "Cannot parse divergence `{}`", val_str))?);
                break;
            }
        }
        // Warn later.
        let Some(val) = opt_val else { continue };
        if val < 0.0 {
            if n_warnings < 10 {
                n_warnings += 1;
                log::warn!("[{}] Cannot use negative divergence values ({} between {} and {})",
                    contigs.tag(), val, split[0], split[5]);
            }
            continue;
        }

        let d: &mut f64 = &mut divergences.get_mut_symmetric(id1.ix(), id2.ix()).expect("id1 != id2");
        if !d.is_nan() && *d != val {
            if n_warnings < 10 {
                n_warnings += 1;
                log::warn!("[{}] Multiple divergence values ({}, {}) for haplotypes {} and {}",
                    contigs.tag(), *d, val, split[0], split[5]);
            }
            continue;
        }
        *d = val;
    }

    let mut missing_count = 0;
    let mut missing_k: Option<usize> = None;
    for (k, val) in divergences.iter_mut().enumerate() {
        if val.is_nan() {
            missing_count += 1;
            missing_k = Some(k);
            *val = repl_missing;
        }
    }
    if missing_count == divergences.linear_len() {
        return Err(error!(InvalidInput, "Divergence missing for all haplotype pairs"));
    }
    if let Some(k) = missing_k {
        let (i, j) = divergences.from_linear_index(k);
        let name1 = contigs.get_name(ContigId::new(i as u16));
        let name2 = contigs.get_name(ContigId::new(j as u16));
        log::warn!("[{}] Divergence missing for {:.1}% haplotypes pairs, for example {} and {}. Replacing with {:.5}",
            contigs.tag(), 100.0 * missing_count as f64 / divergences.linear_len() as f64, name1, name2, repl_missing);
    }
    Ok(divergences)
}

/// Loads `discarded_haplotypes.txt` file, with lines
/// haplotype = haplotype2, haplotype3, ...
fn load_discarded_haplotypes(
    f: impl BufRead,
    contigs: &ContigNames,
) -> io::Result<Vec<(ContigId, Vec<String>)>>
{
    let mut corresp = Vec::new();
    let mut eq_warned = false;
    for line in f.lines() {
        let line = line?;
        let split: Vec<_> = line.split_whitespace().collect();
        let Some(id) = contigs.try_get_id(split[0]) else {
            log::warn!("[{}] Unknown haplotype `{}` among discarded haplotypes",
                contigs.tag(), split[0]);
            continue;
        };
        if split[1] != "=" && !eq_warned {
            eq_warned = true;
            log::warn!("[{}] Input haplotypes were previously pruned. Output newick tree will be incomplete",
                contigs.tag());
        }
        let mut curr_haps = Vec::with_capacity(split.len() - 2);
        for contig in &split[2..] {
            let contig = contig.strip_suffix(',').unwrap_or(contig);
            if contigs.contains(contig) {
                log::warn!("[{}] Haplotype {} is marked as discarded, but present in the haplotypes fasta",
                    contigs.tag(), contig);
                continue;
            }
            curr_haps.push(contig.to_string());
        }
        corresp.push((id, curr_haps));
    }
    Ok(corresp)
}

struct Cluster {
    /// Haplotypes in the cluster and their multiplicity (number of times identical haplotype appeared,
    /// taken from old discarded_haplotypes.txt.gz).
    haps: Vec<(ContigId, u32)>,
    /// Newick representation.
    newick: String,
    /// Divergence from leaf nodes.
    div: f64,
}

impl Cluster {
    fn new(id: ContigId, contigs: &ContigNames) -> Self {
        Self {
            haps: vec![(id, 1)],
            newick: contigs.get_name(id).to_owned(),
            div: 0.0,
        }
    }

    fn add_identical(&mut self, names: &[String]) {
        debug_assert!(self.haps.len() == 1);
        self.haps[0].1 += names.len() as u32;
        self.newick.insert(0, '(');
        write!(self.newick, ":0").unwrap();
        for hap in names {
            write!(self.newick, ",{}:0", hap).unwrap();
        }
        self.newick.push(')');
    }

    /// Merges two clusters and clear contents of the original clusters (for memory efficiency).
    fn merge_and_clear(first: &mut Self, second: &mut Self, div: f64) -> Self {
        // Take old vectors to limit the peak memory usage.
        let mut haps = Vec::with_capacity(first.haps.len() + second.haps.len());
        haps.extend(std::mem::take(&mut first.haps));
        haps.extend(std::mem::take(&mut second.haps));
        let newick = format!("({}:{:.8},{}:{:.8})",
            &first.newick, 0.5 * (div - first.div),
            &second.newick, 0.5 * (div - second.div));
        first.newick.clear();
        second.newick.clear();
        Cluster { haps, newick, div }
    }

    fn select_representative(
        &self,
        divergences: &TriangleMatrix<f64>,
        epsilon: f64,
        power: PowerMean,
        buf: &mut Vec<f64>,
    ) -> ContigId
    {
        buf.clear();
        buf.resize(self.haps.len(), 0.0);
        for (i, &(id1, mult1)) in self.haps.iter().enumerate() {
            if mult1 > 1 {
                buf[i] = power.update_mult(buf[i], epsilon, f64::from(mult1 - 1));
            }
            for (j, &(id2, mult2)) in self.haps[i + 1..].iter().enumerate() {
                let j = i + j + 1;
                let div = epsilon + *divergences.get_symmetric(id1.ix(), id2.ix()).unwrap();
                buf[i] = power.update_mult(buf[i], div, f64::from(mult2));
                buf[j] = power.update_mult(buf[j], div, f64::from(mult1));
            }
        }
        // NOTE: There is no need to finalize the generalized mean,
        // accumulators will be strictly increasing with non-negative power and strictly decreasing otherwise.
        let (k, _) = match power {
            PowerMean::Pow(n) if n < 0 => ext::vec::F64Ext::argmax(&buf),
            _ => ext::vec::F64Ext::argmin(&buf),
        };
        self.haps[k].0
    }

    fn write_discarded_haplotypes(
        &self,
        f: &mut impl Write,
        contigs: &ContigNames,
        repres: ContigId,
    ) -> io::Result<()> {
        assert!(self.haps.len() > 1);
        write!(f, "{} ", contigs.get_name(repres))?;
        let mut separator = '~';
        for &(hap, _) in &self.haps {
            if hap != repres {
                write!(f, "{} {}", separator, contigs.get_name(hap))?;
                separator = ',';
            }
        }
        writeln!(f)?;
        Ok(())
    }
}

/// Select cut threshold.
fn select_cut_threshold(tag: &str, dendrogram: &kodama::Dendrogram<f64>, n_clusters: usize) -> f64 {
    let n_haps = dendrogram.observations();
    // Each step decreases the number of clusters by 1.
    let Some(i) = n_haps.checked_sub(n_clusters + 1) else {
        log::warn!("[{}] There are {} haplotypes, smaller or equal \
            than the required number of pruned haplotypes ({})", tag, n_haps, n_clusters);
        return 0.0;
    };

    let steps = dendrogram.steps();
    // Steps should be always sorted, but since they come from an external library, we cannot guarantee that.
    let thresh = if steps.windows(2).all(|w| w[0].dissimilarity <= w[1].dissimilarity) {
        steps[i].dissimilarity
    } else {
        let mut steps_sorted = steps.to_vec();
        steps_sorted.sort_unstable_by(|a, b| a.dissimilarity.total_cmp(&b.dissimilarity));
        steps_sorted[i].dissimilarity
    };
    log::info!("[{}] Use cut threshold of {:.6}", tag, thresh);
    thresh
}

/// Cluster haplotypes, cut clusters, write new discarded haplotypes & the newick file.
fn cluster_haplotypes(
    contigs: &ContigNames,
    divergences: TriangleMatrix<f64>,
    mut thresh: f64,
    n_clusters: Option<usize>,
    power: PowerMean,
    disc_haps: &[(ContigId, Vec<String>)],
    mut nwk_writer: impl Write,
    mut disc_haps_writer: impl Write,
) -> crate::Result<Vec<ContigId>>
{
    let min_val = divergences.iter().copied().reduce(f64::min).expect("Empty divergence matrix");
    if min_val > thresh {
        log::warn!("[{}] Minimal observed divergence ({}) is larger than the threshold ({}), \
            all haplotypes will be retained", contigs.tag(), min_val, thresh);
    }
    // Value which will replace 0 in generalized mean.
    let epsilon = f64::max(1e-6 * min_val, 1e-12);

    let n_haps = contigs.len();
    let total_clusters = 2 * n_haps - 1;
    // Use Complete method, meaning that we track maximal distance between two points between two clusters.
    // This is done to cut as little as needed.
    let mut div_vec = divergences.linear_data().to_vec();
    let dendrogram = kodama::linkage(&mut div_vec, n_haps, kodama::Method::Complete);
    if let Some(n) = n_clusters {
        thresh = select_cut_threshold(contigs.tag(), &dendrogram, n);
    }

    let mut clusters: Vec<_> = contigs.ids().map(|id| Cluster::new(id, contigs)).collect();
    for (id, haps) in disc_haps {
        clusters[id.ix()].add_identical(&haps);
    }

    let steps = dendrogram.steps();
    let mut buf = Vec::new();
    let mut keep_ids = Vec::new();
    for step in steps {
        let mut curr_clusters = clusters.get_disjoint_mut([step.cluster1, step.cluster2])
            .expect("Two merged clusters must differ");
        if step.dissimilarity > thresh {
            // If this is the first time we exceed divergence threshold,
            // write down discarded haplotypes and clear haplotypes.
            for cluster in curr_clusters.iter_mut() {
                match cluster.haps.len() {
                    0 => continue, // Divergence exceeded the threshold on one the previous iterations.
                    1 => keep_ids.push(cluster.haps[0].0),
                    _ => {
                        let repres = cluster.select_representative(&divergences, epsilon, power, &mut buf);
                        keep_ids.push(repres);
                        cluster.write_discarded_haplotypes(&mut disc_haps_writer, contigs, repres)
                            .map_err(add_path!(!))?;
                    }
                }
                cluster.haps.clear();
            }
        }
        let [cluster1, cluster2] = curr_clusters;
        let new_cluster = Cluster::merge_and_clear(cluster1, cluster2, step.dissimilarity);
        clusters.push(new_cluster);
    }

    assert_eq!(clusters.len(), total_clusters);
    writeln!(nwk_writer, "{};", clusters.last().unwrap().newick).map_err(add_path!(!))?;
    keep_ids.sort_unstable();
    Ok(keep_ids)
}

/// In case where no haplotypes are discarded, simply copy all files from input to output directory.
fn copy_output_files(locus_data: &LocusData) -> crate::Result<bool> {
    let mut all_present = true;
    for &(basename, must_have) in &[
        (paths::LOCUS_FASTA, true),
        (paths::KMERS, true),
        (paths::DISTANCES, false),
        ]
    {
        let in_filename = locus_data.db_dir().join(basename);
        let out_filename = locus_data.out_dir().join(basename);
        if Path::new(&in_filename).exists() {
            std::fs::copy(&in_filename, &out_filename).map_err(add_path!(in_filename, out_filename))?;
        } else if must_have {
            log::error!("[{}] File {} is missing, output will be incomplete",
                locus_data.contig_set().tag(), in_filename.display());
            all_present = false;
        }
    }
    Ok(all_present)
}

fn copy_bed_files(locus_data: &LocusData) -> crate::Result<bool> {
    let ref_bed_filename = OsStr::new("ref.bed");
    let mut has_ref_bed = false;
    for filename in ext::sys::filenames_with_ext(locus_data.db_dir(), "bed")? {
        let basename = filename.file_name().expect("File with BED extension must have basename");
        has_ref_bed |= basename == ref_bed_filename;
        let out_filename = locus_data.out_dir().join(basename);
        std::fs::copy(&filename, &out_filename).map_err(add_path!(filename, out_filename))?;
    }
    if !has_ref_bed {
        log::error!("[{}] File {} is missing, output will be incomplete",
            locus_data.contig_set().tag(), locus_data.db_dir().join(&ref_bed_filename).display());
    }
    Ok(has_ref_bed)
}

/// Thins out files from the input locus directory.
/// Returns true if all necessary input files are present.
fn thin_out_output_files(locus_data: &LocusData, keep_ids: &[ContigId]) -> crate::Result<bool> {
    let mut all_files_present = copy_bed_files(locus_data)?;
    let contig_set = locus_data.contig_set();
    let contigs = contig_set.contigs();
    if keep_ids.len() == contigs.len() {
        all_files_present &= copy_output_files(locus_data)?;
        return Ok(all_files_present);
    }

    let fasta_filename = locus_data.out_dir().join(paths::LOCUS_FASTA);
    let mut fasta_writer = ext::sys::create_gzip(&fasta_filename)?;
    for &id in keep_ids {
        fastx::write_fasta(&mut fasta_writer, contigs.get_name(id).as_bytes(), contig_set.get_seq(id))
            .map_err(add_path!(fasta_filename))?;
    }

    // NOTE, that kmers are also read in ContigSet::load, but this is difficult to remove.
    let in_kmers_filename = locus_data.db_dir().join(paths::KMERS);
    let out_kmers_filename = locus_data.out_dir().join(paths::KMERS);
    let mut kmers_reader = ext::sys::open(&in_kmers_filename)?;
    // Read twice because there are k-mer counts in the KMERS file (off-target & on-target).
    // Currently, only off-target counts are used, but we still need to thin out both.
    let kmer_counts1 = KmerCounts::load(&mut kmers_reader).map_err(add_path!(in_kmers_filename))?;
    let kmer_counts2 = KmerCounts::load(&mut kmers_reader).map_err(add_path!(in_kmers_filename))?;
    if kmer_counts1.validate(contigs).is_ok() && kmer_counts2.validate(contigs).is_ok() {
        let mut kmers_writer = ext::sys::create_lz4_slow(&out_kmers_filename)?;
        for counts in [kmer_counts1, kmer_counts2] {
            counts.thin_out(keep_ids.iter().copied().map(ContigId::ix))
                .save(&mut kmers_writer).map_err(add_path!(out_kmers_filename))?;
        }
    } else {
        log::warn!("[{}] k-mer counts in {} do not match input haplotypes, output will be incomplete",
            contigs.tag(), in_kmers_filename.display());
        all_files_present = false;
    }

    let dist_filename = locus_data.db_dir().join(paths::DISTANCES);
    if dist_filename.exists() {
        let dist_file = ext::sys::open_uncompressed(&dist_filename)?;
        let (k, w, dists) = div::load_divergences(dist_file, &dist_filename, contigs.len())?;
        let subdists = dists.thin_out(keep_ids);
        let new_dist_filename = locus_data.out_dir().join(paths::DISTANCES);
        let new_dist_file = ext::sys::create(&new_dist_filename)?;
        div::write_divergences(new_dist_file, k, w, &subdists, |&d| d).map_err(add_path!(new_dist_filename))?;
    }
    Ok(all_files_present)
}

/// Returns Err if irreversible error occured,
///     Ok(false) if process finished, but output is incomplete, and Ok(true) if everything finished successfully.
fn process_locus(
    locus_data: &LocusData,
    args: &Args,
) -> crate::Result<bool> {
    let contig_set = locus_data.contig_set();
    let contigs = contig_set.contigs();
    if contigs.is_empty() {
        return Err(error!(InvalidData, "No haplotypes found"));
    }
    let paf_filename = args.alignments.replace("{}", contig_set.tag());
    let repl_missing = if args.n_clusters.is_some() { f64::INFINITY } else { 10.0 * args.threshold };
    let divergences = load_divergences(ext::sys::open(&paf_filename)?, contigs, &args.field, repl_missing)?;

    let disc_filename = PathBuf::from(locus_data.db_dir().join(paths::DISCARDED_HAPS));
    let mut disc_haps_data = Vec::new();
    let disc_haplotypes = if args.skip_tree || !disc_filename.exists() {
        Default::default()
    } else {
        ext::sys::open(&disc_filename)?.read_to_end(&mut disc_haps_data).map_err(add_path!(disc_filename))?;
        load_discarded_haplotypes(&disc_haps_data as &[u8], contigs).map_err(add_path!(!))?
    };

    let nwk_writer: Box<dyn Write> = if args.skip_tree {
        Box::new(io::sink())
    } else {
        Box::new(ext::sys::create_gzip(&locus_data.out_dir().join("all_haplotypes.nwk.gz"))?)
    };
    let keep_ids = cluster_haplotypes(contigs, divergences, args.threshold, args.n_clusters, args.power,
        &disc_haplotypes, nwk_writer, &mut disc_haps_data)?;
    if !disc_haps_data.is_empty() {
        let out_disc_haps_filename = locus_data.out_dir().join(paths::DISCARDED_HAPS);
        ext::sys::create(&out_disc_haps_filename)?
            .write_all(&disc_haps_data).map_err(add_path!(out_disc_haps_filename))?;
    }

    if keep_ids.len() < contigs.len() {
        log::info!("[{}] Retained {} / {} haplotypes", contigs.tag(), keep_ids.len(), contigs.len());
    } else {
        log::info!("[{}] Retained all {} haplotypes", contigs.tag(), contigs.len());
    }

    if thin_out_output_files(locus_data, &keep_ids)? {
        super::write_success_file(locus_data.out_dir().join(paths::SUCCESS))?;
        Ok(true)
    } else {
        Ok(false)
    }
}

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();

    let input = args.input.as_ref().expect("Input path must be defined");
    let output = args.output.as_ref().expect("Output path must be defined");
    ext::sys::mkdir(&output)?;
    let loci_dir = output.join(paths::LOCI_DIR);
    ext::sys::mkdir(&loci_dir)?;

    let loci = super::genotype::load_loci(&[input], &output, &args.subset_loci,
        super::Rerun::from_force(args.force))?;
    if loci.is_empty() {
        return Ok(());
    }

    let total = loci.len();
    let mut failed = 0;
    let mut incomplete = 0;
    for locus in loci {
        match process_locus(&locus, &args) {
            Ok(true) => {}
            Ok(false) => incomplete += 1,
            Err(e) => {
                log::error!("Error while analyzing locus {}:\n        {}", locus.contig_set().tag(), e.display());
                failed += 1;
            }
        }
    }

    let succeed = total - failed;
    if succeed == 0 {
        log::error!("Failed to prune all {} loci", failed);
    } else if failed > 0 {
        log::warn!("Successfully pruned {} loci, failed to prune {} loci", succeed, failed);
    } else {
        log::info!("Successfully pruned {} loci", succeed);
    }
    if incomplete > 0 {
        log::warn!("Of the pruned loci, {} output directories were not completely filled", incomplete);
    }
    log::info!("Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
