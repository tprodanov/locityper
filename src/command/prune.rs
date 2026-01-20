use std::{
    path::PathBuf,
    time::Instant,
    io::{self, Write, BufRead},
    fmt::Write as FmtWrite,
};
use colored::Colorize;
use crate::{
    ext::{self, TriangleMatrix},
    algo::HashSet,
    seq::{fastx, ContigId, ContigNames},
    err::{validate_param, error, add_path},
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
    div_field: String,
    div_thresh: f64,
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
            div_field: "dv".to_string(),
            div_thresh: 0.0002,
            force: false,
        }
    }
}

impl Args {
    fn validate(mut self) -> crate::Result<Self> {
        validate_param!(self.input.is_some(), "Input database is not provided (see -i/--input)");
        validate_param!(self.output.is_some(), "Output database is not provided (see -o/--output)");
        validate_param!(!self.div_field.contains(':'), "PAF divergence field ({}) must not contain :",
            self.div_field);
        validate_param!(self.div_thresh >= 0.0, "Divergence threshold ({}) should be non-negative",
            self.div_thresh);

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

    println!("\n{}", "Input arguments:".bold());
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
        "    --div-field".green(), "STR".yellow(), super::fmt_def(&defaults.div_field));
    println!("    {:KEY$} {:VAL$}  Divergence threshold for pruning [{}].",
        "    --div-thresh".green(), "NUM".yellow(), super::fmt_def_f64(defaults.div_thresh));
    println!("    {:KEY$} {:VAL$}  Limit the pruning to loci from this file.",
        "    --subset-loci".green(), "FILE".yellow());
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

            Long("subset-loci") | Long("loci-subset") => {
                for val in parser.values()? {
                    args.subset_loci.insert(val.parse()?);
                }
            }
            Long("skip-tree") => args.skip_tree = true,
            Long("div-field") => args.div_field = parser.value()?.parse()?,
            Long("div-thresh") => args.div_thresh = parser.value()?.parse()?,
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

        let i = usize::min(id1.ix(), id2.ix());
        let j = usize::max(id1.ix(), id2.ix());
        let d = &mut divergences[(i, j)];
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
    let n_nans = divergences.iter().filter(|&d| d.is_nan()).count();
    if n_nans == divergences.linear_len() {
        return Err(error!(InvalidInput, "[{}] Divergence missing for all haplotype pairs", contigs.tag()));
    } else if n_nans > 0 {
        let k = divergences.iter().enumerate().filter(|(_i, &d)| d.is_nan()).next().unwrap().0;
        let (i, j) = divergences.from_linear_index(k);
        log::warn!("[{}] Divergence missing for {}/{} haplotype pairs, for example {} and {}",
            contigs.tag(), n_nans, divergences.linear_len(),
            contigs.get_name(ContigId::new(i as u16)), contigs.get_name(ContigId::new(j as u16)));
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
    /// Haplotypes, where the first one is the representative.
    haps: Vec<ContigId>,
    /// How many haplotypes are in the cluster, including previously discarded haplotypes.
    actual_size: u32,
    /// Newick representation.
    newick: String,
}

impl Cluster {
    fn new(id: ContigId, contigs: &ContigNames) -> Self {
        Self {
            haps: vec![id],
            actual_size: 1,
            newick: contigs.get_name(id).to_owned(),
        }
    }

    fn add_identical(&mut self, names: &[String]) {
        self.actual_size += names.len() as u32;
        self.newick.insert(0, '(');
        write!(self.newick, ":0").unwrap();
        for hap in names {
            write!(self.newick, ",{}:0", hap).unwrap();
        }
        self.newick.push(')');
    }

    /// Merges two clusters and clear contents of the original clusters (for memory efficiency).
    fn merge_and_clear(first: &mut Self, second: &mut Self, distance: f64) -> Self {
        // Larger cluster will be starting.
        if first.actual_size < second.actual_size {
            std::mem::swap(&mut first.haps, &mut second.haps);
        }
        let mut haps = std::mem::take(&mut first.haps);
        haps.extend(std::mem::take(&mut second.haps));
        let newick = format!("({}:{dist:.8},{}:{dist:.8})", &first.newick, &second.newick, dist = 0.5 * distance);
        first.newick.clear();
        second.newick.clear();
        Cluster {
            haps, newick,
            actual_size: first.actual_size + second.actual_size,
        }
    }

    fn write_discarded_haplotypes(&self, f: &mut impl Write, contigs: &ContigNames) -> io::Result<()> {
        if self.haps.len() > 1 {
            write!(f, "{} ~ {}", contigs.get_name(self.haps[0]), contigs.get_name(self.haps[1]))?;
            for &hap in &self.haps[2..] {
                write!(f, ", {}", contigs.get_name(hap))?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

/// Cluster haplotypes, cut clusters, write new discarded haplotypes & the newick file.
fn cluster_haplotypes(
    contigs: &ContigNames,
    divergences: TriangleMatrix<f64>,
    thresh: f64,
    disc_haps: &[(ContigId, Vec<String>)],
    mut nwk_writer: impl Write,
    mut disc_haps_writer: impl Write,
) -> io::Result<Vec<ContigId>> {
    let n = contigs.len();
    let total_clusters = 2 * n - 1;
    // Use Complete method, meaning that we track maximal distance between two points between two clusters.
    // This is done to cut as little as needed.
    let mut divergences = divergences.take_linear();
    let dendrogram = kodama::linkage(&mut divergences, n, kodama::Method::Complete);
    std::mem::drop(divergences);

    let mut clusters: Vec<_> = contigs.ids().map(|id| Cluster::new(id, contigs)).collect();
    for (id, haps) in disc_haps {
        clusters[id.ix()].add_identical(&haps);
    }
    let steps = dendrogram.steps();
    let mut keep_ids = Vec::new();
    for step in steps.iter() {
        if step.dissimilarity > thresh {
            // If this is the first time we exceed divergence threshold,
            // write down discarded haplotypes and clear haplotypes.
            for i in [step.cluster1, step.cluster2] {
                let cluster = &mut clusters[i];
                if cluster.haps.is_empty() {
                    continue;
                }
                keep_ids.push(cluster.haps[0]);
                cluster.write_discarded_haplotypes(&mut disc_haps_writer, contigs)?;
                cluster.haps.clear();
            }
        }
        let [cluster1, cluster2] = clusters.get_disjoint_mut([step.cluster1, step.cluster2])
            .expect("Two merged clusters must differ");
        let new_cluster = Cluster::merge_and_clear(cluster1, cluster2, step.dissimilarity);
        clusters.push(new_cluster);
    }
    assert_eq!(clusters.len(), total_clusters);
    writeln!(nwk_writer, "{};", clusters.last().unwrap().newick)?;
    Ok(keep_ids)
}

fn process_locus(
    locus_data: &LocusData,
    args: &Args,
) -> crate::Result<()> {
    let contig_set = locus_data.contig_set();
    let contigs = contig_set.contigs();
    let paf_filename = args.alignments.replace("{}", contig_set.tag());
    let divergences = load_divergences(ext::sys::open(&paf_filename)?, contigs, &args.div_field)?;

    let disc_filename = PathBuf::from(locus_data.db_locus_dir().join(paths::DISCARDED_HAPS));
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
    let keep_ids = cluster_haplotypes(contigs, divergences, args.div_thresh, &disc_haplotypes,
        nwk_writer, &mut disc_haps_data).map_err(add_path!(!))?;
    if !disc_haps_data.is_empty() {
        let out_disc_haps_filename = locus_data.out_dir().join(paths::DISCARDED_HAPS);
        ext::sys::create(&out_disc_haps_filename)?
            .write_all(&disc_haps_data).map_err(add_path!(out_disc_haps_filename))?;
    }

    let fasta_filename = locus_data.out_dir().join(paths::LOCUS_FASTA);
    let mut fasta_writer = ext::sys::create_gzip(&fasta_filename)?;
    for &id in &keep_ids {
        fastx::write_fasta(&mut fasta_writer, contigs.get_name(id).as_bytes(), contig_set.get_seq(id))
            .map_err(add_path!(fasta_filename))?;
    }
    if keep_ids.len() < contigs.len() {
        log::info!("[{}] Retained {} / {} haplotypes", contigs.tag(), keep_ids.len(), contigs.len());
    } else {
        log::info!("[{}] Retained all {} haplotypes", contigs.tag(), contigs.len());
    }
    Ok(())
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

    let total = loci.len();
    let mut failed = 0;
    for locus in loci {
        if let Err(e) = process_locus(&locus, &args) {
            log::error!("Error while analyzing locus {}:\n        {}", locus.contig_set().tag(), e.display());
            failed += 1;
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
    log::info!("Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
