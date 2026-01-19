use std::{
    path::{Path, PathBuf},
    time::Instant,
    io::{self, Write, BufRead},
    ffi::OsStr,
};
use colored::Colorize;
use crate::{
    ext::{self, TriangleMatrix},
    algo::{HashSet, IntMap},
    seq::{ContigId, ContigNames, ContigSet},
    err::{self, validate_param, error, add_path},
};
use super::genotype::LocusData;

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
            div_thresh: 0.005,
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

        if !self.alignments.contains("{}") {
            // Make path to alignments INPUT/loci/{}/ALIGNMENTS
            let mut new_alignments = self.input.clone().unwrap();
            new_alignments.push(super::paths::LOCI_DIR);
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
        if d.is_nan() && *d != val {
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
    mut f: impl BufRead,
    contigs: &ContigNames,
) -> io::Result<IntMap<ContigId, Vec<String>>>
{
    let mut corresp = IntMap::default();
    let mut eq_warned = false;
    for line in f.lines() {
        let line = line?;
        let split: Vec<_> = line.split_whitespace().collect();
        let Some(id) = contigs.try_get_id(split[0]) else { continue };
        if split[1] != "=" && !eq_warned {
            eq_warned = true;
            log::warn!("[{}] In discarded_haplotypes.txt, some lines contains `{}` separator. \
                This will not be reflected in the newick tree", contigs.tag(), split[1]);
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
        corresp.insert(id, curr_haps);
    }
    Ok(corresp)
}

// fn cluster_haplotypes(
//     mut nwk_writer: impl Write,
//     entries: &[NamedSeq],
//     mut divergences: Vec<f64>,
//     thresh: f64,
// ) -> io::Result<Vec<bool>> {
//     let n = entries.len();
//     let total_clusters = 2 * n - 1;
//     // Use Complete method, meaning that we track maximal distance between two points between two clusters.
//     // This is done to cut as little as needed.
//     let dendrogram = kodama::linkage(&mut divergences, n, kodama::Method::Complete);
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
//             dist = 0.5 * step.dissimilarity));
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
//             if step.dissimilarity <= thresh {
//                 keep_seqs[cluster_repr[i]] = true;
//             } else {
//                 queue.push(step.cluster1);
//                 queue.push(step.cluster2);
//             }
//         }
//     }
//     Ok(keep_seqs)
// }

fn process_locus(
    locus_data: &LocusData,
    args: &Args,
) -> crate::Result<()> {
    let contig_set = locus_data.contig_set();
    let contigs = contig_set.contigs();
    let paf_filename = args.alignments.replace("{}", contig_set.tag());
    let divergences = load_divergences(ext::sys::open(&paf_filename)?, contigs, &args.div_field)?;

    let disc_filename = PathBuf::from(locus_data.out_dir().join(super::paths::DISCARDED_HAPS));
    let disc_haplotypes = if args.skip_tree || !disc_filename.exists() {
        Default::default()
    } else {
        load_discarded_haplotypes(ext::sys::open(&disc_filename)?, contigs).map_err(add_path!(disc_filename))?
    };

    Ok(())
}

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let args = parse_args(argv)?.validate()?;
    super::greet();
    let timer = Instant::now();

    let input = args.input.as_ref().expect("Input path must be defined");
    let output = args.output.as_ref().expect("Output path must be defined");
    let loci = super::genotype::load_loci(&[input], &output, &args.subset_loci,
        super::Rerun::from_force(args.force))?;

    log::info!("Success! Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
