use std::{
    fs,
    io::{self, Write},
    process::{Command, Stdio},
    cmp::{min, max},
    path::{Path, PathBuf},
    time::Instant,
    sync::Arc,
};
use colored::Colorize;
use const_format::{str_repeat, concatcp};
use fnv::FnvHashSet;
use crate::{
    err::{Error, validate_param},
    math::Ln,
    seq::{
        recruit, fastx,
        ContigSet, NamedSeq,
        kmers::Kmer,
    },
    bg::{BgDistr, JsonSer},
    ext,
    model::{
        Params as AssgnParams,
        locs::{self, PrelimAlignments},
        windows::ContigWindows,
        dp_cache::CachedDepthDistrs,
    },
    solvers::scheme,
};
use htslib::bam;
use super::paths;

struct Args {
    input: Vec<PathBuf>,
    database: Option<PathBuf>,
    output: Option<PathBuf>,
    subset_loci: FnvHashSet<String>,
    ploidy: u8,
    solvers: Option<PathBuf>,

    interleaved: bool,
    threads: u16,
    force: bool,
    strobealign: PathBuf,
    samtools: PathBuf,
    seed: Option<u64>,
    debug: bool,

    recr_params: recruit::Params,
    assgn_params: AssgnParams,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: Vec::new(),
            database: None,
            output: None,
            subset_loci: FnvHashSet::default(),
            ploidy: 2,
            solvers: None,

            interleaved: false,
            threads: 8,
            force: false,
            strobealign: PathBuf::from("strobealign"),
            samtools: PathBuf::from("samtools"),
            seed: None,
            debug: false,

            recr_params: Default::default(),
            assgn_params: Default::default(),
        }
    }
}

impl Args {
    /// Validate arguments, modifying some, if needed.
    fn validate(mut self) -> Result<Self, Error> {
        self.threads = max(self.threads, 1);
        let n_input = self.input.len();
        validate_param!(n_input > 0, "Read files are not provided (see -i/--input)");
        validate_param!(n_input != 2 || !self.interleaved,
            "Two read files (-i/--input) are provided, however, --interleaved is specified");
        if !self.interleaved && n_input == 1 {
            log::warn!("Running in single-end mode.");
        }
        validate_param!(self.ploidy > 0 && self.ploidy <= 11, "Ploidy ({}) must be within [1, 10]", self.ploidy);

        validate_param!(self.database.is_some(), "Database directory is not provided (see -d/--database)");
        validate_param!(self.output.is_some(), "Output directory is not provided (see -o/--output)");

        self.strobealign = ext::sys::find_exe(self.strobealign)?;
        self.samtools = ext::sys::find_exe(self.samtools)?;

        self.recr_params.validate()?;
        self.assgn_params.validate()?;
        Ok(self)
    }

    fn is_paired_end(&self) -> bool {
        self.input.len() == 2 || self.interleaved
    }
}

fn print_help() {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Analyze WGS dataset.".yellow());

    println!("\n{} {} analyze -i reads1.fq [reads2.fq] -d db -o out [arguments]",
        "Usage:".bold(), super::PKG_NAME);

    println!("\n{}", "Input/output arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Database directory (initialized with {} & {}).",
        "-d, --db".green(), "DIR".yellow(), concatcp!(super::PKG_NAME, " create").underline(), "add".underline());
    println!("    {:KEY$} {:VAL$}  Output directory   (initialized with {}).",
        "-o, --output".green(), "DIR".yellow(), concatcp!(super::PKG_NAME, " preproc").underline());
    println!("    {:KEY$} {:VAL$}  Input reads are interleaved.",
        "-^, --interleaved".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Optional: only analyze loci with names from this list.",
        "    --subset-loci".green(), "STR+".yellow());

    println!("\n{}", "Read recruitment:".bold());
    println!("    {:KEY$} {:VAL$}  Minimizer k-mer size (no larger than {}) [{}].",
        "-k, --recr-kmer".green(), "INT".yellow(), recruit::Minimizer::MAX_KMER_SIZE, defaults.recr_params.minimizer_k);
    println!("    {:KEY$} {:VAL$}  Take k-mers with smallest hash across {} consecutive k-mers [{}].",
        "-w, --recr-window".green(), "INT".yellow(), "INT".yellow(), defaults.recr_params.minimizer_w);
    println!("    {:KEY$} {:VAL$}  Recruit single-end reads or read pairs with at least this fraction\n\
        {EMPTY}  of minimizers matching one of the targets [{:.1}].",
        "-m, --matches-frac".green(), "FLOAT".yellow(), defaults.recr_params.matches_frac);
    println!("    {:KEY$} {:VAL$}  Recruit reads in chunks of this size [{}].\n\
        {EMPTY}  May impact runtime in multi-threaded read recruitment.",
        "-c, --chunk-size".green(), "INT".yellow(), defaults.recr_params.chunk_size);

    println!("\n{}", "Locus genotyping:".bold());
    println!("    {:KEY$} {:VAL$}  Solution ploidy [{}]. May be very slow for ploidy over 2.",
        "-p, --ploidy".green(), "INT".yellow(), defaults.ploidy);
    println!("    {:KEY$} {:VAL$}  Optional: describe sequence of solvers in a JSON file.\n\
        {EMPTY}  Please see README for information on the file content.",
        "-S, --solvers".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Read depth likelihood contribution relative to\n\
        {EMPTY}  read alignment likelihoods [{:.1}].",
        "-C, --dp-contrib".green(), "FLOAT".yellow(), defaults.assgn_params.depth_contrib);
    println!("    {:KEY$} {:VAL$}  Ignore read alignments that are 10^{} times worse than\n\
        {EMPTY}  the best alignment [{:.1}].",
        "-D, --prob-diff".green(), "FLOAT".yellow(), "FLOAT".yellow(), Ln::to_log10(defaults.assgn_params.prob_diff));
    println!("    {:KEY$} {:VAL$}  Unmapped read mate receives 10^{} penalty [{:.1}].",
        "-U, --unmapped".green(), "FLOAT".yellow(), "FLOAT".yellow(),
        Ln::to_log10(defaults.assgn_params.unmapped_penalty));
    println!("    {:KEY$} {} \n\
        {EMPTY}  Contig windows receive different weight depending on average k-mer\n\
        {EMPTY}  frequency. Windows with values under {} [{}] receive full weight.\n\
        {EMPTY}  Windows with values equal to {} [{}] receive half weight.",
        "    --rare-kmer".green(), "FLOAT FLOAT".yellow(), "FLOAT_1".yellow(), defaults.assgn_params.rare_kmer,
        "FLOAT_2".yellow(), defaults.assgn_params.semicommon_kmer);
    println!("{}", "    Random tweaking:".bold());
    println!("    {:KEY$} {:VAL$}  Randomly move read coordinates by at most {} bp [{}].",
        "    --tweak".green(), "INT".yellow(), "INT".yellow(), defaults.assgn_params.tweak);
    println!("    {:KEY$} {:VAL$}  Average across this number of random tweaking\n\
        {EMPTY}  for each genotype and each solver [{}].",
        "    --tw-count".green(), "INT".yellow(), defaults.assgn_params.tweak_count);
    println!("    {:KEY$} {:VAL$}  Genotype likelihoods for random tweakings are averaged\n\
        {EMPTY}  using generalized mean with this power [{:.0}].",
        "    --tw-aver".green(), "FLOAT".yellow(), defaults.assgn_params.tweak_hoelder_p);

    println!("\n{}", "Execution parameters:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), defaults.threads);
    println!("    {:KEY$} {:VAL$}  Force rewrite output directory.",
        "-F, --force".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Random seed. Ensures reproducibility for the same\n\
        {EMPTY}  input and product version.",
        "-s, --seed".green(), "INT".yellow());
    println!("    {:KEY$} {:VAL$}  Create more files with debug information.",
        "    --debug".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Strobealign executable. [{}].",
        "    --strobealign".green(), "EXE".yellow(), defaults.strobealign.display());
    println!("    {:KEY$} {:VAL$}  Samtools executable [{}].",
        "    --samtools".green(), "EXE".yellow(), defaults.samtools.display());

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
            Short('i') | Long("input") =>
                args.input = parser.values()?.take(2).map(|s| s.parse()).collect::<Result<_, _>>()?,
            Short('d') | Long("database") => args.database = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),
            Long("subset-loci") => {
                for val in parser.values()? {
                    args.subset_loci.insert(val.parse()?);
                }
                if args.subset_loci.is_empty() {
                    return Err(lexopt::Error::MissingValue { option: Some("subset-loci".to_owned()) });
                }
            }

            Short('k') | Long("recr-kmer") => args.recr_params.minimizer_k = parser.value()?.parse()?,
            Short('w') | Long("recr-window") => args.recr_params.minimizer_w = parser.value()?.parse()?,
            Short('m') | Long("matches-frac") | Long("matches-fraction") =>
                args.recr_params.matches_frac = parser.value()?.parse()?,
            Short('c') | Long("chunk") | Long("chunk-size") => args.recr_params.chunk_size = parser.value()?.parse()?,

            Short('S') | Long("solvers") => args.solvers = Some(parser.value()?.parse()?),
            Short('C') | Long("dp-contrib") | Long("depth-contrib") | Long("dp-contribution")
                => args.assgn_params.depth_contrib = parser.value()?.parse()?,
            Short('D') | Long("prob-diff") => args.assgn_params.prob_diff = Ln::from_log10(parser.value()?.parse()?),
            Short('U') | Long("unmapped") =>
                args.assgn_params.unmapped_penalty = Ln::from_log10(parser.value()?.parse()?),
            Long("rare-kmer") => {
                let mut values = parser.values()?;
                args.assgn_params.rare_kmer = values.next()
                    .ok_or_else(|| lexopt::Error::MissingValue { option: Some("rare-kmer".to_owned()) })?
                    .parse()?;
                args.assgn_params.semicommon_kmer = values.next()
                    .ok_or_else(|| lexopt::Error::MissingValue { option: Some("rare-kmer".to_owned()) })?
                    .parse()?;
            }
            Long("tweak") => args.assgn_params.tweak = parser.value()?.parse()?,

            Short('^') | Long("interleaved") => args.interleaved = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Short('F') | Long("force") => args.force = true,
            Short('s') | Long("seed") => args.seed = Some(parser.value()?.parse()?),
            Long("debug") => args.debug = true,
            Long("strobealign") => args.strobealign = parser.value()?.parse()?,
            Long("samtools") => args.samtools = parser.value()?.parse()?,

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

fn locus_name_matches<'a>(path: &'a Path, subset_loci: &FnvHashSet<String>) -> Option<&'a str> {
    match path.file_name().unwrap().to_str() {
        None => {
            log::error!("Skipping directory {:?} - filename is not a valid UTF-8", path);
            None
        }
        Some(name) => {
            if !subset_loci.is_empty() && !subset_loci.contains(name) {
                log::trace!("Skipping locus {} (not in the subset loci)", name);
                None
            } else if name.starts_with(".") {
                log::trace!("Skipping hidden directory {}", name);
                None
            } else {
                Some(name)
            }
        }
    }
}

struct LocusData {
    set: ContigSet,
    db_locus_dir: PathBuf,
    /// Output directory with locus data.
    out_dir: PathBuf,
    /// Temporary file with recruited reads.
    tmp_reads_filename: PathBuf,
    /// Final file with recruited reads.
    reads_filename: PathBuf,
    /// Temporary file with read alignments to the haplotypes.
    tmp_aln_filename: PathBuf,
    /// Final file with read alignments to the haplotypes.
    aln_filename: PathBuf,
    /// File with haplotype-pair likelihoods.
    lik_filename: PathBuf,
}

impl LocusData {
    fn new(set: ContigSet, db_locus_dir: &Path, out_loci_dir: &Path, force: bool) -> io::Result<Self> {
        let out_dir = out_loci_dir.join(set.tag());
        if out_dir.exists() && force {
            fs::remove_dir_all(&out_dir)?;
        }
        ext::sys::mkdir(&out_dir)?;
        Ok(Self {
            db_locus_dir: db_locus_dir.to_owned(),
            tmp_reads_filename: out_dir.join("reads.tmp.fq.gz"),
            reads_filename: out_dir.join("reads.fq.gz"),
            tmp_aln_filename: out_dir.join("aln.tmp.bam"),
            aln_filename: out_dir.join("aln.bam"),
            lik_filename: out_dir.join("lik.csv.gz"),
            set, out_dir,
        })
    }
}

/// Loads all loci from the database. If `subset_loci` is not empty, only loads loci that are contained in it.
/// If `force`, delete old data in the corresponding output directories.
fn load_loci(
    db_path: &Path,
    out_path: &Path,
    subset_loci: &FnvHashSet<String>,
    force: bool
) -> io::Result<Vec<LocusData>>
{
    log::info!("Loading database.");
    let db_loci_dir = db_path.join(paths::LOCI_DIR);
    let out_loci_dir = out_path.join(paths::LOCI_DIR);
    ext::sys::mkdir(&out_loci_dir)?;
    if force {
        log::warn!("Force flag is set: overwriting output directories {}/*", ext::fmt::path(&out_loci_dir));
    }

    let mut loci = Vec::new();
    let mut total_entries = 0;
    for entry in fs::read_dir(&db_loci_dir)? {
        let entry = entry?;
        if !entry.file_type()?.is_dir() {
            continue;
        }

        total_entries += 1;
        let path = entry.path();
        if let Some(name) = locus_name_matches(&path, subset_loci) {
            let fasta_filename = path.join(paths::LOCUS_FASTA);
            let kmers_filename = path.join(paths::KMERS);
            match ContigSet::load(name, &fasta_filename, &kmers_filename, ()) {
                Ok(set) => loci.push(LocusData::new(set, &path, &out_loci_dir, force)?),
                Err(e) => log::error!("Could not load locus information from {}: {:?}", ext::fmt::path(&path), e),
            }
        }
    }
    let n = loci.len();
    if n < total_entries {
        log::info!("Loaded {} loci, skipped {} directories", n, total_entries - n);
    } else {
        log::info!("Loaded {} loci", n);
    }
    Ok(loci)
}

/// Recruits reads to all loci, where neither reads nor alignments are available.
fn recruit_reads(loci: &[LocusData], args: &Args) -> io::Result<()> {
    let filt_loci: Vec<&LocusData> = loci.iter()
        .filter(|locus| !locus.reads_filename.exists() && !locus.aln_filename.exists())
        .collect();
    if filt_loci.is_empty() {
        log::info!("Skipping read recruitment");
        return Ok(());
    }
    let n_filt_loci = filt_loci.len();
    if n_filt_loci < loci.len() {
        log::info!("Skipping read recruitment to {} loci", loci.len() - n_filt_loci);
    }

    log::info!("Generating recruitment targets");
    let mut targets = recruit::Targets::new(&args.recr_params);
    let mut writers = Vec::with_capacity(n_filt_loci);
    let mut total_seqs = 0;
    for locus in filt_loci.iter() {
        let mut fasta_reader = fastx::Reader::from_path(&locus.db_locus_dir.join(paths::LOCUS_FASTA_ALL))?;
        let locus_all_seqs = fasta_reader.read_all()?;
        total_seqs += locus_all_seqs.len();
        targets.add(locus_all_seqs.iter().map(NamedSeq::seq));
        writers.push(ext::sys::create_gzip(&locus.tmp_reads_filename)?);
    }
    log::info!("Collected {} minimizers across {} loci and {} sequences", targets.total_minimizers(),
        n_filt_loci, total_seqs);

    // Cannot put reader into a box, because `FastxRead` has a type parameter.
    if args.input.len() == 1 && !args.interleaved {
        let reader = fastx::Reader::from_path(&args.input[0])?;
        targets.recruit(reader, writers, args.threads, args.recr_params.chunk_size)?;
    } else if args.interleaved {
        let reader = fastx::PairedEndInterleaved::new(fastx::Reader::from_path(&args.input[0])?);
        targets.recruit(reader, writers, args.threads, args.recr_params.chunk_size)?;
    } else {
        let reader = fastx::PairedEndReaders::new(
            fastx::Reader::from_path(&args.input[0])?,
            fastx::Reader::from_path(&args.input[1])?);
        targets.recruit(reader, writers, args.threads, args.recr_params.chunk_size)?;
    }
    for locus in filt_loci.iter() {
        fs::rename(&locus.tmp_reads_filename, &locus.reads_filename)?;
    }
    Ok(())
}

fn map_reads(locus: &LocusData, args: &Args) -> Result<(), Error> {
    if locus.aln_filename.exists() {
        log::info!("    Skipping read mapping");
        return Ok(());
    }

    let in_fasta = locus.db_locus_dir.join(paths::LOCUS_FASTA);
    log::info!("    Mapping reads to {}", ext::fmt::path(&in_fasta));
    // Output at most this number of alignments per read.
    let n_locs = min(30000, locus.set.len() * 30).to_string();
    let start = Instant::now();
    let mut strobealign_cmd = Command::new(&args.strobealign);
    strobealign_cmd.args(&[
        "--no-progress",
        // NOTE: Provide single input FASTQ, and do not use --interleaved option.
        // This is because Strobealign starts to connect read pairs in multiple ways,
        // which produces too large output, and takes a lot of time.
        "-M", &n_locs,   // Try as many secondary locations as possible.
        "-N", &n_locs,   // Output as many secondary alignments as possible.
        "-S", "0.5",     // Try candidate sites with score >= 0.2 * best score.
        "-f", "0",       // Do not discard repetitive minimizers.
        "-k", "15",      // Use smaller minimizers to get more matches.
        // TODO: Use -r MEAN_READ_LEN
        "-t", &args.threads.to_string()])
        .arg(&in_fasta) // Input fasta.
        .arg(&locus.reads_filename) // Input FASTQ.
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = strobealign_cmd.spawn()?;
    let child_stdout = child.stdout.take().unwrap();

    let mut samtools_cmd = Command::new(&args.samtools);
    samtools_cmd.args(&["view", "-b"]) // Output BAM.
        .arg("-o").arg(&locus.tmp_aln_filename)
        .stdin(Stdio::from(child_stdout));

    log::debug!("    {} | {}", ext::fmt::command(&strobealign_cmd), ext::fmt::command(&samtools_cmd));
    let samtools_output = samtools_cmd.output()?;
    let bwa_output = child.wait_with_output()?;
    log::debug!("    Finished in {}", ext::fmt::Duration(start.elapsed()));
    if !bwa_output.status.success() {
        return Err(Error::SubprocessFail(bwa_output));
    } else if !samtools_output.status.success() {
        return Err(Error::SubprocessFail(samtools_output));
    }
    fs::rename(&locus.tmp_aln_filename, &locus.aln_filename)?;
    Ok(())
}

fn analyze_locus(
    locus: &LocusData,
    bg_distr: &BgDistr,
    scheme: &scheme::Scheme,
    cached_distrs: &Arc<CachedDepthDistrs>,
    mut rng: ext::rand::XoshiroRng,
    args: &Args,
) -> Result<(), Error>
{
    log::info!("Analyzing {}", locus.set.tag());
    map_reads(locus, &args)?;

    log::info!("    [{}] Calculating read alignment probabilities", locus.set.tag());
    let bam_reader = bam::Reader::from_path(&locus.aln_filename)?;
    let contigs = locus.set.contigs();

    let err_prof = bg_distr.error_profile();
    let is_paired = args.is_paired_end();
    let mut locs = if args.debug {
        let mut reads_writer = ext::sys::create_gzip(&locus.out_dir.join("reads.csv.gz"))?;
        writeln!(reads_writer, "{}", locs::CSV_HEADER)?;
        PrelimAlignments::load(bam_reader, is_paired, contigs, err_prof, &args.assgn_params, reads_writer)?
    } else {
        PrelimAlignments::load(bam_reader, is_paired, contigs, err_prof, &args.assgn_params, io::sink())?
    };

    let all_alns = locs.identify_locations(bg_distr.insert_distr(), &args.assgn_params);
    let contig_windows = ContigWindows::new_all(&locus.set, bg_distr.depth(), &args.assgn_params);

    let contig_ids: Vec<_> = contigs.ids().collect();
    let tuples = ext::vec::Tuples::repl_combinations(&contig_ids, usize::from(args.ploidy));
    // let mut tuples = ext::vec::Tuples::new(2);
    // tuples.push(&[contigs.get_id("HG00733.1"), contigs.get_id("HG00735.1")]);
    // tuples.push(&[contigs.get_id("HG00673.1"), contigs.get_id("HG02080.1")]);

    let lik_writer = ext::sys::create_gzip(&locus.lik_filename)?;
    let data = scheme::Data {
        scheme: scheme.clone(),
        params: args.assgn_params.clone(),
        contigs: Arc::clone(&contigs),
        cached_distrs: Arc::clone(&cached_distrs),
        debug: args.debug,
        all_alns, contig_windows, tuples,
    };
    scheme::solve(data, lik_writer, &locus.out_dir, &mut rng, args.threads)
}

pub(super) fn run(argv: &[String]) -> Result<(), Error> {
    let timer = Instant::now();
    let args = parse_args(argv)?.validate()?;
    let db_dir = args.database.as_ref().unwrap();
    let out_dir = args.output.as_ref().unwrap();

    let bg_stream = io::BufReader::new(fs::File::open(out_dir.join(paths::BG_DIR).join(paths::BG_DISTR))?);
    let bg_stream = flate2::bufread::GzDecoder::new(bg_stream);
    let bg_distr = BgDistr::load(&json::parse(&io::read_to_string(bg_stream)?)?)?;
    let cached_distrs = Arc::new(CachedDepthDistrs::new(&bg_distr));
    validate_param!(bg_distr.depth().window_padding() <= args.assgn_params.boundary_size,
        "Window padding ({}) must not exceed boundary size ({})", bg_distr.depth().window_padding(),
        args.assgn_params.boundary_size);

    let loci = load_loci(db_dir, out_dir, &args.subset_loci, args.force)?;
    recruit_reads(&loci, &args)?;

    let scheme = match args.solvers.as_ref() {
        Some(filename) => scheme::Scheme::from_json(&ext::sys::load_json(filename)?)?,
        None => scheme::Scheme::default(),
    };

    let mut rng = ext::rand::init_rng(args.seed);
    for locus in loci.iter() {
        let rng_clone = rng.clone();
        // Jump over 2^192 random numbers. This way, all loci have independent random numbers.
        rng.long_jump();
        analyze_locus(locus, &bg_distr, &scheme, &cached_distrs, rng_clone, &args)?;
    }
    log::info!("Success. Total time: {}", ext::fmt::Duration(timer.elapsed()));
    Ok(())
}
