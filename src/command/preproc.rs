//! Preprocess WGS dataset.

use std::{
    fs, thread,
    ffi::OsStr,
    io::{self, Write, BufRead},
    fmt::Write as FmtWrite,
    cmp::{min, max},
    path::{Path, PathBuf},
    process::{Stdio, Command, ChildStdin},
    time::Instant,
    ops::Deref,
    sync::Arc,
    str::FromStr,
};
use bio::io::fasta;
use colored::Colorize;
use const_format::str_repeat;
use htslib::bam;
use regex::{Regex, RegexBuilder};
use crate::{
    ext::{
        self,
        rand::XoshiroRng,
        fmt::{PrettyU32, PrettyU64},
    },
    err::{Error, error, validate_param, add_path},
    seq::{
        Interval,
        kmers::{JfKmerGetter, KmerCounts},
        fastx::{self, FastxRead},
        recruit::{self, RecruitableRecord},
        cigar::{self, Cigar},
        aln::{NamedAlignment, Alignment, ReadEnd, OpCounter},
        contigs::{ContigNames, ContigSet, GenomeVersion},
    },
    bg::{
        self, BgDistr, Technology, SequencingInfo, Windows,
        insertsz::{self, InsertDistr},
        err_prof::{ErrorProfile, SingleEditDistCache},
        depth::ReadDepth,
        ser::{JsonSer, json_get},
    },
    math::implies,
};
use super::{paths, genotype};

fn check_filename(
    filename: &Path,
    should_match: &Regex,
    exp_format: &'static str,
    shouldnt_match: &Regex,
    wrong_format: &'static str,
) -> crate::Result<()> {
    if filename == OsStr::new("-") || filename.starts_with("/dev") {
        log::error!("Stdin and other /dev/* files may not be supported. Continuing for now");
        return Ok(());
    } else if filename == OsStr::new("!") {
        // Ignore this filename.
        return Ok(());
    } else if !filename.exists() {
        log::error!("Input file {} does not exist. Continuing for now", ext::fmt::path(filename));
    }

    if let Some(s) = filename.to_str() {
        if shouldnt_match.is_match(s) {
            return Err(error!(InvalidInput,
                "Incorrect file format for {} (expected {}, found {}). Please check -i and -a arguments",
                ext::fmt::path(filename), exp_format, wrong_format));
        } else if !should_match.is_match(s) {
            log::warn!("Could not guess file format for {} (expected {}). Can lead to problems later",
                ext::fmt::path(filename), exp_format);
        }
    } else {
        log::warn!("Could not guess file format for {} (expected {}). Can lead to problems later",
            ext::fmt::path(filename), exp_format);
    }
    Ok(())
}

/// Collection of input files, same as for `genotype`.
pub struct InputFiles {
    pub reads1: Vec<PathBuf>,
    pub reads2: Vec<PathBuf>,
    pub alns: Vec<PathBuf>,
    pub in_list: Option<PathBuf>,
    pub reference: Option<PathBuf>,
    pub interleaved: bool,
    pub no_index: bool,
}

impl Default for InputFiles {
    fn default() -> Self {
        Self {
            reads1: Vec::new(),
            reads2: Vec::new(),
            alns: Vec::new(),
            in_list: None,
            reference: None,
            interleaved: false,
            no_index: false,
        }
    }
}

impl InputFiles {
    /// Sets input files (except for the reference) based on `in_list` with lines
    /// `<flag>  <file>  [<file2>]`.
    /// Flag is one of `p` (paired-end), `s` (single-end), `pi` (paired-end interleaved),
    /// `a` (alignment), `u` (unmapped BAM/CRAM file) and `ui` (unmapped interleaved BAM/CRAM file).
    ///
    /// Paired-end (`p`) line may contain either two files or one file with `*` inside, which is replaced with 1 and 2.
    pub fn fill_from_inlist(&mut self) -> crate::Result<()> {
        let Some(list_filename) = self.in_list.take() else { return Ok(()) };
        validate_param!(self.reads1.is_empty() && self.reads2.is_empty() && self.alns.is_empty(),
            "Input list (-I) cannot be used together with other input files (-i/-a).");
        validate_param!(!self.interleaved && !self.no_index,
            "Input list (-I) and --interleaved/--no-index cannot be provided together, please see README");
        let dirname = list_filename.parent();

        let mut prev_flag: Option<String> = None;
        for line in ext::sys::open(&list_filename)?.lines() {
            let line = line.map_err(add_path!(!))?;
            if line.starts_with('#') {
                continue;
            }
            let trimmed_line = line.trim_end();
            let split: smallvec::SmallVec::<[&str; 3]> = trimmed_line.split_whitespace().collect();
            validate_param!(split.len() == 2 || (split[0] == "p" && split.len() == 3),
                "Incorrect number of arguments in input line {:?}", trimmed_line);
            let flag = split[0].to_string();
            let filename1 = split[1];
            validate_param!(prev_flag.as_ref().map(|val| val == &flag).unwrap_or(true),
                "All lines in the input list (-I) must contain the same flag ({} != {})", prev_flag.unwrap(), &flag);
            match &flag as &str {
                "s" => {
                    // validate_param!(!self.interleaved,
                    //     "Cannot provide single-end and paired-end interleaved files together");
                    self.reads1.push(ext::sys::add_dir(dirname, filename1));
                }
                "pi" => {
                    // validate_param!(self.interleaved || self.reads1.is_empty(),
                    //     "Cannot provide single-end and paired-end interleaved files together");
                    self.interleaved = true;
                    self.reads1.push(ext::sys::add_dir(dirname, filename1));
                }
                "a" => {
                    // validate_param!(!self.no_index,
                    //     "Cannot provide indexed and non-indexed alignment files together");
                    self.alns.push(ext::sys::add_dir(dirname, filename1));
                }
                "u" => {
                    // validate_param!(self.no_index || self.alns.is_empty(),
                    //     "Cannot provide indexed and non-indexed alignment files together");
                    // validate_param!(!self.interleaved,
                    //     "Cannot provide single-end and paired-end interleaved files together");
                    self.no_index = true;
                    self.alns.push(ext::sys::add_dir(dirname, filename1));
                }
                "ui" => {
                    // validate_param!(self.no_index || self.alns.is_empty(),
                    //     "Cannot provide indexed and non-indexed alignment files together");
                    // validate_param!(self.interleaved || self.alns.is_empty(),
                    //     "Cannot provide single-end and paired-end interleaved files together");
                    self.no_index = true;
                    self.interleaved = true;
                    self.alns.push(ext::sys::add_dir(dirname, filename1));
                }
                "p" => {
                    if split.len() == 3 {
                        self.reads1.push(ext::sys::add_dir(dirname, filename1));
                        self.reads2.push(ext::sys::add_dir(dirname, split[2]));
                    } else {
                        validate_param!(filename1.contains('*'),
                            "Cannot parse line {:?}: paired-end entry requires either two files, or one file with `*`",
                            trimmed_line);
                        self.reads1.push(ext::sys::add_dir(dirname, &filename1.replace('*', "1")));
                        self.reads2.push(ext::sys::add_dir(dirname, &filename1.replace('*', "2")));
                    }
                }
                rem => return Err(error!(ParsingError, "Cannot parse line {:?}: unexpected flag {}", trimmed_line, rem)),
            }
            prev_flag = Some(flag);
        }
        Ok(())
    }

    pub fn validate(&self, ref_required: bool) -> crate::Result<()> {
        let n1 = self.reads1.len();
        let n2 = self.reads2.len();
        let m = self.alns.len();

        validate_param!(n1 == n2 || n2 == 0, "Cannot mix single-end and paired-end reads");
        validate_param!(n1 > 0 || m > 0, "Neither reads (-i) nor alignments (-a) are provided");
        validate_param!(n1 == 0 || m == 0, "Cannot use reads (-i) and alignments (-a) together");
        validate_param!(implies(self.interleaved, n2 == 0),
            "Second end reads are specified together with --interleaved");
            validate_param!(implies(m > 1, self.no_index),
            "Cannot use multilpe indexed BAM/CRAM files (consider --no-index)");

        if ref_required {
            validate_param!(self.reference.is_some(), "Reference file (-r) is not provided");
        } else if self.alns.iter().any(|filename| filename.extension()
                .map(|ext| ext == "cram" || ext == "CRAM").unwrap_or(false)) {
            validate_param!(self.reference.is_some(), "Input CRAM file (-a) requires a reference file (-r)");
        }

        // Should only be run once, so there is no need for lazy static.
        let re_fastx = RegexBuilder::new(r"\.f(ast)?[aq](\.[^.]{1,3})?$").case_insensitive(true).build().unwrap();
        let re_bam = RegexBuilder::new(r"\.(bam|cram)$").case_insensitive(true).build().unwrap();
        let fastx_descr = "fasta/fastq[.gz]";
        let bam_descr = "bam/cram";

        for filename in self.reads1.iter().chain(&self.reads2) {
            check_filename(filename, &re_fastx, fastx_descr, &re_bam, bam_descr)?;
        }
        for filename in self.alns.iter() {
            check_filename(filename, &re_bam, bam_descr, &re_fastx, fastx_descr)?;
        }
        Ok(())
    }

    /// Returns true if need to process indexed alignment file.
    pub fn has_indexed_alignment(&self) -> bool {
        !self.alns.is_empty() && !self.no_index
    }

    /// Returns true if input reads are paired-end. Panics for indexed alignment file.
    pub fn is_paired_end(&self) -> bool {
        assert!(!self.has_indexed_alignment());
        !self.reads2.is_empty() || self.interleaved
    }

    /// Returns sum file size across first and second reads, as well as alignments.
    pub fn sum_file_size(&self) -> crate::Result<u64> {
        let mut s = 0;
        for filename in self.reads1.iter().chain(&self.reads2).chain(&self.alns) {
            s += fs::metadata(filename).map_err(add_path!(filename))?.len();
        }
        Ok(s)
    }
}

struct Args {
    in_files: InputFiles,
    jf_counts: Option<PathBuf>,
    output: Option<PathBuf>,
    similar_dataset: Option<PathBuf>,
    bg_region: Option<String>,

    technology: Technology,
    min_mapq: u8,

    threads: u16,
    recr_threads: f64,
    skip_recruitment: bool,

    subsampling_rate: f64,
    seed: Option<u64>,

    rerun: super::Rerun,
    strobealign: PathBuf,
    minimap: PathBuf,
    samtools: PathBuf,
    jellyfish: PathBuf,
    debug: bool,
    debug_head: Option<u64>,
    /// Was technology explicitely provided?
    explicit_technology: bool,
    /// When calculating insert size distributions and read error profiles,
    /// ignore reads with `clipping > max_clipping * read_len`.
    max_clipping: f64,
    use_file_size: bool,

    bg_params: bg::Params,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            in_files: InputFiles::default(),
            jf_counts: None,
            output: None,
            similar_dataset: None,
            bg_region: None,

            technology: Technology::Illumina,
            min_mapq: 30,

            threads: 8,
            recr_threads: 0.4,
            skip_recruitment: false,

            subsampling_rate: 1.0,
            seed: None,

            rerun: super::Rerun::None,
            debug: false,
            debug_head: None,
            strobealign: PathBuf::from("strobealign"),
            minimap: PathBuf::from("minimap2"),
            samtools: PathBuf::from("samtools"),
            jellyfish: PathBuf::from("jellyfish"),
            explicit_technology: false,
            max_clipping: 0.02,
            bg_params: bg::Params::default(),
            use_file_size: false,
        }
    }
}

impl Args {
    /// Validate arguments, modifying some, if needed.
    fn validate(mut self) -> crate::Result<Self> {
        self.in_files.validate(true)?;
        self.threads = max(self.threads, 1);

        if self.in_files.has_indexed_alignment() {
            validate_param!(self.similar_dataset.is_none(),
                "Similar dataset (-~) cannot be used with indexed alignments:\
                    preprocessing without the similar dataset is faster and more accurate");
        } else {
            let paired_end_allowed = self.technology.paired_end_allowed();
            validate_param!(!self.in_files.is_paired_end() || paired_end_allowed,
                "Paired end reads are not supported by {}", self.technology.long_name());
            if !self.in_files.is_paired_end() && paired_end_allowed {
                log::warn!("Running in single-end mode.");
            }
        }

        if self.debug_head.is_some() {
            self.skip_recruitment = true;
        }
        if !self.skip_recruitment {
            validate_param!(self.recr_threads >= 0.0, "Number of recruitment threads ({}) must be non-negative",
                self.recr_threads);
            validate_param!(self.recr_threads < 1.0 || self.recr_threads.fract() < 1e-8,
                "Number of recruitment threads ({}) must be either integer, or smaller than 1",
                self.recr_threads);
        }

        validate_param!(0.0 < self.subsampling_rate && self.subsampling_rate <= 1.0,
            "Subsampling rate ({}) must be in (0, 1].", self.subsampling_rate);
        if self.in_files.has_indexed_alignment() || self.similar_dataset.is_some() || self.debug_head.is_some() {
            self.subsampling_rate = 1.0;
        }

        validate_param!(self.jf_counts.is_some(), "Jellyfish counts are not provided (see -j/--jf-counts)");
        validate_param!(self.output.is_some(), "Output directory is not provided (see -o/--output)");

        if self.technology == Technology::Illumina {
            self.strobealign = ext::sys::find_exe(self.strobealign)?;
        } else {
            self.minimap = ext::sys::find_exe(self.minimap)?;
        }
        self.samtools = ext::sys::find_exe(self.samtools)?;
        self.jellyfish = ext::sys::find_exe(self.jellyfish)?;

        validate_param!(0.0 <= self.max_clipping && self.max_clipping <= 1.0,
            "Max clipping ({:.5}) must be within [0, 1]", self.max_clipping);
        self.bg_params.validate()?;
        Ok(self)
    }
}

fn print_help(extended: bool) {
    const KEY: usize = 18;
    const VAL: usize = 5;
    const EMPTY: &'static str = str_repeat!(" ", KEY + VAL + 5);

    let defaults = Args::default();
    println!("{}", "Preprocess WGS dataset.".yellow());

    println!("\n{}", "Usage:".bold());
    println!("    {} preproc (-i reads1.fq [reads2.fq] | -a reads.bam [--no-index] | -I in-list) \\", super::PROGRAM);
    println!("        -r reference.fa -j counts.jf -o out [args]");
    println!("    {} preproc -i/-a/-I <input> -~ similar -r reference.fa -j counts.jf -o out [args]", super::PROGRAM);
    if !extended {
        println!("\nThis is a {} help message. Please use {} to see the full help.",
            "short".red(), "-H/--full-help".green());
    }

    println!("\n{}  (please see {} for more information on {}/{}/{} arguments)",
        "Input/output arguments:".bold(), "README".italic(), "-i".green(), "-a".green(), "-I".green());
    println!("    {:KEY$} {:VAL$}  Reads 1 and 2 in FASTA or FASTQ format, optionally gzip compressed.\n\
        {EMPTY}  Reads 1 are required, reads 2 are optional.",
        "-i, --input".green(), "FILE+".yellow());
    println!("    {:KEY$} {:VAL$}  Reads in BAM/CRAM format, mutually exclusive with {}.\n\
        {EMPTY}  By default, mapped, sorted and indexed BAM/CRAM file is expected,\n\
        {EMPTY}  please specify {} otherwise.",
        "-a, --alignment".green(), "FILE".yellow(), "-i/--input".green(), "--no-index".green());
    println!("    {:KEY$} {:VAL$}  File with input filenames (see {}).",
        "-I, --in-list".green(), "FILE".yellow(), "README".italic());
    println!("    {:KEY$} {:VAL$}  Reference FASTA file. Must be indexed with FAIDX.",
        "-r, --reference".green(), "FILE".yellow());
    println!("    {:KEY$} {:VAL$}  Jellyfish k-mer counts (see {}).",
        "-j, --jf-counts".green(), "FILE".yellow(), "README".italic());
    println!("    {:KEY$} {:VAL$}  Output directory.",
        "-o, --output".green(), "DIR".yellow());
    println!("    {:KEY$} {:VAL$}  This dataset is similar to already preprocessed dataset.\n\
        {EMPTY}  {}. Only utilizes difference in the number of reads.",
        "-~, --like".green(), "DIR".yellow(), "Use with care".red());

    println!("\n{}", "Optional arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Interleaved paired-end reads in single input file.",
        "-^, --interleaved".green(), super::flag());
    println!("    {:KEY$} {:VAL$}  Use input full BAM/CRAM file ({}) without index.\n\
        {EMPTY}  Single-end and paired-end interleaved ({}) data is allowed.",
        "    --no-index".green(), super::flag(), "-a".green(), "-^".green());
    println!("    {:KEY$} {:VAL$}  Sequencing technology [{}]:\n\
        {EMPTY}  sr  | illumina : short-read sequencing,\n\
        {EMPTY}    hifi         : PacBio HiFi,\n\
        {EMPTY}  pb  | pacbio   : PacBio CLR,\n\
        {EMPTY}  ont | nanopore : Oxford Nanopore.",
        "-t, --tech".green(), "STR".yellow(), super::fmt_def(defaults.technology));
    println!("    {:KEY$} {:VAL$}  Preprocess WGS data based on this background region,\n\
        {EMPTY}  preferably >3 Mb and without many duplications.\n\
        {EMPTY}  Default regions are defined for CHM13, GRCh38 and GRCh37.",
        "-b, --bg-region".green(), "STR".yellow());

    if extended {
        println!("\n{}", "Read recruitment:".bold());
        println!("    {:KEY$} {:VAL$}  Skip read recruitment before read mapping.\n\
            {EMPTY}  Otherwise, read recruitment is executed with default parameters.",
            "    --skip-recruit".green(), super::flag());
        println!("    {:KEY$} {:VAL$}  Number of threads, used for read recruitment [{}].\n\
            {EMPTY}  Fraction of the total number of threads, if under 1.",
            "    --recr-threads".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.recr_threads));

        println!("\n{}", "Insert size and error profile estimation:".bold());
        println!("    {:KEY$} {:VAL$}  Ignore reads with mapping quality less than {} [{}].",
            "-q, --min-mapq".green(), "INT".yellow(), "INT".yellow(), super::fmt_def(defaults.min_mapq));
        println!("    {:KEY$} {:VAL$}  Ignore reads with soft/hard clipping over {} * read length [{}].",
            "-c, --max-clipping".green(), "FLOAT".yellow(), "FLOAT".yellow(),
            super::fmt_def_f64(defaults.max_clipping));
        println!("    {:KEY$} {:VAL$}\n\
            {EMPTY}  Two p-value thresholds for filtering recruited reads:\n\
            {EMPTY}  on insert size [{}], and on edit distance [{}].",
            "    --pval-thresh".green(), "FLOAT FLOAT".yellow(),
            super::fmt_def_f64(defaults.bg_params.insert_pval),
            super::fmt_def_f64(defaults.bg_params.edit_pval));

        println!("\n{}", "Background read depth estimation:".bold());
        println!("    {:KEY$} {:VAL$}  Specie ploidy [{}].",
            "-p, --ploidy".green(), "INT".yellow(), super::fmt_def(defaults.bg_params.depth.ploidy));
        println!("    {:KEY$} {:VAL$}  Count read depth in windows of this size [{}].\n\
            {EMPTY}  Default: half of the mean read length.",
            "-w, --window".green(), "INT".yellow(), super::fmt_def("auto"));
        println!("    {:KEY$} {:VAL$}  Skip {} bp near the edge of the background region [{}].",
            "    --boundary".green(), "INT".yellow(), "INT".yellow(),
            super::fmt_def(PrettyU32(defaults.bg_params.depth.boundary_size)));
        println!("    {:KEY$} {:VAL$}  Ignore windows, where less than {}% k-mers are unique [{}].",
            "    --kmer-perc".green(), "FLOAT".yellow(), "FLOAT".yellow(),
            super::fmt_def_f64(defaults.bg_params.depth.uniq_kmer_perc));
        println!("    {:KEY$} {:VAL$}  This fraction of all windows is used in LOESS during read depth\n\
            {EMPTY}  estimation [{}]. Smaller values lead to less robust estimates,\n\
            {EMPTY}  larger values - to similar estimates across different GC-contents.",
            "    --frac-windows".green(), "FLOAT".yellow(),
            super::fmt_def_f64(defaults.bg_params.depth.frac_windows));
        println!("    {:KEY$} {}\n\
            {EMPTY}  Read depth estimates are blured for windows with extreme GC-content\n\
            {EMPTY}  (less than {} windows with smaller/larger GC). There, read depth\n\
            {EMPTY}  is set to the last non-extreme depth, while variance is increased\n\
            {EMPTY}  by a {} factor for each addition GC value [{} {}].",
            "    --blur-extreme".green(), "INT FLOAT".yellow(), "INT".yellow(), "FLOAT".yellow(),
            super::fmt_def(defaults.bg_params.depth.min_tail_obs),
            super::fmt_def_f64(defaults.bg_params.depth.tail_var_mult));
        println!("    {:KEY$} {:VAL$}  Estimate read depth by comparing file sizes\n\
            {EMPTY}  with similar dataset ({}). {}.",
            "    --filesize".green(), super::flag(), "-~".green(), "Use with extreme care".red());

        println!("\n{}", "Subsampling:".bold());
        println!("    {:KEY$} {:VAL$}  Subsample input reads by this fraction [{}].\n\
            {EMPTY}  Smaller values increase speed, but impact accuracy.",
            "    --subsample".green(), "FLOAT".yellow(), super::fmt_def_f64(defaults.subsampling_rate));
        println!("    {:KEY$} {:VAL$}  Subsampling seed (optional). Ensures reproducibility\n\
            {EMPTY}  for the same input and program version.",
            "    --seed".green(), "INT".yellow());
    }

    println!("\n{}", "Execution arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Number of threads [{}].",
        "-@, --threads".green(), "INT".yellow(), super::fmt_def(defaults.threads));
    println!("    {:KEY$} {:VAL$}  Rerun mode [{}]. Rerun everything ({}); do not rerun\n\
        {EMPTY}  read mapping ({}); do not rerun ({}).",
        "    --rerun".green(), "STR".yellow(), super::fmt_def(defaults.rerun),
        "all".yellow(), "part".yellow(), "none".yellow());
    println!("    {:KEY$} {:VAL$}  Save debug CSV files.",
        "    --debug".green(), "INT".yellow());
    if extended {
        println!("    {:KEY$} {:VAL$}  Instead of the full preprocessing, map first {} reads to the\n\
            {EMPTY}  reference and extract insert sizes and error profiles from them.\n\
            {EMPTY}  Can be used to validate similar dataset preprocessing ({}).\n\
            {EMPTY}  Recommended to use together with {}.",
            "    --debug-head".green(), "INT".yellow(), "INT".yellow(), "-~".green(), "-q 60".underline());
        println!("    {:KEY$} {:VAL$}  Strobealign executable [{}].",
            "    --strobealign".green(), "EXE".yellow(), super::fmt_def(defaults.strobealign.display()));
        println!("    {:KEY$} {:VAL$}  Minimap2 executable    [{}].",
            "    --minimap".green(), "EXE".yellow(), super::fmt_def(defaults.minimap.display()));
        println!("    {:KEY$} {:VAL$}  Samtools executable    [{}].",
            "    --samtools".green(), "EXE".yellow(), super::fmt_def(defaults.samtools.display()));
        println!("    {:KEY$} {:VAL$}  Jellyfish executable   [{}].",
            "    --jellyfish".green(), "EXE".yellow(), super::fmt_def(defaults.jellyfish.display()));
    }

    println!("\n{}", "Other arguments:".bold());
    println!("    {:KEY$} {:VAL$}  Show short help message.", "-h, --help".green(), "");
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
            Short('i') | Long("input") => {
                let mut values = parser.values()?.take(2);
                args.in_files.reads1.push(values.next().expect("First argument is always present").parse()?);
                if let Some(val) = values.next() {
                    args.in_files.reads2.push(val.parse()?);
                }
            }
            Short('a') | Long("aln") | Long("alns") | Long("alignment") | Long("alignments") =>
                args.in_files.alns.push(parser.value()?.parse()?),
            Short('I') | Long("in-list") | Long("input-list") => args.in_files.in_list = Some(parser.value()?.parse()?),
            Short('r') | Long("reference") => args.in_files.reference = Some(parser.value()?.parse()?),
            Short('o') | Long("output") => args.output = Some(parser.value()?.parse()?),
            Short('t') | Long("tech") | Long("technology") => {
                args.explicit_technology = true;
                args.technology = parser.value()?.parse()?;
            }
            Short('j') | Long("jf-counts") => args.jf_counts = Some(parser.value()?.parse()?),
            Short('~') | Long("like") => args.similar_dataset = Some(parser.value()?.parse()?),
            Short('b') | Long("bg") | Long("bg-region") => args.bg_region = Some(parser.value()?.parse()?),

            Short('q') | Long("min-mapq") | Long("min-mq") => args.min_mapq = parser.value()?.parse()?,
            Short('c') | Long("max-clip") | Long("max-clipping") => args.max_clipping = parser.value()?.parse()?,
            Long("pval-thresh") | Long("pval-threshold") | Long("pvalue-threshold") => {
                args.bg_params.insert_pval = parser.value()?.parse()?;
                args.bg_params.edit_pval = parser.value()?.parse()?;
            }

            Short('p') | Long("ploidy") => args.bg_params.depth.ploidy = parser.value()?.parse()?,
            Short('w') | Long("window") => {
                let val = parser.value()?;
                args.bg_params.depth.window_size = if val == "auto" {
                    None
                } else {
                    Some(val.parse()?)
                };
            }
            Long("boundary") => args.bg_params.depth.boundary_size = parser.value()?.parse::<PrettyU32>()?.get(),
            Long("kmer-perc") | Long("kmer-percentage") | Long("kmer-percentile") =>
                args.bg_params.depth.uniq_kmer_perc = parser.value()?.parse()?,
            Long("frac-windows") | Long("fraction-windows") =>
                args.bg_params.depth.frac_windows = parser.value()?.parse()?,
            Long("blur-extreme") => {
                args.bg_params.depth.min_tail_obs = parser.value()?.parse()?;
                args.bg_params.depth.tail_var_mult = parser.value()?.parse()?;
            }

            Long("subsample") => args.subsampling_rate = parser.value()?.parse()?,
            Long("seed") => args.seed = Some(parser.value()?.parse()?),

            Long("debug-head") => args.debug_head = Some(parser.value()?.parse::<PrettyU64>()?.get()),
            Long("skip-recruit") | Long("skip-recr") | Long("skip-recruitment") => args.skip_recruitment = true,
            Short('^') | Long("interleaved") => args.in_files.interleaved = true,
            Long("no-index") => args.in_files.no_index = true,
            Short('@') | Long("threads") => args.threads = parser.value()?.parse()?,
            Long("recr-threads") | Long("recruit-threads") => args.recr_threads = parser.value()?.parse()?,
            Long("rerun") => args.rerun = parser.value()?.parse()?,
            Long("debug") => args.debug = true,
            Long("filesize") | Long("file-size") => args.use_file_size = true,
            Long("strobealign") => args.strobealign = parser.value()?.parse()?,
            Long("minimap") | Long("minimap2") => args.minimap = parser.value()?.parse()?,
            Long("samtools") => args.samtools = parser.value()?.parse()?,
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

/// Returns default background region.
fn default_region(ver: GenomeVersion) -> &'static str {
    match ver {
        GenomeVersion::Chm13 => "chr17:72950001-77450000",
        GenomeVersion::GRCh38 => "chr17:72062001-76562000",
        GenomeVersion::GRCh37 => "chr17:70060001-74560000",
    }
}

/// Returns the first appropriate interval for the contig names.
/// If `bg_region` is set, try to parse it. Otherwise, iterate over `BG_REGIONS`.
///
/// Returns error if no interval is appropriate (chromosome not in the contig set, or interval is out of bounds).
fn select_bg_interval(
    ref_filename: &Path,
    contigs: &Arc<ContigNames>,
    bg_region: &Option<String>
) -> crate::Result<Interval>
{
    if let Some(s) = bg_region {
        let region = Interval::parse(s, contigs).map_err(|_| error!(InvalidInput,
            "Reference genome {} does not contain region {}", ext::fmt::path(ref_filename), s))?;
        return if contigs.in_bounds(&region) {
            Ok(region)
        } else {
            Err(error!(InvalidInput, "Region {} is out of bounds", s))
        };
    }

    let genome = GenomeVersion::guess(contigs)
        .ok_or_else(|| error!(RuntimeError,
            "Could not recognize reference genome. Please provide background region (-b) explicitely, \
            preferably >3 Mb long and without significant duplications"))?;
    let region = default_region(genome);
    log::info!("Recognized {} reference genome, using background region {}", genome, region);
    // Try to crop `chr` if region cannot be found.
    for &crop in &[0, 3] {
        if let Ok(region) = Interval::parse(&region[crop..], contigs) {
            if contigs.in_bounds(&region) {
                return Ok(region);
            } else {
                return Err(error!(InvalidInput, "Region {} is too long for the reference genome", region));
            }
        }
    }
    Err(error!(InvalidInput, "Reference genome {} does not contain any of the default background regions. \
        Please provide background region (-b) explicitely, preferably >3 Mb long and without significant duplications",
        ext::fmt::path(ref_filename)))
}

fn recruit_reads(
    reader: impl FastxRead<Record = impl RecruitableRecord> + 'static,
    writer: impl io::Write + Send + 'static,
    seq_info: &SequencingInfo,
    bg_region: &BgRegion,
    is_paired_end: bool,
    recr_threads: u16,
    sampling: Option<(f64, XoshiroRng)>,
) -> crate::Result<thread::JoinHandle<Result<Option<u64>, Error>>>
{
    let minimizer_kw = seq_info.technology().default_minim_size();
    let match_frac = seq_info.technology().default_match_frac(is_paired_end);
    let recr_params = recruit::Params::new(minimizer_kw, match_frac, recruit::DEFAULT_MATCH_LEN,
        super::recruit::DEFAULT_KMER_THRESH)?;

    let mut target_builder = recruit::TargetBuilder::new(recr_params);
    let contig_set = ContigSet::new(Arc::new(ContigNames::empty()),
        vec![bg_region.padded_sequence.clone()], bg_region.padded_kmer_counts.clone());
    target_builder.add(&contig_set, seq_info.mean_read_len());
    let mut targets = target_builder.finalize();
    targets.hide_recruited_reads();

    let chunk_size = genotype::calculate_chunk_size(genotype::DEFAULT_CHUNK_LENGTH,
        seq_info.mean_read_len(), is_paired_end);
    Ok(thread::spawn(move ||
        targets.recruit(reader, vec![writer], recr_threads, chunk_size, sampling)
            .map(|stats| Some(stats.processed()))))
}

/// Returns handle, which returns the total number of consumed reads/read pairs after finishing.
fn set_mapping_stdin(
    child_stdin: ChildStdin,
    ref_contigs: &Arc<ContigNames>,
    seq_info: &SequencingInfo,
    bg_region: Option<&BgRegion>,
    recr_threads: u16,
    args: &Args,
) -> crate::Result<thread::JoinHandle<Result<Option<u64>, Error>>>
{
    let sampling = if args.subsampling_rate < 1.0 {
        Some((args.subsampling_rate, ext::rand::init_rng(args.seed)))
    } else { None };

    let mut writer = io::BufWriter::new(child_stdin);
    fastx::process_readers!(args.in_files, Some(ref_contigs); let {mut} reader; {
        if let Some(count) = args.debug_head {
            Ok(thread::spawn(move || reader.copy_first(&mut writer, count)
                .map(|_| None) // Do not return the number of reads.
            ))
        } else if args.skip_recruitment {
            Ok(thread::spawn(move || reader.copy_or_subsample(&mut writer, sampling).map(Some)))
        } else {
            recruit_reads(reader, writer, seq_info, bg_region.expect("Background region must be defined"),
                args.in_files.is_paired_end(), recr_threads, sampling)
        }
    })
}

fn create_mapping_command(args: &Args, seq_info: &SequencingInfo, ref_filename: &Path, threads: u16) -> Command {
    let mut command;
    if seq_info.technology() == Technology::Illumina {
        // Strobealign for short reads.
        command = Command::new(&args.strobealign);
        command.args(&[
            "-N", "0",     // Retain 0 additional alignments,
            "-R", "0",     // Do not rescue reads,
            "-U",          // Do not output unmapped reads,
            "-f", "0.001", // Discard more minimizers to speed up alignment,
            "--eqx",       // Output X/= instead of M operations,
            "-t", &threads.to_string(), // Specify the number of threads,
            "-r", &format!("{:.0}", seq_info.mean_read_len()), // Provide mean read length.
            "--no-progress",
            ]);
        if args.in_files.is_paired_end() {
            command.arg("--interleaved");
        }
    } else {
        // Minimap2 for long reads.
        command = Command::new(&args.minimap);
        command.args(&[
            "-a", // Output SAM format,
            "-x", seq_info.technology().minimap_preset(), // Set mapping preset.
            "-N", "2",     // Examine at most 2 additional alignments,
            "--secondary=no", // Do not output secondary alignments,
            "-f", "0.001", // Discard more minimizers to speed up alignment,
            "--eqx",       // Output X/= instead of M operations,
            "-t", &threads.to_string(), // Specify the number of threads,
            ]);
    }
    // Provide paths to the reference and pipe reads.
    command.arg(&ref_filename).arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    command
}

fn first_step_str(args: &Args, bg_region: Option<&BgRegion>) -> String {
    let mut s = String::new();
    if let Some(count) = args.debug_head {
        write!(s, "_head_ -n {}", PrettyU64(count)).unwrap();
    } else {
        if args.skip_recruitment {
            s.push_str("_interleave_");
        } else {
            write!(s, "_recruit_ --target {}", bg_region.unwrap().interval).unwrap();
        }
        if args.subsampling_rate < 1.0 {
            write!(s, " --subsample {}", crate::math::fmt_signif(args.subsampling_rate, 6)).unwrap();
            if let Some(seed) = args.seed {
                write!(s, " --seed {}", seed).unwrap();
            }
        }
    }

    s.push_str(" -i ");
    if args.in_files.reads1.len() > 1 || args.in_files.alns.len() > 1 {
        s.push_str("...");
    } else if let Some(filename) = args.in_files.alns.first() {
        s.push_str(&ext::fmt::path(filename));
    } else {
        s.push_str(&ext::fmt::path(&args.in_files.reads1[0]));
        if let Some(filename) = args.in_files.reads2.first() {
            s.push(' ');
            s.push_str(&ext::fmt::path(filename));
        }
    }
    s.push_str(" | ");
    s
}

/// Mapping parameters will be stored at `analysis/bg/MAPPING_PARAMS`.
const MAPPING_PARAMS: &'static str = "mapping.json";

struct MappingParams {
    technology: Technology,
    min_mapq: u8,
    subsampling: f64,
    bg_region: String,
}

impl MappingParams {
    fn new(args: &Args, bg_region: Option<&Interval>) -> Self {
        MappingParams {
            technology: args.technology,
            min_mapq: args.min_mapq,
            subsampling: args.subsampling_rate,
            bg_region: bg_region.map(Interval::to_string).unwrap_or_else(|| "None".to_owned()),
        }
    }

    /// Compares mapping parameters against parameters from the previous run.
    /// Returns `Ok` if two parameters match, and there is no reason for remapping.
    /// Otherwise, returns `Err` with explanation.
    fn compare(&self, path: &Path) -> Result<(), String> {
        if !path.exists() {
            log::warn!("Cannot find previous mapping parameters at {}. Continuing anyway", ext::fmt::path(path));
            return Ok(())
        }
        let old = match ext::sys::load_json(&path).map_err(Error::from).and_then(|json| MappingParams::load(&json)) {
            Err(e) => {
                return Err(format!("Cannot load previous mapping parameters: {}", e.display()));
            }
            Ok(val) => val,
        };
        if self.technology != old.technology {
            Err(format!("Sequencing technologies do not match ({} -> {})", old.technology, self.technology))
        } else if (self.subsampling - old.subsampling).abs() > 1e-8 {
            Err(format!("Subsampling rate changed ({} -> {})", old.subsampling, self.subsampling))
        } else if self.min_mapq < old.min_mapq {
            Err(format!("Minimal mapping quality decreased ({} -> {})", old.min_mapq, self.min_mapq))
        } else if self.bg_region != old.bg_region {
            Err(format!("Background region has changed ({} -> {})", old.bg_region, self.bg_region))
        } else {
            Ok(())
        }
    }
}

impl JsonSer for MappingParams {
    fn save(&self) -> json::JsonValue {
        json::object!{
            technology: self.technology.to_str(),
            min_mapq: self.min_mapq,
            subsampling: self.subsampling,
            bg_region: &self.bg_region as &str,
        }
    }

    fn load(obj: &json::JsonValue) -> crate::Result<Self> {
        json_get!(obj => technology (as_str), min_mapq (as_u8), subsampling (as_f64), bg_region (as_str));
        let technology = Technology::from_str(technology).map_err(|e| Error::ParsingError(e))?;
        Ok(Self {
            technology, min_mapq, subsampling,
            bg_region: bg_region.to_owned(),
        })
    }
}

/// Returns true if read mapping can be skipped.
/// Panics if output BAM file exists, but cannot be used.
fn need_mapping(args: &Args, out_dir: &Path, out_bam: &Path, bg_region: Option<&Interval>) -> crate::Result<bool> {
    let params_path = out_dir.join(MAPPING_PARAMS);
    let curr_params = MappingParams::new(args, bg_region);
    if !out_bam.exists() {
        let mut params_file = ext::sys::create_file(&params_path)?;
        curr_params.save().write_pretty(&mut params_file, 4).map_err(add_path!(params_path))?;
        return Ok(true);
    }
    match curr_params.compare(&params_path) {
        Ok(()) => {
            log::warn!("BAM file {} exists, skipping read mapping", ext::fmt::path(out_bam));
            Ok(false)
        }
        Err(s) => {
            log::error!("Problem comparing read mapping parameters against {}", ext::fmt::path(&params_path));
            log::error!("    {}", s);
            log::error!("Please use `--rerun full`, or remove file {} to override this error",
                ext::fmt::path(&params_path));
            std::process::exit(1)
        }
    }
}

/// Split total number of threads into the recruitment and mapping threads.
fn split_n_threads(args: &Args) -> (u16, u16) {
    if args.skip_recruitment {
        let mapping_threads = max(1, args.threads);
        log::info!("Mapping reads to the reference in {} threads", mapping_threads);
        return (0, mapping_threads);
    } else  if args.threads <= 1 {
        return (1, 1);
    }

    let mut recr_threads = if args.recr_threads < 1.0 {
        (f64::from(args.threads) * args.recr_threads).round()
    } else {
        args.recr_threads.round()
    } as u16;
    recr_threads = recr_threads.clamp(1, args.threads - 1);
    let mapping_threads = max(1, args.threads - 1);
    log::info!("Recruiting and mapping reads to the reference in {} and {} threads, respectively",
        recr_threads, mapping_threads);
    (recr_threads, mapping_threads)
}

/// Either load sequencing info from old file, or recompute.
fn prepare_seq_info(
    args: &Args,
    need_mapping: bool,
    ref_contigs: &Arc<ContigNames>,
    out_dir: &Path,
) -> crate::Result<SequencingInfo>
{
    let old_path = out_dir.join(paths::BG_DISTR);
    if !need_mapping && old_path.exists() {
        if let Ok(old_distr) = BgDistr::load_from(&old_path, None) {
            return Ok(old_distr.seq_info().clone());
        }
    }
    let mut seq_info = SequencingInfo::new(read_len_from_reads(args, Some(ref_contigs))?,
        args.technology, args.explicit_technology)?;
    match args.in_files.sum_file_size() {
        Ok(s) => seq_info.set_file_size(s),
        Err(e) => log::warn!("Could not calculate file size, {:?}", e),
    }
    Ok(seq_info)
}

/// Map reads to the whole reference genome, and then take only reads mapped to the corresponding BED file.
/// Subsample reads if the corresponding rate is less than 1.
/// Return the total number of reads/read pairs.
fn run_mapping(
    args: &Args,
    seq_info: &mut SequencingInfo,
    ref_filename: &Path,
    ref_contigs: &Arc<ContigNames>,
    out_dir: &Path,
    out_bam: &Path,
    bg_region: Option<&BgRegion>,
) -> crate::Result<()>
{
    let (recr_threads, mapping_threads) = split_n_threads(args);
    let mut mapping = create_mapping_command(args, seq_info, ref_filename, mapping_threads);
    let mapping_exe = PathBuf::from(mapping.get_program().to_owned());
    let mut mapping_child = mapping.spawn().map_err(add_path!(mapping_exe))?;
    let mapping_stdin = mapping_child.stdin.take();
    let mapping_stdout = mapping_child.stdout.take();
    let mut pipe_guard = ext::sys::PipeGuard::new(mapping_exe, mapping_child);
    let handle = set_mapping_stdin(mapping_stdin.unwrap(), ref_contigs, seq_info, bg_region, recr_threads, args)?;

    let mut tmp_bed = None;
    if let Some(bg_data) = bg_region {
        let filename = out_dir.join("tmp.bed");
        let mut bed_file = ext::sys::create_file(&filename)?;
        writeln!(bed_file, "{}", bg_data.interval.bed_fmt()).map_err(add_path!(filename))?;
        std::mem::drop(bed_file);
        tmp_bed = Some(filename);
    }

    let tmp_bam = out_dir.join("tmp.bam");
    let mut samtools = Command::new(&args.samtools);
    // See SAM flags here: https://broadinstitute.github.io/picard/explain-flags.html.
    samtools.args(&["view",
            "-b", // Output BAM.
            // Ignore reads where any of the mates is unmapped,
            // + ignore secondary & supplementary alignments + ignore failed checks.
            "-F", "3852",
            "-q", &args.min_mapq.to_string(),
            ]);
    if let Some(filename) = &tmp_bed {
        samtools.arg("-L").arg(filename);
    }
    samtools.arg("-o").arg(&tmp_bam)
        .stdin(Stdio::from(mapping_stdout.unwrap()))
        .stdout(Stdio::piped()).stderr(Stdio::piped());
    log::debug!("    {}{} | {}", first_step_str(&args, bg_region),
        ext::fmt::command(&mapping), ext::fmt::command(&samtools));
    let samtools_child = samtools.spawn().map_err(add_path!(args.samtools))?;
    pipe_guard.push(args.samtools.clone(), samtools_child);
    pipe_guard.wait()?;
    if let Some(total_reads) = handle.join().map_err(|e| error!(RuntimeError, "Read mapping failed: {:?}", e))?? {
        seq_info.set_total_reads(total_reads);
    }

    fs::rename(&tmp_bam, out_bam).map_err(add_path!(tmp_bam, out_bam))?;
    if let Some(filename) = tmp_bed {
        fs::remove_file(&filename).map_err(add_path!(filename))?;
    }
    Ok(())
}

/// Loads records and converts them to Alignments.
/// All returned alignments are primary, have MAPQ over the `args.min_mapq`,
/// and clipping under the threshold `bg_params.max_clipping`.
///
/// Ignore reads, for which `cigar_getter` returns None.
///
/// Second returned argument is true if records are paired-end.
///
/// Alignment probabilities are not set.
fn load_alns(
    reader: &mut impl bam::Read,
    cigar_getter: impl Fn(&bam::Record) -> Option<Cigar>,
    contigs: &Arc<ContigNames>,
    args: &Args
) -> crate::Result<(Vec<NamedAlignment>, bool)>
{
    let min_mapq = args.min_mapq;
    let max_clipping = args.max_clipping;
    let mut paired_counts = [0_u64, 0];

    let mut alns = Vec::new();
    let mut record = bam::Record::new();
    let mut ignored = 0;
    let mut wo_cigar = 0;
    while let Some(()) = reader.read(&mut record).transpose()? {
        if record.flags() & 3844 == 0 && record.mapq() >= min_mapq && cigar::clipping_rate(&record) <= max_clipping {
            if let Some(cigar) = cigar_getter(&record) {
                alns.push(NamedAlignment::new(&record, cigar, ReadEnd::from_record(&record),
                    Arc::clone(contigs), f64::NAN));
                paired_counts[usize::from(record.is_paired())] += 1;
            } else {
                wo_cigar += 1;
            }
        } else {
            ignored += 1;
        }
    }
    if paired_counts[0] > 0 && paired_counts[1] > 0 {
        return Err(error!(InvalidData, "BAM file contains both paired and unpaired reads"));
    }
    if alns.is_empty() {
        return Err(error!(InvalidData, "BAM file contains no reads in the target region"));
    }
    log::debug!("    Loaded {} alignments, discarded {}", alns.len(), ignored);
    if wo_cigar > 0 {
        log::warn!("    Could not create extended CIGAR for {} records", wo_cigar);
    }
    Ok((alns, paired_counts[1] > 0))
}

/// Calculate mean read length from existing alignments.
fn read_len_from_alns(alns: &[NamedAlignment]) -> f64 {
    const N_RECORDS: usize = 10000;
    let n = min(alns.len(), N_RECORDS);
    alns[..n].iter()
        .map(|aln| f64::from(aln.cigar().query_len()))
        .sum::<f64>() / n as f64
}

/// Calculate mean read length from input reads.
fn read_len_from_reads(args: &Args, contigs: Option<&ContigNames>) -> crate::Result<f64> {
    const MIN_RECORDS: u64 = 200;
    const MIN_SUM_LEN: u64 = 3_000_000;
    log::info!("Calculating mean read length");
    fastx::process_readers!(args.in_files, contigs; let {mut} reader;
        { fastx::mean_read_len(&mut reader, MIN_RECORDS, MIN_SUM_LEN) })
}

fn get_windows(
    seq_info: &SequencingInfo,
    bg_region: Option<&BgRegion>,
    opt_out_dir: Option<&Path>,
    args: &Args,
) -> crate::Result<(OpCounter, Option<Windows>)>
{
    if let Some(bg_data) = bg_region {
        Ok((
            OpCounter::Bounded(bg_data.interval.clone()),
            Some(Windows::create(&bg_data.interval, bg_data.region_sequence(),
                &bg_data.kmer_counts, &seq_info, &args.bg_params.depth, opt_out_dir)?)
        ))
    } else {
        Ok((OpCounter::Unbounded, None))
    }
}

fn estimate_bg_from_paired(
    alns: Vec<NamedAlignment>,
    seq_info: SequencingInfo,
    opt_out_dir: Option<&Path>,
    bg_region: Option<&BgRegion>,
    args: &Args,
) -> crate::Result<BgDistr>
{
    // Group reads into pairs, and estimate insert size from them.
    let pair_ixs = insertsz::group_mates(&alns)?;
    let insert_distr = InsertDistr::estimate(&alns, &pair_ixs, opt_out_dir)?;
    let conf_lvl = 1.0 - args.bg_params.insert_pval;
    let (min_insert, max_insert) = insert_distr.confidence_interval(conf_lvl);
    log::info!("    Allowed insert size: [{}, {}]  ({}%-confidence interval)",
        min_insert, max_insert, crate::math::fmt_signif(100.0 * conf_lvl, 5));

    // Estimate error profile from read pairs with appropriate insert size.
    let mut errprof_alns = Vec::with_capacity(pair_ixs.len() * 2);
    for &(i, j) in pair_ixs.iter() {
        let first = &alns[i];
        let second = &alns[j];
        let insert = first.insert_size(second.deref());
        if min_insert <= insert && insert <= max_insert {
            errprof_alns.push(first);
            errprof_alns.push(second);
        }
    }

    let (op_counter, opt_windows) = get_windows(&seq_info, bg_region, opt_out_dir, args)?;
    let err_prof = ErrorProfile::estimate(&errprof_alns, &op_counter, opt_windows.as_ref(), opt_out_dir)?;
    let edit_dist_cache = SingleEditDistCache::new(&err_prof, args.bg_params.edit_pval);
    edit_dist_cache.print_log(seq_info.mean_read_len());

    // Estimate backgorund read depth from read pairs with both good probabilities and good insert size prob.
    let mut depth_alns: Vec<&Alignment> = Vec::with_capacity(errprof_alns.len());
    for chunk in errprof_alns.chunks_exact(2) {
        let aln1 = chunk[0];
        let aln2 = chunk[1];
        let dist1 = op_counter.count(aln1).edit_distance();
        let dist2 = op_counter.count(aln2).edit_distance();
        if dist1.edit() <= edit_dist_cache.get(dist1.read_len())
            && dist2.edit() <= edit_dist_cache.get(dist2.read_len())
        {
            depth_alns.push(aln1.deref());
            depth_alns.push(aln2.deref());
        }
    }
    let depth_distr = if let Some(windows) = opt_windows {
        ReadDepth::estimate(&depth_alns, &windows, &args.bg_params.depth,
            args.subsampling_rate, true, &seq_info, opt_out_dir)?
    } else {
        ReadDepth::empty()
    };
    Ok(BgDistr::new(seq_info, insert_distr, err_prof, depth_distr))
}

fn estimate_bg_from_unpaired(
    alns: Vec<NamedAlignment>,
    seq_info: SequencingInfo,
    opt_out_dir: Option<&Path>,
    bg_region: Option<&BgRegion>,
    args: &Args,
) -> crate::Result<BgDistr>
{
    let insert_distr = InsertDistr::undefined();
    let (op_counter, opt_windows) = get_windows(&seq_info, bg_region, opt_out_dir, args)?;
    let err_prof = ErrorProfile::estimate(&alns, &op_counter, opt_windows.as_ref(), opt_out_dir)?;
    let edit_dist_cache = SingleEditDistCache::new(&err_prof, args.bg_params.edit_pval);
    edit_dist_cache.print_log(seq_info.mean_read_len());

    let filt_alns: Vec<&Alignment> = alns.iter()
        .filter(|aln| {
            let dist = op_counter.count(aln).edit_distance();
            dist.edit() <= edit_dist_cache.get(dist.read_len())
        })
        .map(Deref::deref)
        .collect();

    let depth_distr = if let Some(windows) = opt_windows {
        ReadDepth::estimate(&filt_alns, &windows, &args.bg_params.depth,
            args.subsampling_rate, false, &seq_info, opt_out_dir)?
    } else {
        ReadDepth::empty()
    };
    Ok(BgDistr::new(seq_info, insert_distr, err_prof, depth_distr))
}

/// Estimate background distributions from input reads or existing alignments.
fn estimate_bg_distrs(
    args: &Args,
    out_dir: &Path,
    ref_contigs: &Arc<ContigNames>,
    bg_region: Option<&BgRegion>,
) -> crate::Result<BgDistr>
{
    let opt_out_dir = if args.debug { Some(out_dir) } else { None };
    let ref_filename = args.in_files.reference.as_ref().unwrap();
    let aln_is_paired_end: bool;
    let alns: Vec<NamedAlignment>;
    let mut seq_info: SequencingInfo;

    if args.in_files.has_indexed_alignment() {
        assert!(args.in_files.alns.len() == 1);
        let alns_filename = &args.in_files.alns[0];
        log::debug!("Loading mapped reads into memory ({})", ext::fmt::path(alns_filename));
        let mut bam_reader = bam::IndexedReader::from_path(alns_filename)?;
        fastx::set_reference(alns_filename, &mut bam_reader, &args.in_files.reference, Some(ref_contigs))?;

        let bg_data = bg_region.expect("Background region must be defined");
        let interval = &bg_data.interval;

        let padded_start = bg_data.padded_interval.start();
        let padded_seq = &bg_data.padded_sequence;
        bam_reader.fetch((interval.contig_name(), i64::from(interval.start()), i64::from(interval.end())))?;
        (alns, aln_is_paired_end) = load_alns(&mut bam_reader,
            |record| Cigar::infer_ext_cigar(record, padded_seq, padded_start),
            ref_contigs, args)?;
        if aln_is_paired_end && !args.technology.paired_end_allowed() {
            return Err(error!(InvalidInput, "Paired end reads are not supported by {}", args.technology.long_name()));
        }
        seq_info = SequencingInfo::new(read_len_from_alns(&alns), args.technology, args.explicit_technology)?;
    } else {
        let bam_filename = out_dir.join("aln.bam");
        let need_mapping = need_mapping(args, &out_dir, &bam_filename, bg_region.map(|data| &data.interval))?;
        seq_info = prepare_seq_info(args, need_mapping, ref_contigs, &out_dir)?;
        log::info!("Mean read length = {:.1}", seq_info.mean_read_len());
        if need_mapping {
            run_mapping(args, &mut seq_info, &ref_filename, ref_contigs, &out_dir, &bam_filename, bg_region)?;
        }

        log::debug!("Loading mapped reads into memory ({})", ext::fmt::path(&bam_filename));
        let mut bam_reader = bam::Reader::from_path(&bam_filename)?;
        let get_cigar = |record: &bam::Record| -> Option<Cigar> {
            let cigar = Cigar::from_raw(record.raw_cigar());
            assert!(!cigar.has_hard_clipping(), "Cannot process primary alignments with hard clipping");
            Some(cigar)
        };
        (alns, aln_is_paired_end) = load_alns(&mut bam_reader, get_cigar, ref_contigs, args)?;
        if aln_is_paired_end != args.in_files.is_paired_end() {
            return Err(error!(RuntimeError,
                "Input data is {}-end, while alignment file is {}-end. \
                Perhaps preprocessing was run multiple times with different inputs?",
                if args.in_files.is_paired_end() { "paired" } else { "single" },
                if aln_is_paired_end { "paired" } else { "single" },
            ));
        }
    }

    if aln_is_paired_end {
        estimate_bg_from_paired(alns, seq_info, opt_out_dir, bg_region, args)
    } else {
        estimate_bg_from_unpaired(alns, seq_info, opt_out_dir, bg_region, args)
    }
}

/// Calculates rate, by which read depth differs, using difference in the number of reads.
fn read_count_factor(
    args: &Args,
    seq_info: &SequencingInfo,
    similar_seq_info: &SequencingInfo,
) -> crate::Result<f64>
{
    let Some(similar_n_reads) = similar_seq_info.total_reads() else {
        return Err(error!(InvalidInput,
            "Cannot use similar dataset {}: total number of reads was not estimated",
            ext::fmt::path(args.similar_dataset.as_ref().unwrap())))
    };
    let mut total_reads = if !args.in_files.alns.is_empty() {
        log::info!("Counting reads in {}", ext::fmt::paths(&args.in_files.alns));
        args.in_files.alns.iter()
            .map(|filename| fastx::count_reads_bam(filename, &args.samtools, &args.in_files.reference, args.threads))
            .sum::<Result<u64, _>>()?
    } else {
        log::info!("Counting reads in {}", ext::fmt::paths(&args.in_files.reads1));
        args.in_files.reads1.iter()
            .map(|filename| fastx::count_reads_fastx(filename))
            .sum::<Result<u64, _>>()?
    };
    if args.in_files.interleaved {
        total_reads /= 2;
    }
    log::debug!("    Current dataset: {:.0}k reads, previous dataset: {:.0}k reads",
        1e-3 * total_reads as f64, 1e-3 * similar_n_reads as f64);
    // NOTE: total reads are deliberately not provided to `seq_info`, so as not to propagate errors.

    let reads_factor = total_reads as f64 / similar_n_reads as f64;
    let mut lengths_factor = seq_info.mean_read_len() / similar_seq_info.mean_read_len();
    if (lengths_factor - 1.0).abs() > 0.03 {
        log::debug!("    Read length differ by a factor of {:.5}", lengths_factor);
    } else {
        lengths_factor = 1.0;
    }
    Ok(reads_factor * lengths_factor)
}

/// Calculates rate, by which read depth differs, using difference in file sizes.
fn file_size_factor(
    args: &Args,
    seq_info: &SequencingInfo,
    similar_seq_info: &SequencingInfo,
) -> crate::Result<f64>
{
    log::warn!("Calculating read depth difference using file sizes.");
    log::warn!("    This may produce incorrect results.");
    log::warn!("    Please check that the files are very similar,");
    log::warn!("    including similar read names and identical file compression method.");
    let Some(prev_size) = similar_seq_info.file_size() else {
        return Err(error!(InvalidInput,
            "Cannot use similar dataset {}: file size was not estimated",
            ext::fmt::path(args.similar_dataset.as_ref().unwrap())))
    };
    let Some(file_size) = seq_info.file_size() else {
        return Err(error!(RuntimeError, "File size could not be calculated for unknown reason"));
    };
    log::debug!("    Current dataset: {:.0} MiB, previous dataset: {:.0} MiB",
        file_size as f64 / 1048576.0, prev_size as f64 / 1048576.0);
    Ok(file_size as f64 / prev_size as f64)
}

/// Assume that this WGS dataset is similar to another dataset, and use its parameters.
fn estimate_like(
    args: &Args,
    out_dir: &Path,
    like_dir: &Path,
    ref_contigs: &Arc<ContigNames>,
) -> crate::Result<BgDistr>
{
    let similar_path = like_dir.join(paths::BG_DISTR);
    log::info!("Loading distribution parameters from {}", ext::fmt::path(&similar_path));
    let similar_distr = BgDistr::load_from(&similar_path, Some(&like_dir.join(paths::SUCCESS)))?;
    let similar_seq_info = similar_distr.seq_info();

    let seq_info = prepare_seq_info(args, false, ref_contigs, out_dir)?;
    log::info!("Mean read length = {:.1}", seq_info.mean_read_len());
    if similar_seq_info.technology() != seq_info.technology() {
        return Err(error!(InvalidInput,
            "Cannot use similar dataset {}: different sequencing technology ({} and {})",
            ext::fmt::path(like_dir), similar_seq_info.technology(), seq_info.technology()));
    } else if similar_distr.insert_distr().is_paired_end() != args.in_files.is_paired_end() {
        return Err(error!(InvalidInput,
            "Cannot use similar dataset {}: paired-end status does not match", ext::fmt::path(like_dir)));
    } else if !seq_info.technology().is_read_len_similar(seq_info.mean_read_len(), similar_seq_info.mean_read_len()) {
        return Err(error!(InvalidInput,
            "Cannot use similar dataset {}: read lengths are different ({:.0} and {:.0})",
            ext::fmt::path(like_dir), seq_info.mean_read_len(), similar_seq_info.mean_read_len()));
    }

    let factor = if args.use_file_size {
        file_size_factor(args, &seq_info, &similar_seq_info)?
    } else {
        read_count_factor(args, &seq_info, &similar_seq_info)?
    };
    let mut new_distr = similar_distr.clone();
    new_distr.set_seq_info(seq_info);
    if factor < 0.1 {
        return Err(error!(InvalidInput,
            "Read depth changed too much (by a factor of {:.4}), please estimate parameters anew", factor))
    } else if (factor - 1.0).abs() < 0.01 {
        log::debug!("    Almost identical read depth");
    } else {
        log::debug!("    Update read depth by a factor of {:.4}", factor);
        new_distr.depth_mut().mul_depth(factor);
    }
    Ok(new_distr)
}

/// Information about the background region.
struct BgRegion {
    /// Background region.
    interval: Interval,
    /// Background region + padding.
    padded_interval: Interval,
    padded_sequence: Vec<u8>,

    /// k-mer counts on the padded sequence.
    padded_kmer_counts: KmerCounts,
    /// k-mer counts on the unpadded background region.
    kmer_counts: KmerCounts,
}

impl BgRegion {
    fn new(
        args: &Args,
        ref_filename: &Path,
        ref_contigs: &Arc<ContigNames>,
        ref_fasta: &mut fasta::IndexedReader<impl io::Read + io::Seek>,
    ) -> crate::Result<Self>
    {
        let interval = select_bg_interval(&ref_filename, &ref_contigs, &args.bg_region)?;
        const PADDING: u32 = 50_000;
        let padded_interval = interval.add_padding(PADDING);

        let kmer_getter = JfKmerGetter::new(args.jellyfish.clone(), args.jf_counts.clone().unwrap())?;
        let padded_sequence = padded_interval.fetch_seq(ref_fasta)?;

        log::info!("Calculating k-mer counts on the background region");
        let padded_kmer_counts = kmer_getter.fetch([padded_sequence.clone()])?;
        let kmer_counts = padded_kmer_counts.subregion(
            interval.start() - padded_interval.start(), interval.end() - padded_interval.start());
        Ok(Self { interval, padded_interval, padded_sequence, padded_kmer_counts, kmer_counts })
    }

    /// Returns region sequence without padding.
    pub fn region_sequence(&self) -> &[u8] {
        let i = (self.interval.start() - self.padded_interval.start()) as usize;
        let j = (self.interval.end() - self.padded_interval.start()) as usize;
        &self.padded_sequence[i..j]
    }
}

pub(super) fn run(argv: &[String]) -> crate::Result<()> {
    let mut args = parse_args(argv)?;
    args.in_files.fill_from_inlist()?;
    let args = args.validate()?;

    super::greet();
    let timer = Instant::now();
    let out_dir = args.output.as_ref().unwrap();
    if !args.rerun.prepare_dir(&out_dir)? {
        std::process::exit(0);
    }
    ext::sys::mkdir(out_dir)?;

    let ref_filename = args.in_files.reference.as_ref().unwrap();
    let (ref_contigs, mut ref_fasta) = ContigNames::load_indexed_fasta("ref", &ref_filename)?;
    let ref_contigs = Arc::new(ref_contigs);

    let bg_distr = if let Some(similar_dataset) = &args.similar_dataset {
        estimate_like(&args, &out_dir, similar_dataset, &ref_contigs)?
    } else {
        let bg_region = if args.debug_head.is_none() {
            Some(BgRegion::new(&args, &ref_filename, &ref_contigs, &mut ref_fasta)?)
        } else { None };
        estimate_bg_distrs(&args, &out_dir, &ref_contigs, bg_region.as_ref())?
    };

    let distr_filename = out_dir.join(paths::BG_DISTR);
    let mut distr_file = ext::sys::create_gzip(&distr_filename)?;
    bg_distr.save().write_pretty(&mut distr_file, 4).map_err(add_path!(distr_filename))?;
    log::info!("Success. Total time: {}", ext::fmt::Duration(timer.elapsed()));
    super::write_success_file(out_dir.join(paths::SUCCESS))?;
    Ok(())
}
