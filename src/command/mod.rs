pub mod paths;
mod add;
mod preproc;
mod genotype;
#[cfg(feature = "align")]
mod align;
mod recruit;

use std::{
    fs,
    fmt::{self, Write as FmtWrite},
    path::Path,
    str::FromStr,
};
use colored::Colorize;
use crate::{
    ext,
    err::{Error, add_path},
};

pub const PROGRAM: &'static str = env!("CARGO_PKG_NAME");
pub const VERSION: &'static str = env!("CARGO_PKG_VERSION");

pub fn run(argv: &[String]) -> Result<(), Error> {
    if argv.len() <= 1 {
        print_help();
        std::process::exit(1);
    }
    match &argv[1] as &str {
        "a" | "add" => add::run(&argv[2..])?,
        "p" | "preproc" | "preprocess" => preproc::run(&argv[2..])?,
        "g" | "genotype" => genotype::run(&argv[2..])?,
        "r" | "recruit" => recruit::run(&argv[2..])?,

        "aln" | "align" =>
        {
            #[cfg(feature = "align")] {
                align::run(&argv[2..])?;
            }
            #[cfg(not(feature = "align"))] {
                log::error!("{} command unavailable, please enable `align` feature during compilation", "align".bold());
                std::process::exit(1);
            }
        }

        "h" | "help" | "--help" | "-h" => print_help(),
        "V" | "version" | "--version" | "-V" => print_version(),
        "cite" => print_citation(),
        cmd => {
            log::error!("Unknown command {}", cmd.red());
            std::process::exit(1);
        }
    }
    Ok(())
}

fn print_version() {
    println!("{} {}", PROGRAM.underline(), format!("v{}", VERSION).green());
    let authors: Vec<_> = env!("CARGO_PKG_AUTHORS").split(':').collect();
    let n = authors.len();
    if n == 0 {
        return;
    }
    print!("Created by ");
    for (i, author) in authors.iter().enumerate() {
        if i == 0 {
            print!("{}", author.bright_blue());
        } else if i < n - 1 {
            print!(", {}", author.bright_blue());
        } else {
            print!(" and {}", author.bright_blue());
        }
    }
    println!();
}

fn print_citation() {
    println!("{}", "Thank you for using our tool!".bold());
    println!();

    const PREPRINT_LINK: &'static str = "https://doi.org/10.1101/2024.05.03.592358";
    println!("You can find preprint here:");
    println!("    T. Prodanov, E. Plender, G. Seebohm, S. Meuth, E. Eichler, T. Marschall. \
        Locityper: targeted");
    println!("    genotyping of complex polymorphic genes. {} (2024), {}",
        "bioRxiv".italic(), PREPRINT_LINK);
    println!();
    println!("Please check for peer-reviewed publication later!");
}

fn print_help() {
    const WIDTH: usize = 12;
    print_version();
    println!("\n{} {} command [arguments]",
        "Usage:".bold(), PROGRAM);

    println!("\n{}", "[ Loci database ]".bold());
    println!("    {:WIDTH$} Add target locus/loci to the database",
        "a, add".red());

    println!("\n{}", "[ Analysing WGS data ]".bold());
    println!("    {:WIDTH$} Preprocess WGS dataset",
        "p, preproc".red());
    println!("    {:WIDTH$} Genotype target loci",
        "g, genotype".red());

    println!("\n{}", "[ Other utilities ]".bold());
    println!("    {:WIDTH$} Align medium-size sequences to each other",
        "   align".red());
    println!("    {:WIDTH$} Recruit reads to one/multiple loci",
        "r, recruit".red());

    println!("\n{}", "[ General help ]".bold());
    println!("    {:WIDTH$} Show this help message",
        "h, help".red());
    println!("    {:WIDTH$} Show version",
        "V, version".red());
    println!("    {:WIDTH$} Show citation information",
        "   cite".red());
}

fn flag() -> colored::ColoredString {
    "░░░".yellow().dimmed()
}

/// Format default value for the help message.
fn fmt_def(val: impl fmt::Display) -> colored::ColoredString {
    val.to_string().cyan()
}

fn fmt_def_f64(val: f64) -> colored::ColoredString {
    crate::math::fmt_signif(val, 6).cyan()
}

/// Re-runing mode.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Rerun {
    /// Complete rerun analysis.
    All,
    /// Use some existing evaluations.
    Part,
    /// Skip successfully complete analyses.
    None,
}

impl Rerun {
    pub fn to_str(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::Part => "part",
            Self::None => "none",
        }
    }

    /// Converts boolean flag `force` into either `All` or `None` rerun mode.
    fn from_force(force: bool) -> Self {
        if force { Self::All } else { Self::None }
    }

    /// Creates/cleans output directory and returns true if further analysis is needed.
    /// Depending on the rerun mode,
    ///     if none:
    ///         if `success_file` exists, does nothing and returns false.
    ///         otherwise: return true.
    ///     if part:  always returns true, and
    ///         if `success_file` exists, removes it.
    ///     if all: always returns true and always clears output directory.
    ///
    /// If directory already existed, run `clean(dir)`.
    fn prepare_and_clean_dir<F>(self, dir: &Path, clean: F) -> Result<bool, Error>
    where F: FnOnce(&Path) -> Result<(), Error>,
    {
        if !dir.exists() {
            ext::sys::mkdir(dir)?;
            return Ok(true);
        }

        if self == Self::All {
            log::warn!("Clearing directory {}", ext::fmt::path(&dir));
            fs::remove_dir_all(&dir).map_err(add_path!(dir))?;
            ext::sys::mkdir(dir)?;
            return Ok(true);
        }

        let mut need_rerun = true;
        let success_file = dir.join(paths::SUCCESS);
        if success_file.exists() {
            if self == Self::None {
                log::info!("Skipping directory {} (successfully completed)", ext::fmt::path(&dir));
                need_rerun = false;
            } else {
                fs::remove_file(&success_file).map_err(add_path!(success_file))?;
            }
        }
        if need_rerun {
            clean(dir)?;
        }
        Ok(need_rerun)
    }

    fn prepare_dir(self, dir: &Path) -> Result<bool, Error> {
        self.prepare_and_clean_dir(dir, |_| Ok(()))
    }
}

impl FromStr for Rerun {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match &s.to_lowercase() as &str {
            "all" | "full" => Ok(Self::All),
            "part" | "partial" => Ok(Self::Part),
            "none" | "no" => Ok(Self::None),
            _ => Err(format!("Unknown rerun mode {:?} (allowed modes: all, part, none)", s)),
        }
    }
}

impl fmt::Display for Rerun {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(self.to_str())
    }
}

impl fmt::Debug for Rerun {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

/// Debug information at the start of the log.
fn greet() {
    let command = std::env::args().map(ext::fmt::path).collect::<Vec<_>>().join(" ");
    log::debug!("{}", command);
    log::debug!("{} v{} @ {}", PROGRAM, VERSION, chrono::Local::now().format("%Y-%m-%d %H:%M:%S"));
    #[cfg(debug_assertions)] {
        log::warn!("Debug assertions enabled - faster optimization is available");
    }
}

/// How many debug CSV files to save.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DebugLvl {
    None,
    Some,
    Full,
}

impl From<u8> for DebugLvl {
    fn from(val: u8) -> Self {
        match val {
            0 => Self::None,
            1 => Self::Some,
            2 => Self::Full,
            _ => {
                log::warn!("Debug level is too big ({}), setting it to 2", val);
                Self::Full
            }
        }
    }
}

// /// Writes command line arguments, current version and time.
// fn write_command(filename: impl AsRef<Path>) -> Result<(), Error> {
//     let mut s = String::new();
//     for (i, arg) in std::env::args().enumerate() {
//         if i > 0 {
//             s.push(' ');
//         }
//         s.push_str(&ext::fmt::path(&arg));
//     }
//     write!(s, "\n{} v{}\n", PROGRAM, VERSION).unwrap();
//     writeln!(s, "{:?}", chrono::offset::Local::now()).unwrap();
//     fs::write(&filename, s).map_err(add_path!(filename))
// }

fn write_success_file(filename: impl AsRef<Path>) -> Result<(), Error> {
    fs::write(&filename, const_format::formatcp!("v{}\n", VERSION)).map_err(add_path!(filename))
}

/// Describe values for each technology in format
/// "illumina: X, hifi: Y, ..."
pub fn describe_defaults<T, V>(
    options: impl IntoIterator<Item = T>,
    fmt_option: impl Fn(&T) -> String,
    get_default: impl Fn(&T) -> V,
) -> String
where V: fmt::Display + PartialEq,
{
    // Vector of pairs (Default value, all options where it appears).
    let mut tech_vals: Vec<(V, String)> = Vec::new();
    for opt in options {
        let val = get_default(&opt);
        if let Some(i) = tech_vals.iter().position(|(v, _)| v == &val) {
            write!(tech_vals[i].1, ",{}", fmt_option(&opt)).unwrap();
        } else {
            tech_vals.push((val, fmt_option(&opt)));
        }
    }
    tech_vals.into_iter()
        .map(|(val, opts)| format!("{}: {}", opts, fmt_def(val)))
        .collect::<Vec<String>>()
        .join("; ")
}
