mod create;
mod add;
mod preproc;

use colored::Colorize;

use crate::Error;

fn print_version() {
    println!("{} {}", env!("CARGO_PKG_NAME").underline(), format!("v{}", env!("CARGO_PKG_VERSION")).green());
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
    print_version();
    println!();
    println!("{}", "Thank you for using our tool!".bold());
    println!("Publication in progress, please check later.");
}

fn print_help() {
    const WIDTH: usize = 8;
    print_version();
    println!("\n{} {} command [arguments]",
        "Usage:".bold(), env!("CARGO_PKG_NAME"));

    println!("\n{}", "[ Creating database ]".bold());
    println!("    {:WIDTH$}  Create an empty database.", "create".red());
    println!("    {:WIDTH$}  Add complex locus/loci to the database.", "add".red());

    println!("\n{}", "[ Analysing WGS data ]".bold());
    println!("    {:WIDTH$}  Preprocess WGS dataset.", "preproc".red());

    println!("\n{}", "[ General help ]".bold());
    println!("    {:WIDTH$}  Show this help message.", "help".red());
    println!("    {:WIDTH$}  Show version.", "version".red());
    println!("    {:WIDTH$}  Show citation information.", "cite".red());
}

fn flag() -> impl std::fmt::Display {
    "░░░".yellow().dimmed()
}

pub fn run(argv: &[String]) -> Result<(), Error> {
    if argv.len() <= 1 {
        print_help();
        std::process::exit(1);
    }
    match &argv[1] as &str {
        "create" | "c" => create::run(&argv[2..])?,
        "add" | "a" => add::run(&argv[2..])?,
        "preproc" | "preprocess" | "p" => preproc::run(&argv[2..])?,

        "help" | "h" | "--help" | "-h" => print_help(),
        "version" | "V" | "--version" | "-V" => print_version(),
        "cite" => print_citation(),
        cmd => panic!("Unknown command {}", cmd),
    }
    Ok(())
}
