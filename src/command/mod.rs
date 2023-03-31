mod common;
mod create;
mod add;

use colored::Colorize;

use crate::Error;
use common::{find_exe, print_version, file_or_stdin, file_or_stdout};

fn print_citation() {
    print_version();
    println!();
    println!("{}", "Thank you for using our tool!".bold());
    println!("Publication in progress, please check later.");
}

fn print_help() {
    print_version();
    println!("\n{} {} command [arguments]",
        "Usage:".bold(), env!("CARGO_PKG_NAME"));

    println!("\n{}", "[ Creating database ]".bold());
    println!("    {:<7}  Create an empty database.", "create".red());
    println!("    {:<7}  Add complex locus/loci to the database.", "add".red());

    println!("\n{}", "[ General help ]".bold());
    println!("    {:<7}  Show this help message.", "help".red());
    println!("    {:<7}  Show version.", "version".red());
    println!("    {:<7}  Show citation information.", "cite".red());
}

pub fn run(argv: &[String]) -> Result<(), Error> {
    if argv.len() <= 1 {
        print_help();
        std::process::exit(1);
    }
    match &argv[1] as &str {
        "create" => create::run(&argv[2..])?,
        "add" => add::run(&argv[2..])?,
        "help" | "h" | "--help" | "-h" => print_help(),
        "version" | "--version" | "-V" => print_version(),
        "cite" => print_citation(),
        cmd => panic!("Unknown command {}", cmd),
    }
    Ok(())
}
