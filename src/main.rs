pub mod seq;
pub mod algo;
pub mod math;
pub mod bg;
pub mod model;
pub mod solvers;
pub mod ext;
pub mod command;

mod err;
// mod test;
pub use err::{Error, Result};

fn init_logger() {
    use fern::{
        Dispatch,
        colors::{Color, ColoredLevelConfig},
    };
    let colors = ColoredLevelConfig::default()
        .info(Color::Green)
        .debug(Color::Cyan);
    Dispatch::new()
        .format(move |out, message, record| {
            out.finish(format_args!(
                "[{} {:>5}] {}",
                chrono::Local::now().format("%H:%M:%S"),
                colors.color(record.level()),
                message
            ))
        })
        .level(log::LevelFilter::Debug)
        .level_for("highs", log::LevelFilter::Info)
        .chain(std::io::stderr())
        // .chain(fern::log_file("output.log")?)
        .apply()
        .unwrap();
}

fn print_kmers(descr: &str, v: &[(u32, u64)]) {
    println!("=== {} ({}) ===", descr, v.len());
    for (i, kmer) in v {
        println!("   {:3} {:16x}", i, kmer);
    }
}

fn test_minimizers() {
    use crate::seq::kmers::{self, Minimizer};

    // let seq = b"ACGTACGATCGACGTTAGCACAGTTGACTAGCGTATGCTAGTGCTGACCAGTTGTGTACACGTAACGTACGTAGCTGCGGGGATCGTATACGGCTAGAGGTATAGCTGATAGTCGATGCTGAACACGTTGTAGCTGATCGTAGCTACGACTAGATGTTGTGCATCGTATCGTACCGTAGCTAGCTGTAGCGAGCGAGCCATGTCGAT";
    let seq = b"GTTGTGTACACGTACACGTAACACNACACGTGACGTTACAGCTGTAGACGATCTGATGCTACC";
    let k = 10;
    let w = 5;

    const CANON: bool = false;

    let mut all = Vec::new();
    kmers::kmers::<u64, { CANON }>(seq, k, &mut all);
    let all: Vec<_> = all.into_iter().enumerate().map(|(i, kmer)| (i as u32, kmer.fast_hash())).collect();
    print_kmers("Minimizers (all)", &all);

    let mut out0 = Vec::<(u32, u64)>::new();
    kmers::naive_minimizers::<_, _, { CANON }>(seq, k, w, &mut out0);
    print_kmers("Minimizers (naïve)", &out0);

    let mut fast_all = Vec::<(u32, u64)>::new();
    kmers::minimizers::<_, _, { CANON }>(seq, k, 1, &mut fast_all);
    print_kmers("Minimizers (fast all)", &fast_all);

    let mut out1 = Vec::<(u32, u64)>::new();
    kmers::minimizers::<_, _, { CANON }>(seq, k, w, &mut out1);
    print_kmers("Minimizers (fast)", &out1);

    assert_eq!(out0, out1);

    // // let mut out0 = Vec::<(u32, u64)>::new();
    // // kmers::minimizers0::<_, _, false>(seq, k, w, &mut out0);
    // // print_kmers("Minimizers (old)", &out0);

    // let mut out1 = Vec::<(u32, u64)>::new();
    // kmers::minimizers::<_, _, false>(seq, k, w, &mut out1);
    // print_kmers("Minimizers (new)", &out1);
}

fn main() {
    test_minimizers();

    init_logger();
    // let args: Vec<_> = std::env::args().collect();
    // if let Err(e) = command::run(&args) {
    //     log::error!("Finished with an error:\n{}", e.display());
    //     std::process::exit(1);
    // }
}
