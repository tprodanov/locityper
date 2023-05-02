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
pub use err::Error;

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
        .level(log::LevelFilter::Trace)
        .level_for("highs", log::LevelFilter::Info)
        .chain(std::io::stderr())
        // .chain(fern::log_file("output.log")?)
        .apply()
        .unwrap();
}

fn main() {
    init_logger();
    // let args: Vec<_> = std::env::args().collect();
    // command::run(&args).unwrap();

    use seq::wfa;
    let aligner = wfa::Aligner::new(&wfa::Penalties::default());
    println!("{:?}", aligner.align(b"ACGTAGAC", b"ACGTAGACA"));
}
