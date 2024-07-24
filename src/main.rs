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

fn main() {
    init_logger();
    let args: Vec<_> = std::env::args().collect();
    if let Err(e) = command::run(&args) {
        log::error!("Finished with an error:\n{}", e.display());
        std::process::exit(1);
    }
}
