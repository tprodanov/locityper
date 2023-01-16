use simplelog::{TermLogger, ConfigBuilder, LevelPadding, LevelFilter, TerminalMode, ColorChoice};

pub mod seq;
pub mod algo;
pub mod reconstr;
pub mod bg;

mod test;

fn init_logger() {
    TermLogger::init(
        log::LevelFilter::Trace,
        ConfigBuilder::new()
            .set_level_padding(LevelPadding::Right)
            .set_thread_level(LevelFilter::Off)
            .set_target_level(LevelFilter::Off)
            .build(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    ).unwrap();
}

fn main() {
    init_logger();
    test::test();
}
