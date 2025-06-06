use anyhow::Result;
use colored::Colorize;
use std::path::PathBuf;

const VERSION: &str = "0.1.0";

pub fn run(genbank_file: &PathBuf, output_file: &PathBuf, sequence_name: &String) -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;

    log::info!(
        "{}",
        format!("This is {} version {}", "gb-extract".italic(), VERSION)
            .bold()
            .bright_purple()
    );

    Ok(())
}
