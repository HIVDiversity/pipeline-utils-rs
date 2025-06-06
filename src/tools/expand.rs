use crate::utils::fasta_utils::{FastaRecords, load_fasta};
use anyhow::{Context, Result};
use bio::io::fasta;
use colored::Colorize;
use serde_json::from_reader;
use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;
const VERSION: &str = "0.1.1";
type NewToOldNameMapping = HashMap<String, Vec<String>>;

fn uncollapse_and_write_sequences(
    collapsed_seqs: FastaRecords,
    name_mapping: NewToOldNameMapping,
    output_file: &PathBuf,
) -> Result<()> {
    let mut writer = fasta::Writer::to_file(output_file)
        .with_context(|| format!("Trying to write to file {:?}", output_file))?;

    for (collapsed_seq_name, sequence) in collapsed_seqs {
        match name_mapping.get(&collapsed_seq_name) {
            None => log::warn!(
                "The sequence with new name {:?} did not have a corresponding entry in the name mapping",
                &collapsed_seq_name
            ),
            Some(old_seq_names) => {
                for old_seq_name in old_seq_names {
                    writer
                        .write(old_seq_name, None, &sequence)
                        .with_context(|| {
                            format!(
                                "Trying to write sequence {:?} to {:?}",
                                old_seq_name, output_file
                            )
                        })?
                }
            }
        }
    }

    Ok(())
}

pub fn run(input_file: &PathBuf, name_mapping_file: &PathBuf, output_file: &PathBuf) -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;
    log::info!(
        "{}",
        format!("This is {} version {}", "expand".italic(), VERSION)
            .bold()
            .bright_magenta()
    );

    let collapsed_sequences = load_fasta(input_file)
        .with_context(|| format!("Failed to read sequences from {:?}", input_file))?;

    let name_mapping: NewToOldNameMapping = from_reader(File::open(name_mapping_file)?)
        .with_context(|| format!("Failed to read name mapping from {:?}", name_mapping_file))?;

    uncollapse_and_write_sequences(collapsed_sequences, name_mapping, output_file)?;

    Ok(())
}
