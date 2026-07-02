use crate::utils::fasta_utils::{FastaRecords, load_fasta, write_fasta_sequences};
use anyhow::{Context, Result};
use colored::Colorize;
use serde_json::from_reader;
use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;
type NewToOldNameMapping = HashMap<String, Vec<String>>;

pub fn uncollapse_sequences(
    collapsed_seqs: FastaRecords,
    name_mapping: NewToOldNameMapping,
    include_missing_seqs: bool,
) -> Result<FastaRecords> {
    let mut expanded_seqs: FastaRecords = FastaRecords::with_capacity(collapsed_seqs.len());

    for (collapsed_seq_name, sequence) in collapsed_seqs {
        match name_mapping.get(&collapsed_seq_name) {
            None => {
                log::warn!(
                    "The sequence with new name {:?} did not have a corresponding entry in the name mapping",
                    &collapsed_seq_name
                );
                if include_missing_seqs {
                    expanded_seqs.insert(collapsed_seq_name, sequence);
                }
            }
            Some(old_seq_names) => {
                for old_seq_name in old_seq_names {
                    expanded_seqs.insert(old_seq_name.clone(), sequence.clone());
                }
            }
        }
    }

    Ok(expanded_seqs)
}

pub fn run(
    input_file: &PathBuf,
    name_mapping_file: &PathBuf,
    output_file: &PathBuf,
    include_missing_seqs: bool,
) -> Result<()> {
    log::info!(
        "{}",
        format!("This is {} version {}", "expand".italic(), env!("CARGO_PKG_VERSION"))
            .bold()
            .bright_magenta()
    );

    let collapsed_sequences = load_fasta(input_file)
        .with_context(|| format!("Failed to read sequences from {:?}", input_file))?;

    let name_mapping: NewToOldNameMapping = from_reader(File::open(name_mapping_file)?)
        .with_context(|| format!("Failed to read name mapping from {:?}", name_mapping_file))?;

    let expanded_sequences =
        uncollapse_sequences(collapsed_sequences, name_mapping, include_missing_seqs)?;

    write_fasta_sequences(output_file, &expanded_sequences)?;

    Ok(())
}
