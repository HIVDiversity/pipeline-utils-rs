use std::collections::HashMap;
use std::path::PathBuf;
use anyhow::{Context, Result};
use bio::io::fasta;

pub fn write_fasta_sequences(output_file: &PathBuf, sequences: &HashMap<String, Vec<u8>>) -> Result<()>{
    let mut writer = fasta::Writer::to_file(output_file).with_context(|| "Could not open output file")?;

    for (seq_id, seq) in sequences {
        writer.write(&seq_id, None, seq.as_slice())?;
    }

    Ok(())
}