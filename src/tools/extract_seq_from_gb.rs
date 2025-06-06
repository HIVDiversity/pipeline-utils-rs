use anyhow::{Context, Result, anyhow};
use bio::io::fasta;
use colored::Colorize;
use gb_io::reader::parse_file;
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

    log::info!("Reading file {:?}", genbank_file);
    let genbank_contents = parse_file(genbank_file).context("Error parsing genbank file")?;

    // Complex series of steps here.
    // Iterate through the genbank features, looking to see which ones has a feature with the "note"
    // parameter. If it has a note param, then check if the value of that param is set.
    // If the param is set, then check if its value is equal to the name of the sequence we want
    let seq_of_interest = genbank_contents
        .get(0)
        .expect("Genbank file was empty")
        .features
        .to_owned()
        .into_iter()
        .find(|feature| {
            if let Some(note_feature) = feature
                .clone()
                .qualifiers
                .iter()
                .find(|qualifier| qualifier.0 == "note")
            {
                if let Some(note_name) = note_feature.1.as_ref() {
                    note_name == sequence_name
                } else {
                    false
                }
            } else {
                false
            }
        }).with_context(|| anyhow!("We were not able to find a feature in the genbank file that had a 'note' field which matched {}", sequence_name.bold()))?;

    log::debug!("Found sequence of interest! Extracting nucleotide sequence");

    let nt_seq = match seq_of_interest.location.clone().find_bounds() {
        Ok(bounds) => {
            let from_idx = bounds.0 as usize;
            let to_idx = bounds.1 as usize;
            genbank_contents[0].seq[from_idx..to_idx].to_vec()
        }
        Err(e) => {
            anyhow::bail!(
                "Got an error trying to get the bounds of the sequence of interest: {:?}",
                e.to_string()
            );
        }
    };
    log::info!("Successfully extracted nucleotide sequence from main reference.");
    let output_record =
        fasta::Record::with_attrs(sequence_name, None, nt_seq.to_ascii_uppercase().as_slice());

    log::info!("Writing record to {:?}", output_file);
    fasta::Writer::to_file(output_file)
        .with_context(|| anyhow!("Failed to write to file {:?}", output_file))?
        .write_record(&output_record)
        .with_context(|| {
            anyhow!(
                "Could not write record {:?} to file {:?}",
                output_record,
                output_file
            )
        })?;

    Ok(())
}
