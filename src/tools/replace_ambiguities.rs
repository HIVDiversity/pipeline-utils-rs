use crate::utils::translate::AMBIGUOUS_NT_LOOKUP;
use anyhow::Context;
use bio::io::fasta;
use colored::Colorize;
use std::path::PathBuf;
const VERSION: &str = "1.0.0";

fn replace_ambiguities(sequence: &[u8], rng: &mut oorandom::Rand32) -> anyhow::Result<Vec<u8>> {
    let new_sequence: Vec<u8> = sequence
        .iter()
        .cloned()
        .map(|nt| {
            return if AMBIGUOUS_NT_LOOKUP.contains_key(&[nt]) {
                let possible_nts = &AMBIGUOUS_NT_LOOKUP[&[nt]];
                let index = rng.rand_range(0..possible_nts.len() as u32) as usize;
                possible_nts
                    .iter()
                    .nth(index)
                    .with_context(|| format!("Failed to get nucleotide for nt {:?}", nt))
                    .unwrap_or(&&[nt])[0]
            } else {
                nt
            };
        })
        .collect();

    Ok(new_sequence)
}

pub fn run(input_filepath: &PathBuf, output_filepath: &PathBuf, seed: u64) -> anyhow::Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;

    log::info!(
        "{}",
        format!(
            "This is {} version {}",
            "replace-ambiguities".italic(),
            VERSION
        )
        .bold()
        .bright_purple()
    );
    log::info!("Command was run with a random seed = {}", seed);

    let reader = fasta::Reader::from_file(input_filepath).expect("Could not open input file.");
    let mut rng = oorandom::Rand32::new(seed);

    let mut writer =
        fasta::Writer::to_file(output_filepath).with_context(|| "Could not open output file")?;

    log::info!(
        "Reading sequences from {:?} and writing to {:?}.",
        input_filepath,
        output_filepath
    );

    for record_result in reader.records() {
        match record_result {
            Ok(record) => {
                let new_seq = replace_ambiguities(record.seq(), &mut rng)?;
                writer.write(record.id(), None, new_seq.as_slice())?;
            }
            Err(_) => {
                log::error!("Failed to read record from file.");
            }
        }
    }

    log::info!("Done. Exiting.");
    Ok(())
}
