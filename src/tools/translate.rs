use crate::utils::fasta_utils::{FastaRecords, load_fasta, write_fasta_sequences};
use crate::utils::translate::{TranslationOptions, translate};
use anyhow::Result;
use colored::Colorize;
use std::path::PathBuf;

const VERSION: &str = "1.0.0";

pub fn run(
    nt_filepath: &PathBuf,
    output_filepath: &PathBuf,
    translation_options: &TranslationOptions,
) -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;

    log::info!(
        "{}",
        format!("This is {} version {}", "translate".italic(), VERSION)
            .bold()
            .bright_purple()
    );
    log::info!(
        "Command was run with the following options:\n{}",
        translation_options
    );

    log::info!("Reading sequences from {:?}", nt_filepath);
    let nucleotide_sequences = load_fasta(nt_filepath)?;

    log::info!("Translating sequences.");
    let mut translated_sequences: FastaRecords =
        FastaRecords::with_capacity(nucleotide_sequences.capacity());

    for sequence in nucleotide_sequences {
        let translated_seq = translate(sequence.1.as_slice(), translation_options)?;
        translated_sequences.insert(sequence.0.to_string(), translated_seq);
    }

    log::info!("Done. Writing sequences to {:?}", output_filepath);

    write_fasta_sequences(output_filepath, &translated_sequences)?;

    log::info!("Done. Exiting.");
    Ok(())
}
