use crate::utils::fasta_utils::{FastaRecords, load_fasta, write_fasta_sequences};
use crate::utils::translate::translate;
use anyhow::Result;
use colored::Colorize;
use std::path::PathBuf;

// TODO: Replace the NT in the sequence name with AA when translating....
const VERSION: &str = "0.2.0";
fn translate_fasta_file(
    sequences: &FastaRecords,
    strip_gaps: bool,
    ignore_gap_codons: bool,
    drop_incomplete_codons: bool,
    aa_stop_char: Option<char>,
) -> Result<FastaRecords> {
    let mut translated_sequences: FastaRecords = FastaRecords::with_capacity(sequences.capacity());

    for sequence in sequences {
        let translated_seq = translate(
            sequence.1.as_slice(),
            strip_gaps,
            ignore_gap_codons,
            drop_incomplete_codons,
            aa_stop_char,
        )?;
        translated_sequences.insert(sequence.0.to_string(), translated_seq);
    }

    Ok(translated_sequences)
}

pub fn run(
    nt_filepath: &PathBuf,
    output_filepath: &PathBuf,
    strip_gaps: bool,
    ignore_gap_codons: bool,
    drop_incomplete_codons: bool,
    aa_stop_char: Option<char>,
) -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;

    log::info!(
        "{}",
        format!("This is {} version {}", "translate".italic(), VERSION)
            .bold()
            .bright_purple()
    );

    log::info!("Reading sequences from {:?}", nt_filepath);
    let nucleotide_sequences = load_fasta(nt_filepath)?;

    log::info!("Translating sequences.");
    let translated_sequences = translate_fasta_file(
        &nucleotide_sequences,
        strip_gaps,
        ignore_gap_codons,
        drop_incomplete_codons,
        aa_stop_char,
    )?;
    log::info!("Done. Writing sequences to {:?}", output_filepath);

    write_fasta_sequences(output_filepath, &translated_sequences)?;

    log::info!("Done. Exiting.");
    Ok(())
}
