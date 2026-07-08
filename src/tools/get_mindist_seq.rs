use crate::tools::get_consensus::{AmbiguityMode, build_consensus, sequences_to_matrix};
use crate::utils::codon_tables::GAP_CHAR;
use crate::utils::fasta_utils::{FastaRecords, load_fasta, write_fasta_sequences};
use anyhow::{Result, bail};
use clap::ValueEnum;
use colored::Colorize;
use std::path::PathBuf;

#[derive(ValueEnum, Clone, Copy)]
pub enum ComputeMode {
    Exact,
    Heuristic,
}

fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

/// Given an MSA (name -> aligned sequence, all same length), returns the name
/// of the sequence most representative of the alignment: the real sequence
/// closest (by Hamming distance) to the per-column consensus.
///
/// Returns `None` if the map is empty.
///
/// # Panics
/// Panics if sequences aren't all the same length.
pub fn get_most_representative_sequence(
    msa: &FastaRecords,
    ambiguity_mode: AmbiguityMode,
    compute_mode: ComputeMode,
) -> Result<String> {
    assert!(
        msa.len() > 1,
        "There needs to be 2 or more sequences provided."
    );

    let seq_len = match msa.values().next() {
        Some(seq) => seq.len(),
        None => {
            bail!("No sequences have been provided.")
        }
    };

    assert!(
        msa.values().all(|s| s.len() == seq_len),
        "all sequences in the MSA must have the same length"
    );

    let msa_seqs: Vec<Vec<u8>> = msa.values().cloned().collect();
    let msa_matrix = sequences_to_matrix(&msa_seqs)?;
    let consensus = build_consensus(&msa_matrix, ambiguity_mode)?;

    let computed_seq_name = match compute_mode {
        ComputeMode::Exact => msa
            .iter()
            .min_by_key(|(_, seq)| hamming_distance(seq, &consensus))
            .map(|(name, _)| name.clone()),

        ComputeMode::Heuristic => msa
            .iter()
            .min_by_key(|(_, seq_i)| {
                msa.values()
                    .map(|seq_j| hamming_distance(seq_i, seq_j))
                    .sum::<usize>()
            })
            .map(|(name, _)| name.clone()),
    };

    let seq_name = match computed_seq_name {
        Some(name) => name,
        None => unreachable!(
            "This is returned if the msa iterator is empty, but this is not possible since we check for it earlier."
        ),
    };

    Ok(seq_name)
}

pub fn run(
    input_file: &PathBuf,
    output_file: &PathBuf,
    ambiguity_mode: AmbiguityMode,
    compute_mode: ComputeMode,
) -> anyhow::Result<()> {
    log::info!(
        "{}",
        format!(
            "This is 'get-mindist-seq' version {}",
            env!("CARGO_PKG_VERSION")
        )
        .bold()
        .bright_magenta()
    );

    log::info!("Reading input file {:?}", input_file);
    let sequences = load_fasta(input_file)?;
    let representative_seq_name =
        get_most_representative_sequence(&sequences, ambiguity_mode, compute_mode)?;
    log::info!("Most representative sequence: {}", representative_seq_name);

    let mut representative_seq = sequences[&representative_seq_name].clone();
    representative_seq.retain(|&base| base != GAP_CHAR);

    let output_sequences: FastaRecords =
        FastaRecords::from([(representative_seq_name, representative_seq)]);
    write_fasta_sequences(output_file, &output_sequences)?;

    Ok(())
}
