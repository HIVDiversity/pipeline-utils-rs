use crate::utils::codon_tables::AMBIGUOUS_NT_LOOKUP;
use crate::utils::fasta_utils::{load_fasta, write_fasta_sequences, FastaRecords};
use anyhow::Context;
use colored::Colorize;
use itertools::Itertools;
use std::path::PathBuf;

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

pub fn replace_ambiguities_records(
    sequences: FastaRecords,
    seed: u64,
) -> anyhow::Result<FastaRecords> {
    let mut rng = oorandom::Rand32::new(seed);
    let mut new_sequences: FastaRecords = FastaRecords::with_capacity(sequences.capacity());

    // Iterate in a deterministic order (HashMap order is randomized per-process) so the
    // seeded RNG stream is applied to sequences in the same order on every run.
    for seq_id in sequences.keys().sorted().cloned().collect::<Vec<_>>() {
        let sequence = &sequences[&seq_id];
        let new_seq = replace_ambiguities(sequence, &mut rng)?;
        new_sequences.insert(seq_id, new_seq);
    }

    Ok(new_sequences)
}

pub fn run(input_filepath: &PathBuf, output_filepath: &PathBuf, seed: u64) -> anyhow::Result<()> {
    log::info!(
        "{}",
        format!(
            "This is {} version {}",
            "replace-ambiguities".italic(),
            env!("CARGO_PKG_VERSION")
        )
        .bold()
        .bright_purple()
    );
    log::info!("Command was run with a random seed = {}", seed);

    log::info!(
        "Reading sequences from {:?} and writing to {:?}.",
        input_filepath,
        output_filepath
    );

    let sequences = load_fasta(input_filepath).context("Could not open input file.")?;
    let new_sequences = replace_ambiguities_records(sequences, seed)?;
    write_fasta_sequences(output_filepath, &new_sequences)?;

    log::info!("Done. Exiting.");
    Ok(())
}
