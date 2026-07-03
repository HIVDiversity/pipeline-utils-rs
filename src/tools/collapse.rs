use crate::utils::codon_tables::GAP_CHAR;
use crate::utils::fasta_utils::{load_fasta, write_fasta_sequences, FastaRecords};
use anyhow::Result;
use colored::Colorize;
use std::collections::HashMap;
use std::path::PathBuf;

pub(crate) type SeqToNameMapping = HashMap<Vec<u8>, Vec<String>>;

pub(crate) fn collapse_sequences(
    sequences: FastaRecords,
    strip_gaps: bool,
) -> Result<SeqToNameMapping> {
    let mut unique_sequences: SeqToNameMapping =
        SeqToNameMapping::with_capacity(sequences.capacity());

    for fasta_record in sequences {
        let record_id = fasta_record.0;
        let mut record_seq = fasta_record.1;

        if strip_gaps {
            record_seq.retain(|&val| val != GAP_CHAR);
        }

        unique_sequences
            .entry(record_seq)
            .and_modify(|seq_name_vec| seq_name_vec.push(record_id.to_owned()))
            .or_insert(vec![record_id.to_owned()]);
    }

    Ok(unique_sequences)
}

pub(crate) fn build_collapsed_output(
    collapsed_seqs: SeqToNameMapping,
    seq_prefix: &str,
) -> (FastaRecords, HashMap<String, Vec<String>>) {
    let mut collapsed_sequences: FastaRecords = FastaRecords::with_capacity(collapsed_seqs.len());
    let mut name_mapping: HashMap<String, Vec<String>> =
        HashMap::with_capacity(collapsed_seqs.len());

    let mut counter = 0;
    for (sequence, sequence_names) in collapsed_seqs {
        // This will generate a sequence with a unique int for each collapsed seq, and a count
        // for the sequences that make up this collapsed one
        let seq_name = format!(
            "{}_{:0>4}_{:0>4}",
            seq_prefix,
            counter,
            sequence_names.len()
        );

        collapsed_sequences.insert(seq_name.clone(), sequence);
        counter += 1;
        name_mapping.insert(seq_name, sequence_names);
    }

    (collapsed_sequences, name_mapping)
}

fn write_sequences_and_name_mapping(
    collapsed_seqs: SeqToNameMapping,
    output_file: &PathBuf,
    name_mapping_output: &PathBuf,
    seq_prefix: &String,
) -> Result<()> {
    let (collapsed_sequences, name_mapping) = build_collapsed_output(collapsed_seqs, seq_prefix);

    log::info!("Writing unique sequences to file {:?}", output_file);
    write_fasta_sequences(output_file, &collapsed_sequences)?;

    log::info!("Writing name mapping to {:?}", name_mapping_output);
    std::fs::write(
        name_mapping_output,
        serde_json::to_string(&name_mapping).expect("Error serializing the name map."),
    )
    .expect("Error with writing the name map to the disk.");
    Ok(())
}

pub fn run(
    input_file: &PathBuf,
    output_file: &PathBuf,
    namefile_output: &PathBuf,
    seq_name_prefix: &String,
    strip_gaps: bool,
) -> Result<()> {
    log::info!(
        "{}",
        format!("This is 'collapse' version {}", env!("CARGO_PKG_VERSION"))
            .bold()
            .bright_yellow()
    );

    log::info!("Reading input file {:?}", input_file);
    let sequences = load_fasta(input_file)?;
    let collapsed_seqs = collapse_sequences(sequences, strip_gaps)?;

    write_sequences_and_name_mapping(
        collapsed_seqs,
        output_file,
        namefile_output,
        seq_name_prefix,
    )?;

    Ok(())
}
