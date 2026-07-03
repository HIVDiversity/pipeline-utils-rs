use crate::utils::codon_tables::GAP_CHAR;
use crate::utils::fasta_utils::{load_fasta, write_fasta_sequences, FastaRecords};
use anyhow::{bail, Result};
use colored::Colorize;

use itertools::Itertools;
use std::path::PathBuf;

fn transpose_sequences(sequences: Vec<Vec<u8>>) -> Result<Vec<Vec<u8>>> {
    let max_seq_length = match sequences.is_empty() {
        true => {
            bail!("No sequences were provided.")
        }
        false => match sequences.iter().map(|seq| seq.len()).all_equal() {
            true => sequences[0].len(),
            false => bail!("Sequence lengths are not equal."),
        },
    };

    let mut transposed_sequence_columns: Vec<Vec<u8>> = Vec::with_capacity(max_seq_length);
    for col_index in 0..max_seq_length {
        let this_col_vec: Vec<u8> = sequences
            .iter()
            .map(|seq| seq.get(col_index).unwrap_or(&(0u8)).clone())
            .collect();
        transposed_sequence_columns.push(this_col_vec);
    }

    Ok(transposed_sequence_columns)
}

pub(crate) fn strip_gap_columns(
    sequence_records: FastaRecords,
    pct_gap_cols_to_remove: usize,
) -> Result<FastaRecords> {
    let (seq_names, sequences): (Vec<String>, Vec<Vec<u8>>) = sequence_records.into_iter().unzip();
    let num_sequences = sequences.len();
    let transposed_sequences = transpose_sequences(sequences)?;

    let new_columns: Vec<Vec<u8>> = transposed_sequences
        .into_iter()
        .filter(|column| {
            let gap_count = column.iter().filter(|c| **c == GAP_CHAR).count();
            println!("{}", gap_count / num_sequences);
            (((gap_count as f32 / num_sequences as f32) * 100f32) as usize) < pct_gap_cols_to_remove
        })
        .collect();

    let final_sequences = transpose_sequences(new_columns)?;
    let output_sequences = seq_names
        .into_iter()
        .zip(final_sequences.into_iter())
        .collect();

    Ok(output_sequences)
}

pub fn run(input_file: &PathBuf, output_file: &PathBuf, gap_pct_to_remove: usize) -> Result<()> {
    log::info!(
        "{}",
        format!(
            "This is 'strip-gap-cols' version {}",
            env!("CARGO_PKG_VERSION")
        )
        .bold()
        .bright_yellow()
    );

    log::info!("Reading input file {:?}", input_file);
    let sequences = load_fasta(input_file)?;
    let stripped_sequences = strip_gap_columns(sequences, gap_pct_to_remove)?;

    write_fasta_sequences(output_file, &stripped_sequences)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use velcro::hash_map;

    #[test]
    fn basic_test() -> Result<()> {
        let input_seqs: FastaRecords = hash_map!(
            "Test A".to_string(): vec![b'A', b'T', b'-', b'G', b'C', b'C'],
            "Test B".to_string(): vec![b'A', b'T', b'-', b'G', b'-', b'-'],
            "Test C".to_string(): vec![b'A', b'T', b'-', b'G', b'-', b'-'],
            "Test D".to_string(): vec![b'A', b'T', b'-', b'G', b'C', b'-']
        );

        let expected_seqs: FastaRecords = hash_map!(
            "Test A".to_string(): vec![b'A', b'T', b'G', b'C', b'C'],
            "Test B".to_string(): vec![b'A', b'T', b'G', b'-', b'-'],
            "Test C".to_string(): vec![b'A', b'T', b'G', b'-', b'-'],
            "Test D".to_string(): vec![b'A', b'T', b'G', b'C', b'-']
        );

        let obtained_sequences = strip_gap_columns(input_seqs, 100);
        for (seq_name, seq) in obtained_sequences? {
            println!("{}: {}", seq_name, String::from_utf8(seq.clone())?);
            assert!(expected_seqs.get(&seq_name).unwrap().eq(&seq));
        }

        Ok(())
    }

    #[test]
    fn test_transpose_sequences() -> Result<()> {
        let input_seqs: Vec<Vec<u8>> = vec![
            vec![b'A', b'T', b'-', b'C'],
            vec![b'A', b'T', b'-', b'T'],
            vec![b'A', b'T', b'-', b'G'],
        ];

        let expected_seqs: Vec<Vec<u8>> = vec![
            vec![b'A', b'A', b'A'],
            vec![b'T', b'T', b'T'],
            vec![b'-', b'-', b'-'],
            vec![b'C', b'T', b'G'],
        ];

        let obtained_sequences = transpose_sequences(input_seqs)?;

        expected_seqs
            .into_iter()
            .zip(obtained_sequences.into_iter())
            .for_each(|(expected, actual)| {
                assert!(expected.eq(&actual));
            });

        Ok(())
    }

    #[test]
    fn test_unequal_sequences() {
        let input_seqs: FastaRecords = hash_map!(
            "Test A".to_string(): vec![b'A', b'T', b'-', b'G', b'C', b'C'],
            "Test B".to_string(): vec![b'A', b'T', b'-', b'G'],
        );

        assert!(strip_gap_columns(input_seqs, 100).is_err())
    }
}
