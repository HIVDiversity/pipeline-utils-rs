use crate::utils;
use anyhow::{Result, anyhow};
use bio::io::fasta;
use clap::ValueEnum;
use colored::Colorize;
use itertools::Itertools;
use nalgebra::DMatrix;
use rand::seq::IteratorRandom;
use std::collections::HashMap;
use std::path::PathBuf;
use utils::fasta_utils;
use utils::translate::find_ambiguity_code;

const VERSION: &str = "0.2.1";

#[derive(ValueEnum, Clone, Copy)]
pub enum AmbiguityMode {
    UseIUPAC,
    First,
    Random,
    MarkN,
}

fn sequences_to_matrix(sequences: &Vec<Vec<u8>>) -> Result<DMatrix<u8>> {
    // Check if sequences are empty
    if sequences.is_empty() {
        return Err(anyhow!(
            "There are no sequences in the sequence vector passed to the sequence_to_matrix function."
        ));
    }

    // Check that all sequences are the same length (this is an MSA)
    let mut count = 0;
    for seq in sequences {
        count = count + 1;
        if seq.len() != sequences[0].len() {
            return Err(anyhow!(
                "Not all sequences in the MSA have the same length. The length of the 1st seq is {} and the length of the {} seq is {}",
                sequences[0].len(),
                count,
                seq.len()
            ));
        }
    }

    Ok(DMatrix::from_row_slice(
        sequences.len(),
        sequences[0].len(),
        &sequences.concat(),
    ))
}

fn build_consensus(msa: &DMatrix<u8>, ambiguity_mode: AmbiguityMode) -> Result<Vec<u8>> {
    let mut consensus: Vec<u8> = Vec::new();

    for col in msa.column_iter() {
        let mut col_count = HashMap::new();

        for item in col {
            *col_count.entry(item).or_insert(0) += 1;
        }

        // Attempt to get the item in the column with the largest count, or if there
        // are multiple then get the set.
        let largest_items: Vec<&u8> = col_count
            .iter()
            .max_set_by(|a, b| a.1.cmp(&b.1))
            .iter()
            .cloned()
            .map(|(k, _v)| *k)
            .collect();

        if largest_items.len() == 1 {
            consensus.push(*largest_items[0]);
        } else {
            match ambiguity_mode {
                AmbiguityMode::UseIUPAC => {
                    let ambiguity_code = find_ambiguity_code(&largest_items);
                    match ambiguity_code {
                        None => {
                            return Err(anyhow!(
                                "A nucleotide set doesn't have an ambiguity code."
                            ));
                        }
                        Some(code) => {
                            consensus.push(code[0]);
                        }
                    }
                }
                AmbiguityMode::First => {
                    let first_item = largest_items
                        .iter()
                        .sorted()
                        .map(|x| **x)
                        .collect::<Vec<u8>>()
                        .first()
                        .unwrap()
                        .to_owned();

                    consensus.push(first_item);
                }
                AmbiguityMode::Random => {
                    let random_item = largest_items.iter().choose(&mut rand::rng()).unwrap();
                    consensus.push(**random_item);
                }
                AmbiguityMode::MarkN => {
                    consensus.push(b'N');
                }
            }
        }
    }

    Ok(consensus)
}

fn write_consensus(output_file: &PathBuf, seq_name: &str, seq: &Vec<u8>) -> Result<()> {
    let mut writer = fasta::Writer::to_file(output_file)?;
    let mut degapped_seq = seq.clone();
    let gap_char = b'-';
    degapped_seq.retain(|&val| val != gap_char);
    writer.write(seq_name, None, &degapped_seq)?;

    Ok(())
}

pub fn run(
    input_seqs_aligned: &PathBuf,
    output_path: &PathBuf,
    consensus_name: &String,
    ambiguity_mode: AmbiguityMode,
) -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;
    log::info!(
        "{}",
        format!("This is get-consensus version {}", VERSION)
            .bold()
            .bright_green()
    );

    log::info!("Reading input FASTA file: {:?}", input_seqs_aligned);
    let seqs_map = fasta_utils::load_fasta(input_seqs_aligned)?;
    let seqs: Vec<Vec<u8>> = seqs_map.into_iter().map(|(_, seq)| seq).collect();

    log::info!("Successfully read {} sequences into memory.", seqs.len());

    let seq_matrix = sequences_to_matrix(&seqs)?;
    log::info!(
        "Successfully created a {} by {} matrix of sequences.",
        seq_matrix.nrows(),
        seq_matrix.ncols()
    );

    log::info!("Generating consensus.");
    let consensus = build_consensus(&seq_matrix, ambiguity_mode)?;

    log::info!("Writing consensus to {:?}", output_path);
    write_consensus(output_path, consensus_name, &consensus)?;

    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_ambiguities() {
        let input: Vec<Vec<u8>> = vec![vec![b'T', b'T', b'G'], vec![b'A', b'T', b'G']];
        let matrix = sequences_to_matrix(&input).unwrap();
        let consensus_iupac = build_consensus(&matrix, AmbiguityMode::UseIUPAC).unwrap();
        let consensus_first = build_consensus(&matrix, AmbiguityMode::First).unwrap();
        let consensus_markn = build_consensus(&matrix, AmbiguityMode::MarkN).unwrap();

        assert_eq!(
            String::from("WTG"),
            String::from_utf8(consensus_iupac).unwrap()
        );

        assert_eq!(
            String::from("NTG"),
            String::from_utf8(consensus_markn).unwrap()
        );

        assert_eq!(
            String::from("ATG"),
            String::from_utf8(consensus_first).unwrap()
        );
    }
}
