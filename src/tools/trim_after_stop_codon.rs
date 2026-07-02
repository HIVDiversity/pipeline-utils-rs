use crate::utils::fasta_utils::{FastaRecords, load_fasta, write_fasta_sequences};
use crate::utils::codon_tables::STOP_CODONS;
use anyhow::Result;
use colored::Colorize;
use std::collections::HashMap;
use std::path::PathBuf;

fn trim_sequence(sequence: &Vec<u8>, include_stop_codon: bool) -> Result<Vec<u8>> {
    let first_stop_codon_index = sequence
        .chunks(3)
        .into_iter()
        .position(|codon| STOP_CODONS.contains(<&[u8; 3]>::try_from(codon).unwrap_or(&[0; 3])));

    match first_stop_codon_index {
        None => { Ok(sequence.clone()) }
        Some(index) => {
            let trim_index = if include_stop_codon {
                (index + 1) * 3
            } else {
                index * 3
            };
            Ok(sequence.clone()[..trim_index].to_vec())
        }
    }
}

pub(crate) fn process_file(
    sequences: FastaRecords,
    include_stop_codon: bool,
) -> Result<FastaRecords> {
    let mut output_sequences = HashMap::<String, Vec<u8>>::with_capacity(sequences.len());

    for (seq_name, sequence) in sequences {
        let trimmed_sequence = trim_sequence(&sequence, include_stop_codon)?;
        output_sequences.insert(seq_name, trimmed_sequence);
    }

    Ok(output_sequences)
}

pub fn run(input_file: &PathBuf, output_file: &PathBuf, include_stop_codon: bool) -> Result<()> {
    log::info!(
        "{}",
        format!("This is 'trim_after_stop_codon' version {}", env!("CARGO_PKG_VERSION"))
            .bold()
            .bright_yellow()
    );

    log::info!("Reading input file {:?}", input_file);
    let sequences = load_fasta(input_file)?;
    let trimmed_sequences = process_file(sequences, include_stop_codon)?;

    write_fasta_sequences(output_file, &trimmed_sequences)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_no_stop_codon() {
        let input = b"CATAATGTTTGGGCTACACATGCCTGTGTA".to_vec();
        let expected = b"CATAATGTTTGGGCTACACATGCCTGTGTA".to_vec();

        let obtained = trim_sequence(&input, false).unwrap();
        assert_eq!(String::from_utf8(obtained).unwrap(), String::from_utf8(expected).unwrap());
    }

    #[test]
    fn test_trim_normal_include_stop_codon() {
        let input = b"GCTTTGCTATAAGATAGGATTATTGAAGTAGCACAAAGAATC".to_vec();
        let expected = b"GCTTTGCTATAA".to_vec();

        let obtained = trim_sequence(&input, true).unwrap();
        assert_eq!(String::from_utf8(obtained).unwrap(), String::from_utf8(expected).unwrap());
    }

    #[test]
    fn test_trim_at_start_include_stop_codon() {
        let input = b"TAAGCTTTGCTA".to_vec();
        let expected = b"TAA".to_vec();

        let obtained = trim_sequence(&input, true).unwrap();
        assert_eq!(String::from_utf8(obtained).unwrap(), String::from_utf8(expected).unwrap());
    }

    #[test]
    fn test_trim_seq_not_div_3() {
        let input = b"GAAGGAACAGATAGGATTATTGA".to_vec();
        let expected = b"GAAGGAACAGATAGGATTATTGA".to_vec();

        let obtained = trim_sequence(&input, true).unwrap();
        assert_eq!(String::from_utf8(obtained).unwrap(), String::from_utf8(expected).unwrap());
    }

    #[test]
    fn test_trim_normal_exclude_stop_codon() {
        let input = b"GCTTTGCTATAAGATAGGATTATTGAAGTAGCACAAAGAATC".to_vec();
        let expected = b"GCTTTGCTA".to_vec();

        let obtained = trim_sequence(&input, false).unwrap();
        assert_eq!(String::from_utf8(obtained).unwrap(), String::from_utf8(expected).unwrap());
    }

    #[test]
    fn test_trim_at_start_exclude_stop_codon() {
        let input = b"TAAGCTTTGCTA".to_vec();
        let expected = b"".to_vec();

        let obtained = trim_sequence(&input, false).unwrap();
        assert_eq!(String::from_utf8(obtained).unwrap(), String::from_utf8(expected).unwrap());
    }
}
