use anyhow::{anyhow, Context, Result};
use bio::io::fasta;
use colored::Colorize;
use log;
use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;
use std::process::exit;
use crate::utils::fasta_utils::{load_fasta, write_fasta_sequences};
use crate::utils::translate::GAP_CHAR;

type FastaRecords = HashMap<String, Vec<u8>>;
const VERSION: &str = "0.3.0";



pub fn reverse_translate(aa_seq: &Vec<u8>, nt_seq: &Vec<u8>) -> Result<Vec<u8>> {
    let gap_char = "-".as_bytes()[0];
    let mut new_nt_seq = Vec::with_capacity(aa_seq.len() * 3);

    let mut current_nt_idx = 0;

    for amino_acid in aa_seq.iter() {
        if amino_acid == &gap_char {
            new_nt_seq.extend_from_slice(&std::iter::repeat(gap_char).take(3).collect::<Vec<u8>>());
        } else {
            let to_idx = current_nt_idx + 3;

            if to_idx > nt_seq.len() {
                return Err(anyhow!("Failed to grab a codon from {} to {} on the nucleotide sequence. Index out of bounds.", current_nt_idx, to_idx));
            }

            new_nt_seq.extend_from_slice(&nt_seq[current_nt_idx..to_idx]);
            current_nt_idx += 3;
        }
    }

    Ok(new_nt_seq)
}

fn process_sequences(
    aa_sequences: FastaRecords,
    nt_sequences: FastaRecords,
) -> Result<FastaRecords> {

    let mut missing_seqs = 0;
    let mut translation_errors = 0;

    let mut reverse_translated_sequences: FastaRecords = FastaRecords::with_capacity(aa_sequences.capacity());

    for (sequence_id, aa_sequence) in aa_sequences {

        match nt_sequences.get(&sequence_id){
            None => {
                log::error!("The sequence with name {sequence_id} from the amino acid sequences could not be found in the nucleotide sequences");
                missing_seqs += 1;
            },
            Some(nt_sequence) => {
                let mut degapped_nt_seq = nt_sequence.clone();
                degapped_nt_seq.retain(|&base| base != GAP_CHAR);

                match reverse_translate(&aa_sequence, &degapped_nt_seq){
                    Err(e) => {
                        log::error!("Error in reverse-translating the read {}.\n{:?}", sequence_id, e);
                        translation_errors += 1;
                    },
                    Ok(reverse_translated_seq) =>{
                        reverse_translated_sequences.insert(sequence_id, reverse_translated_seq);
                    }
                }
            }
        }
    }

    log::info!("We had {:?} sequences present in the AA file but missing from the NT file.", missing_seqs);
    log::info!("We had {:?} reverse-translation errors.", translation_errors);

    Ok(reverse_translated_sequences)
}

fn load_name_file(name_mapping_path: &PathBuf) -> Result<HashMap<String, Vec<String>>> {
    let reader = std::io::BufReader::new(File::open(name_mapping_path)?);

    let mappings: HashMap<String, Vec<String>> = serde_json::from_reader(reader)?;

    Ok(mappings)
}

pub fn run(
    aa_filepath: &PathBuf,
    nt_filepath: &PathBuf,
    output_file_path: &PathBuf,
) -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;

    let mut amino_acid_sequences: FastaRecords = load_fasta(aa_filepath)?;
    let mut nuc_sequences: FastaRecords = load_fasta(nt_filepath)?;

    let rev_translated_seqs = process_sequences(amino_acid_sequences, nuc_sequences)
        .context("Error occurred while processing the sequences")?;

    write_fasta_sequences(output_file_path, &rev_translated_seqs)
        .with_context(|| format!("Error occurred while trying to write reverse translated sequences to {:?}", output_file_path))?;

    Ok(())
}
