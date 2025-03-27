use std::collections::HashMap;
use std::fs::File;
use std::path::{PathBuf};
use std::process::exit;
use bio::io::fasta;
use clap::Parser;
use log::{info, warn};
use serde::{Deserialize, Serialize};
use anyhow::{Error, Result, Context, anyhow};
use colored::Colorize;

type FastaRecords = HashMap<String, Vec<u8>>;
const VERSION: &str = "0.2.1";


fn load_fasta(file_path: &PathBuf, sequence_store: &mut FastaRecords) {
    let reader = fasta::Reader::from_file(file_path)
        .expect("Could not open file.");

    // let mut parsing_errors = 0;

    for result in reader.records() {
        let record = result.expect("This record is invalid and failed to parse.");
        sequence_store.insert(record.id().to_string(), record.seq().to_vec());
    }

    // if parsing_errors > 0 {
    //     log::warn!("There were {parsing_errors} records that couldn't be parsed.")
    // }
}

fn reverse_translate(aa_seq: &Vec<u8>, nt_seq: &Vec<u8>) -> Result<Vec<u8>> {
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

fn process_sequences(aa_sequences: &FastaRecords,
                     nt_sequences: &FastaRecords,
                     aligned_nt_seqs: &mut FastaRecords,
                     name_mapping: &HashMap<String, Vec<String>>) -> Result<()> {
    let mut missing_seqs = 0;
    let gap_char = "-".as_bytes()[0];
    for (sequence_id, aa_sequence) in aa_sequences {
        if let Some(nt_seq_names) = name_mapping.get(sequence_id) {
            for nt_seq_name in nt_seq_names {
                if let Some(nt_seq) = nt_sequences.get(nt_seq_name) {
                    let mut degapped_nt_seq = nt_seq.clone();
                    degapped_nt_seq.retain(|&val| val != gap_char);

                    let new_seq = reverse_translate(&aa_sequence, &degapped_nt_seq).with_context(|| format!("Error in reverse-translating the read {}", nt_seq_name))?;

                    if aa_sequence.len() * 3 != new_seq.len() {
                        warn!("For seq {sequence_id} -> {nt_seq_name}, the length of the amino acid sequence is different to the translated sequence")
                    } else {
                        aligned_nt_seqs.insert(nt_seq_name.clone(), new_seq);
                    }
                } else {
                    warn!("The sequence with ID {} was not found in the nt file, but it was present in the mapping", &nt_seq_name.bold().cyan());
                    missing_seqs += 1;
                }
            }
        } else {
            warn!("The mapping for the sequence ID {} was not found in the name file, but was present in the AA file", &sequence_id.bold().cyan());
            missing_seqs += 1;
        }
    }
    info!("There were {} sequences which were missing either from the NT file or the mapping.", missing_seqs);
    Ok(())
}

fn load_name_file(name_mapping_path: &PathBuf) -> Result<HashMap<String, Vec<String>>> {
    let reader = std::io::BufReader::new(File::open(name_mapping_path)?);

    let mappings: HashMap<String, Vec<String>> = serde_json::from_reader(reader)?;

    Ok(mappings)
}


pub fn run(aa_filepath: &PathBuf, nt_filepath: &PathBuf, name_mapping: &PathBuf, output_file_path: &PathBuf, check_keys_match: bool) -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;

    let name_mappings = load_name_file(name_mapping)?;

    let mut amino_acid_sequences: FastaRecords = FastaRecords::new();
    let mut nuc_sequences: FastaRecords = FastaRecords::new();
    let mut new_nuc_sequences: FastaRecords = FastaRecords::new();

    load_fasta(&aa_filepath, &mut amino_acid_sequences);
    load_fasta(&nt_filepath, &mut nuc_sequences);
    let mut keys_match: bool = false;

    if check_keys_match {
        keys_match = (amino_acid_sequences.len() == nuc_sequences.len()) && (amino_acid_sequences.keys().all(|key| nuc_sequences.contains_key(key)));
    }

    if check_keys_match && !keys_match {
        log::error!("The keys in the Amino Acid sequence file do not match the keys in the Nucleotide file");
        exit(1);
    } else {
        process_sequences(&amino_acid_sequences, &nuc_sequences, &mut new_nuc_sequences, &name_mappings)?;
    }


    let mut writer = fasta::Writer::to_file(output_file_path)
        .expect("Could not open file for writing.");

    for (seq_id, seq) in new_nuc_sequences {
        writer.write(&seq_id, None, seq.as_slice())?;
    }

    Ok(())
}

