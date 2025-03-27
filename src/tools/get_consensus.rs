use std::collections::HashMap;
use std::fs::File;
use std::path::{PathBuf};
use std::process::exit;
use bio::io::fasta;
use clap::Parser;
use log::{info, warn};
use serde::{Deserialize, Serialize};
use anyhow::{Error, Result, Context, anyhow};
use bio::io::fasta::Record;
use colored::Colorize;
use nalgebra::DMatrix;

fn read_fasta(fasta_file: &PathBuf) -> Result<Vec<Vec<u8>>>{
    let mut reader = fasta::Reader::from_file(fasta_file).expect("Could not open provided FASTA file.");
    let mut seqs: Vec<Vec<u8>> = Vec::new();

    for result in reader.records() {
        let record = result.expect("This record is invalid and failed to parse.");
        seqs.push(record.seq().to_vec());
    }

    Ok(seqs)

}

fn sequences_to_matrix(sequences: &Vec<Vec<u8>>) -> Result<DMatrix<u8>> {
    // Check if sequences are empty
    if sequences.is_empty() {
        return Err(anyhow!("There are no sequences in the sequence vector passed to the sequence_to_matrix function."));
    }

    // Check that all sequences are the same length (this is an MSA)
    let mut count = 0;
    for seq in sequences{
        count = count+1;
        if seq.len() != sequences[0].len(){
            return Err(anyhow!("Not all sequences in the MSA have the same length. The length of the 1st seq is {} and the length of the {} seq is {}",
            sequences[0].len(), count, seq.len()))
        }

    }

    Ok(
        DMatrix::from_row_slice(sequences.len(), sequences[0].len(), &sequences.concat())
    )
}

pub fn run(input_seqs_aligned: &PathBuf, output_path: &PathBuf) -> Result<()>{
    let seqs = read_fasta(input_seqs_aligned)?;
    let seq_matrix = sequences_to_matrix(&seqs);

    println!("{:?}", seq_matrix);

    Ok(())

}