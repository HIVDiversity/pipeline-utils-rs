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

fn read_fasta(fasta_file: &PathBuf) -> Result<()>{
    let mut reader = fasta::Reader::from_file(fasta_file).expect("Could not open provided FASTA file.");
    let mut seqs: Vec<Vec<u8>> = Vec::new();

    for result in reader.records() {
        let record = result.expect("This record is invalid and failed to parse.");
        seqs.push(record.seq().to_vec());
    }

    for seq in seqs {
        for a in seq {
            print!("{}", a)
        }
        println!()
    }

    Ok(())

}

pub fn run(input_seqs_aligned: &PathBuf, output_path: &PathBuf) -> Result<()>{
    read_fasta(input_seqs_aligned)
}