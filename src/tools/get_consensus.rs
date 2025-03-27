use std::collections::HashMap;
use std::path::{PathBuf};
use bio::io::fasta;
use anyhow::{Result, Context, anyhow};
use nalgebra::DMatrix;
use colored::Colorize;
const VERSION: &str = "0.1.0";

fn read_fasta(fasta_file: &PathBuf) -> Result<Vec<Vec<u8>>>{
    let reader = fasta::Reader::from_file(fasta_file).expect("Could not open provided FASTA file.");
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

fn build_consensus(msa: &DMatrix<u8>) -> Result<Vec<u8>>{
    let mut consensus: Vec<u8> = Vec::new();

    for col in msa.column_iter(){
        let mut col_count = HashMap::new();

        for item in col{
            *col_count.entry(item).or_insert(0) += 1;
        }

        let largest_item = col_count
            .iter()
            .max_by(|a, b| a.1.cmp(&b.1))
            .map(|(k,    _v)| k).ok_or(anyhow!("Could not get the most frequent element in this column."))?;

        consensus.push(largest_item.to_owned().to_owned());

    }
    Ok(consensus)
}

pub fn run(input_seqs_aligned: &PathBuf, output_path: &PathBuf) -> Result<()>{
    simple_logger::SimpleLogger::new().env().init()?;
    log::info!("{}" ,format!("This is get-consensus version {}", VERSION).bold().bright_green());

    log::info!("Reading input FASTA file: {:?}", input_seqs_aligned);
    let seqs = read_fasta(input_seqs_aligned)?;
    log::info!("Successfully read {} sequences into memory.", seqs.len());

    log::info!("Converting sequences into a matrix.");
    let seq_matrix = sequences_to_matrix(&seqs)?;
    log::info!("Successfully created a {} by {} matrix of sequences.", seq_matrix.nrows(), seq_matrix.ncols());

    log::info!("Generating consensus.");
    let consensus = build_consensus(&seq_matrix)?;

    let cons_string = String::from_utf8(consensus).with_context(||"Consensus vector is not a valid UTF-8 string.")?;
    log::info!("Done. Consensus:\n{cons_string}");


    Ok(())

}