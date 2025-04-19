use anyhow::{Context, Result};
use bio::io::fasta;
use std::collections::HashMap;
use std::path::PathBuf;

pub type FastaRecords = HashMap<String, Vec<u8>>;

#[derive(Clone, Copy)]
pub enum SequenceType{
    Nucleotide,
    AminoAcid
}
pub fn write_fasta_sequences(
    output_file: &PathBuf,
    sequences: &HashMap<String, Vec<u8>>,
) -> Result<()> {
    let mut writer =
        fasta::Writer::to_file(output_file).with_context(|| "Could not open output file")?;

    for (seq_id, seq) in sequences {
        writer.write(&seq_id, None, seq.as_slice())?;
    }

    Ok(())
}

// TODO: move to a public function
pub fn load_fasta(file_path: &PathBuf) -> Result<FastaRecords> {
    let mut sequences: FastaRecords = FastaRecords::new();
    let reader = fasta::Reader::from_file(file_path).expect("Could not open file.");

    // let mut parsing_errors = 0;

    for result in reader.records() {
        let record = result.expect("This record is invalid and failed to parse.");
        let mut seq = record.seq().to_vec();
        seq.make_ascii_uppercase();
        sequences.insert(record.id().to_string(), seq);
    }

    Ok(sequences)
}
