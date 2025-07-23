use crate::utils::fasta_utils::{FastaRecords, write_fasta_sequences};
use anyhow::{Context, Result};

use bio::bio_types::sequence::SequenceRead;
use colored::Colorize;
use log::warn;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::{bam, bam::Read, bam::Record};
use std::collections::HashMap;
use std::path::PathBuf;

const VERSION: &str = "1.0.0";

fn find_read_pos_from_ref_pos(read: &Record, ref_pos: i64) -> Option<i64> {
    for pair in read.aligned_pairs_full() {
        let current_query_pos = pair[0];
        let current_ref_pos = pair[1];
        if current_ref_pos.is_some_and(|x| x >= ref_pos) {
            if current_query_pos.is_some() {
                return current_query_pos;
            }
        }
    }
    None
}

pub fn run(
    input_file: &PathBuf,
    output_file: &PathBuf,
    trim_from: i64,
    trim_to: i64,
) -> Result<()> {
    // Set up logging with the desired log level
    simple_logger::SimpleLogger::new().env().init()?;

    // Print information about this program
    log::info!(
        "{}",
        format!("This is trim_sam version {}", VERSION)
            .bold()
            .bright_green()
    );

    let mut reader = bam::Reader::from_path(input_file)?;

    let mut output_seqs: FastaRecords = HashMap::new();

    for record in reader.records() {
        let record = record?;

        // We have to subtract 1 from the user-provided idx since those are base 1 and hts-lib works
        // in base 0. We then have to add 1 to the trim_to_seq value since the user provides us with
        // the last base they want INCLUDED
        let trim_from_seq =
            find_read_pos_from_ref_pos(&record, trim_from - 1).unwrap_or_else(|| {
                warn!("Failed to convert the read pos");
                return 0;
            }) as usize;
        let mut trim_to_seq = (find_read_pos_from_ref_pos(&record, trim_to - 1)
            .unwrap_or(record.len() as i64)
            + 1) as usize;

        if trim_to_seq + 1 > record.len() {
            trim_to_seq = record.len();
        }

        // We have to add 1 to the trim_to_seq value since the user provides us with the last base
        // the want INCLUDED
        output_seqs.insert(
            String::from_utf8(record.name().to_vec())?,
            record.seq().as_bytes()[trim_from_seq..trim_to_seq].to_vec(),
        );
    }

    write_fasta_sequences(output_file, &output_seqs)
        .with_context(|| format!("Failed to write output file {:?}", output_file))?;

    Ok(())
}
