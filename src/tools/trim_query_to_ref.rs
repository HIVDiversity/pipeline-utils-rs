use anyhow::{Context, Result};
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use bio::io::fasta;
use std::collections::HashMap;
use std::iter::Iterator;
use std::path::PathBuf;
use colored::Colorize;
use crate::utils::translate::translate;

const VERSION: &str = "0.2.2";


// TODO: Move readfasta to the utils crate
fn read_fasta(fasta_file: &PathBuf) -> Result<Vec<Vec<u8>>> {
    let reader = fasta::Reader::from_file(fasta_file).expect("Could not open provided FASTA file.");
    let mut seqs: Vec<Vec<u8>> = Vec::new();

    for result in reader.records() {
        let record = result.expect("This record is invalid and failed to parse.");
        seqs.push(record.seq().to_vec());
    }

    Ok(seqs)
}

fn write_fasta(output_file: &PathBuf, seq_name: &str, seq: &Vec<u8>) -> Result<()> {
    let mut writer = fasta::Writer::to_file(output_file)?;
    writer.write(seq_name, None, seq)?;

    Ok(())
}

pub fn run(
    reference_file: &PathBuf,
    query_file: &PathBuf,
    output_file: &PathBuf,
    output_seq_name: &str,
    strip_gaps: bool,
    output_type: &String,
) -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;

    log::info!(
        "{}",
        format!("This is pairwise-align-to-ref version {}", VERSION)
            .bold()
            .bright_green()
    );

    let reference_read = read_fasta(reference_file)?;
    let reference = reference_read[0].as_slice();

    let query_read = read_fasta(query_file)?;
    let query = query_read[0].as_slice();

    let ref_aa = translate(reference, true, false, true)
        .context("Error in translating the reference string to amino acids")?;

    let ref_aa_slice = ref_aa.as_slice();

    let mut aligner =
        Aligner::with_capacity(query.len() / 3, ref_aa.len(), -5, -1, bio::scores::blosum62);

    let mut best_score = 0;
    let mut best_frame = 0;
    let mut best_translation: Vec<u8> = Vec::with_capacity(query.len() / 3);
    let mut best_alignment: Alignment = Default::default();

    for frame in 0..3 {
        log::info!("Translating query in frame {:?}", frame + 1);

        // TODO: Add error handling if the cons_aa can't be translated
        let cons_aa = translate(&query[frame..], true, true, true)?;
        let alignment = aligner.local(cons_aa.as_slice(), ref_aa_slice);

        log::info!(
            "Alignment with query in frame {:?} gave a score of {:?}",
            frame,
            alignment.score
        );
        if alignment.score > best_score {
            best_score = alignment.score;
            best_frame = frame;
            best_translation = cons_aa[alignment.xstart..alignment.xend].to_vec();
            best_alignment = alignment.clone();
        }
    }

    log::info!(
        "Choosing translation in frame {:?} with score {:?}:\n{:?}",
        best_frame,
        best_score,
        String::from_utf8(best_translation.clone())?
    );
    log::info!(
        "With this best alignment, we trimmed the query sequence (AA) from position {} to {}",
        best_alignment.xstart,
        best_alignment.xend
    );

    // TODO: Allow outputting both

    if output_type == "AA" {
        log::info!(
            "Writing trimmed query amino acid sequence to {:?}",
            output_file
        );
        write_fasta(output_file, output_seq_name, &best_translation)?;
    } else {
        if output_type != "NT" {
            log::error!("Unrecognized output type, outputting NT sequence")
        }

        let trim_nt_start = (best_alignment.xstart * 3) + best_frame;
        let trim_nt_end = (best_alignment.xend * 3) + best_frame;
        log::info!("Trimming NT from {:?} to {:?}", trim_nt_start, trim_nt_end);
        let trimmed_nt = query[trim_nt_start..trim_nt_end].to_vec();
        write_fasta(output_file, output_seq_name, &trimmed_nt)?;
        log::info!("Outputting NT sequence to {:?}", output_file);
    }

    Ok(())
}
