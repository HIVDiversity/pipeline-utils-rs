use crate::utils::translate::translate;
use anyhow::{Context, Result};
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use bio::io::fasta;
use colored::Colorize;
use std::collections::HashMap;
use std::iter::Iterator;
use std::path::PathBuf;

const VERSION: &str = "0.3.6";

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
    gap_open_penalty: i32,
    gap_extend_penalty: i32,
) -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;

    log::info!(
        "{}",
        format!("This is pairwise-align-to-ref version {}", VERSION)
            .bold()
            .bright_green()
    );

    log::info!(
        "Using a gap open penalty of {} and a gap extend penalty of {}",
        gap_open_penalty,
        gap_extend_penalty
    );

    let reference_read = read_fasta(reference_file)?;
    let reference = reference_read[0].as_slice();

    let query_read = read_fasta(query_file)?;
    let query = query_read[0].as_slice();

    let ref_aa = translate(reference, true, false, true)
        .context("Error in translating the reference string to amino acids")?;

    let ref_aa_slice = ref_aa.as_slice();

    let scoring = Scoring::new(gap_open_penalty, gap_extend_penalty, bio::scores::blosum62)
        .yclip(MIN_SCORE)
        .xclip(-10);

    let mut aligner = Aligner::with_capacity_and_scoring(query.len() / 3, ref_aa.len(), scoring);

    let mut best_functional_score = 0;
    let mut best_score = 0;
    let mut best_functional_frame: Option<usize> = None;
    let mut best_frame: usize = 0;
    let mut best_translation: Vec<u8> = Vec::with_capacity(query.len() / 3);
    let mut best_alignment: Alignment = Default::default();

    for frame in 0..3 {
        log::info!("Translating query in frame {:?}", frame + 1);

        // TODO: Add error handling if the cons_aa can't be translated
        let cons_aa = translate(&query[frame..], true, true, true)?;
        let alignment = aligner.custom(cons_aa.as_slice(), ref_aa_slice);

        log::info!(
            "Alignment with query in frame {:?} gave a score of {:?}",
            frame,
            alignment.score
        );
        log::info!(
            "\n{}",
            alignment.pretty(cons_aa.as_slice(), ref_aa_slice, 120)
        );
        let putative_translation = cons_aa[alignment.xstart..alignment.xend].to_vec();

        let translation_starts_with_m =
            putative_translation.get(0).unwrap_or(&b"?"[0]).eq(&b"M"[0]);
        let alignment_has_better_score = alignment.score > best_functional_score;

        match (translation_starts_with_m, alignment_has_better_score) {
            (true, true) => {
                best_functional_score = alignment.score;
                best_score = alignment.score;
                best_functional_frame = Some(frame);
                best_frame = frame;
                best_translation = putative_translation.clone();
                best_alignment = alignment.clone();
            }
            (true, false) => {}
            (false, true) => {
                log::info!("The translation in frame {:?} had a better score {:?} than the previous best {:?} in frame {:?}. We're still going to keep the previous best, but this could indicate a frameshift in the alignment.",
                frame, alignment.score, best_functional_score,best_functional_frame);
                best_score = alignment.score;
                best_frame = frame;
            }
            (false, false) => {}
        }
    }

    best_frame = match best_functional_frame{
        Some(frame) => frame,
        None =>{
            log::info!("We were unable to find a frame that started with an M. We'll thus use the sequence with the best alignment score, even though it doesn't start with an M.");

            best_translation = translate(&query[best_frame..], true, true, true)?;
            best_alignment = aligner.custom(best_translation.as_slice(), ref_aa_slice);
            best_frame
        }
    };


    log::info!(
        "Choosing translation in frame {:?} with score {:?}:\n{:?}",
        best_frame,
        best_functional_score,
        String::from_utf8(best_translation.clone())?
    );
    log::info!(
        "With this best alignment, we trimmed the query sequence (AA) from position {} to {}",
        best_alignment.xstart,
        best_alignment.xend
    );

    log::info!(
        "Alignment:\n{}",
        best_alignment.pretty(
            translate(&query[best_frame..], true, true, true)?.as_slice(),
            ref_aa_slice,
            120
        )
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
