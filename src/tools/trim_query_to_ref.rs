use crate::utils::translate::translate;
use anyhow::{Context, Result};
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use bio::io::fasta;
use colored::Colorize;
use nalgebra::Vector;
use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::iter::Iterator;
use std::path::PathBuf;

const VERSION: &str = "0.4.1";

#[derive(Clone)]
struct AlignmentResult {
    alignment: Alignment,
    frame: usize,
    score: i32,
    start: usize,
    stop: usize,
    trimmed_query: Vec<u8>,
}

impl Eq for AlignmentResult {}

impl PartialEq<Self> for AlignmentResult {
    fn eq(&self, other: &Self) -> bool {
        self.score.eq(&other.score)
    }
}

impl PartialOrd<Self> for AlignmentResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.score.partial_cmp(&other.score)
    }
}

impl Ord for AlignmentResult {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score.cmp(&other.score)
    }
}

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

fn get_best_translation(
    ref_seq: &[u8],
    query: &[u8],
    scoring_function: Scoring<fn(u8, u8) -> i32>,
) -> AlignmentResult {
    let mut aligner =
        Aligner::with_capacity_and_scoring(query.len() / 3, ref_seq.len(), scoring_function);
    let mut results: Vec<AlignmentResult> = Vec::with_capacity(3);

    for frame in 0..3 {
        let cons_aa = translate(&query[frame..], true, true, true)
            .with_context(|| {
                format!(
                    "Could not translate the query sequence in frame {:?}",
                    frame + 1
                )
            })
            .unwrap();
        let alignment = aligner.custom(cons_aa.as_slice(), ref_seq);

        log::info!(
            "Alignment with query in frame {:?} gave a score of {:?}",
            frame,
            alignment.score
        );

        results.push(AlignmentResult {
            alignment: alignment.clone(),
            frame,
            score: alignment.score,
            start: alignment.xstart,
            stop: alignment.xend,
            trimmed_query: cons_aa[alignment.xstart..alignment.xend].to_vec(),
        })
    }

    results.sort();
    for (idx, result) in results.iter().enumerate() {
        if result.trimmed_query.starts_with(b"M") {
            if idx == 0 {
                log::info!("The best alignment started with an 'M':");
            } else {
                log::warn!("The number {:?} best alignment starts with an 'M'. Returning this instead of the first best.", idx+1);
            }

            log::info!("Score: {:?}", result.score);
            log::info!("Frame: {:?}", result.frame);
            log::info!(
                "Alignment:\n{}",
                result.alignment.pretty(
                    translate(&query[result.frame..], true, true, true)
                        .unwrap()
                        .as_slice(),
                    ref_seq,
                    120
                )
            );

            return result.clone();
        } else {
            log::warn!(
                "Alignment number {:?} with score {:?} in frame {:?} did not start with 'M'",
                idx + 1,
                result.score,
                result.frame + 1
            );
        }
    }

    log::warn!(
        "None of the alignments started with an 'M'. Returning the one with the highest score"
    );
    results[0].clone()
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

    let scoring = Scoring::new(
        gap_open_penalty,
        gap_extend_penalty,
        bio::scores::blosum62 as fn(u8, u8) -> i32,
    )
        .yclip(MIN_SCORE)
        .xclip(-10);

    let trimmed_alignment = get_best_translation(ref_aa_slice, query, scoring);

    let trim_nt_start = (trimmed_alignment.start * 3) + trimmed_alignment.frame;
    let trim_nt_end = (trimmed_alignment.stop * 3) + trimmed_alignment.frame;

    log::info!("Trimming NT from {:?} to {:?}", trim_nt_start, trim_nt_end);

    let trimmed_nt = query[trim_nt_start..trim_nt_end].to_vec();
    write_fasta(output_file, output_seq_name, &trimmed_nt)?;
    log::info!("Outputting NT sequence to {:?}", output_file);

    Ok(())
}
