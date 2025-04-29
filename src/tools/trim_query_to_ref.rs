use crate::utils::translate::translate;
use anyhow::{Context, Result};
use bio::alignment::Alignment;
use bio::alignment::pairwise::*;
use bio::alignment::sparse::{find_kmer_matches, lcskpp};
use bio::io::fasta;
use bio::io::fasta::Record;
use colored::Colorize;
use nalgebra::Vector;
use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::f64::MIN;
use std::iter::Iterator;
use std::path::PathBuf;
use std::process::exit;
use clap::ValueEnum;


const VERSION: &str = "0.4.2";

#[derive(ValueEnum, Copy, Clone)]
pub enum AlignmentMode {
    Standard,
    Sparse,
}

#[derive(Clone)]
struct AlignmentResult {
    alignment: Option<Alignment>,
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

fn get_alignment_in_three_frames(ref_seq: &[u8],
                                 query: &[u8],
                                 scoring_function: Scoring<fn(u8, u8) -> i32>,
                                 alignment_mode: AlignmentMode,
) -> Vec<AlignmentResult> {
    let ref_seq_aa = translate(ref_seq, true, true, true).unwrap();
    let mut aligner =
        Aligner::with_capacity_and_scoring(query.len() / 3, ref_seq_aa.len(), scoring_function);
    let mut results: Vec<AlignmentResult> = Vec::with_capacity(3);

    println!("{:?}", String::from_utf8(ref_seq_aa.clone()));

    for frame in 0..3 {
        let query_aa = translate(&query[frame..], true, true, true)
            .with_context(|| {
                format!(
                    "Could not translate the query sequence in frame {:?}",
                    frame + 1
                )
            })
            .unwrap();

        let alignment = aligner.local(query_aa.as_slice(), ref_seq_aa.as_slice());

        let result = AlignmentResult {
            alignment: Some(alignment.clone()),
            frame,
            score: alignment.score,
            start: alignment.xstart,
            stop: alignment.xend,
            trimmed_query: query_aa[alignment.xstart..alignment.xend].to_vec(),
        };


        log::info!(
            "Alignment with query in frame {:?} gave a score of {:?}",
            frame,
            alignment.score
        );

        results.push(result)
    }

    results
}

fn get_alignment_sparse(ref_seq: &[u8], query: &[u8], num_kmers: usize) -> AlignmentResult {
    let matches = find_kmer_matches(query, ref_seq, num_kmers);
    let sparse_aln = lcskpp(&matches, num_kmers);
    let match_path: Vec<(u32, u32)> = sparse_aln.path.iter().map(|i| matches[*i]).collect();

    if match_path.len() < 0 {
        log::error!("No match found!");
        exit(1);
    }

    let start = (match_path[0].1 / 3) as usize;
    let stop = (match_path[match_path.len() - 1].1 / 3) as usize;
    let translated_query = translate(&query[start..stop], false, false, true).unwrap();
    println!(
        "{}",
        String::from_utf8(translated_query.clone()).unwrap()
    );

    AlignmentResult {
        alignment: None,
        frame: 0,
        score: sparse_aln.score as i32,
        start,
        stop,
        trimmed_query: translated_query,
    }
}

fn get_best_translation(
    ref_seq: &[u8],
    query: &[u8],
    scoring_function: Scoring<fn(u8, u8) -> i32>,
    num_kmers: i32,
    alignment_mode: AlignmentMode,
) -> AlignmentResult {
    let mut results;
    match alignment_mode {
        AlignmentMode::Standard => { results = get_alignment_in_three_frames(ref_seq, query, scoring_function, alignment_mode); }
        AlignmentMode::Sparse => { results = vec![get_alignment_sparse(ref_seq, query, num_kmers as usize)] }
    }


    results.sort();
    for (idx, result) in results.iter().enumerate() {
        if result.trimmed_query.starts_with(b"M") {
            if idx == 0 {
                log::info!("The best alignment started with an 'M':");
            } else {
                log::warn!(
                    "The number {:?} best alignment starts with an 'M'. Returning this instead of the first best.",
                    idx + 1
                );
            }

            log::info!("Score: {:?}", result.score);
            log::info!("Frame: {:?}", result.frame);

            match &result.alignment {
                None => {}
                Some(result_aln) => {
                    log::info!(
                        "Alignment:\n{}",
                        result_aln.pretty(
                            translate(&query[result.frame..], true, true, true)
                                .unwrap()
                                .as_slice(),
                            translate(ref_seq, true, true, true).unwrap().as_slice(),
                            120
                        )
                    );
                }
            }

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
    output_type: &String,
    gap_open_penalty: i32,
    gap_extend_penalty: i32,
    num_kmers: i32,
    alignment_mode: AlignmentMode,
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


    let scoring = Scoring::new(
        gap_open_penalty,
        gap_extend_penalty,
        bio::scores::blosum62 as fn(u8, u8) -> i32,
    )
        .yclip(MIN_SCORE)
        .xclip(-10);

    let trimmed_alignment = get_best_translation(reference, query, scoring, num_kmers, alignment_mode);

    let trim_nt_start = (trimmed_alignment.start * 3) + trimmed_alignment.frame;
    let trim_nt_end = (trimmed_alignment.stop * 3) + trimmed_alignment.frame;

    log::info!("Trimming NT from {:?} to {:?}", trim_nt_start, trim_nt_end);

    let trimmed_nt = query[trim_nt_start..trim_nt_end].to_vec();
    write_fasta(output_file, output_seq_name, &trimmed_nt)?;
    log::info!("Outputting NT sequence to {:?}", output_file);

    let k = 8;
    let matches = find_kmer_matches(reference, query, k);
    let sparse_aln = lcskpp(&matches, k);
    let match_path: Vec<(u32, u32)> = sparse_aln.path.iter().map(|i| matches[*i]).collect();
    if match_path.len() > 0 {
        let new_query =
            &query[match_path[0].1 as usize..(match_path[match_path.len() - 1].1) as usize];
        let new_seq = Record::with_attrs("zoop", None, new_query);
    }


    Ok(())
}
