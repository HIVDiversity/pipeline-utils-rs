use crate::utils;
use anyhow::{Context, Result};
use bio::alignment::Alignment;
use bio::pattern_matching::myers::Myers;
use fasta_utils::FastaRecords;
use std::iter::Iterator;
use std::path::PathBuf;
use colored::Colorize;
use utils::fasta_utils;
use utils::translate;

const VERSION: &str = "0.1.2";

fn find_best_alignment(pattern: &[u8], query: &[u8], max_distance: u8) -> Option<Alignment> {
    let mut pattern = Myers::<u64>::new(pattern);
    let mut matches = pattern.find_all_lazy(query, max_distance);

    // TODO: What happens if we have multiple acceptable matches?
    let (best_match_end_idx, dist) = matches.by_ref().min_by_key(|&(_, dist)| dist)?;

    log::info!(
        "Best match found ending at {} with distance {}",
        best_match_end_idx,
        dist
    );

    let mut alignment = Alignment::default();
    matches.alignment_at(best_match_end_idx, &mut alignment);

    Some(alignment)
}

fn process_sequence(
    consensus_start_kmer: &[u8],
    consensus_end_kmer: &[u8],
    query: &[u8],
    max_align_distance: u8,
    output_type: &String,
) -> Result<Vec<u8>> {
    // Note - the end kmer is assumed to be reversed already!
    let query_reversed = query.iter().rev().cloned().collect::<Vec<u8>>();
    let start_aln = find_best_alignment(consensus_start_kmer, query, max_align_distance)
        .with_context(|| "No best alignment found.")?;
    let end_aln = find_best_alignment(consensus_end_kmer, query_reversed.as_slice(), max_align_distance)
        .with_context(|| "No best alignment found")?;

    log::info!("Found an alignment for the start k-mer from {} to {} (dist {}) and alignment for the end k-mer from {} to {} (dist {})",
        start_aln.ystart, start_aln.yend, start_aln.score,
        end_aln.ystart, end_aln.yend, end_aln.score
    );

    let start_trim = start_aln.ystart;
    let end_trim = query.len() - end_aln.ystart;
    let trimmed_query = &query[start_trim..end_trim].to_owned();

    if output_type == "AA" {
        let translated_query = translate::translate(trimmed_query, false, false, false)?;
        Ok(translated_query)
    } else {
        if output_type != "NT" {
            log::warn!("Output type not recognized, outputting NT sequence.")
        }

        Ok(trimmed_query.to_vec())
    }
}

fn process_file(
    query_file: &PathBuf,
    consensus: &[u8],
    kmer_size: i32,
    max_align_distance: u8,
    output_type: &String,
) -> Result<FastaRecords> {
    let final_index = consensus.len() as i32 - kmer_size;
    let start_query = &consensus[0..kmer_size as usize];
    let end_query = consensus[final_index as usize..]
        .iter()
        .rev()
        .cloned()
        .collect::<Vec<u8>>();

    let query_sequences = fasta_utils::load_fasta(query_file)?;
    let mut trimmed_sequences: FastaRecords = FastaRecords::new();

    for query_sequence in query_sequences {
        let trimmed_sequence = process_sequence(
            start_query,
            end_query.as_slice(),
            query_sequence.1.as_slice(),
            max_align_distance,
            output_type,
        )?;
        trimmed_sequences.insert(query_sequence.0, trimmed_sequence);
    }

    Ok(trimmed_sequences)
}

pub fn run(
    input_file: &PathBuf,
    consensus_file: &PathBuf,
    output_file: &PathBuf,
    kmer_size: i32,
    output_type: &String,
    max_align_distance: i32
) -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;

    log::info!(
        "{}",
        format!("This is align-and-trim version {}", VERSION)
            .bold()
            .bright_green()
    );

    let consensus_seq = fasta_utils::load_fasta(consensus_file)?;
    let consensus = consensus_seq
        .values()
        .next()
        .with_context(|| "Consensus file contained no sequences.")?
        .as_slice();

    let output_seqs = process_file(input_file, consensus, kmer_size, max_align_distance as u8, output_type)?;

    fasta_utils::write_fasta_sequences(output_file, &output_seqs)?;

    Ok(())
}
