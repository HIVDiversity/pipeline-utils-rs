use crate::utils;
use crate::utils::fasta_utils::SequenceType;
use crate::utils::translate::{STOP_CHAR, translate};
use anyhow::{Context, Result};
use bio::alignment::Alignment;
use bio::pattern_matching::myers::Myers;
use clap::ValueEnum;
use colored::Colorize;
use fasta_utils::FastaRecords;
use std::iter::Iterator;
use std::path::PathBuf;
use utils::fasta_utils;
use utils::translate;

const VERSION: &str = "0.3.0";

#[derive(ValueEnum, Clone, Copy)]
pub enum OperatingMode {
    DoubleMatch,
    SingleMatch,
}

fn find_best_alignment(pattern: &[u8], query: &[u8], max_distance: u8) -> Option<Alignment> {
    let mut pattern = Myers::<u64>::new(pattern);
    let mut matches = pattern.find_all_lazy(query, max_distance);

    // TODO: What happens if we have multiple acceptable matches?
    let (best_match_end_idx, dist) = matches.by_ref().min_by_key(|&(_, dist)| dist)?;

    let mut alignment = Alignment::default();
    matches.alignment_at(best_match_end_idx, &mut alignment);

    Some(alignment)
}

/// Given a query sequence and a consensus sequence, use the first k nucleotides of the consensus
/// to determine the reading frame of the resulting coding sequence. Translate the nt sequence into
/// amino acids. Then trim the sequence to the first available stop codon. If there is no stop
/// codon, return the whole sequence.
fn process_sequence_single_match(
    consensus_start_kmer: &[u8],
    query: &[u8],
    max_align_distance: u8,
    output_type: SequenceType,
) -> Result<Vec<u8>> {
    let start_aln = find_best_alignment(consensus_start_kmer, query, max_align_distance)
        .with_context(|| "No best alignment found.")?;

    let new_nt_seq = &query[start_aln.ystart..].to_owned();
    let new_aa_seq = translate(new_nt_seq, false, false, true)?;

    // Find the first stop codon, or set it to the length of the string
    let first_stop_codon = new_aa_seq
        .iter()
        .position(|&amino_acid| amino_acid == STOP_CHAR)
        .unwrap_or(new_aa_seq.len() - 1);

    match output_type {
        SequenceType::Nucleotide => {
            // If we return nucleotides, then we convert aa_idx to nt_idx
            let nt_end_idx = ((first_stop_codon + 1) * 3);
            Ok(new_nt_seq[..nt_end_idx].to_vec())
        }
        SequenceType::AminoAcid => Ok(new_aa_seq[..first_stop_codon].to_vec()),
    }
}

fn process_sequence_double_match(
    consensus_start_kmer: &[u8],
    consensus_end_kmer: &[u8],
    query: &[u8],
    seq_name: &String,
    max_align_distance: u8,
    output_type: SequenceType,
) -> Result<Vec<u8>> {
    let query_reversed = query.iter().rev().cloned().collect::<Vec<u8>>();

    let Some(start_aln) = find_best_alignment(consensus_start_kmer, query, max_align_distance)
    else {
        log::warn!("No best start alignment found for {:?}.", seq_name);
        return Ok(query.to_vec());
    };

    // Note - the end kmer is assumed to be reversed already!
    let Some(end_aln) = find_best_alignment(
        consensus_end_kmer,
        query_reversed.as_slice(),
        max_align_distance,
    ) else {
        log::warn!("No best end alignment found for {:?}.", seq_name);

        // If we don't find the end alignment, we just return the protein trimmed from the start to the whole alignment
        // But we need to make sure trimming is viable
        return Ok(query
            .get(start_aln.ystart..)
            .unwrap_or_else(|| {
                log::warn!(
                    "Trimming the sequence {:?} failed. Tried to trim from {:?} to the end",
                    seq_name,
                    start_aln.ystart
                );
                return query;
            })
            .to_vec());
    };

    log::info!(
        "<{:?}> Found an alignment:\n - start k-mer from {} to {} (dist {})\n - end k-mer from {} to {} (dist {})",
        seq_name,
        start_aln.ystart,
        start_aln.yend,
        start_aln.score,
        end_aln.ystart,
        end_aln.yend,
        end_aln.score
    );

    let start_trim = start_aln.ystart;
    let mut end_trim = query.len() - end_aln.ystart;
    let remainder = (end_trim - start_trim) % 3;

    if remainder > 0 {
        if end_trim + remainder > query.len() {
            end_trim = end_trim - remainder;
        } else {
            end_trim = end_trim + (3 - remainder);
        }
    }

    let trimmed_query = query.get(start_trim..end_trim).unwrap_or_else(|| {
        log::warn!(
            "Trimming the sequence {:?} failed. Tried to trim from {:?} to {:?}",
            seq_name,
            start_trim,
            end_trim
        );
        query
    });

    match output_type {
        SequenceType::Nucleotide => Ok(trimmed_query.to_vec()),
        SequenceType::AminoAcid => {
            let translated_query = translate::translate(trimmed_query, false, false, false)?;
            Ok(translated_query)
        }
    }
}

fn process_file(
    query_file: &PathBuf,
    consensus: &[u8],
    kmer_size: i32,
    max_align_distance: u8,
    output_type: SequenceType,
    operating_mode: OperatingMode,
) -> Result<FastaRecords> {
    // No matter which mode we operate in, we need a start kmer
    let start_query = &consensus[0..kmer_size as usize];
    let query_sequences = fasta_utils::load_fasta(query_file)?;
    let mut trimmed_sequences: FastaRecords = FastaRecords::new();

    match operating_mode {
        OperatingMode::DoubleMatch => {
            // If we're operating in double match mode, then we need both a start and end kmer
            let final_index = consensus.len() as i32 - kmer_size;

            // The end query also needs to be reversed
            let end_query = consensus[final_index as usize..]
                .iter()
                .rev()
                .cloned()
                .collect::<Vec<u8>>();

            for (seq_id, seq) in query_sequences {
                trimmed_sequences.insert(
                    seq_id.clone(),
                    process_sequence_double_match(
                        start_query,
                        end_query.as_slice(),
                        seq.as_slice(),
                        &seq_id,
                        max_align_distance,
                        output_type,
                    )?,
                );
            }
        }
        OperatingMode::SingleMatch => {
            for (seq_id, seq) in query_sequences {
                trimmed_sequences.insert(
                    seq_id,
                    process_sequence_single_match(
                        start_query,
                        seq.as_slice(),
                        max_align_distance,
                        output_type,
                    )?,
                );
            }
        }
    }

    Ok(trimmed_sequences)
}

/// Consumes full-length sequences and trims them to the region covered by the reference.
///
/// Given a set of unaligned, untranslated sequences S, and a query sequence Q, we will trim
/// all the sequences in S to the same region as covered by Q.
pub fn run(
    input_file: &PathBuf,
    consensus_file: &PathBuf,
    output_file: &PathBuf,
    kmer_size: i32,
    output_type_str: &String,
    max_align_distance: i32,
    mode: OperatingMode,
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

    let mut output_type: SequenceType;

    if output_type_str == "AA" {
        output_type = SequenceType::AminoAcid;
    } else if output_type_str == "NT" {
        output_type = SequenceType::Nucleotide;
    } else {
        log::warn!("Uknown Sequence Type. Defaulting to NT");
        output_type = SequenceType::Nucleotide;
    }

    let output_seqs = process_file(
        input_file,
        consensus,
        kmer_size,
        max_align_distance as u8,
        output_type,
        mode,
    )?;

    fasta_utils::write_fasta_sequences(output_file, &output_seqs)?;

    Ok(())
}
