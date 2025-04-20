mod tools;
mod utils;

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;
use crate::tools::trim_seqs_to_query::OperatingMode;

#[derive(clap::ValueEnum, Clone)]
enum SequenceOutputType {
    AA,
    NT,
}

#[derive(Parser)]
#[command(name = "ap-utils")]
#[command(about = "A collection of CLI utilities for the alignment pipeline")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    ReverseTranslate {
        /// Path to the aligned amino acid FASTA file
        #[arg(short = 'i', long)]
        aa_filepath: PathBuf,

        /// Path to the unaligned FASTA file containing nucleotide sequences
        #[arg(short = 'n', long)]
        nt_filepath: PathBuf,

        /// Where to write the translated, aligned nt FASTA file
        #[arg(short, long)]
        output_file_path: PathBuf,

    },
    GetConsensus {
        /// Path to the input MSA FASTA file
        #[arg(short = 'i', long)]
        input_msa: PathBuf,

        /// Path to the consensus sequence as a FASTA file
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        ///What to name the consensus sequence in the FASTA file
        #[arg(short = 'n', long)]
        consensus_name: String,
    },
    /// Given a long consensus sequence containing a shorter reference sequence, extract the shorter
    /// reference sequence by aligning the ref seq to the cons seq and trimming.
    AlignConsensus {
        /// Path to the FASTA file containing the reference seq. Note that only the first sequence in the file is used if multiple are present.
        #[arg(short = 'r', long)]
        reference_file: PathBuf,

        /// Path to the FASTA file containing the query seq. Note that only the first sequence in the file is used if multiple are present.
        #[arg(short = 'q', long)]
        query_file: PathBuf,

        /// Path to write the resulting trimmed sequence
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        /// Output sequence name
        #[arg(short = 'n', long)]
        output_seq_name: String,

        /// Strip gaps from both the reference and the query before translating and aligning
        #[arg(short = 's', long, default_value_t = true)]
        strip_gaps: bool,

        /// The gap open penalty. Do not use a negative number here - enter a positive number and it
        /// will be converted to a negative one. e.g 10 becomes -10
        #[arg(long, default_value_t = 5)]
        gap_open_penalty: i32,

        /// The gap extension penalty. Do not use a negative number here - enter a positive number and it
        /// will be converted to a negative one. e.g 10 becomes -10
        #[arg(long, default_value_t = 1)]
        gap_extension_penalty: i32,

        /// What type of sequence to write, either AA or NT
        #[arg(short='t', long, default_value_t = String::from("AA"))]
        output_type: String,
    },
    AlignAndTrim {
        /// The FASTA file containing the untrimmed nucleotide sequences
        #[arg(short = 'q', long)]
        query_sequences: PathBuf,

        /// The file with the NT consensus sequence.
        #[arg(short = 'c', long)]
        consensus_sequence: PathBuf,

        /// The output file
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        /// The size of the kmer to use to compare the consensus to the query. Increasing this value will increase runtime. Decreasing it will reduce specificity.
        #[arg(short='k', long, default_value_t = 10 as i32)]
        kmer_size: i32,

        /// The maximum distance that the kmer can be from where it matches. Higher numbers here mean less specific matches
        #[arg(short='m', long, default_value_t = 2 as i32)]
        max_dist: i32,

        /// What type of sequence to write, either AA or NT
        #[arg(short='t', long, default_value_t = String::from("AA"))]
        output_type: String,

        /// Operating mode
        #[clap(short='d', long, value_enum, default_value_t = OperatingMode::DoubleMatch)]
        operating_mode: OperatingMode
    },
    Translate {
        /// The FASTA-formatted file containing the nucleotide sequences to translate
        #[arg(short = 'i', long)]
        input_file: PathBuf,

        /// The output file to write the translated amino acid sequences to.
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        /// Remove all the gaps from the sequence before translating
        #[arg(short = 's', long, default_value_t = false)]
        strip_gaps: bool,

        /// Don't strip gaps, but if a 'codon' contains only gaps, do not include it in the final
        /// translation. Using this flag will preserve frameshifts as 'X' but avoid '-' in the
        /// translation
        #[arg(short = 'g', long, default_value_t = false)]
        ignore_gap_codons: bool,

        /// The total number of nucleotides in a sequence may not be a multiple of three. The last
        /// codon thus may not be a complete codon. Setting this flag drops that amino acid.
        /// Alternatively, leaving this flag off will replace codons with insufficient nucleotides
        /// with a special character.
        #[arg(short = 'd', long, default_value_t = false)]
        drop_incomplete_codons: bool,
    },
    Collapse{
        /// The FASTA-formatted file containing the uncollapsed sequences
        #[arg(short = 'i', long)]
        input_file: PathBuf,

        /// The output file to write the collapsed sequences to
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        /// The file to write the name mapping to (JSON)
        #[arg(short = 'n', long)]
        name_output_file: PathBuf,

        /// If set, sequences are collapsed purely on nucleotide/amino acid identity, not taking
        /// into account the gap patterns. By default, this is off, and two sequences with identical
        /// sequence information but different gap patterns are considered different.
        #[arg(short='s', long, default_value_t = false)]
        strip_gaps: bool,

        /// The prefix to append to the new sequences when they are collapsed. By default, a unique
        /// integer will be assigned to each sequence, but we can add a string before it
        #[arg(short = 'p', long)]
        sequence_prefix: String


    },

    Expand{
        /// The FASTA-formatted file containing the collapsed sequences
        #[arg(short = 'i', long)]
        input_file: PathBuf,

        /// The JSON file containing a mapping between the current names and the old names
        #[arg(short = 'n', long)]
        name_input_file: PathBuf,

        /// The output file to write the un-collapsed sequences to
        #[arg(short = 'o', long)]
        output_file: PathBuf,



    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match &cli.command {
        Commands::ReverseTranslate {
            aa_filepath,
            nt_filepath,
            output_file_path,
        } => tools::reverse_translate::run(
            aa_filepath,
            nt_filepath,
            output_file_path,
        )?,
        Commands::GetConsensus {
            input_msa,
            output_file,
            consensus_name,
        } => tools::get_consensus::run(input_msa, output_file, consensus_name)?,
        Commands::AlignConsensus {
            reference_file,
            query_file,
            output_file,
            output_seq_name,
            strip_gaps,
            gap_open_penalty,
            gap_extension_penalty,
            output_type,
        } => {
            tools::trim_query_to_ref::run(
                reference_file,
                query_file,
                output_file,
                output_seq_name,
                *strip_gaps,
                output_type,
                (*gap_open_penalty) * -1,
                (*gap_extension_penalty) * -1
            )?;
        }
        Commands::AlignAndTrim {
            query_sequences,
            consensus_sequence,
            output_file,
            kmer_size,
            max_dist,
            output_type,
            operating_mode
        } => {
            tools::trim_seqs_to_query::run(
                query_sequences,
                consensus_sequence,
                output_file,
                *kmer_size,
                output_type,
                *max_dist,
                operating_mode.clone()
            )?;
        }
        Commands::Translate {
            input_file,
            output_file,
            strip_gaps,
            ignore_gap_codons,
            drop_incomplete_codons,
        } => {
            tools::translate::run(
                input_file,
                output_file,
                *strip_gaps,
                *ignore_gap_codons,
                *drop_incomplete_codons,
            )?;
        }
        Commands::Collapse {
            input_file,
            output_file,
            name_output_file,
            strip_gaps,
            sequence_prefix
        } => {
            tools::collapse::run(input_file, output_file, name_output_file, sequence_prefix, *strip_gaps)?;
        }
        Commands::Expand {input_file, name_input_file, output_file} =>{
            tools::expand::run(input_file, name_input_file, output_file)?;
        }

    }
    Ok(())
}
