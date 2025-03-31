mod tools;
mod utils;

use std::path::PathBuf;
use std::process::Output;
use clap::{Subcommand, Parser};
use anyhow::{Error, Result, Context, anyhow};

#[derive(Parser)]
#[command(name = "ap-utils")]
#[command(about = "A collection of CLI utilities for the alignment pipeline")]
struct Cli {
    #[command(subcommand)]
    command: Commands
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

        /// Optionally provide a mapping between the names in the AA fasta file, and the ones in the NT file
        #[arg(short = 'm', long)]
        name_mapping: PathBuf,

        /// Where to write the translated, aligned nt FASTA file
        #[arg(short, long)]
        output_file_path: PathBuf,

        /// Check if the keys in the two files match.
        #[arg(short, long)]
        check_keys_match: bool,
    },
    GetConsensus {
        /// Path to the input MSA FASTA file
        #[arg(short = 'i', long)]
        input_msa: PathBuf,

        /// Path to the consensus sequence as a FASTA file
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        ///What to name the consensus sequence in the FASTA file
        #[arg(short='n', long)]
        consensus_name: String
    },
    AlignConsensus {
        /// Path to the FASTA file containing the reference seq. Note that only the first sequence in the file is used if multiple are present.
        #[arg(short='r', long)]
        reference_file: PathBuf,

        /// Path to the FASTA file containing the query seq. Note that only the first sequence in the file is used if multiple are present.
        #[arg(short='q', long)]
        query_file: PathBuf,

        /// Path to write the resulting trimmed sequence
        #[arg(short='o', long)]
        output_file: PathBuf,

        /// Output sequence name
        #[arg(short='n', long)]
        output_seq_name: String,

        /// Strip gaps from both the reference and the query before translating and aligning
        #[arg(short='s', long, default_value_t = true)]
        strip_gaps: bool,

        /// What type of sequence to write, either AA or NT
        #[arg(short='t', long, default_value_t = String::from("AA"))]
        output_type: String



    },
    AlignAndTrim{

    }
}

fn main() -> Result<()>{
    let cli = Cli::parse();

    match &cli.command {
        Commands::ReverseTranslate { aa_filepath, nt_filepath, name_mapping, output_file_path, check_keys_match} => {
            tools::reverse_translate::run(aa_filepath, nt_filepath, name_mapping, output_file_path, *check_keys_match)?
        },
        Commands::GetConsensus { input_msa, output_file, consensus_name} => {
            tools::get_consensus::run(input_msa, output_file, consensus_name)?
        },
        Commands::AlignConsensus { reference_file, query_file, output_file, output_seq_name, strip_gaps, output_type } => {
            tools::pairwise_align_to_ref::run(reference_file, query_file, output_file, output_seq_name, *strip_gaps, output_type)?;
        },
        Commands::AlignAndTrim {} =>{
            tools::align_and_trim::run();
        }
    }
    Ok(())
}
