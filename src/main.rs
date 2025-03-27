mod tools;

use std::path::PathBuf;
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
    SecondTool {
        // second tool specific arguments
    },
    ThirdTool {
        // third tool specific arguments
    }
}

fn main() -> Result<()>{
    let cli = Cli::parse();

    match &cli.command {
        Commands::ReverseTranslate { aa_filepath, nt_filepath, name_mapping, output_file_path, check_keys_match} => {
            tools::reverse_translate::run(aa_filepath, nt_filepath, name_mapping, output_file_path, *check_keys_match)?
        },
        Commands::SecondTool { /* args */ } => {
            // Call second tool's implementation
        },
        Commands::ThirdTool { /* args */ } => {
            // Call third tool's implementation
        }
    }
    Ok(())
}
