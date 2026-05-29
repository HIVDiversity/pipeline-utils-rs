mod cli;
mod tools;
mod utils;

use crate::tools::get_consensus::AmbiguityMode;
use crate::utils::translate::TranslationOptions;
use anyhow::Result;
use clap::builder::styling;
use clap::{Args, Parser, Subcommand};
use log::LevelFilter;
use std::path::PathBuf;

#[derive(clap::ValueEnum, Clone)]
enum SequenceOutputType {
    AA,
    NT,
}

const STYLES: styling::Styles = styling::Styles::styled()
    .header(styling::AnsiColor::Green.on_default().bold())
    .usage(styling::AnsiColor::Green.on_default().bold())
    .literal(styling::AnsiColor::Blue.on_default().bold())
    .placeholder(styling::AnsiColor::Cyan.on_default());

#[derive(Parser)]
#[command(name = "pipeline-utils-rs")]
#[command(about = "A collection of CLI utilities for manipulating sequencing files.")]
#[command(styles = STYLES)]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Args)]
#[group(required = false, multiple = true)]
struct TranslateCliOptions {
    #[arg(long, default_value_t = TranslationOptions::default().unknown_aa as char)]
    unknown_aa: char,
    #[arg(long, default_value_t = TranslationOptions::default().stop_aa as char)]
    stop_aa: char,
    #[arg(long, default_value_t = TranslationOptions::default().incomplete_aa as char)]
    incomplete_aa: char,
    #[arg(long, default_value_t = TranslationOptions::default().frameshift_aa as char)]
    frameshift_aa: char,
    #[arg(long, default_value_t = TranslationOptions::default().reading_frame)]
    reading_frame: usize,
    #[arg(long, default_value_t = TranslationOptions::default().allow_ambiguities)]
    allow_ambiguities: bool,
    #[arg(long, default_value_t = TranslationOptions::default().strip_gaps)]
    strip_gaps: bool,
    #[arg(long, default_value_t = TranslationOptions::default().ignore_gap_codons)]
    ignore_gap_codons: bool,
    #[arg(long, default_value_t = TranslationOptions::default().drop_incomplete_codons)]
    drop_incomplete_codons: bool,
}

impl Into<TranslationOptions> for &TranslateCliOptions {
    fn into(self) -> TranslationOptions {
        TranslationOptions {
            unknown_aa: self.unknown_aa as u8,
            stop_aa: self.stop_aa as u8,
            incomplete_aa: self.incomplete_aa as u8,
            frameshift_aa: self.frameshift_aa as u8,
            reading_frame: self.reading_frame,
            allow_ambiguities: self.allow_ambiguities,
            strip_gaps: self.strip_gaps,
            ignore_gap_codons: self.ignore_gap_codons,
            drop_incomplete_codons: self.drop_incomplete_codons,
        }
    }
}

#[derive(Subcommand)]
/// General purpose utities for the nf-codon-align pipeline, specifically relating to pre-processing data
enum Commands {
    /// Reverse Translate a multiple sequence alignment.
    /// Converts an amino acid alignment back into nucleotides, using the unaligned nucleotide
    /// sequences as a guide. Ensures the original codons are used as the codons in the output.
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
    /// Get the consensus sequence of a Multiple Sequence Alignment.
    /// Produces a single sequence representing all the sequences in the input file, where each
    /// nucleotide in the output sequence is the most common nucleotide at that position.
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

        /// How to handle ambiguous characters
        #[arg(short = 'a', long)]
        ambiguity_mode: AmbiguityMode,
    },

    /// Translate sequences from nucleotides into amino acids.
    Translate {
        /// The FASTA-formatted file containing the nucleotide sequences to translate
        #[arg(short = 'i', long)]
        input_file: PathBuf,

        /// The output file to write the translated amino acid sequences to.
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        /// Options to use when translating
        #[command(flatten)]
        translation_options: TranslateCliOptions,
    },
    /// Remove non-unique sequences in a file. Resulting file contains only unique sequences.
    Collapse {
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
        #[arg(short = 's', long, default_value_t = false)]
        strip_gaps: bool,

        /// The prefix to append to the new sequences when they are collapsed. By default, a unique
        /// integer will be assigned to each sequence, but we can add a string before it
        #[arg(short = 'p', long)]
        sequence_prefix: String,
    },
    /// Re-introduce the duplicate sequences that were removed from the collapse function.
    Expand {
        /// The FASTA-formatted file containing the collapsed sequences
        #[arg(short = 'i', long)]
        input_file: PathBuf,

        /// The JSON file containing a mapping between the current names and the old names
        #[arg(short = 'n', long)]
        name_input_file: PathBuf,

        /// The output file to write the un-collapsed sequences to
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        /// If this flag is set, when a sequence in the input file isn't present in the name file
        /// it will still be included in the output file
        #[arg(short = 'm', long, default_value_t = false)]
        include_missing: bool,
    },
    /// Extract a feature from a genbank file and write it to a FASTA file.
    GbExtract {
        /// The input Genbank File
        #[arg(short = 'i', long)]
        input_file: PathBuf,

        /// The output file to write the sequence to
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        /// The name of the sequence to extract from the genbank file
        #[arg(short = 'n', long)]
        seq_name: String,
    },
    /// Trim a SAM file using coordinates on the reference sequence.
    TrimSam {
        /// The input SAM file
        #[arg(short = 'i', long)]
        input_file: PathBuf,

        /// The output file to write the trimmed sequences to as a FASTA
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        /// The reference position to trim from inclusive, and base 1
        #[arg(short = 'f', long)]
        trim_from: i64,

        /// The reference position to trim to inclusive and base 1
        #[arg(short = 't', long)]
        trim_to: i64,
    },
    /// Convert IUPAC ambiguity codes to one of their possible nucleotides randomly
    ReplaceAmbiguities {
        /// The input FASTA file
        #[arg(short = 'i', long)]
        input_file: PathBuf,

        /// The output file to write the new sequences to as a FASTA
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        /// A seed for the random number generator
        #[arg(short = 's', long, default_value_t = 42)]
        seed: u64,
    },
    /// Given the TSV output from a Diamond Blastx search, return the trimmed templates from a FASTA
    /// file.
    ProcessMiniprot {
        /// The input FASTA file containing nucleotide sequences
        #[arg(short = 'i', long)]
        input_file: PathBuf,

        /// The TSV file output from the blastx command.
        #[arg(short = 'p', long)]
        paf_file: PathBuf,

        /// Optionally, a comma-separated list of strings to prepend onto the filename
        #[arg(long)]
        prepend: Option<String>,

        /// The output directory to write the resulting files to
        #[arg(short = 'o', long)]
        output_dir: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    // simple_logger::SimpleLogger::new().env().init()?;
    // log::info!("This is pipeline-utils-rs version??");

    match &cli.command {
        Commands::ReverseTranslate {
            aa_filepath,
            nt_filepath,
            output_file_path,
        } => tools::reverse_translate::run(aa_filepath, nt_filepath, output_file_path)?,
        Commands::GetConsensus {
            input_msa,
            output_file,
            consensus_name,
            ambiguity_mode,
        } => tools::get_consensus::run(input_msa, output_file, consensus_name, *ambiguity_mode)?,

        Commands::Translate {
            input_file,
            output_file,
            translation_options,
        } => {
            tools::translate::run(input_file, output_file, &(translation_options.into()))?;
        }
        Commands::Collapse {
            input_file,
            output_file,
            name_output_file,
            strip_gaps,
            sequence_prefix,
        } => {
            tools::collapse::run(
                input_file,
                output_file,
                name_output_file,
                sequence_prefix,
                *strip_gaps,
            )?;
        }
        Commands::Expand {
            input_file,
            name_input_file,
            output_file,
            include_missing,
        } => {
            tools::expand::run(input_file, name_input_file, output_file, *include_missing)?;
        }
        Commands::GbExtract {
            input_file,
            output_file,
            seq_name,
        } => {
            tools::extract_seq_from_gb::run(input_file, output_file, seq_name)?;
        }
        Commands::TrimSam {
            input_file,
            output_file,
            trim_from,
            trim_to,
        } => tools::trim_sam::run(input_file, output_file, *trim_from, *trim_to)?,
        Commands::ReplaceAmbiguities {
            input_file,
            output_file,
            seed,
        } => tools::replace_ambiguities::run(input_file, output_file, *seed)?,
        Commands::ProcessMiniprot {
            input_file,
            paf_file,
            prepend,
            output_dir,
        } => tools::process_miniprot::run(input_file, paf_file, prepend, output_dir)?,
    }
    Ok(())
}
