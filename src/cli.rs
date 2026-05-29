use crate::tools::get_consensus::AmbiguityMode;
use crate::utils::translate::TranslationOptions;
use clap::builder::styling;
use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

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
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(clap::ValueEnum, Clone)]
pub enum SequenceOutputType {
    AA,
    NT,
}

#[derive(Args)]
#[group(required = false, multiple = true)]
pub struct TranslateCliOptions {
    #[arg(long, default_value_t = TranslationOptions::default().unknown_aa as char)]
    pub unknown_aa: char,
    #[arg(long, default_value_t = TranslationOptions::default().stop_aa as char)]
    pub stop_aa: char,
    #[arg(long, default_value_t = TranslationOptions::default().incomplete_aa as char)]
    pub incomplete_aa: char,
    #[arg(long, default_value_t = TranslationOptions::default().frameshift_aa as char)]
    pub frameshift_aa: char,
    #[arg(long, default_value_t = TranslationOptions::default().reading_frame)]
    pub reading_frame: usize,
    #[arg(long, default_value_t = TranslationOptions::default().allow_ambiguities)]
    pub allow_ambiguities: bool,
    #[arg(long, default_value_t = TranslationOptions::default().strip_gaps)]
    pub strip_gaps: bool,
    #[arg(long, default_value_t = TranslationOptions::default().ignore_gap_codons)]
    pub ignore_gap_codons: bool,
    #[arg(long, default_value_t = TranslationOptions::default().drop_incomplete_codons)]
    pub drop_incomplete_codons: bool,
}

impl From<&TranslateCliOptions> for TranslationOptions {
    fn from(opts: &TranslateCliOptions) -> Self {
        TranslationOptions {
            unknown_aa: opts.unknown_aa as u8,
            stop_aa: opts.stop_aa as u8,
            incomplete_aa: opts.incomplete_aa as u8,
            frameshift_aa: opts.frameshift_aa as u8,
            reading_frame: opts.reading_frame,
            allow_ambiguities: opts.allow_ambiguities,
            strip_gaps: opts.strip_gaps,
            ignore_gap_codons: opts.ignore_gap_codons,
            drop_incomplete_codons: opts.drop_incomplete_codons,
        }
    }
}

#[derive(Subcommand)]
pub enum Commands {
    /// Reverse translate a multiple sequence alignment.
    /// Converts an amino acid alignment back into nucleotides, using the unaligned nucleotide
    /// sequences as a guide. Ensures the original codons are used in the output.
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

    /// Get the consensus sequence of a multiple sequence alignment.
    /// Produces a single sequence where each position is the most common nucleotide.
    GetConsensus {
        /// Path to the input MSA FASTA file
        #[arg(short = 'i', long)]
        input_msa: PathBuf,
        /// Path to write the consensus sequence as a FASTA file
        #[arg(short = 'o', long)]
        output_file: PathBuf,
        /// Name for the consensus sequence in the FASTA file
        #[arg(short = 'n', long)]
        consensus_name: String,
        /// How to handle ambiguous characters
        #[arg(short = 'a', long)]
        ambiguity_mode: AmbiguityMode,
    },

    /// Translate sequences from nucleotides into amino acids.
    Translate {
        /// The FASTA file containing nucleotide sequences to translate
        #[arg(short = 'i', long)]
        input_file: PathBuf,
        /// The output file to write the translated amino acid sequences to
        #[arg(short = 'o', long)]
        output_file: PathBuf,
        #[command(flatten)]
        translation_options: TranslateCliOptions,
    },

    /// Remove non-unique sequences. Output contains only unique sequences.
    Collapse {
        /// The input FASTA file containing uncollapsed sequences
        #[arg(short = 'i', long)]
        input_file: PathBuf,
        /// The output file to write collapsed sequences to
        #[arg(short = 'o', long)]
        output_file: PathBuf,
        /// The file to write the name mapping to (JSON)
        #[arg(short = 'n', long)]
        name_output_file: PathBuf,
        /// Collapse on sequence identity only, ignoring gap patterns
        #[arg(short = 's', long, default_value_t = false)]
        strip_gaps: bool,
        /// Prefix to prepend to new sequence names after collapsing
        #[arg(short = 'p', long)]
        sequence_prefix: String,
    },

    /// Re-introduce duplicate sequences removed by the collapse command.
    Expand {
        /// The FASTA file containing collapsed sequences
        #[arg(short = 'i', long)]
        input_file: PathBuf,
        /// The JSON file mapping current names to original names
        #[arg(short = 'n', long)]
        name_input_file: PathBuf,
        /// The output file to write the expanded sequences to
        #[arg(short = 'o', long)]
        output_file: PathBuf,
        /// Include sequences not present in the name mapping file
        #[arg(short = 'm', long, default_value_t = false)]
        include_missing: bool,
    },

    /// Extract a feature from a GenBank file and write it to a FASTA file.
    GbExtract {
        /// The input GenBank file
        #[arg(short = 'i', long)]
        input_file: PathBuf,
        /// The output file to write the sequence to
        #[arg(short = 'o', long)]
        output_file: PathBuf,
        /// The name of the sequence to extract
        #[arg(short = 'n', long)]
        seq_name: String,
    },

    /// Trim a SAM file using coordinates on the reference sequence.
    TrimSam {
        /// The input SAM file
        #[arg(short = 'i', long)]
        input_file: PathBuf,
        /// The output FASTA file to write the trimmed sequences to
        #[arg(short = 'o', long)]
        output_file: PathBuf,
        /// The reference position to trim from (inclusive, 1-based)
        #[arg(short = 'f', long)]
        trim_from: i64,
        /// The reference position to trim to (inclusive, 1-based)
        #[arg(short = 't', long)]
        trim_to: i64,
    },

    /// Convert IUPAC ambiguity codes to one of their possible nucleotides randomly.
    ReplaceAmbiguities {
        /// The input FASTA file
        #[arg(short = 'i', long)]
        input_file: PathBuf,
        /// The output FASTA file to write the resolved sequences to
        #[arg(short = 'o', long)]
        output_file: PathBuf,
        /// Seed for the random number generator
        #[arg(short = 's', long, default_value_t = 42)]
        seed: u64,
    },

    /// Given PAF output from miniprot, return trimmed templates from a FASTA file.
    ProcessMiniprot {
        /// The input FASTA file containing nucleotide sequences
        #[arg(short = 'i', long)]
        input_file: PathBuf,
        /// The PAF file output from miniprot
        #[arg(short = 'p', long)]
        paf_file: PathBuf,
        /// Optional comma-separated strings to prepend onto output filenames
        #[arg(long)]
        prepend: Option<String>,
        /// The output directory to write the resulting files to
        #[arg(short = 'o', long)]
        output_dir: PathBuf,
    },

    /// Trims the nucleotides after the first stop codon in a sequence
    TrimAfterStop {
        /// The input FASTA file
        #[arg(short = 'i', long)]
        input_file: PathBuf,
        /// The output FASTA file to write the trimmed sequences to
        #[arg(short = 'o', long)]
        output_file: PathBuf,
        /// Include the stop codon in the output
        #[arg(long, default_value_t = true)]
        include_stop: bool,
    },
}
