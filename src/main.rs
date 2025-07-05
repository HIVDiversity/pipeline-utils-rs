mod tools;
mod utils;

use crate::tools::kmer_trim::OperatingMode;
use crate::tools::pairwise_align_trim::AlignmentMode;
use crate::utils::translate::{DEFAULT_STOP_CHAR, TranslationOptions};
use anyhow::Result;
use clap::builder::styling;
use clap::{Args, Parser, Subcommand};
use log::{Level, LevelFilter};
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
#[command(name = "ap-utils")]
#[command(about = "A collection of CLI utilities for the alignment pipeline")]
#[command(styles = STYLES)]
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
    },
    /// Align and trim sequences to a reference sequence.
    /// Given a long consensus sequence containing a shorter reference sequence, extract the shorter
    /// reference sequence by aligning the ref seq to the cons seq and trimming.
    AlignTrim {
        /// Path to the FASTA file containing the reference seq. Note that only the first sequence in the file is used if multiple are present.
        #[arg(short = 'r', long)]
        reference_file: PathBuf,

        /// Path to the FASTA file containing the query seq. Note that all the sequences in the file will be processed.
        #[arg(short = 'i', long)]
        query_file: PathBuf,

        /// Path to write the resulting trimmed sequence
        #[arg(short = 'o', long)]
        output_file: PathBuf,

        /// The gap open penalty. Do not use a negative number here - enter a positive number and it
        /// will be converted to a negative one. e.g 10 becomes -10
        #[arg(long, default_value_t = 5)]
        gap_open_penalty: i32,

        /// The gap extension penalty. Do not use a negative number here - enter a positive number and it
        /// will be converted to a negative one. e.g 1 becomes -1
        #[arg(long, default_value_t = 1)]
        gap_extension_penalty: i32,

        /// What algorithm to use under the hood. Custom uses a semi-global approach, while local is a simple local alignment algorithm
        #[arg(short = 'a', long, value_enum, default_value_t = AlignmentMode::Local)]
        alignment_mode: AlignmentMode,

        /// Options to use when translating
        #[command(flatten)]
        translation_options: TranslateCliOptions,

        /// Number of threads to use. Set to 0 to use all threads.
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,

        /// Turns on verbose logging. Not recommended for multiple sequences.
        #[arg(short = 'v', long, default_value_t = false)]
        verbose: bool,

        /// Turns off logging except for errors. Will override the verbose setting.
        #[arg(short = 'q', long, default_value_t = false)]
        quiet: bool,
    },
    /// Trim sequences to a reference sequence using a k-mer matching approach.
    KmerTrim {
        /// The FASTA file containing the untrimmed nucleotide sequences
        #[arg(short = 'i', long)]
        query_seq_file: PathBuf,

        /// The file with the NT consensus sequence.
        #[arg(short = 'r', long)]
        ref_seq_file: PathBuf,

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
        operating_mode: OperatingMode,

        /// Options to use when translating
        #[command(flatten)]
        translation_options: TranslateCliOptions,
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
    /// Remove repeated sequences in a file. Resulting file contains only unique sequences.
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
}

fn main() -> Result<()> {
    let cli = Cli::parse();

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
        } => tools::get_consensus::run(input_msa, output_file, consensus_name)?,
        Commands::AlignTrim {
            reference_file,
            query_file,
            output_file,
            gap_open_penalty,
            gap_extension_penalty,
            alignment_mode,
            translation_options,
            threads,
            verbose,
            quiet,
        } => {
            let log_level = match (verbose, quiet) {
                (true, true) => LevelFilter::Error,
                (false, true) => LevelFilter::Error,
                (true, false) => LevelFilter::Debug,
                (false, false) => LevelFilter::Info,
            };
            tools::pairwise_align_trim::run(
                reference_file,
                query_file,
                output_file,
                (*gap_open_penalty) * -1,
                (*gap_extension_penalty) * -1,
                *alignment_mode,
                &(translation_options.into()),
                *threads,
                log_level,
            )?;
        }
        Commands::KmerTrim {
            query_seq_file,
            ref_seq_file,
            output_file,
            kmer_size,
            max_dist,
            output_type,
            operating_mode,
            translation_options,
        } => {
            tools::kmer_trim::run(
                query_seq_file,
                ref_seq_file,
                output_file,
                *kmer_size,
                output_type,
                *max_dist,
                operating_mode.clone(),
                &(translation_options.into()),
            )?;
        }
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
    }
    Ok(())
}
