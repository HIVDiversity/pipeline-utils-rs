use anyhow::Result;
use clap::Parser;
use purs::cli;
use purs::cli::Commands;
use purs::tools;

fn main() -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;

    let cli = cli::Cli::parse();

    match cli.command {
        Commands::ReverseTranslate {
            aa_filepath,
            nt_filepath,
            output_file_path,
        } => {
            tools::reverse_translate::run(&aa_filepath, &nt_filepath, &output_file_path)?;
        }
        Commands::GetConsensus {
            input_msa,
            output_file,
            consensus_name,
            ambiguity_mode,
        } => {
            tools::get_consensus::run(&input_msa, &output_file, &consensus_name, ambiguity_mode)?;
        }
        Commands::Translate {
            input_file,
            output_file,
            translation_options,
        } => {
            tools::translate::run(&input_file, &output_file, &(&translation_options).into())?;
        }
        Commands::Collapse {
            input_file,
            output_file,
            name_output_file,
            strip_gaps,
            sequence_prefix,
        } => {
            tools::collapse::run(
                &input_file,
                &output_file,
                &name_output_file,
                &sequence_prefix,
                strip_gaps,
            )?;
        }
        Commands::Expand {
            input_file,
            name_input_file,
            output_file,
            include_missing,
        } => {
            tools::expand::run(&input_file, &name_input_file, &output_file, include_missing)?;
        }
        Commands::FilterByLength {
            input_file,
            output_file,
            report_file,
            rejected_seq_output,
            threshold,
        } => {
            tools::filter_by_length::run(
                &input_file,
                &output_file,
                report_file.as_ref(),
                rejected_seq_output.as_ref(),
                (&threshold).into(),
            )?;
        }
        Commands::GbExtract {
            input_file,
            output_file,
            seq_name,
        } => {
            tools::gb_extract::run(&input_file, &output_file, &seq_name)?;
        }
        #[cfg(feature = "trim-sam")]
        Commands::TrimSam {
            input_file,
            output_file,
            trim_from,
            trim_to,
        } => {
            tools::trim_sam::run(&input_file, &output_file, trim_from, trim_to)?;
        }
        Commands::ReplaceAmbiguities {
            input_file,
            output_file,
            seed,
        } => {
            tools::replace_ambiguities::run(&input_file, &output_file, seed)?;
        }
        #[cfg(feature = "process-miniprot")]
        Commands::ProcessMiniprot {
            input_file,
            paf_file,
            prepend,
            output_dir,
        } => {
            tools::process_miniprot::run(&input_file, &paf_file, &prepend, &output_dir)?;
        }
        Commands::TrimAfterStop {
            input_file,
            output_file,
            include_stop,
        } => {
            tools::trim_after_stop_codon::run(&input_file, &output_file, include_stop)?;
        }
        Commands::StripGapCols {
            input_file,
            output_file,
            min_gap_pct,
        } => {
            tools::strip_gap_cols::run(&input_file, &output_file, min_gap_pct)?;
        }
    }

    Ok(())
}
