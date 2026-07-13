use crate::utils::fasta_utils::{FastaRecords, load_fasta, write_fasta_sequences};
use anyhow::{Result, bail};
use colored::Colorize;

use std::path::PathBuf;

use regex::Regex;

pub(crate) fn filter_by_name(
    sequences: FastaRecords,
    pattern: Regex,
    exclude: bool,
) -> Result<(FastaRecords, FastaRecords)> {
    if sequences.is_empty() {
        bail!("No sequences were provided.")
    }

    let mut kept_sequences = FastaRecords::with_capacity(sequences.len());
    let mut rejected_sequences = FastaRecords::new();

    for (seq_name, seq) in sequences {
        match pattern.is_match(seq_name.as_str()) {
            true => {
                if exclude {
                    rejected_sequences.insert(seq_name, seq);
                } else {
                    kept_sequences.insert(seq_name, seq);
                }
            }
            false => {
                if exclude {
                    kept_sequences.insert(seq_name, seq);
                } else {
                    rejected_sequences.insert(seq_name, seq);
                }
            }
        }
    }
    Ok((kept_sequences, rejected_sequences))
}

pub fn run(
    input_file: &PathBuf,
    output_file: &PathBuf,
    rejected_seq_output: Option<&PathBuf>,
    pattern_string: String,
    exclude: bool,
) -> Result<()> {
    log::info!(
        "{}",
        format!(
            "This is 'filter-by-name' version {}",
            env!("CARGO_PKG_VERSION")
        )
        .bold()
        .bright_white()
    );

    log::info!("Reading input file {:?}", input_file);
    let sequences = load_fasta(input_file)?;
    let pattern = Regex::new(pattern_string.as_str())?;
    let (kept_sequences, rejected_sequences) = filter_by_name(sequences, pattern, exclude)?;

    write_fasta_sequences(output_file, &kept_sequences)?;

    if let Some(rejected_seq_output) = rejected_seq_output {
        log::info!("Writing rejected sequences to {:?}", rejected_seq_output);
        write_fasta_sequences(rejected_seq_output, &rejected_sequences)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn to_fasta_records(names: Vec<&str>, sequences: Vec<&str>) -> FastaRecords {
        names
            .iter()
            .zip(sequences.iter())
            .map(|(name, seq)| (name.to_owned().to_owned(), seq.as_bytes().to_vec()))
            .collect::<HashMap<String, Vec<u8>>>()
    }

    #[test]
    fn test_simple_exclude_filter() -> Result<()> {
        let names = vec![
            "this_is_a_name",
            "_here_is_a_nother",
            "thisNameHasNoUnderscores",
        ];
        let sequences = vec!["xyz", "abc", "def"];
        let records = to_fasta_records(names, sequences);

        let pattern = Regex::new("^_")?;

        let (kept, rejected) = filter_by_name(records, pattern, true)?;
        assert!(rejected.contains_key("_here_is_a_nother"));
        assert!(kept.contains_key("this_is_a_name"));
        assert!(kept.contains_key("thisNameHasNoUnderscores"));
        Ok(())
    }

    #[test]
    fn test_simple_include_filter() -> Result<()> {
        let names = vec![
            "this_is_a_name",
            "_here_is_a_nother",
            "thisNameHasNoUnderscores",
        ];
        let sequences = vec!["xyz", "abc", "def"];
        let records = to_fasta_records(names, sequences);

        let pattern = Regex::new("^_")?;

        let (kept, rejected) = filter_by_name(records, pattern, false)?;
        assert!(kept.contains_key("_here_is_a_nother"));
        assert!(rejected.contains_key("this_is_a_name"));
        assert!(rejected.contains_key("thisNameHasNoUnderscores"));
        Ok(())
    }

    #[test]
    fn test_exclude_alphanumeric_regex() -> Result<()> {
        let names = vec![
            "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
            "C002_CAP177_1110_env_pb-GAGGGGAT_p-fs22-ma095",
            "C002_CAP177_1110_env_pb-AGATAGTA_p-fs41-ma090",
        ];
        let sequences = vec!["xyz", "abc", "def"];
        let records = to_fasta_records(names, sequences);

        let pattern = Regex::new("HXB2")?;

        let (kept, rejected) = filter_by_name(records, pattern, true)?;
        assert!(rejected.contains_key("B.FR.83.HXB2_LAI_IIIB_BRU.K03455"));
        assert!(kept.contains_key("C002_CAP177_1110_env_pb-GAGGGGAT_p-fs22-ma095"));
        assert!(kept.contains_key("C002_CAP177_1110_env_pb-AGATAGTA_p-fs41-ma090"));
        Ok(())
    }
    #[test]
    fn test_exclude_complex_regex() -> Result<()> {
        let names = vec![
            "C002_CAP177_1110_env_pb-AGATAGTA_p-fs41-ma090",
            "C002_CAP007_2000_env_pb-GAGGGGAT_p-fs22-ma095",
            "C002_CAP177_3500_env_pb-GAGGGGAT_p-fs22-ma095",
        ];
        let sequences = vec!["xyz", "abc", "def"];
        let records = to_fasta_records(names, sequences);

        let pattern = Regex::new("CAP.+_[1,2].+")?;

        let (kept, rejected) = filter_by_name(records, pattern, false)?;
        assert!(rejected.contains_key("C002_CAP177_3500_env_pb-GAGGGGAT_p-fs22-ma095"));
        assert!(kept.contains_key("C002_CAP177_1110_env_pb-AGATAGTA_p-fs41-ma090"));
        assert!(kept.contains_key("C002_CAP007_2000_env_pb-GAGGGGAT_p-fs22-ma095"));
        Ok(())
    }
}
