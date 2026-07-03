use crate::utils::fasta_utils::{load_fasta, write_fasta_sequences, FastaRecords};
use anyhow::{bail, Result};
use colored::Colorize;
use std::path::PathBuf;

pub enum LengthThreshold {
    Fixed(usize),
    Median,
    Mean,
}

pub(crate) struct FilterReportRow {
    pub(crate) seq_name: String,
    pub(crate) length: usize,
    pub(crate) kept: bool,
}

fn threshold_value(lengths: &[usize], threshold: &LengthThreshold) -> f64 {
    match threshold {
        LengthThreshold::Fixed(n) => *n as f64,
        LengthThreshold::Median => {
            let mut sorted_lengths = lengths.to_vec();
            sorted_lengths.sort_unstable();
            let mid = sorted_lengths.len() / 2;
            if sorted_lengths.len() % 2 == 0 {
                (sorted_lengths[mid - 1] + sorted_lengths[mid]) as f64 / 2f64
            } else {
                sorted_lengths[mid] as f64
            }
        }
        LengthThreshold::Mean => {
            let sum: usize = lengths.iter().sum();
            sum as f64 / lengths.len() as f64
        }
    }
}

pub(crate) fn filter_by_length(
    sequences: FastaRecords,
    threshold: LengthThreshold,
) -> Result<(FastaRecords, Vec<FilterReportRow>)> {
    if sequences.is_empty() {
        bail!("No sequences were provided.")
    }

    let lengths: Vec<usize> = sequences.values().map(|seq| seq.len()).collect();
    let threshold_value = threshold_value(&lengths, &threshold);

    let mut output_sequences = FastaRecords::with_capacity(sequences.len());
    let mut report_rows = Vec::with_capacity(sequences.len());

    for (seq_name, seq) in sequences {
        let length = seq.len();
        let kept = length as f64 >= threshold_value;

        report_rows.push(FilterReportRow {
            seq_name: seq_name.clone(),
            length,
            kept,
        });

        if kept {
            output_sequences.insert(seq_name, seq);
        }
    }

    report_rows.sort_unstable_by(|a, b| a.seq_name.cmp(&b.seq_name));

    Ok((output_sequences, report_rows))
}

fn write_report(report_file: &PathBuf, rows: &[FilterReportRow]) -> Result<()> {
    let mut writer = csv::Writer::from_path(report_file)?;
    writer.write_record(["seq_name", "length", "filter_result"])?;

    for row in rows {
        writer.write_record([
            row.seq_name.as_str(),
            row.length.to_string().as_str(),
            row.kept.to_string().as_str(),
        ])?;
    }

    writer.flush()?;
    Ok(())
}

pub fn run(
    input_file: &PathBuf,
    output_file: &PathBuf,
    report_file: Option<&PathBuf>,
    threshold: LengthThreshold,
) -> Result<()> {
    log::info!(
        "{}",
        format!(
            "This is 'filter-by-length' version {}",
            env!("CARGO_PKG_VERSION")
        )
        .bold()
        .bright_yellow()
    );

    log::info!("Reading input file {:?}", input_file);
    let sequences = load_fasta(input_file)?;
    let (filtered_sequences, report_rows) = filter_by_length(sequences, threshold)?;

    write_fasta_sequences(output_file, &filtered_sequences)?;

    if let Some(report_file) = report_file {
        log::info!("Writing filter report to {:?}", report_file);
        write_report(report_file, &report_rows)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use velcro::hash_map;

    #[test]
    fn test_fixed_threshold() -> Result<()> {
        let input_seqs: FastaRecords = hash_map!(
            "A".to_string(): vec![b'A'; 5],
            "B".to_string(): vec![b'A'; 10],
            "C".to_string(): vec![b'A'; 15],
        );

        let (output, report) = filter_by_length(input_seqs, LengthThreshold::Fixed(10))?;

        assert_eq!(output.len(), 2);
        assert!(output.contains_key("B"));
        assert!(output.contains_key("C"));
        assert!(!output.contains_key("A"));

        assert_eq!(report.len(), 3);
        assert_eq!(report[0].seq_name, "A");
        assert!(!report[0].kept);
        assert_eq!(report[1].seq_name, "B");
        assert!(report[1].kept);
        assert_eq!(report[2].seq_name, "C");
        assert!(report[2].kept);

        Ok(())
    }

    #[test]
    fn test_median_threshold_odd_count() -> Result<()> {
        let input_seqs: FastaRecords = hash_map!(
            "A".to_string(): vec![b'A'; 5],
            "B".to_string(): vec![b'A'; 10],
            "C".to_string(): vec![b'A'; 15],
        );

        // Median length is 10.
        let (output, _) = filter_by_length(input_seqs, LengthThreshold::Median)?;

        assert_eq!(output.len(), 2);
        assert!(output.contains_key("B"));
        assert!(output.contains_key("C"));

        Ok(())
    }

    #[test]
    fn test_median_threshold_even_count() -> Result<()> {
        let input_seqs: FastaRecords = hash_map!(
            "A".to_string(): vec![b'A'; 5],
            "B".to_string(): vec![b'A'; 10],
            "C".to_string(): vec![b'A'; 20],
            "D".to_string(): vec![b'A'; 25],
        );

        // Median length is (10 + 20) / 2 = 15.
        let (output, _) = filter_by_length(input_seqs, LengthThreshold::Median)?;

        assert_eq!(output.len(), 2);
        assert!(output.contains_key("C"));
        assert!(output.contains_key("D"));

        Ok(())
    }

    #[test]
    fn test_mean_threshold() -> Result<()> {
        let input_seqs: FastaRecords = hash_map!(
            "A".to_string(): vec![b'A'; 5],
            "B".to_string(): vec![b'A'; 10],
            "C".to_string(): vec![b'A'; 15],
        );

        // Mean length is 10.
        let (output, _) = filter_by_length(input_seqs, LengthThreshold::Mean)?;

        assert_eq!(output.len(), 2);
        assert!(output.contains_key("B"));
        assert!(output.contains_key("C"));

        Ok(())
    }

    #[test]
    fn test_empty_input() {
        let input_seqs: FastaRecords = FastaRecords::new();
        assert!(filter_by_length(input_seqs, LengthThreshold::Fixed(10)).is_err());
    }
}
