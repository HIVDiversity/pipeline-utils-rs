use crate::utils::fasta_utils::{load_fasta, write_fasta_sequences, FastaRecords};
use anyhow::{bail, Result};
use colored::Colorize;
use std::fmt;
use std::path::PathBuf;
use std::str::FromStr;
use crate::utils::codon_tables::GAP_CHAR;

pub enum LengthThreshold {
    Fixed(usize),
    Median,
    Mean,
}

/// A margin around a center length, expressed either as an absolute number of
/// bases or as a percentage of the center value (e.g. "20" vs "20%").
#[derive(Debug, Clone, Copy)]
pub enum Tolerance {
    Absolute(f64),
    Percent(f64),
}

impl Tolerance {
    /// Resolve this tolerance to an absolute number of bases, given the center
    /// length it is relative to.
    fn resolve(&self, center: f64) -> f64 {
        match self {
            Tolerance::Absolute(bases) => *bases,
            Tolerance::Percent(pct) => center * pct / 100.0,
        }
    }
}

#[derive(Debug)]
pub struct ToleranceParseError(String);

impl fmt::Display for ToleranceParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for ToleranceParseError {}

impl FromStr for Tolerance {
    type Err = ToleranceParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.strip_suffix('%') {
            Some(pct) => pct
                .parse::<f64>()
                .map(Tolerance::Percent)
                .map_err(|e| ToleranceParseError(format!("invalid percent tolerance {s:?}: {e}"))),
            None => s
                .parse::<f64>()
                .map(Tolerance::Absolute)
                .map_err(|e| ToleranceParseError(format!("invalid tolerance {s:?}: {e}"))),
        }
    }
}

/// A length filter: a center (fixed/median/mean), optionally widened below and/or
/// above by a tolerance, producing an inclusive `[min, max]` acceptance range.
pub struct LengthRange {
    pub center: LengthThreshold,
    pub min_tolerance: Option<Tolerance>,
    pub max_tolerance: Option<Tolerance>,
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

pub(crate) struct FilterReportRow {
    pub(crate) seq_name: String,
    pub(crate) length: usize,
    pub(crate) kept: bool,
}

pub(crate) fn filter_by_length(
    sequences: FastaRecords,
    range: LengthRange,
    exclude_gaps: bool,
) -> Result<(FastaRecords, FastaRecords, Vec<FilterReportRow>)> {
    if sequences.is_empty() {
        bail!("No sequences were provided.")
    }

    let lengths: Vec<usize> = sequences.values().map(|seq| seq.len()).collect();
    let center_value = threshold_value(&lengths, &range.center);
    let lower_bound = range
        .min_tolerance
        .map_or(center_value, |t| (center_value - t.resolve(center_value)).max(0.0));
    let upper_bound = range.max_tolerance.map(|t| center_value + t.resolve(center_value));

    log::info!(
        "Center length: {center_value}, acceptance range: [{lower_bound}, {}]",
        upper_bound.map_or("inf".to_string(), |u| u.to_string())
    );

    let mut kept_sequences = FastaRecords::with_capacity(sequences.len());
    let mut rejected_sequences = FastaRecords::new();
    let mut report_rows = Vec::with_capacity(sequences.len());

    for (seq_name, seq) in sequences {
        let length = match exclude_gaps {
            true => { seq.iter().filter(|x| { **x != GAP_CHAR }).count() }
            false => { seq.len() }
        };

        let length_f = length as f64;
        let kept = length_f >= lower_bound && upper_bound.is_none_or(|u| length_f <= u);

        report_rows.push(FilterReportRow {
            seq_name: seq_name.clone(),
            length,
            kept,
        });

        if kept {
            kept_sequences.insert(seq_name, seq);
        } else {
            rejected_sequences.insert(seq_name, seq);
        }
    }

    report_rows.sort_unstable_by(|a, b| a.seq_name.cmp(&b.seq_name));

    Ok((kept_sequences, rejected_sequences, report_rows))
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
    rejected_seq_output: Option<&PathBuf>,
    range: LengthRange,
    exclude_gaps: bool,
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
    let (kept_sequences, rejected_sequences, report_rows) = filter_by_length(sequences, range, exclude_gaps)?;

    write_fasta_sequences(output_file, &kept_sequences)?;

    if let Some(rejected_seq_output) = rejected_seq_output {
        log::info!("Writing rejected sequences to {:?}", rejected_seq_output);
        write_fasta_sequences(rejected_seq_output, &rejected_sequences)?;
    }

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

    fn center_only(center: LengthThreshold) -> LengthRange {
        LengthRange {
            center,
            min_tolerance: None,
            max_tolerance: None,
        }
    }

    #[test]
    fn test_fixed_threshold() -> Result<()> {
        let input_seqs: FastaRecords = hash_map!(
            "A".to_string(): vec![b'A'; 5],
            "B".to_string(): vec![b'A'; 10],
            "C".to_string(): vec![b'A'; 15],
        );

        let (output, rejected, report) =
            filter_by_length(input_seqs, center_only(LengthThreshold::Fixed(10)), false)?;

        assert_eq!(output.len(), 2);
        assert!(output.contains_key("B"));
        assert!(output.contains_key("C"));
        assert!(!output.contains_key("A"));

        assert_eq!(rejected.len(), 1);
        assert!(rejected.contains_key("A"));

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
    fn test_gap_exclusion() -> Result<()> {
        let input_seqs: FastaRecords = hash_map!(
            "A".to_string(): vec![b'A', b'T', b'-', b'-', b'G'],
            "B".to_string(): vec![b'A'; 10],
        );

        let (output, _, _) =
            filter_by_length(input_seqs.clone(), center_only(LengthThreshold::Fixed(4)), false)?;

        assert_eq!(output.len(), 2);
        assert!(output.contains_key("B"));
        assert!(output.contains_key("A"));

        let (output, _, _) =
            filter_by_length(input_seqs, center_only(LengthThreshold::Fixed(4)), true)?;
        assert_eq!(output.len(), 1);
        assert!(output.contains_key("B"));
        assert!(!output.contains_key("A"));


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
        let (output, _, _) = filter_by_length(input_seqs, center_only(LengthThreshold::Median), false)?;

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
        let (output, _, _) = filter_by_length(input_seqs, center_only(LengthThreshold::Median), false)?;

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
        let (output, _, _) = filter_by_length(input_seqs, center_only(LengthThreshold::Mean), false)?;

        assert_eq!(output.len(), 2);
        assert!(output.contains_key("B"));
        assert!(output.contains_key("C"));

        Ok(())
    }

    #[test]
    fn test_empty_input() {
        let input_seqs: FastaRecords = FastaRecords::new();
        assert!(filter_by_length(input_seqs, center_only(LengthThreshold::Fixed(10)), false).is_err());
    }

    #[test]
    fn test_min_tolerance_absolute() -> Result<()> {
        let input_seqs: FastaRecords = hash_map!(
            "A".to_string(): vec![b'A'; 75],
            "B".to_string(): vec![b'A'; 80],
            "C".to_string(): vec![b'A'; 100],
        );

        // length 100, min-tolerance 20 -> keep [80, inf)
        let (output, rejected, _) = filter_by_length(
            input_seqs,
            LengthRange {
                center: LengthThreshold::Fixed(100),
                min_tolerance: Some(Tolerance::Absolute(20.0)),
                max_tolerance: None,
            },
            false,
        )?;

        assert_eq!(output.len(), 2);
        assert!(output.contains_key("B"));
        assert!(output.contains_key("C"));
        assert_eq!(rejected.len(), 1);
        assert!(rejected.contains_key("A"));

        Ok(())
    }

    #[test]
    fn test_symmetric_percent_tolerance_around_median() -> Result<()> {
        let input_seqs: FastaRecords = hash_map!(
            "A".to_string(): vec![b'A'; 50],
            "B".to_string(): vec![b'A'; 100],
            "C".to_string(): vec![b'A'; 105],
            "D".to_string(): vec![b'A'; 150],
        );

        // Median length is (100 + 105) / 2 = 102.5, 10% tolerance -> keep [92.25, 112.75]
        let (output, rejected, _) = filter_by_length(
            input_seqs,
            LengthRange {
                center: LengthThreshold::Median,
                min_tolerance: Some(Tolerance::Percent(10.0)),
                max_tolerance: Some(Tolerance::Percent(10.0)),
            }, false,
        )?;

        assert_eq!(output.len(), 2);
        assert!(output.contains_key("B"));
        assert!(output.contains_key("C"));
        assert_eq!(rejected.len(), 2);
        assert!(rejected.contains_key("A"));
        assert!(rejected.contains_key("D"));

        Ok(())
    }
}
