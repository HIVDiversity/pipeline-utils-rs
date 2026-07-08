use crate::utils::codon_tables::AMBIGUOUS_NT_LOOKUP;
use crate::utils::fasta_utils::{load_fasta, write_fasta_sequences, FastaRecords};
use anyhow::{bail, Result};
use colored::Colorize;
use std::path::PathBuf;

/// Expand a single base to the set of concrete bases it can represent (a singleton set for
/// a concrete A/C/G/T, or the IUPAC ambiguity expansion for an ambiguity code).
fn expand_base(base: u8) -> Vec<u8> {
    match AMBIGUOUS_NT_LOOKUP.get(&[base]) {
        Some(set) => set.iter().map(|b| b[0]).collect(),
        None => vec![base],
    }
}

/// Two bases are compatible if the sets of concrete bases they can represent intersect, so
/// an ambiguity code in either the query k-mer or the sequence matches any base it represents.
pub(crate) fn bases_compatible(query: u8, seq: u8) -> bool {
    let query_set = expand_base(query);
    let seq_set = expand_base(seq);
    query_set.iter().any(|q| seq_set.contains(q))
}

pub(crate) fn matches_kmer_at_start(seq: &[u8], kmer: &[u8]) -> bool {
    seq.len() >= kmer.len()
        && seq
            .iter()
            .zip(kmer.iter())
            .all(|(&s, &k)| bases_compatible(k, s))
}

pub(crate) fn matches_kmer_at_end(seq: &[u8], kmer: &[u8]) -> bool {
    seq.len() >= kmer.len()
        && seq[seq.len() - kmer.len()..]
            .iter()
            .zip(kmer.iter())
            .all(|(&s, &k)| bases_compatible(k, s))
}

pub(crate) struct FilterReportRow {
    pub(crate) seq_name: String,
    pub(crate) start_match: Option<bool>,
    pub(crate) end_match: Option<bool>,
    pub(crate) kept: bool,
}

pub(crate) fn filter_by_kmer(
    sequences: FastaRecords,
    start_kmers: Option<&[Vec<u8>]>,
    end_kmers: Option<&[Vec<u8>]>,
) -> Result<(FastaRecords, FastaRecords, Vec<FilterReportRow>)> {
    if sequences.is_empty() {
        bail!("No sequences were provided.")
    }

    let mut kept_sequences = FastaRecords::with_capacity(sequences.len());
    let mut rejected_sequences = FastaRecords::new();
    let mut report_rows = Vec::with_capacity(sequences.len());

    for (seq_name, seq) in sequences {
        let start_match =
            start_kmers.map(|kmers| kmers.iter().any(|k| matches_kmer_at_start(&seq, k)));
        let end_match =
            end_kmers.map(|kmers| kmers.iter().any(|k| matches_kmer_at_end(&seq, k)));

        let kept = start_match.unwrap_or(true) && end_match.unwrap_or(true);

        report_rows.push(FilterReportRow {
            seq_name: seq_name.clone(),
            start_match,
            end_match,
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

fn fmt_match(m: Option<bool>) -> String {
    match m {
        Some(b) => b.to_string(),
        None => "n/a".to_string(),
    }
}

fn write_report(report_file: &PathBuf, rows: &[FilterReportRow]) -> Result<()> {
    let mut writer = csv::Writer::from_path(report_file)?;
    writer.write_record(["seq_name", "start_match", "end_match", "kept"])?;

    for row in rows {
        writer.write_record([
            row.seq_name.as_str(),
            fmt_match(row.start_match).as_str(),
            fmt_match(row.end_match).as_str(),
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
    start_kmers: Option<&[Vec<u8>]>,
    end_kmers: Option<&[Vec<u8>]>,
) -> Result<()> {
    log::info!(
        "{}",
        format!(
            "This is 'filter-by-kmer' version {}",
            env!("CARGO_PKG_VERSION")
        )
        .bold()
        .bright_yellow()
    );

    log::info!("Reading input file {:?}", input_file);
    let sequences = load_fasta(input_file)?;
    let (kept_sequences, rejected_sequences, report_rows) =
        filter_by_kmer(sequences, start_kmers, end_kmers)?;

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

    #[test]
    fn test_bases_compatible_exact_match() {
        assert!(bases_compatible(b'A', b'A'));
        assert!(!bases_compatible(b'A', b'C'));
    }

    #[test]
    fn test_bases_compatible_ambiguity_in_query() {
        // N in the query k-mer should match any concrete sequence base.
        assert!(bases_compatible(b'N', b'A'));
        assert!(bases_compatible(b'N', b'T'));
        // R (A or G) should match A and G but not C or T.
        assert!(bases_compatible(b'R', b'A'));
        assert!(bases_compatible(b'R', b'G'));
        assert!(!bases_compatible(b'R', b'C'));
    }

    #[test]
    fn test_bases_compatible_ambiguity_in_sequence() {
        // An ambiguity code in the sequence should match a concrete query base it represents.
        assert!(bases_compatible(b'A', b'N'));
        assert!(bases_compatible(b'A', b'R'));
        assert!(!bases_compatible(b'C', b'R'));
    }

    #[test]
    fn test_bases_compatible_two_ambiguity_codes() {
        // R = {A, G}, S = {C, G} -> overlap at G.
        assert!(bases_compatible(b'R', b'S'));
        // R = {A, G}, Y = {C, T} -> no overlap.
        assert!(!bases_compatible(b'R', b'Y'));
    }

    #[test]
    fn test_matches_kmer_at_start() {
        assert!(matches_kmer_at_start(b"ATGACG", b"ATG"));
        assert!(!matches_kmer_at_start(b"GTGACG", b"ATG"));
        // Sequence shorter than the k-mer always fails.
        assert!(!matches_kmer_at_start(b"AT", b"ATG"));
        // Ambiguity code in the sequence matches.
        assert!(matches_kmer_at_start(b"NTGACG", b"ATG"));
    }

    #[test]
    fn test_matches_kmer_at_end() {
        assert!(matches_kmer_at_end(b"GACTAA", b"TAA"));
        assert!(!matches_kmer_at_end(b"GACTAC", b"TAA"));
        assert!(!matches_kmer_at_end(b"AA", b"TAA"));
        assert!(matches_kmer_at_end(b"GACTAN", b"TAA"));
    }

    #[test]
    fn test_filter_by_kmer_start_only() -> Result<()> {
        let sequences: FastaRecords = hash_map!(
            "A".to_string(): b"ATGACGT".to_vec(),
            "B".to_string(): b"GTGACGT".to_vec(),
        );

        let start_kmers = vec![b"ATG".to_vec()];
        let (kept, rejected, report) =
            filter_by_kmer(sequences, Some(&start_kmers), None)?;

        assert_eq!(kept.len(), 1);
        assert!(kept.contains_key("A"));
        assert_eq!(rejected.len(), 1);
        assert!(rejected.contains_key("B"));

        assert_eq!(report.len(), 2);
        for row in &report {
            assert!(row.end_match.is_none());
        }

        Ok(())
    }

    #[test]
    fn test_filter_by_kmer_end_only() -> Result<()> {
        let sequences: FastaRecords = hash_map!(
            "A".to_string(): b"ATGACGTAA".to_vec(),
            "B".to_string(): b"ATGACGTAC".to_vec(),
        );

        let end_kmers = vec![b"TAA".to_vec(), b"TAG".to_vec(), b"TGA".to_vec()];
        let (kept, rejected, _) = filter_by_kmer(sequences, None, Some(&end_kmers))?;

        assert_eq!(kept.len(), 1);
        assert!(kept.contains_key("A"));
        assert_eq!(rejected.len(), 1);
        assert!(rejected.contains_key("B"));

        Ok(())
    }

    #[test]
    fn test_filter_by_kmer_start_and_end() -> Result<()> {
        let sequences: FastaRecords = hash_map!(
            // Passes both checks.
            "A".to_string(): b"ATGACGTAA".to_vec(),
            // Fails end check only.
            "B".to_string(): b"ATGACGTAC".to_vec(),
            // Fails start check only.
            "C".to_string(): b"GTGACGTAA".to_vec(),
        );

        let start_kmers = vec![b"ATG".to_vec()];
        let end_kmers = vec![b"TAA".to_vec(), b"TAG".to_vec(), b"TGA".to_vec()];
        let (kept, rejected, _) =
            filter_by_kmer(sequences, Some(&start_kmers), Some(&end_kmers))?;

        assert_eq!(kept.len(), 1);
        assert!(kept.contains_key("A"));
        assert_eq!(rejected.len(), 2);
        assert!(rejected.contains_key("B"));
        assert!(rejected.contains_key("C"));

        Ok(())
    }

    #[test]
    fn test_filter_by_kmer_fails_all_candidates() -> Result<()> {
        let sequences: FastaRecords = hash_map!(
            "A".to_string(): b"ATGACGTCC".to_vec(),
        );

        let end_kmers = vec![b"TAA".to_vec(), b"TAG".to_vec(), b"TGA".to_vec()];
        let (kept, rejected, _) = filter_by_kmer(sequences, None, Some(&end_kmers))?;

        assert_eq!(kept.len(), 0);
        assert_eq!(rejected.len(), 1);

        Ok(())
    }

    #[test]
    fn test_filter_by_kmer_empty_input() {
        let sequences: FastaRecords = FastaRecords::new();
        let start_kmers = vec![b"ATG".to_vec()];
        assert!(filter_by_kmer(sequences, Some(&start_kmers), None).is_err());
    }
}
