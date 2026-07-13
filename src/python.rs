use crate::tools;
use crate::tools::get_consensus::AmbiguityMode;
use crate::utils::fasta_utils::FastaRecords;
use crate::utils::translate::TranslationOptions;

fn to_pyerr(e: anyhow::Error) -> pyo3::PyErr {
    pyo3::exceptions::PyRuntimeError::new_err(e.to_string())
}

fn dict_to_records(seqs: std::collections::HashMap<String, String>) -> FastaRecords {
    seqs.into_iter()
        .map(|(name, seq)| {
            let mut seq = seq.into_bytes();
            seq.make_ascii_uppercase();
            (name, seq)
        })
        .collect()
}

fn records_to_dict(
    records: FastaRecords,
) -> pyo3::PyResult<std::collections::HashMap<String, String>> {
    records
        .into_iter()
        .map(|(name, seq)| {
            String::from_utf8(seq)
                .map(|s| (name, s))
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
        })
        .collect()
}

#[pyo3::pymodule]
mod purs {
    use super::*;
    use pyo3::prelude::*;
    use std::collections::HashMap;

    #[pyfunction]
    fn get_consensus(seqs: Vec<String>, ambiguity_mode: String) -> PyResult<String> {
        let msa: Vec<Vec<u8>> = seqs.into_iter().map(String::into_bytes).collect();
        let mode = match ambiguity_mode.as_str() {
            "IUPAC" => AmbiguityMode::UseIUPAC,
            "First" => AmbiguityMode::First,
            "Random" => AmbiguityMode::Random,
            "MarkN" => AmbiguityMode::MarkN,
            other => {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "Unknown ambiguity mode: {other:?}. Expected one of \"IUPAC\", \"First\", \"Random\", \"MarkN\"."
                )));
            }
        };

        let matrix = tools::get_consensus::sequences_to_matrix(&msa).map_err(to_pyerr)?;
        let consensus = tools::get_consensus::build_consensus(&matrix, mode).map_err(to_pyerr)?;

        String::from_utf8(consensus)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    #[pyfunction]
    #[pyo3(signature = (
        seqs,
        unknown_aa='X',
        stop_aa='*',
        incomplete_aa='?',
        frameshift_aa='X',
        reading_frame=0,
        allow_ambiguities=true,
        strip_gaps=false,
        ignore_gap_codons=false,
        drop_incomplete_codons=true,
    ))]
    fn translate(
        seqs: HashMap<String, String>,
        unknown_aa: char,
        stop_aa: char,
        incomplete_aa: char,
        frameshift_aa: char,
        reading_frame: usize,
        allow_ambiguities: bool,
        strip_gaps: bool,
        ignore_gap_codons: bool,
        drop_incomplete_codons: bool,
    ) -> PyResult<HashMap<String, String>> {
        let options = TranslationOptions {
            unknown_aa: unknown_aa as u8,
            stop_aa: stop_aa as u8,
            incomplete_aa: incomplete_aa as u8,
            frameshift_aa: frameshift_aa as u8,
            reading_frame,
            allow_ambiguities,
            strip_gaps,
            ignore_gap_codons,
            drop_incomplete_codons,
        };

        let translated = tools::translate::translate_records(dict_to_records(seqs), &options)
            .map_err(to_pyerr)?;
        records_to_dict(translated)
    }

    #[pyfunction]
    fn reverse_translate(
        aa_seqs: HashMap<String, String>,
        nt_seqs: HashMap<String, String>,
    ) -> PyResult<HashMap<String, String>> {
        let result = tools::reverse_translate::process_sequences(
            dict_to_records(aa_seqs),
            dict_to_records(nt_seqs),
        )
            .map_err(to_pyerr)?;
        records_to_dict(result)
    }

    #[pyfunction]
    fn replace_ambiguities(
        seqs: HashMap<String, String>,
        seed: u64,
    ) -> PyResult<HashMap<String, String>> {
        let result =
            tools::replace_ambiguities::replace_ambiguities_records(dict_to_records(seqs), seed)
                .map_err(to_pyerr)?;
        records_to_dict(result)
    }

    #[pyfunction]
    #[pyo3(signature = (seqs, include_stop_codon=true))]
    fn trim_after_stop_codon(
        seqs: HashMap<String, String>,
        include_stop_codon: bool,
    ) -> PyResult<HashMap<String, String>> {
        let result =
            tools::trim_after_stop_codon::process_file(dict_to_records(seqs), include_stop_codon)
                .map_err(to_pyerr)?;
        records_to_dict(result)
    }

    #[pyfunction]
    #[pyo3(signature = (seqs, length=None, median=false, mean=false, exclude_gaps=true))]
    fn filter_by_length(
        seqs: HashMap<String, String>,
        length: Option<usize>,
        median: bool,
        mean: bool,
        exclude_gaps: bool,
    ) -> PyResult<(
        HashMap<String, String>,
        HashMap<String, String>,
        Vec<(String, usize, bool)>,
    )> {
        let threshold = match (length, median, mean) {
            (Some(l), false, false) => tools::filter_by_length::LengthThreshold::Fixed(l),
            (None, true, false) => tools::filter_by_length::LengthThreshold::Median,
            (None, false, true) => tools::filter_by_length::LengthThreshold::Mean,
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "Exactly one of length, median, or mean must be specified.",
                ));
            }
        };

        let range = tools::filter_by_length::LengthRange {
            center: threshold,
            min_tolerance: None,
            max_tolerance: None,
        };
        let (records, rejected, report) =
            tools::filter_by_length::filter_by_length(dict_to_records(seqs), range, exclude_gaps)
                .map_err(to_pyerr)?;
        let report_rows = report
            .into_iter()
            .map(|r| (r.seq_name, r.length, r.kept))
            .collect();
        Ok((
            records_to_dict(records)?,
            records_to_dict(rejected)?,
            report_rows,
        ))
    }

    #[pyfunction]
    #[pyo3(signature = (seqs, seq_prefix, strip_gaps=false))]
    fn collapse(
        seqs: HashMap<String, String>,
        seq_prefix: String,
        strip_gaps: bool,
    ) -> PyResult<(HashMap<String, String>, HashMap<String, Vec<String>>)> {
        let collapsed = tools::collapse::collapse_sequences(dict_to_records(seqs), strip_gaps)
            .map_err(to_pyerr)?;
        let (records, name_mapping) =
            tools::collapse::build_collapsed_output(collapsed, &seq_prefix);
        Ok((records_to_dict(records)?, name_mapping))
    }

    #[pyfunction]
    #[pyo3(signature = (seqs, name_mapping, include_missing=false))]
    fn expand(
        seqs: HashMap<String, String>,
        name_mapping: HashMap<String, Vec<String>>,
        include_missing: bool,
    ) -> PyResult<HashMap<String, String>> {
        let expanded = tools::expand::uncollapse_sequences(
            dict_to_records(seqs),
            name_mapping,
            include_missing,
        )
            .map_err(to_pyerr)?;
        records_to_dict(expanded)
    }

    #[pyfunction]
    #[pyo3(signature = (seqs, start_kmers=None, end_kmers=None))]
    fn filter_by_kmer(
        seqs: HashMap<String, String>,
        start_kmers: Option<Vec<String>>,
        end_kmers: Option<Vec<String>>,
    ) -> PyResult<(
        HashMap<String, String>,
        HashMap<String, String>,
        Vec<(String, Option<bool>, Option<bool>, bool)>,
    )> {
        if start_kmers.is_none() && end_kmers.is_none() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "At least one of start_kmers or end_kmers must be specified.",
            ));
        }

        fn to_bytes(kmers: Option<Vec<String>>) -> Option<Vec<Vec<u8>>> {
            kmers.map(|list| {
                list.into_iter()
                    .map(|k| k.to_ascii_uppercase().into_bytes())
                    .collect()
            })
        }

        let start_kmers = to_bytes(start_kmers);
        let end_kmers = to_bytes(end_kmers);

        let (records, rejected, report) = tools::filter_by_kmer::filter_by_kmer(
            dict_to_records(seqs),
            start_kmers.as_deref(),
            end_kmers.as_deref(),
        )
            .map_err(to_pyerr)?;

        let report_rows = report
            .into_iter()
            .map(|r| (r.seq_name, r.start_match, r.end_match, r.kept))
            .collect();

        Ok((
            records_to_dict(records)?,
            records_to_dict(rejected)?,
            report_rows,
        ))
    }
}
