use crate::utils::fasta_utils::{load_fasta, write_fasta_sequences, FastaRecords};
use anyhow::Result;

use colored::Colorize;
use polars::prelude::LazyFrame;
use polars::prelude::*;
use std::collections::HashMap;
use std::path::PathBuf;

fn read_fasta_into_lazyframe(fasta_file: &PathBuf) -> Result<LazyFrame> {
    let records = load_fasta(fasta_file)?;
    let names = records.keys().cloned().collect::<Vec<String>>();
    let sequences = records
        .values()
        .cloned()
        .map(|seq_vec| String::from_utf8(seq_vec).unwrap())
        .collect::<Vec<String>>();
    let fasta_df = df![
        "seq_name" => names,
        "seq_record" => sequences
    ]?;

    Ok(fasta_df.lazy())
}

fn write_dataframe_to_fasta(seq_df: DataFrame, output_file: &PathBuf) -> Result<()> {
    let names = seq_df["query"].clone().take_materialized_series();
    let sequences = seq_df["new_seq_rec"].clone().take_materialized_series();
    let mut fasta_seqs: FastaRecords = HashMap::with_capacity(names.len());

    names
        .iter()
        .zip(sequences.iter())
        .for_each(|(name, sequence)| {
            fasta_seqs.insert(
                name.get_str().unwrap().to_string(),
                sequence.get_str().unwrap().as_bytes().to_vec(),
            );
        });

    write_fasta_sequences(output_file, &fasta_seqs)?;

    Ok(())
}

pub fn run(
    input_file: &PathBuf,
    paf_file: &PathBuf,
    prepend: &Option<String>,
    output_dir: &PathBuf,
) -> Result<()> {
    log::info!(
        "{}",
        format!(
            "This is process-miniprot version {}",
            env!("CARGO_PKG_VERSION")
        )
        .bold()
        .bright_green()
    );

    let paf_schema: Schema = Schema::from_iter(vec![
        Field::new("ref_name".into(), DataType::String),
        Field::new("ref_len".into(), DataType::Int32),
        Field::new("ref_start".into(), DataType::Int32),
        Field::new("ref_end".into(), DataType::Int32),
        Field::new("strand".into(), DataType::String),
        Field::new("query".into(), DataType::String),
        Field::new("query_len".into(), DataType::Int32),
        Field::new("query_start".into(), DataType::Int32),
        Field::new("query_end".into(), DataType::Int32),
        Field::new("matches_nt".into(), DataType::Int32),
        Field::new("nt_excl_introns".into(), DataType::Int32),
        Field::new("qual".into(), DataType::Int32),
        Field::new("AS".into(), DataType::String),
        Field::new("ms".into(), DataType::String),
        Field::new("np".into(), DataType::String),
        Field::new("fs".into(), DataType::String),
        Field::new("st".into(), DataType::String),
        Field::new("da".into(), DataType::String),
        Field::new("do".into(), DataType::String),
        Field::new("cg".into(), DataType::String),
        Field::new("cs".into(), DataType::String),
    ]);

    let paf_df: LazyFrame = LazyCsvReader::new(paf_file.to_str().unwrap().into())
        .with_schema(Some(Arc::from(paf_schema)))
        .with_comment_prefix(Some("#".into()))
        .with_separator("\t".as_bytes()[0])
        .with_has_header(false)
        .finish()?;

    let seq_df = read_fasta_into_lazyframe(&input_file)?;

    let query_col = col("query");
    let query_start_col = col("query_start");
    let query_end_col = col("query_end");
    let seq_name_col = col("seq_name");
    let seq_record_col = col("seq_record");
    let new_seq_col = col("new_seq_rec");

    let trimmed_seq_df = paf_df
        .select([
            query_col.clone(),
            query_start_col.clone(),
            query_end_col.clone(),
        ])
        .group_by([query_col.clone()])
        .agg([
            query_start_col.clone().mode(false).first(),
            query_end_col.clone().mode(false).first(),
        ])
        .join(
            seq_df,
            [query_col.clone()],
            [seq_name_col.clone()],
            JoinArgs::new(JoinType::Inner),
        )
        .with_columns([seq_record_col
            .str()
            .slice(
                query_start_col.clone(),
                query_end_col.clone() - query_start_col.clone(),
            )
            .alias("new_seq_rec")]);
    write_dataframe_to_fasta(trimmed_seq_df.collect()?, &output_dir)?;

    Ok(())
}
