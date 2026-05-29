use crate::utils::fasta_utils::{FastaRecords, load_fasta};
use crate::utils::translate::GAP_CHAR;
use anyhow::{Context, Result};
use bio::io::fasta;
use colored::Colorize;
use std::collections::HashMap;
use std::path::PathBuf;
