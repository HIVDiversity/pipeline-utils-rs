use std::collections::HashMap;
use std::iter::Iterator;
use std::path::{PathBuf};
use std::process::Output;
use bio::io::fasta;
use anyhow::{Result, Context, anyhow};
use bio::alignment::Alignment;
use nalgebra::DMatrix;
use colored::Colorize;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;
use bio::scores::blosum62;
use log::log;
use crate::tools;

const VERSION: &str = "0.1.0";


const GAP_CHAR: u8 = b"-"[0];
const FRAMESHIFT_CHAR: u8 = b"X"[0];
 const UNKNOWN_AA_CHAR: u8 = b"?"[0];

fn translate(dna_seq: &[u8], strip_gaps: bool, ignore_gap_codons: bool, codon_table: &HashMap<&[u8; 3], &[u8; 1]>) -> Result<Vec<u8>>{
    let mut new_seq = dna_seq.to_vec();
    if strip_gaps{
        new_seq = new_seq.iter().copied().filter(|character| *character != GAP_CHAR).collect();
    }

    let mut amino_acids = Vec::with_capacity(new_seq.len()/3);
    for codon in new_seq.chunks(3) {

        // If the codon is not a multiple of 3, we will always want to replace it with an incomplete amino acid, so we don't need to
        // check anything else.

        if codon.len() != 3 {
            log::warn!("The codon {:?} had a length of {} so we're adding a {:?}", String::from_utf8(codon.to_vec())?, codon.len(), UNKNOWN_AA_CHAR as char);
            amino_acids.push(FRAMESHIFT_CHAR);
        } else {
            let nt_triplet: [u8; 3] = codon.try_into().expect("The codon should always be a triplet vector since we've checked for it.");
            let amino_acid = codon_table
                .get(&nt_triplet)
                .with_context(|| format!("Couldn't find amino acid {:?}", String::from_utf8(nt_triplet.to_vec())))
                .unwrap_or(&&[UNKNOWN_AA_CHAR,]);

            amino_acids.push(amino_acid.clone()[0]);
        }
    }


    Ok(amino_acids)
}

// TODO: Move readfasta to the utils crate
fn read_fasta(fasta_file: &PathBuf) -> Result<Vec<Vec<u8>>>{
    let reader = fasta::Reader::from_file(fasta_file).expect("Could not open provided FASTA file.");
    let mut seqs: Vec<Vec<u8>> = Vec::new();

    for result in reader.records() {
        let record = result.expect("This record is invalid and failed to parse.");
        seqs.push(record.seq().to_vec());
    }

    Ok(seqs)

}

fn write_fasta(output_file: &PathBuf, seq_name: &str, seq: &Vec<u8>) -> Result<()>{
    let mut writer = fasta::Writer::to_file(output_file)?;
    writer.write(seq_name, None, seq)?;

    Ok(())
}

pub fn run(reference_file: &PathBuf,
           query_file: &PathBuf,
           output_file: &PathBuf,
           output_seq_name: &str,
           strip_gaps: bool,
           output_type: &String) -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;
    let codon_table: HashMap<&[u8; 3], &[u8; 1]> = HashMap::from([
        (b"TTT", b"F"),
        (b"TTC", b"F"),
        (b"TTA", b"L"),
        (b"TTG", b"L"),
        (b"CTT", b"L"),
        (b"CTC", b"L"),
        (b"CTA", b"L"),
        (b"CTG", b"L"),
        (b"ATT", b"I"),
        (b"ATC", b"I"),
        (b"ATA", b"I"),
        (b"ATG", b"M"),
        (b"GTT", b"V"),
        (b"GTC", b"V"),
        (b"GTA", b"V"),
        (b"GTG", b"V"),
        (b"TCT", b"S"),
        (b"TCC", b"S"),
        (b"TCA", b"S"),
        (b"TCG", b"S"),
        (b"CCT", b"P"),
        (b"CCC", b"P"),
        (b"CCA", b"P"),
        (b"CCG", b"P"),
        (b"ACT", b"T"),
        (b"ACC", b"T"),
        (b"ACA", b"T"),
        (b"ACG", b"T"),
        (b"GCT", b"A"),
        (b"GCC", b"A"),
        (b"GCA", b"A"),
        (b"GCG", b"A"),
        (b"TAT", b"Y"),
        (b"TAC", b"Y"),
        (b"CAT", b"H"),
        (b"CAC", b"H"),
        (b"CAA", b"Q"),
        (b"CAG", b"Q"),
        (b"AAT", b"N"),
        (b"AAC", b"N"),
        (b"AAA", b"K"),
        (b"AAG", b"K"),
        (b"GAT", b"D"),
        (b"GAC", b"D"),
        (b"GAA", b"E"),
        (b"GAG", b"E"),
        (b"TGT", b"C"),
        (b"TGC", b"C"),
        (b"TGG", b"W"),
        (b"CGT", b"R"),
        (b"CGC", b"R"),
        (b"CGA", b"R"),
        (b"CGG", b"R"),
        (b"AGT", b"S"),
        (b"AGC", b"S"),
        (b"AGA", b"R"),
        (b"AGG", b"R"),
        (b"GGT", b"G"),
        (b"GGC", b"G"),
        (b"GGA", b"G"),
        (b"GGG", b"G"),
        (b"TAA", b"*"),
        (b"TAG", b"*"),
        (b"TGA", b"*")]);

    let reference_read = read_fasta(reference_file)?;
    let reference = reference_read[0].as_slice();

    let query_read = read_fasta(query_file)?;
    let query = query_read[0].as_slice();

    let ref_aa = translate(reference, false, false, &codon_table)?;
    let ref_aa_slice = ref_aa.as_slice();

    let mut aligner = Aligner::with_capacity(query.len() / 3, ref_aa.len(), -5, -1, bio::scores::blosum62);
    let mut best_score = 0;
    let mut best_frame = 0;
    let mut best_translation: Vec<u8> = Vec::with_capacity(query.len() / 3);
    let mut best_alignment: Alignment = Default::default();

    for frame in 0..3 {
        log::info!("Translating query in frame {:?}", frame+1);
        let cons_aa = translate(&query[frame..], true, true, &codon_table)?;
        let alignment = aligner.local(cons_aa.as_slice(), ref_aa_slice);

        log::info!("Alignment with query in frame {:?} gave a score of {:?}", frame, alignment.score);
        if alignment.score > best_score {
            best_score = alignment.score;
            best_frame = frame;
            best_translation = cons_aa[alignment.xstart..alignment.xend].to_vec();
            best_alignment = alignment;
        }
    }

    log::info!("Choosing translation in frame {:?} with score {:?}:\n{:?}", best_frame, best_score, String::from_utf8(best_translation.clone())? );
    log::info!("With this best alignment, we trimmed the query sequence (AA) from position {} to {}",best_alignment.xstart, best_alignment.xend);

    if output_type == "AA"{
        log::info!("Writing trimmed query amino acid sequence to {:?}", output_file);
        write_fasta(output_file, output_seq_name, &best_translation)?;

    }else{
        if output_type!= "NT"{
            log::error!("Unrecognized output type, outputting NT sequence")
        }

        let trim_nt_start = best_alignment.xstart * 3;
        let trim_nt_end = best_alignment.xend * 3;
        let trimmed_nt = query[trim_nt_start..trim_nt_end].to_vec();
        write_fasta(output_file, output_seq_name, &trimmed_nt)?;
        log::info!("Outputting NT sequence to {:?}", output_file);
    }

//
    //
    // println!("{}", alignment.pretty(cons_aa.as_slice(), ref_aa.as_slice(), 160));
    // println!("{:?}", alignment.score);

    Ok(())
}