use anyhow::{Context, Result};
use phf::{phf_map, phf_set};
use std::convert::TryInto;
use std::io::repeat;

pub const GAP_CHAR: u8 = b"-"[0];
const FRAMESHIFT_CHAR: u8 = b"X"[0];
const UNKNOWN_AA_CHAR: u8 = b"X"[0];
const INCOMPLETE_AA_CHAR: u8 = b"~"[0];

pub const DEFAULT_STOP_CHAR: u8 = b"*"[0];

static CODON_TABLE: phf::Map<&[u8; 3], &[u8; 1]> = phf_map! {
        b"TTT" => b"F",
        b"TTC" => b"F",
        b"TTA" => b"L",
        b"TTG" => b"L",
        b"CTT" => b"L",
        b"CTC" => b"L",
        b"CTA" => b"L",
        b"CTG" => b"L",
        b"ATT" => b"I",
        b"ATC" => b"I",
        b"ATA" => b"I",
        b"ATG" => b"M",
        b"GTT" => b"V",
        b"GTC" => b"V",
        b"GTA" => b"V",
        b"GTG" => b"V",
        b"TCT" => b"S",
        b"TCC" => b"S",
        b"TCA" => b"S",
        b"TCG" => b"S",
        b"CCT" => b"P",
        b"CCC" => b"P",
        b"CCA" => b"P",
        b"CCG" => b"P",
        b"ACT" => b"T",
        b"ACC" => b"T",
        b"ACA" => b"T",
        b"ACG" => b"T",
        b"GCT" => b"A",
        b"GCC" => b"A",
        b"GCA" => b"A",
        b"GCG" => b"A",
        b"TAT" => b"Y",
        b"TAC" => b"Y",
        b"CAT" => b"H",
        b"CAC" => b"H",
        b"CAA" => b"Q",
        b"CAG" => b"Q",
        b"AAT" => b"N",
        b"AAC" => b"N",
        b"AAA" => b"K",
        b"AAG" => b"K",
        b"GAT" => b"D",
        b"GAC" => b"D",
        b"GAA" => b"E",
        b"GAG" => b"E",
        b"TGT" => b"C",
        b"TGC" => b"C",
        b"TGG" => b"W",
        b"CGT" => b"R",
        b"CGC" => b"R",
        b"CGA" => b"R",
        b"CGG" => b"R",
        b"AGT" => b"S",
        b"AGC" => b"S",
        b"AGA" => b"R",
        b"AGG" => b"R",
        b"GGT" => b"G",
        b"GGC" => b"G",
        b"GGA" => b"G",
        b"GGG" => b"G",
        b"---" => b"-",
};

static STOP_CODONS: phf::Set<&[u8; 3]> = phf_set! {b"TAA", b"TAG", b"TGA"};

// Thanks https://cran.r-project.org/web/packages/MLMOI/vignettes/StandardAmbiguityCodes.html
static AMBIGUOUS_CODON_TABLE: phf::Map<&[u8; 3], &[u8; 1]> = phf_map! {
    b"GCN" =>  b"A",
    b"TGY"=> b"C",
    b"GAY" => b"D",
    b"GAR" => b"E",
    b"TTY" => b"F",
    b"GGN" => b"G",
    b"CAY" => b"H",
    b"ATH" => b"I",
    b"AAR" => b"K",
    b"YTR" => b"L",
    b"CTN" => b"L",
    b"AAY" => b"N",
    b"CCN" => b"P",
    b"CAR" => b"Q",
    b"CGN" => b"R",
    b"MGR" => b"R",
    b"TCN" => b"S",
    b"AGY" => b"S",
    b"ACN" => b"T",
    b"GTN" => b"V",
    b"TAY" => b"Y"


};

// https://en.wikipedia.org/wiki/International_Union_of_Pure_and_Applied_Chemistry#Amino_acid_and_nucleotide_base_codes
static AMBIGUOUS_CODON_AND_AA_TABLE: phf::Map<&[u8; 3], &[u8; 1]> = phf_map! {
    b"RAY" => b"B",
    b"SAR" => b"Z"
};

pub fn translate(
    dna_seq: &[u8],
    allow_ambiguities: bool,
    strip_gaps: bool,
    ignore_gap_codons: bool,
    drop_incomplete_codons: bool,
    alternate_stop_codon: Option<char>,
    alternate_frameshift_character: Option<char>,
    alternate_unknown_character: Option<char>,
    alternate_incomplete_character: Option<char>,
) -> Result<Vec<u8>> {
    let stop_codon_aa_char: u8 = match alternate_stop_codon {
        None => DEFAULT_STOP_CHAR,
        Some(character) => character as u8,
    };

    let frameshift_aa_char: u8 = match alternate_frameshift_character {
        None => FRAMESHIFT_CHAR,
        Some(character) => character as u8,
    };

    let unknown_aa_char: u8 = match alternate_unknown_character {
        None => UNKNOWN_AA_CHAR,
        Some(character) => character as u8,
    };

    let incomplete_aa_char: u8 = match alternate_incomplete_character {
        None => INCOMPLETE_AA_CHAR,
        Some(character) => character as u8,
    };

    let mut new_seq = dna_seq.to_vec();
    if strip_gaps {
        new_seq = new_seq
            .iter()
            .copied()
            .filter(|character| *character != GAP_CHAR)
            .collect();
    }

    let mut amino_acids = Vec::with_capacity(new_seq.len() / 3);
    for codon in new_seq.chunks(3) {
        // If the codon is not a multiple of 3, we will always want to replace it with an incomplete amino acid, so we don't need to
        // check anything else.

        if codon.len() != 3 {
            if !drop_incomplete_codons {
                log::debug!(
                    "The codon {:?} had a length of {} so we're adding a {:?}",
                    String::from_utf8(codon.to_vec())?,
                    codon.len(),
                    incomplete_aa_char as char
                );
                amino_acids.push(incomplete_aa_char);
            }
            continue;
        }
        let nt_triplet: [u8; 3] = codon
            .try_into()
            .expect("The codon should always be a triplet vector since we've checked for it.");

        if !strip_gaps {
            let num_gaps = nt_triplet.iter().filter(|char| **char == GAP_CHAR).count();
            if (num_gaps == 1) | (num_gaps == 2) {
                amino_acids.push(frameshift_aa_char);
                continue;
            }
        }
        let amino_acid;

        if CODON_TABLE.contains_key(&nt_triplet) {
            amino_acid = &CODON_TABLE[&nt_triplet][0];
        } else if allow_ambiguities && AMBIGUOUS_CODON_TABLE.contains_key(&nt_triplet) {
            amino_acid = &AMBIGUOUS_CODON_TABLE[&nt_triplet][0];
        } else if allow_ambiguities && AMBIGUOUS_CODON_AND_AA_TABLE.contains_key(&nt_triplet) {
            amino_acid = &AMBIGUOUS_CODON_AND_AA_TABLE[&nt_triplet][0];
        } else if STOP_CODONS.contains(&nt_triplet) {
            amino_acid = &stop_codon_aa_char;
        } else {
            log::debug!(
                "Could not find a suitable character for the codon {:?}",
                String::from_utf8(nt_triplet.to_vec())
            );
            amino_acid = &unknown_aa_char;
        }

        if ignore_gap_codons & (amino_acid.eq(&GAP_CHAR)) {
            continue;
        } else {
            amino_acids.push(amino_acid.clone());
        }
    }

    Ok(amino_acids)
}

pub fn main() {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_test() {
        let dna_seq = "ATGTTATAA";
        let expected_translation = "ML*";
        let translation = translate(
            dna_seq.as_bytes(),
            true,
            false,
            false,
            false,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        assert_eq!(expected_translation.as_bytes(), translation.as_slice());
    }

    #[test]
    fn strip_gaps_true() {
        let dna_seq = "ATGTTA-TAA";
        let expected_translation = "ML*";
        let translation = translate(
            dna_seq.as_bytes(),
            true,
            true,
            false,
            false,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        assert_eq!(expected_translation.as_bytes(), translation.as_slice());
    }

    #[test]
    fn strip_gaps_false() {
        let dna_seq = "ATGTTA-TAA";
        let expected_translation = "MLX~";
        let translation = translate(
            dna_seq.as_bytes(),
            true,
            false,
            false,
            false,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        assert_eq!(expected_translation.as_bytes(), translation.as_slice());
    }

    #[test]
    fn ignore_gap_codons() {
        let dna_seq = "ATGTTA---TAA";
        let expected_translation = "ML*";
        let translation = translate(
            dna_seq.as_bytes(),
            true,
            false,
            true,
            false,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        assert_eq!(expected_translation.as_bytes(), translation.as_slice());
    }

    #[test]
    fn test_ambiguity() {
        let test_cases = vec!["ATGTTACTNTAA", "NNNATGGGG", "ATGRAY---GTA"];
        let expected_outputs = vec!["MLL*", "?MG", "MB-V"];

        for (idx, test_case) in test_cases.iter().enumerate() {
            let expected_translation = expected_outputs[idx].as_bytes();
            let translation = translate(
                test_case.as_bytes(),
                true,
                false,
                false,
                false,
                None,
                None,
                None,
                None,
            )
            .unwrap();

            assert_eq!(expected_translation, translation.as_slice());
        }
    }

    #[test]
    fn test_alternate_stop_codon_char() {
        let translation_standard = translate(
            "ATGTTACTNTAA".as_bytes(),
            true,
            false,
            false,
            false,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        let translation_custom = translate(
            "ATGTTACTNTAA".as_bytes(),
            true,
            false,
            false,
            false,
            Some('X'),
            None,
            None,
            None,
        )
        .unwrap();

        assert_eq!("MLL*".as_bytes(), translation_standard.as_slice());
        assert_eq!("MLLX".as_bytes(), translation_custom.as_slice());
    }

    // TODO: Add more tests lol
}
