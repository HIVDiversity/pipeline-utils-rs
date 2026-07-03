use crate::utils::codon_tables::{
    AMBIGUOUS_CODON_AND_AA_TABLE, AMBIGUOUS_CODON_TABLE, AMBIGUOUS_NT_LOOKUP, CODON_TABLE,
    GAP_CHAR, STOP_CODONS,
};
use anyhow::Result;
use itertools::Itertools;
use std::collections::HashSet;
use std::convert::TryInto;
use std::fmt;

#[derive(Clone, Copy)]
pub struct TranslationOptions {
    pub unknown_aa: u8,
    pub stop_aa: u8,
    pub incomplete_aa: u8,
    pub frameshift_aa: u8,
    pub reading_frame: usize,
    pub allow_ambiguities: bool,
    pub strip_gaps: bool,
    pub ignore_gap_codons: bool,
    pub drop_incomplete_codons: bool,
}

impl Default for TranslationOptions {
    fn default() -> Self {
        Self {
            unknown_aa: b'X',
            frameshift_aa: b'X',
            stop_aa: b'*',
            incomplete_aa: b'?',
            reading_frame: 0,
            allow_ambiguities: true,
            strip_gaps: false,
            ignore_gap_codons: false,
            drop_incomplete_codons: true,
        }
    }
}

impl fmt::Display for TranslationOptions {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Translation Options: {{\n\t")?;
        write!(f, "unknown_aa_char: {:?}\n\t", self.unknown_aa as char)?;
        write!(
            f,
            "frameshift_aa_char: {:?}\n\t",
            self.frameshift_aa as char
        )?;
        write!(f, "stop_aa_char: {:?}\n\t", self.stop_aa as char)?;
        write!(f, "reading_frame: {:?}\n\t", self.reading_frame)?;
        write!(f, "allow_ambiguities: {:?}\n\t", self.allow_ambiguities)?;
        write!(f, "strip_gaps: {:?}\n\t", self.strip_gaps)?;
        write!(f, "ignore_gap_codons: {:?}\n\t", self.ignore_gap_codons)?;
        write!(
            f,
            "drop_incomplete_codons: {:?}\n",
            self.drop_incomplete_codons
        )?;
        write!(f, "}}")
    }
}

pub fn find_ambiguity_code(nts: &Vec<&u8>) -> Option<&'static [u8; 1]> {
    let query_set: HashSet<&u8> = nts.iter().copied().sorted().collect();

    for (code, nt_set) in AMBIGUOUS_NT_LOOKUP.entries() {
        let code_set: HashSet<&u8> = nt_set.iter().map(|ambig_char| &ambig_char[0]).collect();
        if query_set == code_set {
            return Some(code);
        }
    }
    None
}

pub fn translate(dna_seq: &[u8], options: &TranslationOptions) -> Result<Vec<u8>> {
    let mut new_seq = dna_seq[options.reading_frame..].to_vec();
    if options.strip_gaps {
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
            if !options.drop_incomplete_codons {
                log::debug!(
                    "The codon {:?} had a length of {} so we're adding a {:?}",
                    String::from_utf8(codon.to_vec())?,
                    codon.len(),
                    options.incomplete_aa as char
                );
                amino_acids.push(options.incomplete_aa);
            }
            continue;
        }
        let nt_triplet: [u8; 3] = codon
            .try_into()
            .expect("The codon should always be a triplet vector since we've checked for it.");

        if !options.strip_gaps {
            let num_gaps = nt_triplet.iter().filter(|char| **char == GAP_CHAR).count();
            if (num_gaps == 1) | (num_gaps == 2) {
                amino_acids.push(options.frameshift_aa as u8);
                continue;
            }
        }
        let amino_acid;

        if CODON_TABLE.contains_key(&nt_triplet) {
            amino_acid = &CODON_TABLE[&nt_triplet][0];
        } else if options.allow_ambiguities && AMBIGUOUS_CODON_TABLE.contains_key(&nt_triplet) {
            amino_acid = &AMBIGUOUS_CODON_TABLE[&nt_triplet][0];
        } else if options.allow_ambiguities
            && AMBIGUOUS_CODON_AND_AA_TABLE.contains_key(&nt_triplet)
        {
            amino_acid = &AMBIGUOUS_CODON_AND_AA_TABLE[&nt_triplet][0];
        } else if STOP_CODONS.contains(&nt_triplet) {
            amino_acid = &options.stop_aa;
        } else {
            log::debug!(
                "Could not find a suitable character for the codon {:?}",
                String::from_utf8(nt_triplet.to_vec())
            );
            amino_acid = &options.unknown_aa;
        }

        if options.ignore_gap_codons & (amino_acid.eq(&GAP_CHAR)) {
            continue;
        } else {
            amino_acids.push(amino_acid.clone());
        }
    }

    Ok(amino_acids)
}

#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn basic_test() -> Result<()> {
        let dna_seq = "ATGTTATAA";
        let expected_translation = "ML*";
        let translation = translate(dna_seq.as_bytes(), &TranslationOptions::default())?;

        assert_eq!(expected_translation.to_owned(), String::from_utf8(translation)?);

        Ok(())
    }

    #[test]
    fn strip_gaps_true() -> Result<()> {
        let dna_seq = "ATGTTA-TAA";
        let expected_translation = "ML*";
        let translation = translate(dna_seq.as_bytes(), &TranslationOptions {
            strip_gaps: true,
            stop_aa: b'*',
            ..TranslationOptions::default()
        }, )?;

        assert_eq!(expected_translation.to_owned(), String::from_utf8(translation)?);

        Ok(())
    }

    #[test]
    fn strip_gaps_false() -> Result<()> {
        let dna_seq = "ATGTTA-TAA";
        let expected_translation = "MLX~";
        let translation = translate(dna_seq.as_bytes(), &TranslationOptions {
            strip_gaps: false,
            drop_incomplete_codons: false,
            incomplete_aa: b'~',
            ..TranslationOptions::default()
        }, )?;

        assert_eq!(expected_translation.to_owned(), String::from_utf8(translation)?);

        Ok(())
    }

    #[test]
    fn ignore_gap_codons() {
        let dna_seq = "ATGTTA---TAA";
        let expected_translation = "ML*";
        let translation = translate(
            dna_seq.as_bytes(),
            &TranslationOptions {
                ignore_gap_codons: true,
                ..TranslationOptions::default()
            },
        )
            .unwrap();

        assert_eq!(expected_translation.as_bytes(), translation.as_slice());
    }

    #[test]
    fn test_ambiguity() -> Result<()> {
        let test_cases = vec!["ATGTTACTNTAA", "NNNATGGGG", "ATGRAY---GTA"];
        let expected_outputs = vec!["MLL*", "?MG", "MB-V"];

        for (idx, test_case) in test_cases.iter().enumerate() {
            let expected_translation = expected_outputs[idx].as_bytes();
            let translation =
                translate(test_case.as_bytes(), &TranslationOptions {
                    unknown_aa: b'?',
                    ..TranslationOptions::default()
                }, )?;

            assert_eq!(String::from_utf8(expected_translation.to_vec())?, String::from_utf8(translation)?);
        }

        Ok(())
    }

    #[test]
    fn test_alternate_stop_codon_char() -> Result<()> {
        let translation_standard =
            translate("ATGTTACTNTAA".as_bytes(), &TranslationOptions::default()).unwrap();

        let translation_custom =
            translate("ATGTTACTNTAA".as_bytes(), &TranslationOptions {
                stop_aa: b';',
                ..TranslationOptions::default()
            })?;

        assert_eq!("MLL*".to_owned(), String::from_utf8(translation_standard)?);
        assert_eq!("MLL;".to_owned(), String::from_utf8(translation_custom)?);
        Ok(())
    }

    // TODO: Add more tests lol
}
