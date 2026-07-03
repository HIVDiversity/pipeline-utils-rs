use phf::{phf_map, phf_set};

pub const GAP_CHAR: u8 = b"-"[0];
pub const DEFAULT_STOP_CHAR: u8 = b"*"[0];

pub(crate) static CODON_TABLE: phf::Map<&[u8; 3], &[u8; 1]> = phf_map! {
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

pub static STOP_CODONS: phf::Set<&[u8; 3]> = phf_set! {b"TAA", b"TAG", b"TGA"};

// Thanks https://cran.r-project.org/web/packages/MLMOI/vignettes/StandardAmbiguityCodes.html
pub(crate) static AMBIGUOUS_CODON_TABLE: phf::Map<&[u8; 3], &[u8; 1]> = phf_map! {
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
pub(crate) static AMBIGUOUS_CODON_AND_AA_TABLE: phf::Map<&[u8; 3], &[u8; 1]> = phf_map! {
    b"RAY" => b"B",
    b"SAR" => b"Z"
};

// https://www.bioinformatics.org/sms/iupac.html or https://www.hiv.lanl.gov/content/sequence/HelpDocs/IUPAC.html
pub static AMBIGUOUS_NT_LOOKUP: phf::Map<&[u8; 1], phf::Set<&[u8; 1]>> = phf_map! {
    b"R" => phf_set!(b"A", b"G"),
    b"Y" => phf_set!(b"C", b"T"),
    b"S" => phf_set!(b"C", b"G"),
    b"W" => phf_set!(b"A", b"T"),
    b"K" => phf_set!(b"G", b"T"),
    b"M" => phf_set!(b"A", b"C"),
    b"B" => phf_set!(b"T", b"C", b"G"),
    b"H" => phf_set!(b"T", b"C", b"A"),
    b"D" => phf_set!(b"T", b"A", b"G"),
    b"V" => phf_set!(b"C", b"A", b"G"),
    b"N" => phf_set!(b"T", b"A", b"G", b"C"),
    b"X" => phf_set!(b"T", b"A", b"G", b"C"),
};
