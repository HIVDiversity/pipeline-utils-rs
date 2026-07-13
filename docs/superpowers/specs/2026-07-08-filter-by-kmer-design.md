# Design: `filter-by-kmer`

## Purpose

Add a new subcommand that filters FASTA sequences by whether they start and/or end with a
user-provided k-mer (or one of several allowed k-mers), e.g. requiring a start codon (`ATG`)
and/or any stop codon (`TAA`/`TAG`/`TGA`). Matching is IUPAC-ambiguity-aware: an ambiguity code
in either the query k-mer or the sequence matches any base it can represent.

## CLI

New `Commands::FilterByKmer` variant in `src/cli.rs`:

```
purs filter-by-kmer -i <input_file> -o <output_file>
    [--start-kmers ATG]
    [--end-kmers TAA,TAG,TGA]
    [--rejected-seq-output <path>]
    [--report-file <path>]
```

- `-i, --input-file <PathBuf>`: input FASTA (required).
- `-o, --output-file <PathBuf>`: FASTA to write sequences that pass all requested checks
  (required).
- `--start-kmers <STR>`: comma-separated list of allowed k-mers to match against the start of
  each sequence. A sequence passes the start check if it matches *any* one of them.
- `--end-kmers <STR>`: comma-separated list of allowed k-mers to match against the end of each
  sequence, same semantics as `--start-kmers`.
- At least one of `--start-kmers` / `--end-kmers` must be provided (enforced via a clap
  `ArgGroup` with `required = true, multiple = true`). Providing only one runs only that check;
  the other is simply not evaluated (not "n/a-and-fail").
- `--rejected-seq-output <PathBuf>` (optional): FASTA file to write sequences that failed at
  least one requested check.
- `--report-file <PathBuf>` (optional): CSV report, one row per input sequence.

K-mers within a single `--start-kmers`/`--end-kmers` list may have different lengths (e.g. this
is not needed for the stop-codon case but is not disallowed).

## Matching semantics

- Two individual bases are *compatible* if their IUPAC-expanded base sets intersect. Concrete
  bases (A/C/G/T) expand to a singleton set containing themselves. Ambiguity codes (R, Y, S, W,
  K, M, B, D, H, V, N, X) expand via the existing `AMBIGUOUS_NT_LOOKUP` table in
  `src/utils/codon_tables.rs`. This applies symmetrically: an ambiguity code in the query k-mer
  matches any concrete base it represents in the sequence, and vice versa (e.g. sequence base
  `N` is compatible with query base `A`).
- A sequence matches a given k-mer at the **start** if the sequence is at least as long as the
  k-mer and every position of the k-mer is compatible with the corresponding position of the
  sequence's prefix. Matching at the **end** is the mirrored check against the sequence's
  suffix. A sequence shorter than a candidate k-mer automatically fails that candidate.
- A sequence passes the **start check** if it is compatible with at least one k-mer in
  `--start-kmers` (when that flag was given); analogously for the **end check**.
- A sequence is **kept** if it passes every check that was actually requested. If only one of
  `--start-kmers`/`--end-kmers` is given, only that check determines keep/reject.

## Module structure

New `src/tools/filter_by_kmer.rs`, following the shape of `src/tools/filter_by_length.rs`:

- `pub(crate) fn bases_compatible(query: u8, seq: u8) -> bool` — the IUPAC compatibility check
  described above.
- `pub(crate) fn matches_kmer_at_start(seq: &[u8], kmer: &[u8]) -> bool` and
  `matches_kmer_at_end(seq: &[u8], kmer: &[u8]) -> bool`.
- `pub(crate) struct FilterReportRow { seq_name: String, start_match: Option<bool>, end_match:
  Option<bool>, kept: bool }` — `None` when that check wasn't requested.
- `pub(crate) fn filter_by_kmer(sequences: FastaRecords, start_kmers: Option<&[Vec<u8>]>,
  end_kmers: Option<&[Vec<u8>]>) -> Result<(FastaRecords, FastaRecords, Vec<FilterReportRow>)>`
  — core logic, mirroring `filter_by_length`'s kept/rejected/report-rows return shape. Bails if
  `sequences` is empty, consistent with `filter_by_length`.
- `fn write_report(...)` — CSV writer, columns: `seq_name,start_match,end_match,kept`, values
  `true`/`false`/`n/a`.
- `pub fn run(input_file: &PathBuf, output_file: &PathBuf, report_file: Option<&PathBuf>,
  rejected_seq_output: Option<&PathBuf>, start_kmers: Option<&[Vec<u8>]>, end_kmers:
  Option<&[Vec<u8>]>) -> Result<()>` — loads FASTA, calls `filter_by_kmer`, writes kept output,
  optionally writes rejected output and report, matching `filter_by_length::run`'s structure and
  logging style (version banner via `log::info!` + `colored`).

K-mer strings from the CLI are parsed (split on `,`, uppercased to match `load_fasta`'s
uppercasing behavior) into `Vec<Vec<u8>>` in `cli.rs`, analogous to how `LengthThresholdArgs`
converts into `LengthThreshold`.

## Wiring

- `src/tools/mod.rs`: add `pub mod filter_by_kmer;`.
- `src/cli.rs`: add `Commands::FilterByKmer { ... }` variant with doc comments matching the
  style of neighboring variants.
- `src/main.rs`: add a dispatch arm calling `tools::filter_by_kmer::run(...)`.
- `src/python.rs`: add a `#[pyfunction]` (e.g. `filter_by_kmer`) registered in the `#[pymodule]
  mod purs` block, converting Python dict/Vec args via the existing `dict_to_records`/
  `records_to_dict` helpers, following the pattern used for `collapse`/`get_consensus`.

## Errors

- Empty input FASTA is an error (`bail!`), consistent with `filter_by_length`.
- Neither `--start-kmers` nor `--end-kmers` given is a clap usage error (`ArgGroup`), not a
  runtime error.

## Testing

Unit tests in `src/tools/filter_by_kmer.rs` covering:
- `bases_compatible`: exact match, ambiguity code in query matching concrete base, ambiguity
  code in sequence matching concrete query base, two incompatible concrete bases, two
  ambiguity codes with overlapping/non-overlapping expansions.
- `matches_kmer_at_start` / `matches_kmer_at_end`: exact match, ambiguity match, sequence
  shorter than k-mer (fails), k-mer longer than remaining sequence at a non-zero offset case
  not applicable (single k-mer vs single end).
- `filter_by_kmer`: start-only check, end-only check, both checks, multiple candidate k-mers
  (stop-codon-like case, one candidate matches), sequence failing all candidates, empty input
  errors.
