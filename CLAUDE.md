# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Pipeline Utils RS (PURS) is a Rust CLI providing sequence-manipulation utilities for the DeepLEAP alignment
pipeline (collapse/expand FASTA by identity, consensus generation, GenBank extraction, ambiguity resolution,
translation, SAM trimming, reverse translation, miniprot PAF processing, stop-codon trimming). Each subcommand
is also usable standalone.

Note: the crate is split into a lib (`src/lib.rs`, `[lib] name = "purs"`, `crate-type = ["cdylib", "rlib"]`) and
a bin (`src/main.rs`) target. The lib is built as a `cdylib` via `pyo3` (see `src/python.rs`, `#[pymodule] mod
purs`) so this codebase is also importable as a Python extension module in addition to being a CLI binary; the
bin depends on the lib crate (`purs::cli`, `purs::tools`, ...) and stays a thin dispatcher.

## Build & development commands

This project is normally built and run inside a Docker container (image `purs-build`) via `just`, not with a
bare local `cargo`. This is required, not just convenient: `rust-htslib`'s `bindgen` feature builds `hts-sys`
from the vendored htslib C sources at compile time, and the Clang/libclang version installed on this system
does not produce a working build of that C code — the container pins a Clang version that does. Check whether
the `purs-build` image/toolchain is available before assuming plain `cargo` commands work:

```sh
just build                 # cargo build inside the purs-build docker container
just run <args>             # cargo run -- <args> inside the container
just cargo <args>           # arbitrary cargo subcommand inside the container
just build-docker           # build the release docker image (tagged with latest git tag)
just push-docker             # push the release image to dockerhub
just docker                  # build-docker + push-docker
```

If working with a local Rust toolchain directly (`rust-version = 1.93.0`, edition 2024), the usual `cargo build`,
`cargo test`, `cargo check`, and `cargo run -- <subcommand> <args>` apply. Run a single test with
`cargo test <test_name>`.

CLI usage follows `pipeline-utils-rs <subcommand> [options]`; run `cargo run -- --help` or
`cargo run -- <subcommand> --help` to see options for any tool.

## Architecture

- `src/lib.rs` — declares the library's public modules (`cli`, `python`, `tools`, `utils`) exposed as crate
  `purs`.
- `src/main.rs` — the bin target; parses CLI args via `cli::Cli::parse()` (imported from `purs::cli`) and
  dispatches each `Commands` variant to the corresponding `tools::<name>::run(...)` function. This file should
  stay a thin dispatcher; command-specific logic belongs in `src/tools/`.
- `src/python.rs` — the `pyo3` bindings (`#[pymodule] mod purs`). Each exposed `#[pyfunction]` wraps a
  `tools::<name>` function, converting Python dict/Vec args to/from the internal `FastaRecords` representation
  via the local `dict_to_records`/`records_to_dict` helpers. Currently exposes: `get_consensus`, `translate`,
  `reverse_translate`, `replace_ambiguities`, `trim_after_stop_codon`, `collapse`, `expand`. Adding Python
  support for a tool means adding a `#[pyfunction]` here *and* registering it in the `mod purs` block — pyo3
  auto-registers functions defined inside the `#[pymodule] mod purs { ... }` block.
- `src/cli.rs` — all `clap` argument/subcommand definitions (the `Cli` struct and `Commands` enum). Adding a new
  subcommand means adding a variant here *and* a matching arm in `main.rs` *and* a new module in `src/tools/`.
- `src/tools/` — one module per subcommand, each exposing a `pub fn run(...) -> Result<()>` for CLI use; several
  also expose lower-level functions consumed directly by `src/python.rs` (e.g.
  `tools::get_consensus::{sequences_to_matrix, build_consensus}`, `tools::collapse::{collapse_sequences,
  build_collapsed_output}`). `mod.rs` re-exports each tool module. Note `extract_seq_from_gb.rs` was renamed to
  `gb_extract.rs`.
- `src/utils/` — shared helpers used across tools, notably `fasta_utils.rs` (FASTA I/O helpers, `FastaRecords`
  type), `translate.rs` (codon translation logic and `TranslationOptions`, which `cli::TranslateCliOptions`
  converts into), and `codon_tables.rs`.
- `rust-htslib` is pulled from crates.io (`bindgen` feature, vendored htslib build) rather than via the local
  git submodules `lib/hts-sys`/`lib/rust-htslib`, which have been removed from this working tree — see the
  Docker note above for why this dependency needs the container's Clang toolchain to build.
- `new_test_data/` contains fixture input/output FASTA/PAF files used by the `just test-*` recipes (e.g.
  `test-align-trim`, `test-kmer-trim`) — note these `just` recipes reference subcommands (`align-trim`,
  `kmer-trim`) that are not currently defined in `cli.rs`/`src/tools`, so they may be stale relative to the
  current command set in `cli.rs`.
