# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Pipeline Utils RS (PURS) is a Rust CLI providing sequence-manipulation utilities for the DeepLEAP alignment
pipeline (collapse/expand FASTA by identity, consensus generation, GenBank extraction, ambiguity resolution,
translation, SAM trimming, reverse translation, miniprot PAF processing, stop-codon trimming). Each subcommand
is also usable standalone.

Note: the crate is built as a `cdylib` (see `[lib]` in `Cargo.toml`, name `pygetcons`) via `pyo3`, so this
codebase is also importable as a Python extension module in addition to being a CLI binary.

## Build & development commands

This project is normally built and run inside a Docker container (image `purs-build`) via `just`, not with a
bare local `cargo` — check whether that image/toolchain is available before assuming plain `cargo` commands work:

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

- `src/main.rs` — parses CLI args via `cli::Cli::parse()` and dispatches each `Commands` variant to the
  corresponding `tools::<name>::run(...)` function. This file should stay a thin dispatcher; command-specific
  logic belongs in `src/tools/`.
- `src/cli.rs` — all `clap` argument/subcommand definitions (the `Cli` struct and `Commands` enum). Adding a new
  subcommand means adding a variant here *and* a matching arm in `main.rs` *and* a new module in `src/tools/`.
- `src/tools/` — one module per subcommand, each exposing a `pub fn run(...) -> Result<()>` that implements the
  actual behavior. `mod.rs` re-exports each tool module.
- `src/utils/` — shared helpers used across tools, notably `fasta_utils.rs` (FASTA I/O helpers) and
  `translate.rs` (codon translation logic and `TranslationOptions`, which `cli::TranslateCliOptions` converts
  into).
- Depends on local git submodules `lib/hts-sys` and `lib/rust-htslib` (vendored/patched versions of the htslib
  bindings); the `Cargo.toml` currently pulls `rust-htslib` from crates.io by default, with the local submodule
  path commented out — check which is active before assuming submodule changes take effect.
- `new_test_data/` contains fixture input/output FASTA/PAF files used by the `just test-*` recipes (e.g.
  `test-align-trim`, `test-kmer-trim`) — note these `just` recipes reference subcommands (`align-trim`,
  `kmer-trim`) that are not currently defined in `cli.rs`/`src/tools`, so they may be stale relative to the
  current command set in `cli.rs`.
