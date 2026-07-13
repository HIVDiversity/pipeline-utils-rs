#!/usr/bin/env just --justfile

latest-tag := `git describe --tags --abbrev=0`
image-name := "dlejeune/pipeline-utils-rs"

default:
    just --list

# Builds with cargo
[group('general')]
build:
    cargo build

# Runs this project with cargo and desired arguments
[group('general')]
run *args="":
    cargo run -- {{ args }}

# Cleans the project
[group('general')]
clean:
    cargo clean

# Formats this project
[group('quality-control')]
format:
    cargo fmt

# Runs the test suite
[group('quality-control')]
test:
    cargo test

# Run clippy and surface warnings
[group('quality-control')]
lint:
    cargo clippy -- -D warnings

[group('versioning')]
release level:
    cargo release {{ level }} --execute

# Builds a docker image with the most recent git tag
[group('docker')]
build-docker:
    sudo docker build -t {{ image-name }}:{{ latest-tag }} .

# Pushes the docker image with the most recent git tag to dockerhub
[group('docker')]
push-docker:
    sudo docker push {{ image-name }}:{{ latest-tag }}

# Runs an interactive docker container with the current wd mounted at /data
[group('docker')]
run-docker-it tag=latest-tag:
    sudo docker run --rm -it -v ./:/data {{ image-name }}:{{ tag }} bash

# Builds and pushed the most recently tagged branch in a docker container
[group('docker')]
docker: build-docker push-docker

# Builds this project in a docker container with relevant clang libraries
[group('dockerized')]
build-in-docker:
    just cargo-in-docker build

# Runs the project via cargo run -- <args> in a docker container with necessary clang dependencies
[group('dockerized')]
run-in-docker *args="":
    just cargo-in-docker run -- {{ args }}

# Runs an arbitrary cargo command in a docker container with necessary clang dependencies
[group('dockerized')]
cargo-in-docker *args:
    docker run --rm  -u $(id -u):$(id -g) -i -v "$HOME:$HOME" -w "/home/dlejeune/Projects/pipeline-utils-rs"  purs-build cargo {{ args }}
