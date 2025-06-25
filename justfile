#!/usr/bin/env just --justfile

latest-tag := `git describe --tags --abbrev=0`
image-name := "dlejeune/pipeline-utils-rs"

default:
    just --list

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
    sudo docker run --rm -it -v ./:/data dlejeune/{{ image-name }}:{{ tag }} bash

# Builds and pushed the most recently tagged branch in a docker container
[group('docker')]
docker: build-docker push-docker

build:
    cargo build

run *args="":
    cargo run {{ args }}

test-align-trim *args:
    just run align-trim -r new_test_data/align-trim/ref.fasta -i new_test_data/align-trim/query.fasta -o new_test_data/align-trim/output.fasta {{ args }}
