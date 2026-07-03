FROM rust:1.93.0 AS planner
RUN cargo install cargo-chef

WORKDIR /app
# Copy the whole project
COPY . .
# Prepare a build plan ("recipe")
RUN cargo chef prepare --recipe-path recipe.json

FROM rust:1.93.0 AS builder
RUN cargo install cargo-chef

# Copy the build plan from the previous Docker stage
COPY --from=planner /app/recipe.json recipe.json

RUN apt-get update && apt-get install -y lsb-release gnupg\
                                                      clang \
                                                      libclang-dev \
                                                      llvm \
                                                      cmake \
                                                      build-essential \
                                                      && rm -rf /var/lib/apt/lists/*


# Build dependencies - this layer is cached as long as `recipe.json`
# doesn't change.
RUN cargo chef cook --release --features trim-sam,process-miniprot --recipe-path recipe.json

# Build the whole project
COPY . .
RUN cargo build --release --features trim-sam,process-miniprot

FROM debian:trixie AS release

RUN apt-get update && apt-get install -y procps

COPY --from=builder ./target/release/pipeline-utils-rs /usr/local/bin/pipeline-utils-rs