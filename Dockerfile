# ---- Rust Build Stage ----
FROM rust:1.86 AS chef

RUN apt-get update && apt-get install -y clang llvm-dev
RUN cargo install cargo-chef
WORKDIR /app

FROM chef AS planner
COPY tie_break/Cargo.* .
COPY tie_break/src ./src
RUN cargo chef prepare --recipe-path recipe.json

FROM chef AS builder
COPY --from=planner /app/recipe.json recipe.json
RUN cargo chef cook --release --recipe-path recipe.json

# Copy the source code and build
COPY tie_break/Cargo.* .
COPY tie_break/src ./src
COPY tie_break/test_data ./test_data

ARG TESTING=false
RUN if [ "$TESTING" = "true" ]; then \
    apt-get update && \
    apt-get install -y samtools && \
    cargo test --release; \
fi

RUN cargo build --release

# ---- Conda Build Stage ----
FROM continuumio/miniconda3 AS conda_builder
WORKDIR /app

COPY env.yml /app/env.yml
RUN conda env create --file env.yml


# Pack conda installation
RUN conda install -c conda-forge conda-pack
RUN conda-pack -n competitive_mapping -o /app/conda_env.tar.gz


# ---- Conda Build Stage ----
FROM debian:stable-slim AS runtime
WORKDIR /app

RUN apt-get update && apt-get install -y git libssl-dev procps && rm -rf /var/lib/apt/lists/*

COPY --from=conda_builder /app/conda_env.tar.gz /app/conda_env.tar.gz
RUN mkdir -p /opt/conda
RUN tar -xzf /app/conda_env.tar.gz -C /opt/conda && rm /app/conda_env.tar.gz

ENV PATH="/opt/conda/bin:$PATH"


# Install Python code
COPY ./src /app/src
COPY ./pyproject.toml /app/pyproject.toml

# Install pytest if TESTING is true (default false)
ARG TESTING=false
RUN if [ "$TESTING" = "true" ]; then \
        pip install .[dev]; \
    else \
        pip install .; \
    fi

# Copy the tie_break rust build
COPY --from=builder /app/target/release/tie_break /usr/local/bin/tie_break
