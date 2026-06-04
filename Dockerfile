FROM anaconda/miniconda:26.3.2
ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=true

# Set the working directory within the container
WORKDIR /app

# Needed by nextflow
RUN apt-get update && \
    apt-get install -y procps && \
    rm -rf /var/lib/apt/lists/*

# Set up bioconda
COPY env.yml /app/env.yml
RUN conda env update -n base --file env.yml && conda clean -afy

# Add Python files from repo to Docker image
COPY ./src /app/src
COPY ./pyproject.toml /app/pyproject.toml

ARG TESTING=false
RUN if [ "$TESTING" = "true" ]; then \
        pip install .[dev]; \
    else \
        pip install .; \
    fi
