# Use minconda base image 
FROM continuumio/miniconda3

# Set the working directory within the container
WORKDIR /app

# Conda install dependencies
COPY env.yml /app/env.yml
RUN conda env update -n base --file env.yml

# Install Python code for processing output of minimap and samtools
COPY ./src /app/src
COPY ./pyproject.toml /app/pyproject.toml
RUN pip install .
