# Use minconda base image 
FROM continuumio/miniconda3

# Set the working directory within the container
WORKDIR /app

# Set up bioconda
RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge

# Install minimap2
RUN conda install minimap2~=2.26

# Install samtools
RUN conda install samtools~=1.17

# Install csvkit
RUN conda install csvkit~=1.1.1

# Install jq
RUN conda install jq=1.6

# Run jq (In practice this command will be overriden by NextFlow)
CMD ["jq"]
