# Use minconda base image 
FROM continuumio/miniconda3

# Set the working directory within the container
WORKDIR /app

# Set up bioconda
RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict

# Install minimap2
RUN conda install minimap2~=2.26

# Install samtools
RUN conda install samtools~=1.17

# Include temporary dummy data in container
COPY lib/test_data/cov_sorted_h37rv_100k.json /app/cov_sorted_h37rv_100k.json

# Run minimap2 (In practice this command will be overriden by NextFlow)
CMD ["minimap2"]