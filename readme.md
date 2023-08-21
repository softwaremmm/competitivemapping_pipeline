
# Competitive Mapping Pipeline

Competitive Mapping is an algorithm that compares the sample reads with the references in the manifest and makes a positive selection of the reads matching a specific `rname`. 

The pipeline for competitive mapping that takes a pair of FASTQ files and outputs the positive filtering of the h37_rv reads, a report with the mapping rank and an error report.

## Dependencies
* Docker
* 1G hard-disk 


## Run Locally

Clone the project

```bash
  git clone git@github.com:GlobalPathogenAnalysisService/competitivemapping_pipeline.git
```

Go to the project directory

```bash
  cd competitivemapping_pipeline

```

Run the pipeline

```bash
  nextflow run . -input_dir $PATH -manifest $PATH_TO_MANIFEST_FILE
```

where $PATH is the path to a folder that contains a pair of FAST.GZ files following a *{1,2}.f*q.gz regex convention, and
$PATH_TO_MANIFEST_FILE is the path to the FASTA file putting together the list of target sequences.


## Running Tests
The tests are executed using [nf-test](https://github.com/askimed/nf-test). 

To run tests, run the following command

```bash
 nf-test test tests/nextflow/*.test
```


## Screenshots

Expected output:

![App Screenshot](/home/marcela/Desktop/Screenshot from 2023-08-21 10-45-19.png)



## Integrating to a pipeline

If you want to use the Competitive Mapping Pipeline as subworkflow in your pipeline, use the competitive_mapping named workflow. 
Parameters: 

Channel with a tuple of _val_, _path_, _path_ 
_path_ to tha manifest file


