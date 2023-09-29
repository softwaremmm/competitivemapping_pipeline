
# Competitive Mapping Pipeline

Competitive Mapping is an algorithm that compares the sample reads with the references in the manifest and makes a positive selection of the reads matching a specific `rname`. 

The pipeline for competitive mapping takes a pair of FASTQ files and outputs the positive filtering of the h37_rv reads (i.e. those reads that are judged to map to the *Mycobacterium tuberculosis* H37RV reference genome `h37rv_1.fastq.gz` and `h37rv_2.fastq.gz`), a report with the mapping rank (a list of species in the manifest to which reads have mapped `competitivemapping_report.json`) and an error report (the concatenated standard error output from the minimap and samtools tools `competitivemapping_error.json` - the contents are not in JSON format).

## Dependencies
* Docker
* 1G hard-disk 
* Nextflow


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
  nextflow run . --input_dir $PATH --manifest $PATH_TO_MANIFEST_FILE
```

where $PATH is the path to a folder that contains a pair of FAST.GZ files following a *{1,2}.f*q.gz regex convention, and
$PATH_TO_MANIFEST_FILE is the path to the FASTA file putting together the list of target sequences. Those paths do not need to be absolute paths.


## Running Tests
The tests are executed using [nf-test](https://github.com/askimed/nf-test). 

Before running the tests, you will need to ensure appropriate test data is available by setting up the following symlinks in
the root of the repository:

* `data/manifest/target_101_new.fasta` - a manifest (reference data)
* `data/WTCHG_885333_73205296_1` containing `WTCHG_885333_73205296_1.fastq.gz` and `WTCHG_885333_73205296_1.fastq.gz` - these should be TB FASTQs
* `data/abscessus` containing `file_R1.fastq.gz` and `file_R2.fastq.gz` - these should be NTM FASTQs

To run tests, run the following command

```bash
 nf-test test tests/nextflow/*.test
```


## Screenshots

Expected output:
![output](https://github.com/GlobalPathogenAnalysisService/competitivemapping_pipeline/assets/65816841/b43b60f9-87c9-417a-b11c-67ef9c491f2d)


## Integrating to a pipeline

If you want to use the Competitive Mapping Pipeline as subworkflow in your pipeline, use the competitive_mapping named workflow. 
Parameters: 

* Channel with a tuple of _val_, _path_, _path_ 

* _path_ to tha manifest file


