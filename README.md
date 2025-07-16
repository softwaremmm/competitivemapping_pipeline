# Competitive Mapping Pipeline

Competitive Mapping is an algorithm that compares the sample reads with the references in the manifest and makes a positive selection of the reads matching a specific `rname`.

The pipeline for competitive mapping takes a pair of FASTQ files and outputs the positive filtering of the h37_rv reads (i.e. those reads that are judged to map to the *Mycobacterium tuberculosis* H37RV reference genome) along with unmapped reads, a report with the mapping rank (a list of species in the manifest to which reads have mapped `competitivemapping_report.json`).

### Dependencies
* Docker
* Nextflow

### Needed data
* you will need to download manifest to: `$projectDir/data/manifest/manifest_20231001`. See [extra readme](data/manifest/README.md) for details about manifest.
* species list is provided at `$projectDir/test_data/species_list_manifest_20250324.csv`
* (for sylph) should have GTDB representative genomes at path: `$projectDir/data/sylph/gtdb_genomes_reps_r220`. Can be found [here](https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/genomic_files_reps/)

These can all be found in the (dev) knowledge bucket.

## Overview
Competitive mapping uses minimap2 to map reads against manifest (multifasta of reference genomes) and then analyse the resulting bam file.
A species list file is required to match contigs (`rnames`) in the multifasta manifest to a reference as some references have many contigs.

This is co-ordinated by three scripts:

### manifest_builder [In development]
This takes a sylph query report and builds a manifest from it.
This is used for the metagenomic pipeline in development. Myco has a fixed manifest.

parameters:
- sylph_report: Path to the sylph query/profile output file.
- genome_dirs: path to directory containing all the reference genomes. Must contain a genome_paths.tsv matching what you find in gtdb_genomes_reps (see note on data needed earlier)
- metadata_files: file with taxonomy mapping for each reference
- include_whole_genus: if this is set then for each species sylph finds the whole genus will be added to the manifest

### manifest_mapper
Map reads against the manifest. It will select minimap2 settings based on the seq platform provided.
This step produces a bam file.

parameters:
- manifest: path to multifasta
- reads: path to input fastqs
- seq_platform: `ont` or `illumina`
- filter_secondary/filter_supplementary: Will exclude secondary/supplementary alignments from resulting bam.
- sort_by_name: this sorts the bam file by name rather than position. (Used for development)
- equality_in_cigar: add more equality information to CIGAR strings in bam file. (Used for development)

### competitive_mapping
This takes the bam file and calculates various coverage stats for the references in the manifest. It:
- Run samtools coverage and aggregate results with `process_coverage.py`
- Calculate alignment stats from bam file use `process_aln_stats.py`
- Produce overall csv and report json

parameters:
- input_bam: bam from previous step
- contigs: the species list file
- seq_platform: `ont` or `illumina`
- ref_for_fastq: If set then the reads which mapped to the provided reference will be extracted from the bam file. Used in myco to select TB reads.
- reference_name (optional): Set this to the name of the reference you want to filter reads for e.g. `M.tuberculosis` to output filtered _M. tuberculosis_ reads. Can also take a comma seperated list.
Note: currently unmapped reads will also be extracted by default

## Notes on Bam to Fastq
One step in competitive mapping is to filter the bam file (created by mapping against manifest) for tb reads and unmapped reads, and extracting these to a fastq file. This is complex for paired reads!! And so leads to seeming discrepencies with the `species_comparison_report.json`

Steps:
* filter bam file for reads mapping to h37rv or "\*".
  - "\*" is used for unmapped reads, but only when both reads in a pair are unmapped (as far as Matthew can see from minimap2 outputs).
  - Reads which are unmapped but have a mapped pair will list rname to match the pair, but have the sam flag set for being unmapped.
* `samtools fastq` is used with:
  - `--excl-flags 0x100` which excludes secondary reads
  - `-s /dev/null` to exclude singleton reads, so that the resulting fastq files are properly paired.

The result of this is that reads are only converted to fastq if
1. Both in a pair are unmapped
2. Both in a pair have a primary or supplementary mapping to h37rv

In the future we could change this to out put a read pair as long as **either** read in a pair map to h37rv.

## Running Nextflow
The workflow takes the following inputs
- input_dir. Path to directory containing input fastq files
- seq_platform. `illumina` or `ont`.
- manifest. Path to manifest file
- species_list. Patht to species list file which has the contig to genome mapping

When running locally can save outputs by using `--publish_dir`.

Example using test data:
```bash
nextflow run . \
		--seq_platform illumina \
		--input_dir test_data/chloro_10k \
		--manifest data/manifest/manifest_20231001 \
		--species_list test_data/species_list_manifest_20250324.csv \
    --publish_dir results
```


By default it will look for files in the input directory based on the following params:
- `params.input_paired_suffix = "*_{1,2}.fastq.gz"`
- `params.input_single_suffix = "*.fastq.gz"`

but these can be overriden. e.g.
```
nextflow run ... --input_paired_suffix "tb_sample*_{1,2}.fna.gz"
```


### Running Tests
The tests are executed using [nf-test](https://github.com/askimed/nf-test).

Before running the tests check that you have [needed data](#needed-data).

To run tests, run the following command

```bash
nf-test test tests/nextflow/*.test
```

If you have made changes to the python code, you may need to build a local test container:
```bash
docker build -t test_container_cm .
nf-test test tests/nextflow/*.test --profile local_docker
```

## Python

A Python package that processes the output from command line tools orchestrated by NextFlow is included in this repository. This gets installed in the docker image used by nextflow.

### Installation

Clone the repo as described above. Create a virtual environment and install the package using `pip install -e .[dev]`. Set up pre-commit with `pre-commit install`.


### Execution

Each python module can be run individually. Check arguments with `--help`, e.g. `process_coverage --help` .

### Testing
The tests take 3 minutes as they include running the full process of mapping reads as well.
```
pytest tests/
```

### Outputs

The output from the Python CLI is
validated against a [JSON Schema](src/competitivemapping/competitivemapping.schema.json). [Documentation for the schema](schema_doc.md)
can be built / updated using:

```
generate-schema-doc src/competitivemapping/competitivemapping.schema.json --config template_name=md
```

## Integrating to a pipeline

If you want to use the Competitive Mapping Pipeline as subworkflow in your pipeline, use the competitive_mapping named workflow.


## Manifest remarks
There is a reference (AP018410.1 - M.pseudoshottsii) which is in the manifest species list, but not in the manifest itself.
This will be flagged with a warning when running the code.


## Tags, Releases, and Committing
Use conventional commits. This is enforced with commitizen validate action and pre-commit hooks:
```bash
pre-commit install --hook-type commit-msg
```

This repo uses a standard gitflow approach, but with some changes to deal with docker containers in nextflow:
- There is a pyproject version which should be updated to match semantic version releases (from main/release branches)
- There is an `active_version` controlled by version_bumper which allows develop to have commit hash based versions.
- Every push to develop will cause an action to run `bumper bump <commit-hash> --no-tag --active`. This:
    - bumps the `active_version` in pyproject.toml
    - bumps the container tag used by nextflow processes
    - leaves the pyproject version as is

  Another workflow then builds and pushes the new container.\
  **Important: To deploy this you will need to use the hash of the bump commit, not the hash of the merge commit/container.**

- In a release branch you can create a release candidate with `bumper bump a.b.c-rcX`. This also updates the pyproject version. Pushing the changes and new tag (automatically created) will trigger a build action.
- When release branch is ready for main run `bumper bump a.b.c --no-tag`. Push these changes to main and make a release there to build the container.
