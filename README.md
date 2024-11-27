# Competitive Mapping Pipeline

Competitive Mapping is an algorithm that compares the sample reads with the references in the manifest and makes a positive selection of the reads matching a specific `rname`.

The pipeline for competitive mapping takes a pair of FASTQ files and outputs the positive filtering of the h37_rv reads (i.e. those reads that are judged to map to the *Mycobacterium tuberculosis* H37RV reference genome) along with unmapped reads, a report with the mapping rank (a list of species in the manifest to which reads have mapped `competitivemapping_report.json`).

## Overview
Competitive mapping uses minimap2 to map reads against manifest (multifasta of reference genomes). A species list file is required to match contigs to `rname`'s. This is coordinated by `competitive_mapping.py` which does:
- (Optional) create manifest if using sylph report as input
- Map reads against manifest with minimap2
- Run samtools coverage and aggregate results with `process_coverage.py`
- Calculate alignment stats from bam file use `process_aln_stats.py`
- Produce overall csv and report json



### Dependencies
* Docker
* Nextflow

### Needed data
* manifest should be at path: `$projectDir/data/manifest/manifest_20231001`
* species list is provided at `$projectDir/test_data/species_list_manifest_20240710.csv`
* (for sylph approach) should have GTDB representative genomes at path: `$projectDir/data/sylph/gtdb_genomes_reps_r220`. Can be found [here](https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/genomic_files_reps/)

These can all be found in the (dev) knowledge bucket.

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

## Nextflow

### Run Locally

Run the pipeline

```bash
nextflow run . --input_dir $PATH --manifest $PATH_TO_MANIFEST_FILE --species_list $PATH_TO_SPECIES_LIST_FILE --seq_platform $SEQ_PLATFORM
```

where $PATH is the path to a folder that contains a pair of FASTQ.GZ files following a *{1,2}.f*q.gz regex convention,
$PATH_TO_MANIFEST_FILE is the path to a manifest file containing a list of target contigs and $PATH_TO_SPECIES_LIST_FILE is the path to a species list file where contig rnames are mapped to genomes. Those paths do not need to be absolute paths.
$SEQ_PLATFORM should be 'ont' or 'illumina' depending on platform used.


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

Clone the repo as described above. Create a virtual environment and install the package using `pip install -e .[dev]`. Set up pre-commit
with `pre-commit install`.


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

# Conventional Commits
Use [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) when developing for this repo.
You should install the pre-commit hooks to check your commit messages. This can be done using the tool `pre-commit` which is a dev dependency in the `pyproject.toml`.
You can also use `commitizen` (another dev dependency) to help with writing conventional commits.

To install hooks run
```bash
pre-commit install --hook-type commit-msg
pre-commit install # to get other hooks for formatting etc
```

To make commit with commitizen run
```bash
cz c
```

## Tags and Releases

[Commitizen](https://commitizen-tools.github.io/commitizen/) is used to manage versioning of releases. This tool
can be used to make commits to this repository. Regardless, [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/)
are required to ensure correct version numbering and changelog population.

**Do not add tags by hand.**

On merging a Pull Request a [GitHub action will run](.github/workflows/version.yaml), causing Commitizen to:
* Determine the new [semver](https://semver.org/) based on conventional commits.
* Replace the previous semver in [pyproject.toml](pyproject.toml) and other files as specified therein.
* Update the [CHANGELOG](CHANGELOG.md) based on commit messages.
* Commit these changes to the `main` branch.
* Create a tag for this commit with the tag name of the newly determined semver.
* Build a docker container for this new tag
* Create a new release from this tag.

There is a workflow for manually triggering a docker build action.
This shouldn't be required unless something has gone wrong with the commitizen action.
