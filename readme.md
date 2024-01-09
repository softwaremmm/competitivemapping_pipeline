
# Competitive Mapping Pipeline

Competitive Mapping is an algorithm that compares the sample reads with the references in the manifest and makes a positive selection of the reads matching a specific `rname`. 

The pipeline for competitive mapping takes a pair of FASTQ files and outputs the positive filtering of the h37_rv reads (i.e. those reads that are judged to map to the *Mycobacterium tuberculosis* H37RV reference genome) along with unmapped reads, a report with the mapping rank (a list of species in the manifest to which reads have mapped `competitivemapping_report.json`) and an error report (the concatenated standard error output from the minimap and samtools tools `competitivemapping_error.json` - the contents are not in JSON format).

## Commits

[Commitizen](https://commitizen-tools.github.io/commitizen/) is used to manage versioning of releases. This tool
can be used to make commits to this repository. Regardless, [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) 
are required to ensure correct version numbering and changelog population.

## Tags and Releases

**Do not add tags by hand.**

On merging a Pull Request a [GitHub action will run](.github/workflows/version.yaml), causing Commitizen to:
* Determine the new [semver](https://semver.org/) based on conventional commits.
* Replace the previous semver in [pyporject.toml](pyproject.toml) and other files as specified therein.
* Update the [CHANGELOG](CHANGELOG.md) based on commit messages.
* Commit these changes to the `main` branch.
* Create a tag for this commit with the tag name of the newly determined semver.

If you wish to release this version and make it available for use in the product, **do this by hand** e.g.
by navigating to the repository on GitHub, clicking "Tags", clicking the desired tag, clicking "Generate
Release Notes", then "Create Release from Tag".

## Nextflow

### Dependencies
* Docker
* 1G hard-disk 
* Nextflow

### Run Locally

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
  nextflow run . --input_dir $PATH --manifest $PATH_TO_MANIFEST_FILE --species_list $PATH_TO_SPECIES_LIST_FILE
```

where $PATH is the path to a folder that contains a pair of FAST.GZ files following a *{1,2}.f*q.gz regex convention, 
$PATH_TO_MANIFEST_FILE is the path to a manifest file containing a list of target contigs and $PATH_TO_SPECIES_LIST_FILE is the path to a species list file where contig rnames are mapped to genomes. Those paths do not need to be absolute paths.

 Manifest and species list can be found in a [bucket on OCI](https://cloud.oracle.com/object-storage/buckets/lrbvkel2wjot/dev-relatedness/objects?region=uk-london-1).

### Running Tests
The tests are executed using [nf-test](https://github.com/askimed/nf-test). 

Before running the tests, you will need to ensure appropriate test data is available by setting up the following symlinks in
the root of the repository:

* `data/manifest/manifest_20231001` - a manifest (reference data)
* `data/WTCHG_885333_73205296_1` containing `WTCHG_885333_73205296_1.fastq.gz` and `WTCHG_885333_73205296_1.fastq.gz` - these should be TB FASTQs
* `data/abscessus` containing `file_R1.fastq.gz` and `file_R2.fastq.gz` - these should be NTM FASTQs

To run tests, run the following command

```bash
 nf-test test tests/nextflow/*.test
```

_A copy of a species list file is [included in this repository](test_data/species_list_manifest_20231001.csv) for convenience when running the tests._

### Screenshots

Expected output:
![output](https://github.com/GlobalPathogenAnalysisService/competitivemapping_pipeline/assets/65816841/b43b60f9-87c9-417a-b11c-67ef9c491f2d)


## Python

A Python package that processes the output from command line tools orchestrated by NextFlow is included in this repository. For normal
use it doesn't need to be installed by the user, instructions here are for testing and development. 

### Installation

Clone the repo as described above. Create a virtual environment and install the package using `pip install -e .[dev]`. Set up pre-commit
with `pre-commit install`.

### Execution

Use `process_mapping --help` to get CLI options for running the Python code seperately from NextFlow. Files are included in the repo for
testing e.g.

```
process_mapping --coverage test_data/cov_WTCHG_885333_73205296.tsv --species_list test_data/species_list_manifest_20231001.csv
```

### Outputs

The output from the Python CLI is
validated against a [JSON Schema](src/competitivemapping/competitivemapping.schema.json). [Documentation for the schema](schema_doc.md) 
can be built / updated using:

```
generate-schema-doc src/competitivemapping/competitivemapping.schema.json --config template_name=md
```

### Tests

Tests are written in `pytest`. New code should be covered by tests.

## Integrating to a pipeline

If you want to use the Competitive Mapping Pipeline as subworkflow in your pipeline, use the competitive_mapping named workflow. 
Parameters: 

* Channel with a tuple of _val_, _path_, _path_ 

* _path_ to the manifest file


