# Competitive Mapping Pipeline

Competitive Mapping is an algorithm that compares the sample reads with the references in the manifest and makes a positive selection of the reads matching a specific `rname`.

The pipeline for competitive mapping takes a pair of FASTQ files and outputs the positive filtering of the h37_rv reads (i.e. those reads that are judged to map to the *Mycobacterium tuberculosis* H37RV reference genome) along with unmapped reads, a report with the mapping rank (a list of species in the manifest to which reads have mapped `competitivemapping_report.json`).

### Dependencies
* Docker
* Nextflow

To install the for development:
```bash
conda env create -f env.yml
conda activate competitive_mapping
pip install -e .[dev]
```

### Needed data
* All data for testing is included in `$projectDir/test_data`.
* Myco manifest can be found in the knowledge bucket under `manifest`. See [extra readme](data/manifest/README.md) for details about manifest.
* Flu manifest is under `influenza_virus/manifest`.
* (for sylph) reference databases are found in the knowledge bucket under `sylph`.


## Overview
Competitive mapping uses minimap2 to map reads against manifest (multifasta of reference genomes) and then analyse the resulting bam file.
A species list file is required to match contigs (`rnames`) in the multifasta manifest to a reference as some references have many contigs.

This is co-ordinated by three scripts:

### manifest_builder [In development]
This takes a sylph query report and builds a manifest from it.
This is used for the metagenomic pipeline in development. Myco has a fixed manifest.

parameters:
- sylph_report: Path to the sylph query/profile output file.
- genome_dirs: path to directory containing all the reference genomes. Must contain a `genome_paths.tsv`. Look in knowledge bucket for examples.
- metadata_files: Either the taxonomy file tsv or the metadata csv if ani_groups should be fixed.
- include_whole_genus: if this is set then for each species sylph finds the whole genus will be added to the manifest.

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

### Tie break
This is a rust based replacement to the competitive mapping python script.
See its [readme](tie_break/README.md) for more info.


## Running Nextflow
The workflow takes the following inputs
- input_dir. Path to directory containing input fastq files
- seq_platform. `illumina` or `ont`.
- manifest. Path to manifest file
- species_list. Path to csv which has the contig to genome mapping

When running locally can save outputs by using `--publish_dir`.
And to use locally built container add `-profile local_docker`.

Example using test data:
```bash
nextflow run . \
  --seq_platform illumina \
  --input_dir test_data/samples/illumina/chloro \
  --manifest test_data/myco_manifest/manifest.fasta.gz \
  --species_list test_data/myco_manifest/contigs.csv \
  --publish_dir results \
  -resume
```


By default it will look for files in the input directory based on the following params:
- `params.input_paired_suffix = "*_{1,2}.fastq.gz"`
- `params.input_single_suffix = "*.fastq.gz"`

but these can be overriden. e.g.
```bash
nextflow run ... --input_paired_suffix "tb_sample*_{1,2}.fna.gz"
```

## Testing
There a tests for tie_break, python and nextflow using [nf-test](https://github.com/askimed/nf-test).
All test data is included in the repo.

They are all triggered by running:
```bash
make test
# Or for local container building
make test_local
```

### Outputs

The output from the Python CLI is
validated against a [JSON Schema](src/competitivemapping/competitivemapping.schema.json). [Documentation for the schema](schema_doc.md)
can be built / updated using:

```bash
generate-schema-doc src/competitivemapping/competitivemapping.schema.json --config template_name=md
```

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

### Using a custom tag
If working on a branch you can bump to a custom tag. Once the tag is pushed it will build a container.
```bash
bumper bump -a new_feature_1.0.0
git push
git push --tag
```
Note: only the active version needs to be updated.
