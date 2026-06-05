# Competitive Mapping Pipeline

Competitive Mapping uses Minimap2 to map fastq reads against multiple references at the same time, in order to assign reads to their correct reference.
Competitive Mapping will produce a summary csv and json file detailing how many reads each reference got as well as the resulting coverage.

Within the Myco pipeline the reads assigned to TB (The H37RV reference) and unmapped reads get put into a fastq file to be used for later assembly.

There are two workflows:
1. Standard: This takes a fixed manifest (multi fasta file) of references and a csv mapping the contigs (`rname`) to the actual reference names.
2. Dynamic: This takes a larger database and uses sylph (profile) to select likely references before mapping. The param `fixed_refs` can be used to provide a comma separated list of accessions to always include.

### Dependencies
* Docker
* Nextflow

To install the for development:
```bash
conda env create -f env.yml
conda activate competitive_mapping
pip install -e .[dev]
```

### Running Nextflow
The [Makefile](Makefile) includes several examples of how to run the nextflow.

When running locally can save outputs by using `--publish_dir`.
And to use locally built container add `-profile local_docker`.


By default it will look for files in the input directory based on the following params:
- `params.input_paired_suffix = "*_{1,2}.fastq.gz"`
- `params.input_single_suffix = "*.fastq.gz"`

but these can be overriden. e.g.
```bash
nextflow run ... --input_paired_suffix "tb_sample*_{1,2}.fna.gz"
```

The dynamic workflow can take lists of databases. When running locally this means that you can set `--ref_genome_dirs`, `--sylph_dbs` and `--taxonomy_files` to a comma separated list of file paths. Alternatively, you can set the params to a list within a `nextflow.config`.

### Needed data
* All data for testing is included in `$projectDir/test_data`.
* (Myco) Myco manifest can be found in the knowledge bucket under `manifest`. See [extra readme](data/manifest/README.md) for details about manifest.
* (Flu) Flu manifest is under `influenza_virus/manifest`.
* (Dynamic) reference databases are found in the knowledge bucket under `sylph`. For more details and how to create them see [reference_set_manager](https://github.com/softwaremmm/reference_set_manager).


## Overview
1. (Dynamic) Run sylph against database to produce a sylph profile and query report
2. (Dynamic) Use `manifest_builder` to create a manifest fasta and contigs csv from a sylph report
3. Use `manifest_mapper` to run minimap2 with the desired settings
4. Use `competitive_mapping` to investigate the bam file and calculate the output stats

### manifest_builder
Builds the manifest from a sylph report.

parameters:
- sylph_report: Path to the sylph query/profile output file.
- genome_dirs: path to directory containing all the reference genomes. Must contain a `genome_paths.tsv`. Look in knowledge bucket for examples.
- taxonomy_files: csv/tsv with taxonomy info for each reference used. Can optionally have an "ani_group" column.
- ani_threshold: if set, will use this threshold for grouping similar references. Will not be used if taxonomy file has "ani_group" column already.
- fixed_refs: optional comma separated list of references to always include.

### manifest_mapper
Map reads against the manifest using Minimap2.

parameters:
- manifest: path to multifasta
- reads: path to input fastqs
- seq_platform: `ont` or `illumina`
- filter_secondary/filter_supplementary: Will exclude secondary/supplementary alignments from resulting bam.

### competitive_mapping
This takes the bam file and calculates various coverage stats for the references in the manifest. It:
- Run samtools coverage and aggregate results with `process_coverage.py`
- Calculate alignment stats from bam file use `process_aln_stats.py`
- Produce overall csv and report json

parameters:
- input_bam: bam from previous step
- contigs: file with mapping from contigs to reference. Often called species_list.
- seq_platform: `ont` or `illumina`
- ref_for_fastq: If set then the reads which mapped to the provided reference will be extracted from the bam file. Used in myco to select TB reads.
Note: currently unmapped reads will also be extracted by default


## Testing
There a tests for python and nextflow using [nf-test](https://github.com/askimed/nf-test).
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

In the future we could change this to output a read pair as long as **either** read in a pair map to h37rv.


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
