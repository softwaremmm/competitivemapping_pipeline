## New

- Split the minimap2 and the analysis steps into separate python scripts
- Breaking: Make reference name a parameter so reads mapping to references other than _M. tuberculosis_ can be output if needed
- Tests for reference name as a parameter, and for flu specifically

## 2.0.0 (2025-02-12)

### Feat
- Dynamic competitive mapping. This system uses sylph output report to select references (gzipped fastas) from a directory to create a manifest on-the-fly
- Calculate coverage and meandepth both without and without secondary reads
- Helper scripts for downloading from NCBI and making manifest manually
- Produce a csv as well as json comparison report

### Refactor
- Move competitive mapping into python controlled code
- Move to gitflow system, this includes a system for building commit hash based releases on PR to develop
- Use nextflow language server style
- improve local running options by adding publish_dir, simpler configs and examples

### Fix
- remove cm_error channel as not used
- Output details for all references in manifest, even if no reads mapped to them
- Change cpu from 8 to 4. Make memory a steep ramp up starting from 12G
- Run pytest within docker container
- better clean up for CI test workflow
- added standard pre-commits
- use a conda env.yml
- rename test container and param to "test_container_cm"

### Other
- Added folder with metadata, readme, and scripts about the myco manifest
- added notes on the bam to fastq step as it has some confusing details


## 1.3.1 (2024-08-23)

### Fix

- Don't throw errors for missing values
- Return false on read threshold error
- check_read_count CLI params

### Refactor

- More concise Python
- Drop jq from Docker
- Check read count using Python, not jq

## 1.3.0 (2024-07-15)

### Feat

- also push to gpasltd ocr

## 1.2.6 (2024-07-10)

### Fix

- update species list

## 1.2.5 (2024-07-09)

### Fix

- make secondary cov optional and remove from pipeline

## 1.2.4 (2024-07-08)

### Fix

- try using sam file for coverage

## 1.2.3 (2024-07-04)

### Fix

- use -o instead of >

## 1.2.2 (2024-07-03)

### Fix

- test using knowledge bucket
- leave test mem at 16
- try higher memory

## 1.2.1 (2024-07-01)

### Fix

- ensure empty secondary reads is valid

## 1.2.0 (2024-06-27)

### Feat

- add overall summary and secondary coverage stats

## 1.1.2 (2024-06-03)

### Fix

- Limit RAM when testing

## 1.1.1 (2024-05-08)

### Fix

- order competitive mapping output by meandepth

## 1.1.0 (2024-04-09)

### Feat

- adds more information about the alignments vs reads to the summary json

### Fix

- better read breakdown
- do not update numreads
- account for secondary reads in alignment summary
- add read summary function

## 1.0.2 (2024-01-18)

### Fix

- spelling

## 1.0.1 (2024-01-18)

### Fix

- include template in minimal files

## 1.0.0 (2024-01-18)

### Feat

- support for ont reads through seq_platform parameter

### Fix

- add trap to process template
- seq_platform now a string

## 0.7.0 (2024-01-16)

### Feat

- commitizen for auto build and release

### Fix

- update version files

## 0.6.3 (2024-01-16)

### Fix

- remove integrate subworkflows action

## 0.6.2 (2024-01-05)

- Manual bump to re-sync commitizen

## v0.6.2 (2023-12-13)

### Fix

- **s3fs**: ensure buckets are unmounted on exit

## 0.6.1 (2023-11-20)

### Fix

- Increase resources

## 0.6.0 (2023-10-30)

### Feat

- Added pod labels to processes

## 0.5.3 (2023-10-17)

### Fix

- Use latest checkout action
- Don't skip CI

## 0.5.2 (2023-10-17)

### Fix

- Rename reads output

## 0.5.1 (2023-10-16)

### Fix

- Unmount buckets when done

## 0.5.0 (2023-10-10)

### Feat

- Use Python to process output
- Accept species_list file

### Fix

- Allocate more resources
- Use correct variable name

## 0.4.1 (2023-10-10)

### Fix

- commitizen ci

## 0.4.0 (2023-10-06)

### Feat

- Add schema checking of output
- Add cli entrypoint

### Fix

- Sort output
- pre-commit and mypy
- Add JSON Schema

## 0.3.0 (2023-10-06)

### Feat

- Enable CLI use
- Function to determine overall coverage
- Contig coverage aggreation function
- QA checks on species list
- Function to add references to samtools out
- Empty CLI function
- Add cli entrypoint
- Add cli args class
- Accept species_list file

### Fix

- meandepth calculation
- Match previous field names
- Better naming
- Add fields to output
- Add species list for testing
- Remove symlink

## 0.2.2 (2023-10-06)

### Fix

- Typo for the sake of commitizen

## 0.2.1 (2023-10-05)

### Fix

- Silent change to trigger commitzen bump in CI

## 0.2.0 (2023-10-05)

### Fix

- Better starting version
- Docker action
- Reference correct file
- Correct YAML
- Use versioned container for has_enough_reads
- Docker action exactly as per summary pipeline
- Docker action as summary pipeline
- Use initial version 0.1.0 in nextflow
- Start with version 0.1.0
- Enable commitizen version bumps
- Establish Python project
- Remove redundant readme

## v0.1.12 (2023-10-02)

### Fix

- Extract unmapped reads, merge with TB reads

## v0.1.11 (2023-09-20)

### Fix

- fix typo in container image

## v0.1.10 (2023-09-20)

## v0.1.9 (2023-09-06)

### Fix

- Increase memory request (fix k8s runs)

## v0.1.8 (2023-08-30)

## v0.1.7 (2023-08-21)

## v0.1.6 (2023-08-16)

## v0.1.5 (2023-08-16)

## v0.1.4 (2023-08-15)

## v0.1.3 (2023-08-14)

## v0.1.2 (2023-08-08)

## v0.1.1 (2023-08-02)

## 0.1.0 (2023-07-14)
