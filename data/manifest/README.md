# Manifest

The metadata for the manifest is provided in manifest_metadata_20231001.csv

Download files with:
```
python3 download_manifest.py manifest_metadata_20231001.csv reference_genomes
```

Can then join all the references to make a manifest
```
cat reference_genomes/* > new_manifest.fasta
```

Note: manifest is not gzipped for some reason of history

## Dependencies
- ncbi download datasets tool: https://www.ncbi.nlm.nih.gov/datasets/docs/v1/download-and-install/
- pandas
- unzip

## metadata

## Errors

In version 20231001 there was a naming issue:
- saopaulense had the assembly/strains for basiliense
- saopaolens has the assembly/strains for saopaulense. Indeed saopaolens was likely a typo or older name.
This was fixed for 20250324
