# Manifest

The metadata for the manifest is provided in manifest_metadata_20250324.csv

Download files with:
```
python3 download_manifest.py manifest_metadata_20250324.csv reference_genomes
```

Can then join all the references to make a manifest
```
cat reference_genomes/* > all_refs.fasta
```

Lastly need to remove plasmids:
```
python3 remove_plasmids.py all_refs.fasta new_manifest.fasta
```

Beware that this will still be different to the manifest in the knowledge bucket.


Note: manifest is not gzipped for some reason of history

## Dependencies
- ncbi download datasets tool: https://www.ncbi.nlm.nih.gov/datasets/docs/v1/download-and-install/
- pandas
- unzip
- pyfastx

## metadata

Each reference in the manifest has a ncbi assembly accession and the strain it was sequenced from.

Each species in Mykrobe has a list of assembly accessions used for it.

`combined_metadata.csv` list all the references used by manifest, Mykrobe and Tim Walker strains.
Note that some names have been combined where sources use slightly different variants.
The NCBI names (on the right) have been used.
- arupense/arupensis
- heraklionense/heraklionensis
- koreense/koreensis
- kumamotonense/kumamotonensis
- longobardum/longobardus
- minnesotense/minnesotensis
- nonchromogenicum/nonchromogenicus
- parakoreense/parakoreensis
- triviale/trivialis
- virginiense/virginiensis


## Errors

In version 20231001 there was a naming issue:
- saopaulense had the assembly/strains for basiliense
- saopaolens has the assembly/strains for saopaulense. Indeed saopaolens was likely a typo or older name.
This was fixed for 20250324
