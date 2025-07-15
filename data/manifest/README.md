# Manifest

## Dependencies
- ncbi download datasets tool: https://www.ncbi.nlm.nih.gov/datasets/docs/v1/download-and-install/
`conda install -c conda-forge ncbi-datasets-cli`
- pandas
- unzip
- pyfastx

Note: original manifest is not gzipped for some reason of history.

## Downloading manifest from metadata

The metadata for the manifest is provided in manifest_metadata files.
Note: Beware that this may still be different to the manifest in the knowledge bucket.

Download files with:
```bash
python3 download_manifest.py manifest_metadata_20250709.csv reference_genomes
```

If wanted can remove plasmids:
```bash
python3 remove_plasmids.py reference_genomes/*
```

Can then join all the references to make a manifest using
```bash
# Need to gzip all first
pigz reference_genomes/*
manifest_builder_manual reference_genomes manifest_metadata_20250709.csv \
    --output_root new --cpus 20 --ani_threshold 97
```

## metadata

Each reference in the manifest has a ncbi assembly accession and the strain it was sequenced from.
The `reference` name is a simplified name to avoid genus confusion.

Each species in Mykrobe has a list of assembly accessions used for it.

## Updating Manifest (WIP)
Current manifest seems to be a fair bit behind that available from RefSeq.
Particularly there is a improved gordonae reference out.

From NCBI/RefSeq we need:
- The list of reference genomes [assembly_summary_refseq.tsv](https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/).
- The taxonomy [taxdump.tar.gz](https://ftp.ncbi.nih.gov/pub/taxonomy/). This should be extracted to a folder.
Both of these are updated regularly, so note when they were collected and aim to be at the same time.

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt

mkdir ncbi_taxonomy
cd ncbi_taxonomy
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf taxdump.tar.gz
rm taxdump.tar.gz
cd ..
```

Also need to install the following to add taxonomy data in python:
```bash
pip install git+https://github.com/GlobalPathogenAnalysisService/kraken2_analysis.git
```

Finally run code:
```bash
python3 get_ncbi_myco_references.py ncbi_taxonomy assembly_summary_refseq.txt
```

This produces `refs_ncbi.csv` with all the RefSeq Myco references.
Any extra references can be added now. The required columns is `reference` and `assembly_accession`.

#### Additions

Tim has requested:
```
M.intracellulare_yongonense,GCF_000418535.1,1767,1138871,Mycobacterium intracellulare subsp. yongonense,strain=05-1390 (KCTC_19555),,,,,,,
M.intracellulare_chimaera,GCF_002219285.1,1767,222805,Mycobacterium intracellulare subsp. chimaera,strain=DSM_44623,,,,,,,
```


### Comparing
To compare this to other sources of reference (mykrobe and tim walker) run:

```bash
python3 compare_refs.py ncbi_taxonomy refs_ncbi.csv other_sources/tim_walker_samples.csv other_sources/mykrobe_accessions_202309.tsv
```

### Mykrobe
This seems to have used GTDB at some point (hence the _A and _B species), so does not match NCBI always.
It also seems to have different species names in places, which get corrected by the code:

e.g.
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

The list of genomes (mykrobe_accessions_202309.tsv) was compiled by Tim Peto.

Note: it contains many non Myco species.


### Tim Walker
The tim walker collection comes from strain collections.
All the reads are on ENA under project PRJNA1169681 (see tim_walker_samples.csv).

All samples have the strain id, the name assigned by Tim Walker, and the NCBI name.

Note: some samples like NCTC_00525 seem to have issues as Walker called it smegmatis, but looking online it seems NCTC_00525 should be phlei.

Note: There are multiple copies of some species. For instance of Mycobacterium avium subsp. avium or smegmatis.


## Previous Errors

In version 20231001 there was a naming issue:
- saopaulense had the assembly/strains for basiliense
- saopaolens has the assembly/strains for saopaulense. Indeed saopaolens was likely a typo or older name.
This was fixed for 20250324

## ANI similarity in manifest
For finding ANI similarity between references in the manifest use the following.
```
skani triangle -t 100 -s 90 --medium reference_genomes/*.fasta -E > manifest_edge_list.tsv

python3 group_similar.py manifest_edge_list.tsv manifest_metadata_20250709.csv manifest_groups.csv --threshold 97
```

can adjust threshold to group at different cutoffs.
