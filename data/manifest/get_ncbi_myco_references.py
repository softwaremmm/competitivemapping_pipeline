"""Get references from NCBI, Mykrobe, and Walker databases."""

import argparse

import pandas as pd
from k2_analysis.kraken_tree import Tree, ncbi_taxonomy_to_tree

PROVIDED_COLS = [
    "assembly_accession",
    "bioproject",
    "biosample",
    "wgs_master",
    "refseq_category",
    "taxid",
    "species_taxid",
    "organism_name",
    "infraspecific_name",
    "isolate",
    "version_status",
    "assembly_level",
    "release_type",
    "genome_rep",
    "seq_rel_date",
    "asm_name",
    "asm_submitter",
    "gbrs_paired_asm",
    "paired_asm_comp",
    "ftp_path",
    "excluded_from_refseq",
    "relation_to_type_material",
    "asm_not_live_date",
    "assembly_type",
    "group",
    "genome_size",
    "genome_size_ungapped",
    "gc_percent",
    "replicon_count",
    "scaffold_count",
    "contig_count",
    "annotation_provider",
    "annotation_name",
    "annotation_date",
    "total_gene_count",
    "protein_coding_gene_count",
    "non_coding_gene_count",
    "pubmed_id",
]

SELECTED_COLS = [
    "assembly_accession",
    "bioproject",
    "biosample",
    "refseq_category",
    "taxid",
    "species_taxid",
    "group",
    "organism_name",
    "infraspecific_name",
    "assembly_level",
    "genome_size",
    "contig_count",
]

LEVEL_CODES = {
    "superkingdom": "d",
    "phylum": "p",
    "class": "c",
    "order": "o",
    "family": "f",
    "genus": "g",
    "species": "s",
}


def get_taxonomy(node: Tree | None) -> list[str]:
    """Produce list with taxonomy heirarchy."""
    if node is None:
        return ["unknown_taxa"]

    if node.level_code in ["R"] or node.name == "cellular organisms":
        return []

    level_code = LEVEL_CODES.get(node.level_code, "x")

    taxa_string = level_code + "__" + str(node.name)
    return get_taxonomy(node.parent) + [taxa_string]


def get_simple_name(node: Tree | None) -> str:
    """Simplify the name by removing the genus and keeping only the species."""
    if node is None:
        return "unknown_species"

    name = str(node.name)
    name = name.replace("Candidatus", "").strip()

    name = name.split(" ", maxsplit=1)[-1]
    name = f"M.{name}"

    return name


def main():
    """Summarise NCBI, Mykrobe, and Walker references."""
    parser = argparse.ArgumentParser(
        description="Select reference genomes from refseq metadata"
    )
    parser.add_argument(
        "ncbi_taxonomy_dir",
        type=str,
        help="Path to directory with names.dmp and nodes.dmp",
    )
    parser.add_argument("ncbi_metadata", type=str, help="Path to the metadata file")
    args = parser.parse_args()

    # Load NCBI taxonomy
    taxonomy = ncbi_taxonomy_to_tree(
        f"{args.ncbi_taxonomy_dir}/nodes.dmp", f"{args.ncbi_taxonomy_dir}/names.dmp"
    )
    taxa_to_node = {taxon.taxid: taxon for taxon in taxonomy.to_list()}
    print("Taxonomy loaded")

    # Load ncbi references
    ncbi_df = pd.read_csv(
        args.ncbi_metadata, sep="\t", comment="#", names=PROVIDED_COLS
    )[SELECTED_COLS]
    ncbi_df = ncbi_df[ncbi_df["group"].isin(["bacteria"])]
    ncbi_df = ncbi_df[ncbi_df["refseq_category"] == "reference genome"]

    ncbi_df["taxonomy"] = ncbi_df["taxid"].apply(
        lambda x: ";".join(get_taxonomy(taxa_to_node.get(x, None)))
    )
    ncbi_df = ncbi_df[ncbi_df["taxonomy"] != "unknown_taxa"]

    # Filter for Mycobacteriaceae family
    ncbi_df = ncbi_df[ncbi_df["taxonomy"].str.contains("f__Mycobacteriaceae")]

    ncbi_df["reference"] = ncbi_df["species_taxid"].apply(
        lambda x: get_simple_name(taxa_to_node.get(x, None))
    )
    ncbi_df.sort_values(by="reference", inplace=True)

    ncbi_df.rename(
        columns={
            "infraspecific_name": "strain",
        }
    )

    final_cols = [
        "reference",
        "assembly_accession",
        "species_taxid",
        "taxid",
        "organism_name",
        "strain",
        "bioproject",
        "biosample",
        "assembly_level",
        "genome_size",
        "contig_count",
        "taxonomy",
        "download",
    ]
    ncbi_df = ncbi_df[final_cols]

    ncbi_df.to_csv("refs_ncbi.csv", index=False)


if __name__ == "__main__":
    main()
