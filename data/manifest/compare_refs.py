"""Get references from NCBI, Mykrobe, and Walker databases."""

import argparse

import pandas as pd
from k2_analysis.kraken_tree import Tree, ncbi_taxonomy_to_tree

LEVEL_CODES = {
    "superkingdom": "d",
    "phylum": "p",
    "class": "c",
    "order": "o",
    "family": "f",
    "genus": "g",
    "species": "s",
}


def get_species_taxid(node: Tree | None) -> int:
    """Get the species taxid from the node."""
    if node is None:
        return 0

    if node.level_code != "species":
        return get_species_taxid(node.parent)

    return int(node.taxid)


mykrobe_renaming_map = {
    "algericum": "algericus",
    "arupense": "arupensis",
    "bovis": "tuberculosis",
    # "orygis": "tuberculosis",
    # "mungi": "tuberculosis",
    "heraklionense": "heraklionensis",
    "koreense": "koreensis",
    "kumamotonense": "kumamotonensis",
    "longobardum": "longobardus",
    "minnesotense": "minnesotensis",
    "nonchromogenicum": "nonchromogenicus",
    "parakoreense": "parakoreensis",
    "triviale": "trivialis",
    "virginiense": "virginiensis",
    "canettii": "canetti",
    "senuense": "senuensis",
    "sinense": "sinensis",
}


def rename_mykrobe_species_name(name: str) -> str:
    """Rename Mycobacterium species names to a simpler format."""
    name = name.replace("Mycobacterium_", "")
    name = mykrobe_renaming_map.get(name, name)

    name = name.split("_")[0]
    name = mykrobe_renaming_map.get(name, name)

    return f"M.{name}"


def get_mykrobe_references(mykrobe_metadata: str) -> pd.DataFrame:
    """Get summary of Mykrobe references."""
    mykrobe = pd.read_csv(mykrobe_metadata, sep="\t")
    mykrobe["reference"] = mykrobe["Species"].apply(rename_mykrobe_species_name)
    mykrobe = (
        mykrobe.groupby("reference")["Accessions"]
        .apply(",".join)
        .reset_index(name="mykrobe_accessions")
    )
    return mykrobe


def get_simple_name(node: Tree | None) -> str:
    """Simplify the name by removing the genus and keeping only the species."""
    if node is None:
        return "unknown_species"

    name = str(node.name)
    name = name.replace("Candidatus", "").strip()

    name = name.split(" ", maxsplit=1)[-1]
    name = f"M.{name}"

    return name


def get_walker_references(walker_metadata: str, taxa_to_node: dict) -> pd.DataFrame:
    """Get summary of Walker references."""
    walker_df = pd.read_csv(walker_metadata).rename(
        columns={
            "tax_id": "taxid",
        }
    )

    walker_df["species_taxid"] = walker_df["taxid"].apply(
        lambda x: get_species_taxid(taxa_to_node.get(x, None))
    )

    walker_df["reference"] = walker_df["species_taxid"].apply(
        lambda x: get_simple_name(taxa_to_node.get(x, None))
    )

    walker_df = walker_df[
        [
            "reference",
            "species_taxid",
            "taxid",
            "walker_taxon",
            "strain",
            "run_accession",
        ]
    ].rename(
        columns={
            "taxid": "walker_taxid",
            "strain": "walker_strain",
            "run_accession": "walker_run_accession",
        }
    )
    return walker_df


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
    parser.add_argument("walker_metadata", type=str, help="Path to the metadata file")
    parser.add_argument("mykrobe_metadata", type=str, help="Path to the metadata file")
    args = parser.parse_args()

    # Load NCBI taxonomy
    taxonomy = ncbi_taxonomy_to_tree(
        f"{args.ncbi_taxonomy_dir}/nodes.dmp", f"{args.ncbi_taxonomy_dir}/names.dmp"
    )
    taxa_to_node = {taxon.taxid: taxon for taxon in taxonomy.to_list()}
    print("Taxonomy loaded")

    # Load ncbi references
    ncbi_df = pd.read_csv(args.ncbi_metadata)

    mykrobe_df = get_mykrobe_references(args.mykrobe_metadata)
    mykrobe_df.to_csv("refs_mykrobe.csv", index=False)

    walker_df = get_walker_references(args.walker_metadata, taxa_to_node)
    walker_df.to_csv("refs_walker.csv", index=False)

    ## Combine
    index_cols = [
        "reference",
    ]

    ncbi_df = ncbi_df[index_cols]
    ncbi_df["ncbi"] = True

    mykrobe_df = mykrobe_df[index_cols]
    mykrobe_df["mykrobe"] = True

    walker_df = walker_df[index_cols]
    walker_df["walker"] = True

    df = ncbi_df.merge(mykrobe_df, on=index_cols, how="outer").merge(
        walker_df, on=index_cols, how="outer"
    )
    df = df.fillna(False)
    df = df.sort_values(by="reference").reset_index(drop=True).drop_duplicates()
    df.to_csv("refs_sources_comparison.csv", index=False)


if __name__ == "__main__":
    main()
