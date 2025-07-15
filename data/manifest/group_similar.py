"""Find reference genomes groups based on ANI threshold."""
import argparse

import networkx as nx
import pandas as pd


def find_groups(edge_list, ref_file, ani_threshold, output):
    """Find groups of similar genomes based on ANI threshold."""
    df = pd.read_csv(edge_list, sep="\t")
    refs = pd.read_csv(ref_file)

    df = df[df["ANI"] >= ani_threshold]

    df["Ref_file"] = df["Ref_file"].str.split("/").str[-1].str.replace(".fasta.gz", "")
    df["Query_file"] = (
        df["Query_file"].str.split("/").str[-1].str.replace(".fasta.gz", "")
    )

    name_lookup = refs.set_index("assembly_accession")["reference"].to_dict()

    df.sort_values(["Ref_file", "Query_file"], inplace=True)

    # build graph from filtered data
    graph = nx.from_pandas_edgelist(df, "Ref_file", "Query_file")

    groups = []
    for i, c in enumerate(nx.connected_components(graph)):
        for genome in c:
            groups.append([genome, name_lookup[genome], i])
    groups_df = pd.DataFrame(groups, columns=["genome", "organism_name", "group"])
    groups_df.to_csv(output, index=False)


def main():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Group similar genomes")
    parser.add_argument("edge_list", help="tsv file with edge list")
    parser.add_argument("references", help="csv file with reference genomes")
    parser.add_argument("output", help="output csv file with groups")
    parser.add_argument(
        "--threshold", type=float, default=95, help="ANI threshold for grouping"
    )
    args = parser.parse_args()
    find_groups(args.edge_list, args.references, args.threshold, args.output)


if __name__ == "__main__":
    main()
