"""Script to aid in building manual manifest files."""

import argparse
import gzip
import os

import pandas as pd
from Bio import SeqIO


def make_manifest(
    genomes_df: pd.DataFrame, output_root: str
) -> tuple[str, pd.DataFrame]:
    """Make a manifest from a list of genome fasta files.

    Args:
        genomes_df (pd.DataFrame): Dataframe with columns 'reference' and 'genome_file'
        output_root (str): Path to output root

    Returns:
        tuple[str, pd.DataFrame]: Path to manifest file and dataframe of contigs
    """

    manifest_file = f"{output_root}manifest.fasta.gz"

    # check that all files are unique
    if len(genomes_df["genome_file"].unique()) != len(genomes_df):
        raise ValueError("Duplicate genome files found in manifest")

    # Cat all genomes into a single file
    with open(manifest_file, "wb") as outfile:
        for genome_file in genomes_df["genome_file"]:
            with open(genome_file, "rb") as infile:
                outfile.write(infile.read())

    # Get IDs of genomes
    contigs = []
    for ref, genome_file in zip(genomes_df["reference"], genomes_df["genome_file"]):
        with gzip.open(genome_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                contigs.append(
                    {
                        "reference": ref,
                        "rname": record.id,
                        "length": len(record.seq),
                    }
                )
    contigs_df = pd.DataFrame(contigs)
    # calculate total length per reference
    contigs_df["totallength"] = contigs_df.groupby("reference")["length"].transform(
        "sum"
    )

    contigs_df.to_csv(f"{output_root}contigs.csv", index=False)
    return manifest_file, contigs_df


def make_list_from_dir(genomes_dir: str) -> pd.DataFrame:
    """Make a list of genomes from a directory of gzipped fasta files"""
    genomes = os.listdir(genomes_dir)
    genomes = [g for g in genomes if g.endswith(".gz")]

    def remove_extensions(g):
        result = g
        for ext in [".gz", ".fasta", ".fa", "_genomic", ".fna"]:
            result = result.replace(ext, "")
        return result

    return pd.DataFrame(
        [
            {
                "reference": remove_extensions(g),
                "genome_file": os.path.join(genomes_dir, g),
            }
            for g in genomes
        ]
    )


def main():
    """CLI entry point"""
    parser = argparse.ArgumentParser(
        description="Build a manifest from collection of gzipped fasta genomes."
    )
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "-o", "--output_root", type=str, help="Path to output root"
    )

    subparsers = parser.add_subparsers(dest="subcommand", required=True)
    dir_subparser = subparsers.add_parser(
        "dir", parents=[common_parser], help="Build manifest from directory"
    )
    dir_subparser.add_argument(
        "genomes_dir", type=str, help="Path to directory of genomes."
    )

    list_subparser = subparsers.add_parser(
        "list", parents=[common_parser], help="Build manifest from list of genomes"
    )
    list_subparser.add_argument(
        "list",
        type=str,
        help="Path to list of genomes. Should be a tsv or csv with columns 'reference' and 'file'",
    )
    args = parser.parse_args()

    if args.subcommand == "dir":
        genomes_list_df = make_list_from_dir(args.genomes_dir)
    elif args.subcommand == "list":
        genomes_list_df = pd.read_csv(
            args.list, sep="\t" if args.list.endswith("tsv") else ","
        )
    else:
        raise ValueError(f"Invalid subcommand {args.subcommand}")

    make_manifest(genomes_list_df, args.output_root)


if __name__ == "__main__":
    main()
