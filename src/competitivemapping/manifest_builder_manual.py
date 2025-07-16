# pylint: disable=too-many-locals
"""Build manifest from directory of genomes.
Used when making manifest manually, e.g. for Myco"""

import argparse
import os
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

from competitivemapping.manifest_builder import (
    Config,
    assign_ani_groups,
    get_ani_distances,
    read_contigs,
)


def dir_to_df(genome_directory: str) -> pd.DataFrame:
    """Make a list of genomes from a directory of (gzipped) fasta files."""
    genomes = os.listdir(genome_directory)
    df = pd.DataFrame({"path": genomes})
    df["accession"] = (
        df["path"]
        .str.replace(".gz", "", regex=False)
        .str.replace(".fasta", "", regex=False)
    )
    df["path"] = df["path"].apply(lambda x: os.path.join(genome_directory, x))
    return df


def cli_entry_point():
    """Entry point for the CLI"""
    parser = argparse.ArgumentParser(
        description="Build manifest from directory of genomes."
    )
    # Need to provide parent directory to work with nextflow symlinks
    parser.add_argument(
        "genome_directory",
        help="Path to directory of genomes. Should contain gzipped fasta files.",
        type=str,
    )
    parser.add_argument(
        "metadata",
        help="path to file with columns assembly_accession, and reference",
        type=str,
    )
    parser.add_argument(
        "--ani_threshold",
        help="ANI threshold for grouping genomes",
        default=97.0,
        type=float,
    )
    parser.add_argument(
        "-c", "--cpus", help="Number of CPUs to use", default=4, type=int
    )
    parser.add_argument(
        "-o", "--output_root", required=True, help="Path to the output files"
    )

    args = parser.parse_args()

    config = Config(
        cpus=args.cpus,
        include_whole_genus=False,
        ani_threshold=args.ani_threshold,
        output_root=args.output_root,
    )

    refs_df = dir_to_df(args.genome_directory)
    with open(f"{config.output_root}manifest.fasta.gz", "wb") as outfile:
        for filepath in refs_df["path"]:
            if not filepath.endswith(".gz"):
                raise ValueError(
                    f"Expected gzipped files, but found {filepath} without .gz extension."
                )
            with open(filepath, "rb") as infile:
                outfile.write(infile.read())

    # Read contigs in parallel
    with ProcessPoolExecutor(max_workers=config.cpus) as executor:
        results = list(
            executor.map(read_contigs, zip(refs_df["accession"], refs_df["path"]))
        )

    contigs_df = pd.DataFrame([contig for result in results for contig in result])
    contigs_df["totallength"] = contigs_df.groupby("reference")["length"].transform(
        "sum"
    )

    # add ani information
    ani_df = get_ani_distances(refs_df, config)
    contigs_df = assign_ani_groups(contigs_df, ani_df, config.ani_threshold)

    contigs_df.rename(columns={"reference": "assembly_accession"}, inplace=True)

    meta = pd.read_csv(args.metadata)[["assembly_accession", "reference"]]

    contigs_df = contigs_df.merge(meta, on="assembly_accession", how="left")
    contigs_df.sort_values(by=["reference", "rname"], inplace=True)
    contigs_df.to_csv(f"{config.output_root}contigs.csv", index=False)


if __name__ == "__main__":
    cli_entry_point()
