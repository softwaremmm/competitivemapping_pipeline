# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments
# pylint: disable=too-many-locals
# pylint: disable=duplicate-code
"""Run Competitive Mapping and aggregate the results"""

import argparse
import logging
import subprocess

import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    format="%(asctime)s — %(relativeCreated)d — %(levelname)s — %(funcName)s:%(lineno)d — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S%z",
    level=logging.DEBUG,
)


def get_megahit_contig_stats(contigs: list[str]) -> pd.DataFrame:
    """Reads metadata from contig fasta header lines"""
    contig_stats = []
    for contig in contigs:
        with open(contig, "rt", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                attributes = record.description.split(" ")
                stats = {"qseqid": record.id}
                for attr in attributes:
                    if "multi=" in attr:
                        stats["multi"] = float(attr.split("=")[1])
                    elif "len=" in attr:
                        stats["length"] = int(attr.split("=")[1])
                contig_stats.append(stats)

    df = pd.DataFrame(contig_stats)
    return df


def blastn_contigs(
    manifest: str, contigs: list[str], cpus: int, output_root: str
) -> pd.DataFrame:
    """Run blastn to map contigs against manifest, returning sorted bam.

    Args:
        manifest (str): Path to the manifest file
        contigs (list[str]): List of paths to the contigs
        cpus (int): number of cores to use
        output_root (str): Path to the output root

    Returns:
        pd.DataFrame: dataframe of blast results
    """
    logging.info("Making blast db")
    blast_db = f"{output_root}blast_db"
    # if manifest is gzipped then unzip first
    if manifest.endswith(".gz"):
        command = (
            f"gzip -dc {manifest} | makeblastdb -in - -dbtype nucl -parse_seqids"
            + f" -out {blast_db} -title manifest_db"
        )
    else:
        command = f"makeblastdb -in {manifest} -dbtype nucl -parse_seqids -out {blast_db} -title manifest_db"
    subprocess.run(
        command,
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    logging.info("Running blastn")
    blast_mapping = f"{output_root}blast_mapping.tsv"
    command = f"blastn -db {blast_db} -query {' '.join(contigs)} -outfmt 6 -num_threads {cpus} -out {blast_mapping}"
    subprocess.run(
        command,
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    blast_df = pd.read_csv(
        blast_mapping,
        sep="\t",
        names=[
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ],
    )
    blast_df.to_csv(f"{output_root}blast_mapping.tsv", sep="\t", index=False)
    # Only keep best hit per contig, use highest bitscore
    blast_df.sort_values("bitscore", ascending=False, ignore_index=True, inplace=True)

    blast_df = blast_df.drop_duplicates(subset=["qseqid"])

    return blast_df


def competitive_map_contigs(
    mainfest: str,
    contigs_df: pd.DataFrame,
    contigs: list[str],
    cpus: int,
    output_root: str,
) -> pd.DataFrame:
    """Map contigs to manifest and aggregate results

    Args:
        mainfest (str): Path to the manifest file
        contigs_df (pd.DataFrame): DataFrame of contigs metadata
        contigs (list[str]): Paths to contig fasta files
        cpus (int): Number of cores to use
        output_root (str): Path to the output root

    Returns:
        pd.DataFrame: DataFrame of aggregated results
    """
    contig_multiplicity = get_megahit_contig_stats(contigs).rename(
        columns={"length": "contig_length", "multi": "depth"}
    )
    blast_df = blastn_contigs(mainfest, contigs, cpus, output_root).rename(
        columns={"sseqid": "rname"}
    )[["qseqid", "rname"]]
    blast_df = blast_df.merge(contig_multiplicity, on="qseqid", how="inner")

    df = contigs_df.merge(blast_df, on="rname", how="inner")
    aggregated = (
        df.groupby("reference")
        .apply(
            lambda ref: pd.Series(
                {
                    "meandepth": (ref.contig_length * ref.depth).sum()
                    / ref.totallength.iloc[0],
                }
            ),
            include_groups=False,
        )
        .reset_index()
    )

    result_df = contigs_df[["species", "reference", "totallength"]].drop_duplicates()
    result_df = result_df.merge(aggregated, on="reference", how="left").fillna(0)
    result_df.rename(columns={"reference": "genome_name"}, inplace=True)
    result_df.sort_values("meandepth", ascending=False, inplace=True)
    return result_df


def cli_entry_point():
    """Entry point for the CLI"""
    parser = argparse.ArgumentParser(description="Process sylph report and genomes.")

    parser.add_argument("--manifest", required=True, help="Path to the manifest file")
    parser.add_argument(
        "--manifest_contigs", required=True, help="Path to the contigs file"
    )
    parser.add_argument(
        "--contigs", required=True, help="Fasta files of contigs to map", nargs="+"
    )
    parser.add_argument("--cpus", help="Number of CPUs to use", default=4, type=int)
    parser.add_argument("--output_root", required=True, help="Path to the output files")

    args = parser.parse_args()

    manifest_contigs_df = pd.read_csv(args.manifest_contigs)

    result_df = competitive_map_contigs(
        args.manifest, manifest_contigs_df, args.contigs, args.cpus, args.output_root
    )

    result_df.to_csv(f"{args.output_root}species_comparison.csv", index=False)


if __name__ == "__main__":
    cli_entry_point()
