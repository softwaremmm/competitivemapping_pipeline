# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments
"""
Extract reads from fastq file based on alignment csv.
"""

import argparse
import subprocess
import tempfile

import pandas as pd


def extract_reads(
    input_fastqs: list[str],
    alns_csv: str,
    contigs_csv: str,
    reference: str,
    require_both_in_pair: bool,
    output_root: str,
):
    """Extract reads from provided reference"""
    df = pd.merge(
        pd.read_csv(alns_csv),
        pd.read_csv(contigs_csv)[["target_id", "reference"]],
        on="target_id",
    )
    df = df[df["reference"] == reference]

    if require_both_in_pair:
        reads = set(df[df["is_second_in_pair"] == "false"]["query_name"]).intersection(
            set(df[df["is_second_in_pair"] == "true"]["query_name"])
        )
        df = df[df["query_name"].isin(reads)]

    df = df.drop_duplicates(subset=["query_name"])
    read_ids = df["query_name"].astype(str)
    # Get read IDs with "/1" and "/2" suffixes
    read_ids = pd.concat([read_ids, read_ids + "/1", read_ids + "/2"])

    with tempfile.NamedTemporaryFile(mode="w") as temp_file:
        read_ids.to_csv(temp_file.name, index=False, header=False)

        for index, input_fq in enumerate(input_fastqs, 1):
            if len(input_fastqs) == 1:
                output_fq = output_root + ".fastq"
            else:
                output_fq = f"{output_root}_{index}.fastq"

            command = f"seqtk subseq {input_fq} {temp_file.name} > {output_fq}"
            print(f"Running command: {command}")
            subprocess.run(command, shell=True, check=False)


def cli_entry_point():
    """Entry point for the command line interface."""
    parser = argparse.ArgumentParser(
        description="Selected reads from fastq files based on query ids."
    )
    parser.add_argument(
        "-f",
        "--fastq",
        type=str,
        required=True,
        help="Path to the input fastq(s).",
        nargs="+",
    )
    parser.add_argument(
        "-a",
        "--alns_csv",
        type=str,
        required=True,
        help="Path to the input alignments CSV file.",
    )
    parser.add_argument(
        "-c",
        "--contigs",
        type=str,
        required=True,
        help="Path to csv mapping contigs to references.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=str,
        required=True,
        help="Reference name to filter reads for.",
    )
    parser.add_argument(
        "-o",
        "--output_root",
        type=str,
        required=True,
        help="Path to the output fastq root (before _1.fastq)",
    )
    parser.add_argument(
        "--require_both_in_pair",
        action="store_true",
        help="Require both reads in a pair to be present in the output.",
        default=False,
    )
    args = parser.parse_args()

    extract_reads(
        args.fastq,
        args.alns_csv,
        args.contigs,
        args.reference,
        args.require_both_in_pair,
        args.output_root,
    )


if __name__ == "__main__":
    cli_entry_point()
