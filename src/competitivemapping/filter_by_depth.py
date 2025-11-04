import argparse
from pathlib import Path

import pandas as pd


def filter_by_depth(
    tie_break_report: Path, species_list: Path, min_depth: float
) -> pd.DataFrame:
    """Filter references by mean depth from a tie break report. Output references and accessions."""
    tie_break = pd.read_csv(tie_break_report)
    filtered = tie_break[
        (tie_break["mean_depth"] >= min_depth) & (tie_break["depth_type"] == "unique")
    ]
    # Print statents for debugging prototype, remove later
    pd.set_option("display.max_columns", None)
    print(filtered.dtypes)
    print(filtered)
    contig_to_reference = pd.read_csv(species_list)
    contig_to_reference = contig_to_reference[
        ["reference", "assembly_accession", "ref_for_assembly"]
    ].drop_duplicates()
    print(contig_to_reference.dtypes)
    print(contig_to_reference)
    references = filtered.merge(contig_to_reference, on="reference", how="left")

    references = references[["reference", "assembly_accession", "ref_for_assembly"]]

    return references


def write_list(items: list[str], output_file: Path):
    with open(output_file, "w") as f:
        for item in items:
            f.write(f"{item}\n")


def cli_entry_point():
    """Entry point for the command line interface."""
    parser = argparse.ArgumentParser(
        description="Selected reads from fastq files based on query ids."
    )
    parser.add_argument(
        "-t",
        "--tie_break_report",
        type=str,
        required=True,
        help="Path to the input tie break report CSV file.",
    )
    parser.add_argument(
        "-s",
        "--species_list",
        type=str,
        required=False,
        help="Path to a list of species to consider.",
        default=None,
    )
    parser.add_argument(
        "-m",
        "--min_depth",
        type=float,
        required=True,
        help="Minimum mean depth to retain a reference.",
    )
    parser.add_argument(
        "-o", "--output", type=str, help="Output prefix", default="high_depth"
    )

    args = parser.parse_args()

    references = filter_by_depth(
        args.tie_break_report,
        args.species_list,
        args.min_depth,
    )

    references.to_csv(Path(args.output + "_refs.txt"), header=False, index=False)


if __name__ == "__main__":
    cli_entry_point()
