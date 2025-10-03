import argparse
from pathlib import Path

import pandas as pd


def filter_by_depth(tie_break_report: Path, min_depth: float) -> list[str]:
    tie_break = pd.read_csv(tie_break_report)
    filtered = tie_break[
        (tie_break["mean_depth"] >= min_depth) & (tie_break["depth_type"] == "unique")
    ]["reference"].to_list()
    return filtered

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
        "-m",
        "--min_depth",
        type=float,
        required=True,
        help="Minimum mean depth to retain a reference.",
    )
    parser.add_argument("-o", "--output", type=str, help="Output file", default="high_depth_refs.txt")

    args = parser.parse_args()

    high_depth_refs = filter_by_depth(
        args.tie_break_report,
        args.min_depth,
    )

    write_list(high_depth_refs, Path(args.output))


if __name__ == "__main__":
    cli_entry_point()
