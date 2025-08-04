import argparse

import pandas as pd


def summarise_tie_break(aln_file, output_file):
    # Load the stats file
    alns = pd.read_csv(aln_file)
    alns.query("depth_type == 'final'", inplace=True)

    alns["cov_exp_fold_diff"] = alns["coverage"] / alns["simple_expected_coverage"]

    alns["spurious"] = (alns["coverage"] < 50) & (alns["cov_exp_fold_diff"] < 0.5)

    alns.to_csv(output_file, index=False)


def cli_entry_point():
    parser = argparse.ArgumentParser(description="Summarise tie break results.")
    parser.add_argument(
        "-a", "--aln", type=str, required=True, help="Path to the aln CSV file"
    )
    # parser.add_argument(
    #     "--stats", type=str, required=True, help="Path to the stats yaml file"
    # )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="tie_break_summary.csv",
        help="Output CSV file",
    )
    args = parser.parse_args()
    summarise_tie_break(args.aln, args.output)


if __name__ == "__main__":
    cli_entry_point()
