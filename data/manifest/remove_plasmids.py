"""Script to remove all plasmid contigs from multifasta"""

import argparse
import shutil

import pyfastx


def main():
    """Remove plasmid contigs from a multifasta file."""
    parser = argparse.ArgumentParser(description="Remove plasmids from a FASTA file.")
    parser.add_argument(
        "input_fasta",
        help="Input FASTA file(s) containing plasmids potentially.",
        nargs="+",
    )

    args = parser.parse_args()

    for input_fasta in args.input_fasta:
        output_fasta = f"{input_fasta}.filtered"

        with open(output_fasta, "wt", encoding="utf-8") as output_file:
            for name, seq in pyfastx.Fasta(
                input_fasta, build_index=False, full_name=True
            ):
                if "plasmid" in name or "Plasmid" in name:
                    print(f"Skipping {name} because it contains 'plasmid' or 'Plasmid'")
                    continue
                output_file.write(f">{name}\n")
                output_file.write(f"{seq}\n")

        shutil.move(output_fasta, input_fasta)


if __name__ == "__main__":
    main()
