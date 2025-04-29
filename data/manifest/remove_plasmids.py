""" Script to remove all plasmid contigs from multifasta """
import argparse

import pyfastx


def main():
    """Remove plasmid contigs from a multifasta file."""
    parser = argparse.ArgumentParser(description="Remove plasmids from a FASTA file.")
    parser.add_argument("input_fasta", help="Input FASTA file containing plasmids.")
    parser.add_argument("output_fasta", help="Output FASTA file without plasmids.")

    args = parser.parse_args()

    with open(args.output_fasta, "wt", encoding="utf-8") as output_file:
        for name, seq in pyfastx.Fasta(
            args.input_fasta, build_index=False, full_name=True
        ):
            if "plasmid" in name or "Plasmid" in name:
                print(f"Skipping {name} because it contains 'plasmid' or 'Plasmid'")
                continue
            output_file.write(f">{name}\n")
            output_file.write(f"{seq}\n")


if __name__ == "__main__":
    main()
