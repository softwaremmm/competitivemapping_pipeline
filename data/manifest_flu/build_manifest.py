"""Given a multi-fasta of all RefSeq influenza A sequences, build the manifest."""
import argparse
import re
import pandas as pd

class Fasta:
    """Hold a fasta sequence and its header."""
    def __init__(self, header: str, sequence: str):
        self.header = header
        self.rname = header.split("|")[0].strip()
        self.sequence = sequence

def parse_multi_fasta(filename: str) -> list[Fasta]:
    """Parse a multi-fasta file and return a list of Fasta objects.

    Args:
        filename (str): Path to the multi-fasta file.

    Returns:
        list[Fasta]: List of Fasta objects containing headers and sequences.
    """
    fastas = []
    with open(filename) as f:
        current_header = ""
        current_sequence = ""
        for line in f:
            line = line.strip()
            if line:
                if line[0] == ">":
                    if current_header:
                        fastas.append(Fasta(current_header, current_sequence))
                        current_header = line[1:]
                        current_sequence = ""
                    else:
                        current_header = line[1:]
                else:
                    current_sequence += line
        # Append the last sequence if it exists
        if current_header:
            fastas.append(Fasta(current_header, current_sequence))
    return fastas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_in", type=str, required=True, help="Input multi-fasta file containing all RefSeq influenza A sequences.")
    parser.add_argument("--manifest_out", type=str, required=True, help="Output manifest filepath.")
    args = parser.parse_args()

    fasta_in = parse_multi_fasta(args.fasta_in)
    subtype_regex = re.compile(r"H(.)N(.)")
    manifest = []
    for fasta in fasta_in:
        subtype_match = subtype_regex.search(fasta.header)
        H, N = subtype_match.groups()
        segment = None

        if "segment 4" in fasta.header or "hemagglutinin" in fasta.header.lower():
            # We have an H segment
            segment = "H" + H
        elif "segment 6" in fasta.header or "neuraminidase" in fasta.header.lower():
            # We have an N segment
            segment = "N" + N
        else:
            segment = fasta.header.split("|")[1].strip() 
        
        manifest.append([fasta.rname, len(fasta.sequence), segment])

    manifest = pd.DataFrame(sorted(manifest, key=lambda x: x[2]), columns=["rname", "totallength", "reference"])
    manifest.to_csv(args.manifest_out, index=False)


