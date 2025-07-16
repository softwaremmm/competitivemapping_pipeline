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
    parser.add_argument("--fasta_out", type=str, required=True, help="Output multi-fasta file for filtered sequences.")
    parser.add_argument("--manifest_out", type=str, required=True, help="Output manifest filepath.")
    parser.add_argument("--ani_duplicates", type=str, required=False, default=None, help="Line separated list of fasta rnames to skip due to being ANI duplicates. If not provided, no sequences will be skipped.")
    args = parser.parse_args()

    if args.ani_duplicates:
        with open(args.ani_duplicates, "r") as f:
            ani_duplicates = set(line.strip() for line in f if line.strip())
    else:
        ani_duplicates = set()

    fasta_in = parse_multi_fasta(args.fasta_in)
    subtype_regex = re.compile(r"H(\d+)N(\d+)")


    manifest = []
    ani_seen = []
    seen_segments = set()
    duplicates = set()
    for fasta in fasta_in:
        if fasta.rname in ani_duplicates:
            # Skip sequences that are in the ANI duplicates list
            continue
        subtype_match = subtype_regex.search(fasta.header)
        if subtype_match:
            H, N = subtype_match.groups()
        segment = None

        if "segment 4" in fasta.header or "hemagglutinin" in fasta.header.lower():
            # We have an H segment
            segment = "H" + H + f" ({fasta.rname})"
            if "H" + H not in ani_seen:
                ani_group = len(ani_seen)
                ani_seen.append("H" + H)
            else:
                ani_group = ani_seen.index("H" + H)
        elif "segment 6" in fasta.header or "neuraminidase" in fasta.header.lower():
            # We have an N segment
            segment = "N" + N + f" ({fasta.rname})"
            if "N" + N not in ani_seen:
                ani_group = len(ani_seen)
                ani_seen.append("N" + N)
            else:
                ani_group = ani_seen.index("N" + N)
        else:
            segment = fasta.header.split("|")[1].strip() 
            if segment in ani_seen:
                ani_group = ani_seen.index(segment)
            else:
                ani_group = len(ani_seen)
                ani_seen.append(segment)
        
        if segment in seen_segments:
            # If we have already seen this segment, skip it
            duplicates.add(fasta.rname)
            print("Duplicate segment with rname: ", fasta.rname)
            continue
        else:
            seen_segments.add(segment)
        
        manifest.append([fasta.rname, len(fasta.sequence), segment, ani_group])

    manifest = pd.DataFrame(sorted(manifest, key=lambda x: x[2]), columns=["rname", "totallength", "reference", "ani_group"])
    manifest.to_csv(args.manifest_out, index=False)

    with open(args.fasta_out, "w") as f_out:
        for fasta in fasta_in:
            if fasta.rname not in duplicates and fasta.rname not in ani_duplicates:
                f_out.write(f">{fasta.header}\n{fasta.sequence}\n")
    



