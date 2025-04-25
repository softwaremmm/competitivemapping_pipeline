"""Build a multi-fasta manifest based on the metadata file."""

import argparse
import pandas as pd
import requests

from requests.adapters import HTTPAdapter, Retry

session = requests.Session()

retries = Retry(total=5, backoff_factor=0.1, status_forcelist=[500, 502, 503, 504])

session.mount("http://", HTTPAdapter(max_retries=retries))


def main(manifest_metadata: str):
    """
    Build a multi-fasta manifest based on the metadata file.

    Parameters:
        manifest_metadata (str): Path to the manifest metadata file.
    """
    # Read the manifest metadata file
    df = pd.read_csv(manifest_metadata)
    # Get the references
    sequences = df["reference"].tolist()

    # Download each assembly
    multifasta = ""
    for sequence in sequences:
        url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={sequence}&db=nuccore&report=fasta"
        multifasta = multifasta + session.get(url).text

    # Save the multifasta file
    with open("multifasta.fasta", "w", encoding="UTF-8") as f:
        f.write(multifasta)


if __name__ == "__main__":
    # use argparse to get the manifest file and output directory
    parser = argparse.ArgumentParser(description="Download assemblies from NCBI.")
    parser.add_argument(
        "manifest_metadata", type=str, help="Path to the manifest metadata file."
    )
    args = parser.parse_args()

    main(args.manifest_metadata)
