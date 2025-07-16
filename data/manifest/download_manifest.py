"""
Script for downloading assemblies from NCBI given a manifest file.

Requires ncbi datasets tool to be downloaded: https://www.ncbi.nlm.nih.gov/datasets/docs/v1/download-and-install/
Can be installed with conda: conda install -c conda-forge ncbi-datasets-cli
"""

import argparse
import os
import shutil
import subprocess
import zipfile

import pandas as pd


def download_assembly(accession, output_dir):
    """
    Downloads an assembly from NCBI given its accession ID.

    Parameters:
        accession (str): The NCBI assembly accession.
        output_dir (str): Directory to save the downloaded assembly.
    """
    # check if output file is already downloaded
    if os.path.exists(f"{output_dir}/{accession}.fasta"):
        print(f"{accession} already downloaded.")
        return

    try:
        # Download the assembly
        print(f"Downloading {accession}...")
        subprocess.run(
            [
                "datasets",
                "download",
                "genome",
                "accession",
                accession,
                "--filename",
                f"tmp_{accession}.zip",
                "--include",
                "genome",
            ],
            check=True,
        )

        with zipfile.ZipFile(f"tmp_{accession}.zip", "r") as zip_ref:
            zip_ref.extractall(f"{accession}_dir")

        # Move the assembly to the output directory. Makes use of glob since not actually sure what filename is
        subprocess.run(
            f"mv {accession}_dir/ncbi_dataset/data/{accession}/* {output_dir}/{accession}.fasta",
            check=True,
            shell=True,
        )

        # Clean up
        shutil.rmtree(f"{accession}_dir", ignore_errors=True)
        os.remove(f"tmp_{accession}.zip")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading {accession}: {e}")


def main(manifest_metadata, output_dir):
    """
    Downloads assemblies from NCBI given a manifest file.

    Parameters:
        manifest_metadata (str): Path to the manifest metadata file.
        output_dir (str): Directory to save the downloaded assemblies.
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(manifest_metadata)

    # Remove the excluded assemblies
    if "not_used_in_manifest" in df.columns:
        df = df[df["not_used_in_manifest"] != "y"]

    # Check for duplicate assemblies
    duplicated_assemblies = df[df.duplicated(subset="assembly_accession", keep=False)]
    if not duplicated_assemblies.empty:
        print("Duplicated assemblies:")
        print(duplicated_assemblies)
        return

    for accession in df["assembly_accession"]:
        download_assembly(accession.strip(), output_dir)


if __name__ == "__main__":
    # use argparse to get the manifest file and output directory
    parser = argparse.ArgumentParser(description="Download assemblies from NCBI.")
    parser.add_argument(
        "manifest_metadata", type=str, help="Path to the manifest metadata file."
    )
    parser.add_argument(
        "output_dir", type=str, help="Directory to save the downloaded assemblies."
    )
    args = parser.parse_args()

    main(args.manifest_metadata, args.output_dir)
