"""Script to aid in downloading genomes from NCBI.
Genomes have accessions like GCF_000005845.2
You need to get a table from ncbi including the ftp_path"""

import argparse
import os
import subprocess

import pandas as pd


def download_ncbi(df: pd.DataFrame, outdir):
    """Download genomes from NCBI using ftp_path"""

    def ftp_path_to_url(ftp_path: str) -> str:
        assembly_code = ftp_path.split("/")[-1]
        url = ftp_path + f"/{assembly_code}_genomic.fna.gz"
        return url

    df["download_url"] = df["ftp_path"].apply(ftp_path_to_url)

    os.makedirs(outdir, exist_ok=True)
    for url, accession in zip(df["download_url"], df["assembly_accession"]):
        print(accession, url)
        subprocess.run(
            ["wget", url, "-O", f"{outdir}/{accession}_genomic.fna.gz"],
            check=True,
        )


def main():
    """CLI entry point"""
    parser = argparse.ArgumentParser(
        description="""Download genome fastas from NCBI.
            Genomes have accessions like GCF_000005845.2.
            You need to get a table from ncbi including the ftp_path."""
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="Path to input csv/tsv. Requires column 'assembly_accession' and 'ftp_path' if refseq table not provided.",
    )
    parser.add_argument(
        "-r",
        "--refseq_table",
        type=str,
        help="""Get ftp_path from assembly_summary_refseq.txt.
            Available from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt""",
    )
    parser.add_argument("-o", "--outdir", type=str, help="Path to output directory")

    args = parser.parse_args()

    input_df = pd.read_csv(args.input, sep="\t" if args.input.endswith("tsv") else ",")
    if args.refseq_table:
        refseq_df = pd.read_csv(
            args.refseq_table,
            sep="\t",
            comment="#",
            usecols=["assembly_accession", "ftp_path"],
        )
        ftp_lookup = refseq_df.set_index("assembly_accession")["ftp_path"].to_dict()
        input_df["ftp_path"] = input_df["assembly_accession"].map(ftp_lookup)

        # report on missing ftp_paths and continue without them
        missing_df = input_df[input_df["ftp_path"].isna()]
        if not missing_df.empty:
            print(f"Missing ftp_path for {missing_df['assembly_accession'].tolist()}")
        input_df = input_df.dropna(subset=["ftp_path"])

    print(input_df)
    download_ncbi(input_df, args.outdir)


if __name__ == "__main__":
    main()
