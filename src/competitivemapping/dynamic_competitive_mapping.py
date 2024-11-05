import argparse
import gzip
import os
import subprocess

import pandas as pd
from Bio import SeqIO

from competitivemapping.process_aln_stats import get_alignment_stats
from competitivemapping.process_mapping import process_mapping


def make_manifest(report_path: str, genomes_path: str) -> tuple[str, pd.DataFrame]:
    manifest_file = "manifest.fasta.gz"
    df = pd.read_csv(report_path, sep="\t")

    genomes_files = df["Genome_file"].unique().tolist()

    # Cat all genomes into a single file
    with open(manifest_file, "wb") as outfile:
        for f in genomes_files:
            with open(os.path.join(genomes_path, f), "rb") as infile:
                outfile.write(infile.read())

    # Get IDs of genomes
    contigs = []
    for genome in genomes_files:
        with gzip.open(os.path.join(genomes_path, genome), "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                contigs.append(
                    {
                        "reference": genome,
                        "rname": record.id,
                        "length": len(record.seq),
                    }
                )
    contigs_df = pd.DataFrame(contigs)
    # calculate total length per reference
    contigs_df["totallength"] = contigs_df.groupby("reference")["length"].transform("sum")
    print(contigs_df)

    return manifest_file, contigs_df


def map_reads(manifest, reads, seq_platform: str, cpus: int) -> str:
    # use subprocess to run minimap2
    if seq_platform == "ont":
        command = (
            f"minimap2 -ax map-ont -t {cpus} --secondary yes -N 1000 {manifest} {' '.join(reads)} > alignment_raw.sam"
        )
    else:
        command = f"minimap2 -ax sr -t {cpus} --secondary yes -N 1000 {manifest} {' '.join(reads)} > alignment_raw.sam"
    subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE)

    # sort with samtools
    subprocess.run(
        "samtools sort -@ {cpus} -o alignment.bam alignment_raw.sam",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    return "alignment.bam"


def get_aln_stats(aln_bam: str, contigs_df: pd.DataFrame) -> pd.DataFrame:
    """Get stats from bam files using pysam"""
    name_mapping = contigs_df.set_index("rname")["reference"].to_dict()
    mapping_stats, aln_stats_df = get_alignment_stats(aln_bam, name_mapping)
    print(mapping_stats)
    aln_stats_df.to_csv("aln_stats.csv", index=False)

    return aln_stats_df


def get_coverage_stats(aln_bam: str, contigs_df: pd.DataFrame) -> pd.DataFrame:
    """Get stats using samtools coverage"""

    # samtools coverage
    subprocess.run(
        f"samtools coverage {aln_bam} -o primary_coverage.tsv",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )
    subprocess.run(
        f"samtools coverage {aln_bam} --excl-flags 1540 -o full_coverage.tsv",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    coverage_df = process_mapping("primary_coverage.tsv", "full_coverage.tsv", contigs_df)

    return coverage_df


def output_fastqs():
    pass


def run_competitive_mapping(
    manifest: str,
    contigs: pd.DataFrame,
    reads: list[str],
    seq_platform: str,
    cpus: int,
    output: str,
):
    aln_bam = map_reads(manifest, reads, seq_platform, cpus)

    coverage_df = get_coverage_stats(aln_bam, contigs)
    print(coverage_df)

    aln_stats = get_aln_stats(aln_bam, contigs)
    print(aln_stats)

    df = pd.merge(coverage_df, aln_stats, on="genome_name")
    df.to_csv(output, index=False)


def run_dynamic_competitive_mapping(
    sylph_report: str,
    genomes: str,
    reads: list[str],
    seq_platform: str,
    cpus: int,
    output: str,
):
    manifest, contigs = make_manifest(sylph_report, genomes)
    run_competitive_mapping(manifest, contigs, reads, seq_platform, cpus, output)


def cli_entry_point():
    parser = argparse.ArgumentParser(description="Process sylph report and genomes.")

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("--reads", required=True, help="Path to the reads", nargs="+")
    common_parser.add_argument("--seq_platform", help="Sequencing platform", default="illumina")
    common_parser.add_argument("--cpus", help="Number of CPUs to use", default=4)
    common_parser.add_argument("--output", required=True, help="Path to the output file")

    # add subcommand "with_manifest" to run_dynamic_competitive_mapping
    subparsers = parser.add_subparsers(dest="subcommand", required=True)
    manifest_parser = subparsers.add_parser(
        "manifest",
        parents=[common_parser],
        help="Run competitive mapping with set manifest and contigs list",
    )
    manifest_parser.add_argument("--manifest", required=True, help="Path to the manifest file")
    manifest_parser.add_argument("--contigs", required=True, help="Path to the contigs file")

    sylph_parser = subparsers.add_parser(
        "sylph",
        parents=[common_parser],
        help="Run competitive mapping dynamically based on sylph report",
    )
    sylph_parser.add_argument("--sylph_report", required=True, help="Path to the sylph report TSV file")
    sylph_parser.add_argument("--genomes", required=True, help="Path to the genomes")

    args = parser.parse_args()

    if args.subcommand == "manifest":
        run_competitive_mapping(
            args.manifest,
            pd.read_csv(args.contigs),
            args.reads,
            args.seq_platform,
            args.cpus,
            args.output,
        )
    elif args.subcommand == "sylph":
        run_dynamic_competitive_mapping(
            args.sylph_report,
            args.genomes,
            args.reads,
            args.seq_platform,
            args.cpus,
            args.output,
        )


if __name__ == "__main__":
    cli_entry_point()
