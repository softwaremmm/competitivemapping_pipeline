# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments
# pylint: disable=too-many-locals
"""Run Competitive Mapping and aggregate the results"""

import argparse
import gzip
import json
import logging
import os
import subprocess
import typing

import pandas as pd
from Bio import SeqIO

from competitivemapping.process_aln_stats import get_alignment_stats
from competitivemapping.process_coverage import process_coverage

FINAL_COLUMNS = [
    "genome_name",
    "length",
    "numreads",
    "coverage",
    "meandepth",
    "coverage_including_secondary",
    "meandepth_including_secondary",
    "total_reads",
    "total_alns",
    "exclusive_reads",
    "primary_reads",
    "secondary_reads",
    "supplementary_reads",
    "supplementary_alns",
]

COLUMN_EXPLANATIONS = {
    "length": "Length of all contigs in reference",
    "numreads": "Number of reads which map well to this reference (primary + supplementary)",
    "coverage": "Percentage of reference covered by reads (not including secondary alignments)",
    "meandepth": "Mean depth of coverage of reference (not including secondary alignments)",
    "coverage_including_secondary": "Percentage of reference covered by reads (including secondary alignments)",
    "meandepth_including_secondary": "Mean depth of coverage of reference (including secondary alignments)",
    "total_reads": "Total number of reads which map (in any way) to chrom",
    "total_alns": "Total number of alignments, so will count supplementary separately",
    "exclusive_reads": "reads which only maps to this chrom",
    "primary_reads": "reads which map best to this chrom",
    "secondary_reads": "reads which map to this chrom but not with primary",
    "supplementary_reads": "reads which have supplementary alignments but not primary",
    "supplementary_alns": "number of supplementary alignments",
}


def make_manifest(
    report_path: str, genomes_path: str, output_root: str
) -> tuple[str, pd.DataFrame]:
    """Produce a multifasta manifest and contig df from a sylph report and genomes folder"""
    manifest_file = f"{output_root}manifest.fasta.gz"
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
    contigs_df["totallength"] = contigs_df.groupby("reference")["length"].transform(
        "sum"
    )

    return manifest_file, contigs_df


def map_reads(manifest, reads, seq_platform: str, cpus: int, output_root: str) -> str:
    """Run minimap2 to map reads to a manifest, and sort to bam file"""
    aln_bam = f"{output_root}alignment.bam"

    command = f"minimap2 -t {cpus} --secondary yes -N 1000"
    if seq_platform == "ont":
        command += " -ax map-ont"
    else:
        command += " -ax sr"
    command += f" {manifest} {' '.join(reads)}"

    command += f"| samtools sort -@ {cpus} -o {aln_bam}"

    subprocess.run(
        command,
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    return aln_bam


def get_aln_stats(aln_bam: str, contigs_df: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    """Get stats from bam files using pysam"""
    name_mapping = contigs_df.set_index("rname")["reference"].to_dict()
    overall_stats, aln_stats_df = get_alignment_stats(aln_bam, name_mapping)

    return overall_stats, aln_stats_df


def get_coverage_stats(
    aln_bam: str, contigs_df: pd.DataFrame, output_root: str
) -> pd.DataFrame:
    """Get coverage stats using samtools coverage"""

    primary_coverage = f"{output_root}coverage_primary.tsv"
    full_coverage = f"{output_root}coverage_full.tsv"

    # samtools coverage
    subprocess.run(
        f"samtools coverage {aln_bam} -o {primary_coverage}",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )
    subprocess.run(
        f"samtools coverage {aln_bam} --excl-flags 1540 -o {full_coverage}",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    coverage_df = process_coverage(primary_coverage, full_coverage, contigs_df)

    return coverage_df


def output_fastqs(
    aln_bam: str,
    reference: str,
    contigs_df: pd.DataFrame,
    include_unmapped: bool,
    seq_platform: str,
    cpus: int,
    output_root: str,
):
    """Extract reads from BAM file and output to FASTQ"""
    rnames = contigs_df[contigs_df["reference"] == reference]["rname"].tolist()
    if include_unmapped:
        rnames.append('"*"')

    sorted_ref_bam = f"{output_root}output_aln.bam"

    subprocess.run(
        f"samtools index {aln_bam}",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )
    # In future could add flag -P to always include read pairs
    # But note that this fails to fetch the pair for supplementary alignments
    subprocess.run(
        f"samtools view {aln_bam} -u {' '.join(rnames)} | samtools sort -@ {cpus} -o {sorted_ref_bam}",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    # Convert BAM output to FASTQ
    # Default excl-flag is 0x900 which is secondary (0x100) and supplementary (0x800)
    # So we need to exclude secondary alignments (0x100) only
    # Note that singletons will be discarded;
    # this is where only one of the read pair is in the bam or passes the flags

    # So for paired reads, only output if both are unmapped, or both map
    # (maybe with supplementary) to the reference
    command = f"samtools fastq --excl-flags 0x100 -@ {cpus} "
    if seq_platform == "ont":
        command += f"-0 {output_root}reads.fastq.gz {sorted_ref_bam}"
    else:
        command += (
            f"-1 {output_root}reads_1.fastq.gz -2 {output_root}reads_2.fastq.gz"
            f" -0 /dev/null -s /dev/null {sorted_ref_bam}"
        )
    subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE)


def produce_empty_outputs(output_root: str):
    """Create empty output files"""
    report = {
        "total_read_counts": {
            "mapped_reads": 0,
            "unmapped_reads": 0,
        },
        "references": [],
        "definitions": COLUMN_EXPLANATIONS,
    }
    with open(f"{output_root}species_comparison.json", "w", encoding="utf-8") as file:
        json.dump(report, file, indent=4)

    with open(f"{output_root}species_comparison.csv", "w", encoding="utf-8") as file:
        file.write(",".join(FINAL_COLUMNS) + "\n")


def run_competitive_mapping(
    manifest: str,
    contigs: pd.DataFrame,
    reads: list[str],
    ref_for_fastq: str | None,
    seq_platform: str,
    cpus: int,
    output_root: str,
):
    """Map reads to a manifest and produce a comparison of coverage and alignment stats"""
    aln_bam = map_reads(manifest, reads, seq_platform, cpus, output_root)

    coverage_df = get_coverage_stats(aln_bam, contigs, output_root)

    overall_stats, aln_stats = get_aln_stats(aln_bam, contigs)
    df = pd.merge(coverage_df, aln_stats, on="genome_name")

    # Can update numreads to actually reflect reads
    # As samtools coverage actually counts alignments
    df["numreads"] = df["primary_reads"] + df["supplementary_reads"]

    for col in (
        "coverage",
        "meandepth",
        "coverage_including_secondary",
        "meandepth_including_secondary",
    ):
        if col in df.columns:
            df[col] = df[col].apply(lambda x: float(f"{x:.4g}"))

    df = df[FINAL_COLUMNS].copy()

    # add final row with unmapped read count
    new_row: dict[str, typing.Any] = {col: 0 for col in FINAL_COLUMNS}
    new_row["genome_name"] = "unmapped"
    for col in "numreads", "primary_reads", "total_reads":
        new_row[col] = overall_stats["unmapped_reads"]
    df.loc[-1] = new_row

    # incorporate species information if available
    if "species" in contigs.columns:
        species_lookup = contigs.set_index("reference")["species"].to_dict()
        species_lookup["unmapped"] = "unmapped"
        df["species"] = df["genome_name"].map(species_lookup)
        df = df[["species"] + FINAL_COLUMNS]

    df.to_csv(f"{output_root}species_comparison.csv", index=False)
    # remove the last row with unmapped read count
    df = df.iloc[:-1]

    output = {
        "total_read_counts": overall_stats,
        "references": df.to_dict(orient="records"),
        "definitions": COLUMN_EXPLANATIONS,
    }
    with open(f"{output_root}species_comparison.json", "w", encoding="utf-8") as file:
        json.dump(output, file, indent=4)

    if ref_for_fastq:
        output_fastqs(
            aln_bam,
            ref_for_fastq,
            contigs,
            include_unmapped=True,
            seq_platform=seq_platform,
            cpus=cpus,
            output_root=output_root,
        )


def run_dynamic_competitive_mapping(
    sylph_report: str,
    genomes: str,
    reads: list[str],
    ref_for_fastq: str | None,
    db_metadata: str | None,
    seq_platform: str,
    cpus: int,
    output_root: str,
):
    """Create manifest then competitive mapping"""

    # check if sylph_report is empty
    if os.stat(sylph_report).st_size == 0 or pd.read_csv(sylph_report).empty:
        logging.warning("Sylph report is empty")
        produce_empty_outputs(output_root)
        return

    manifest, contigs = make_manifest(sylph_report, genomes, output_root)
    if db_metadata:
        metadata = pd.read_csv(
            db_metadata, sep="\t", header=None, names=["assembly", "taxonomy"]
        )
        metadata["species"] = (
            metadata["taxonomy"].str.split(";").str[-1].str.replace("s__", "")
        )
        species_lookup = metadata.set_index("assembly")["species"].to_dict()

        contigs["reference"] = (
            contigs["reference"].str.split("/").str[-1].str.split("_genomic").str[0]
        )
        contigs["species"] = contigs["reference"].map(species_lookup)

    run_competitive_mapping(
        manifest, contigs, reads, ref_for_fastq, seq_platform, cpus, output_root
    )


def cli_entry_point():
    """Entry point for the CLI"""
    parser = argparse.ArgumentParser(description="Process sylph report and genomes.")

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "--reads", required=True, help="Path to the reads", nargs="+"
    )
    common_parser.add_argument(
        "--seq_platform", help="Sequencing platform", default="illumina"
    )
    common_parser.add_argument("--cpus", help="Number of CPUs to use", default=4)
    common_parser.add_argument(
        "--ref_for_fastq", help="Reference to extract reads for", default=None
    )
    common_parser.add_argument(
        "--output_root", required=True, help="Path to the output files"
    )

    # add subcommand "with_manifest" to run_dynamic_competitive_mapping
    subparsers = parser.add_subparsers(dest="subcommand", required=True)
    manifest_parser = subparsers.add_parser(
        "manifest",
        parents=[common_parser],
        help="Run competitive mapping with set manifest and contigs list",
    )
    manifest_parser.add_argument(
        "--manifest", required=True, help="Path to the manifest file"
    )
    manifest_parser.add_argument(
        "--contigs", required=True, help="Path to the contigs file"
    )

    sylph_parser = subparsers.add_parser(
        "sylph",
        parents=[common_parser],
        help="Run competitive mapping dynamically based on sylph report",
    )
    sylph_parser.add_argument(
        "--sylph_report", required=True, help="Path to the sylph report TSV file"
    )
    sylph_parser.add_argument("--genomes", required=True, help="Path to the genomes")
    sylph_parser.add_argument(
        "--db_metadata",
        required=False,
        help="Path to the db metadata with assembly to species mapping",
    )

    args = parser.parse_args()

    if args.subcommand == "manifest":
        run_competitive_mapping(
            args.manifest,
            pd.read_csv(args.contigs),
            args.reads,
            args.ref_for_fastq,
            args.seq_platform,
            args.cpus,
            args.output_root,
        )
    elif args.subcommand == "sylph":
        run_dynamic_competitive_mapping(
            args.sylph_report,
            args.genomes,
            args.reads,
            args.ref_for_fastq,
            args.db_metadata,
            args.seq_platform,
            args.cpus,
            args.output_root,
        )


if __name__ == "__main__":
    cli_entry_point()
