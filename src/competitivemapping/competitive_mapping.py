# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments
# pylint: disable=too-many-locals
# pylint: disable=logging-fstring-interpolation
"""Run Competitive Mapping against manifest and aggregate the results"""

import argparse
import json
import logging
import subprocess
import typing
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

from competitivemapping.process_aln_stats import get_alignment_stats
from competitivemapping.process_coverage import process_coverage

logging.basicConfig(
    format="%(asctime)s — %(relativeCreated)d — %(levelname)s — %(funcName)s:%(lineno)d — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S%z",
    level=logging.DEBUG,
)


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


def map_reads(
    manifest: str, reads: list[str], seq_platform: str, cpus: int, output_root: str
) -> str:
    """Run minimap2 to map reads against manifest, returning sorted bam.

    Args:
        manifest (str): Path to the manifest file
        reads (list[str]): List of paths to the read fastqs
        seq_platform (str): Sequencing platform
        cpus (int): number of cores to use
        output_root (str): Path to the output root

    Returns:
        str: path to the alignment bam file
    """
    logging.info("Mapping reads")
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
    """Get stats from bam files using pysam.

    Args:
        aln_bam (str): path to the alignment bam file
        contigs_df (pd.DataFrame): dataframe of contigs

    Returns:
        tuple[dict, pd.DataFrame]: dict with overall stats and
            dataframe with alignment stats per reference
    """
    logging.info("Getting alignment stats")
    name_mapping = contigs_df.set_index("rname")["reference"].to_dict()
    overall_stats, aln_stats_df = get_alignment_stats(aln_bam, name_mapping)

    return overall_stats, aln_stats_df


def run_samtools_coverage(args: tuple[str, str, str]) -> str:
    """Run samtools coverage on a bam file

    Args:
        args (tuple[str, str, str]): (bam file, flags, output file)

    Returns:
        str: output file
    """
    aln_bam, flags, output_file = args
    subprocess.run(
        f"samtools coverage {aln_bam} {flags} -o {output_file}",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )
    return output_file


def get_coverage_stats(
    aln_bam: str, contigs_df: pd.DataFrame, cpus: int, output_root: str
) -> pd.DataFrame:
    """Get coverage stats using samtools coverage

    Args:
        aln_bam (str): path to the alignment bam file
        contigs_df (pd.DataFrame): dataframe of contigs
        cpus (int): number of cores to use (uses 2 max)
        output_root (str): Path to the output root

    Returns:
        pd.DataFrame: summary of coverage metrics for each reference
    """
    logging.info("Getting coverage stats")

    primary_coverage = f"{output_root}coverage_primary.tsv"
    full_coverage = f"{output_root}coverage_full.tsv"
    args = [
        (aln_bam, "", primary_coverage),
        (aln_bam, "--excl-flags 1540", full_coverage),
    ]

    with ProcessPoolExecutor(max_workers=cpus) as executor:
        _results = list(executor.map(run_samtools_coverage, args))

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
    """Extract reads from BAM file and output to FASTQ

    Args:
        aln_bam (str): path to the alignment bam file
        reference (str): desired reference to extract reads for
        contigs_df (pd.DataFrame): dataframe of contigs
        include_unmapped (bool): whether to include unmapped reads
        seq_platform (str): sequencing platform
        cpus (int): number of cores to use
        output_root (str): Path to the output root
    """
    logging.info("Outputting FASTQs")
    rnames = contigs_df[contigs_df["reference"] == reference]["rname"].tolist()
    if include_unmapped:
        rnames.append('"*"')
    logging.info("Extracting reads for reference {reference} using rnames: {rnames}")

    sorted_ref_bam = f"{output_root}output_aln.bam"

    subprocess.run(
        f"samtools index {aln_bam}",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    for i, rname in enumerate(rnames):
        # In future could add flag -P to always include read pairs
        # But note that this fails to fetch the pair for supplementary alignments
        command = f"samtools view -h {aln_bam} -u {rname} | samtools sort -n -@ {cpus} -o {output_root}.{i}.bam"
        logging.info("Running command: %s", command)
        subprocess.run(
            command,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
        )

    command = f"samtools merge -fo {sorted_ref_bam} {' '.join([f'{output_root}.{i}.bam' for i in range(len(rnames))])}"
    logging.info("Running command: %s", command)
    subprocess.run(
        command,
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
    logging.info("Running command: %s", command)
    subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE)


def produce_empty_outputs(output_root: str):
    """Create empty output files for species comparison json/csv"""
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
    """Map reads to a manifest and produce a comparison of coverage and alignment stats

    Args:
        manifest (str): Path to the manifest file
        contigs (pd.DataFrame): Dataframe of contigs
        reads (list[str]): List of paths to the read fastqs
        ref_for_fastq (str | None): Reference to extract reads for
        seq_platform (str): Sequencing platform
        cpus (int): Number of cores to use
        output_root (str): Path to the output root
    """
    logging.info("Running competitive mapping")
    aln_bam = map_reads(manifest, reads, seq_platform, cpus, output_root)

    coverage_df = get_coverage_stats(aln_bam, contigs, cpus, output_root)

    overall_stats, aln_stats = get_aln_stats(aln_bam, contigs)
    df = pd.merge(coverage_df, aln_stats, on="genome_name")

    logging.info("Writing output files")

    # Can update numreads to actually reflect reads
    # As samtools coverage actually counts alignments
    df["numreads"] = df["primary_reads"] + df["supplementary_reads"]

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

    logging.info("Finished competitive mapping")


def cli_entry_point():
    """Entry point for the CLI"""
    parser = argparse.ArgumentParser(description="Process sylph report and genomes.")

    parser.add_argument("--manifest", required=True, help="Path to the manifest file")
    parser.add_argument("--contigs", required=True, help="Path to the contigs file")
    parser.add_argument("--reads", required=True, help="Path to the reads", nargs="+")
    parser.add_argument(
        "--seq_platform", help="Sequencing platform", default="illumina"
    )
    parser.add_argument("--cpus", help="Number of CPUs to use", default=4, type=int)
    parser.add_argument(
        "--ref_for_fastq", help="Reference to extract reads for", default=None
    )
    parser.add_argument("--output_root", required=True, help="Path to the output files")

    args = parser.parse_args()

    contigs = pd.read_csv(args.contigs)

    if contigs.empty:
        logging.warning("No contigs found in input file. Producing empty results.")
        produce_empty_outputs(args.output_root)
        return

    run_competitive_mapping(
        args.manifest,
        contigs,
        args.reads,
        args.ref_for_fastq,
        args.seq_platform,
        args.cpus,
        args.output_root,
    )


if __name__ == "__main__":
    cli_entry_point()
