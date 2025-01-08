# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments
# pylint: disable=too-many-locals
"""Run Competitive Mapping against manifest and aggregate the results"""

import argparse
import json
import logging
import subprocess
import typing
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from pysam import AlignmentFile  # pylint: disable = no-name-in-module

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
    manifest: str, reads: list[str], platform: str, cpus: int, output_root: str
) -> str:
    """Run minimap2 to map reads against manifest, returning sorted bam.

    Args:
        manifest (str): Path to the manifest file
        reads (list[str]): List of paths to the read fastqs
        platform (str): Sequencing platform, or "fasta" for assembly contigs
        cpus (int): number of cores to use
        output_root (str): Path to the output root

    Returns:
        str: path to the alignment bam file
    """
    logging.info("Mapping reads")
    aln_bam = f"{output_root}alignment.bam"

    command = f"minimap2 -t {cpus} --secondary yes -N 1000"
    if platform == "ont":
        command += " -ax map-ont"
    elif platform == "illumina":
        command += " -ax sr"
    else:
        command += " -a"
    command += f" {manifest} {' '.join(reads)}"

    command += f"| samtools sort -@ {cpus} -o {aln_bam}"

    subprocess.run(
        command,
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    return aln_bam


def get_aln_stats(
    aln_bam: str, contigs_df: pd.DataFrame, weighting_df: pd.DataFrame | None = None
) -> tuple[dict, pd.DataFrame]:
    """Get stats from bam files using pysam.

    Args:
        aln_bam (str): path to the alignment bam file
        contigs_df (pd.DataFrame): dataframe of contigs
        weighting_df (pd.DataFrame): dataframe of weights for reads (optional)

    Returns:
        tuple[dict, pd.DataFrame]: dict with overall stats and
            dataframe with alignment stats per reference
    """
    logging.info("Getting alignment stats")
    name_mapping = contigs_df.set_index("rname")["reference"].to_dict()
    overall_stats, aln_stats_df = get_alignment_stats(
        aln_bam, name_mapping, weighting_df=weighting_df
    )

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


def estimate_mean_depth(
    bam_file: str, manifest_contigs: pd.DataFrame, query_contig_stats: pd.DataFrame
) -> pd.DataFrame:
    """Use pysam to estimate mean depth of coverage when mapping assembled contigs

    Args:
        bam_file (str): path to the alignment bam file
        manifest_contigs (pd.DataFrame): dataframe of manifest contigs
        query_contig_stats (pd.DataFrame): stats about assembled contigs

    Returns:
        pd.DataFrame: Dataframe with columns genome_name, meandepth, meandepth_including_secondary
    """
    alignments = []
    with AlignmentFile(bam_file, "rb") as bam:  # ignore: no-member
        for read in bam:
            if read.is_qcfail or read.is_duplicate or read.is_unmapped:
                continue

            alignments.append(
                {
                    "contig_name": read.query_name,
                    "rname": read.reference_name,
                    "secondary": read.is_secondary,
                }
            )

    df = pd.DataFrame(alignments)
    df = pd.merge(df, manifest_contigs, on="rname").drop(columns=["rname"])

    # So df now has contig_name, secondary, reference, totallength
    # But only let contig appear once for each reference (favouring primary alignments)
    df.sort_values(
        by=["contig_name", "reference", "secondary"],
        inplace=True,
        ascending=[True, True, True],
    )
    df.drop_duplicates(subset=["contig_name", "reference"], inplace=True)

    df = pd.merge(
        df, query_contig_stats[["contig_name", "length", "meandepth"]], on="contig_name"
    )
    df["depth_contribution"] = df["length"] * df["meandepth"] / df["totallength"]

    all_aggregated = df.groupby("reference", as_index=False).agg(
        meandepth_including_secondary=("depth_contribution", "sum")
    )
    primary_aggregated = (
        df[~df["secondary"]]
        .groupby("reference", as_index=False)
        .agg(meandepth=("depth_contribution", "sum"))
    )
    aggregated = (
        manifest_contigs[["reference"]]
        .drop_duplicates()
        .merge(primary_aggregated, on="reference", how="left")
        .merge(all_aggregated, on="reference", how="left")
        .fillna(0)
        .rename(columns={"reference": "genome_name"})
    )
    return aggregated


def output_fastqs(
    aln_bam: str,
    reference: str,
    contigs_df: pd.DataFrame,
    include_unmapped: bool,
    platform: str,
    cpus: int,
    output_root: str,
):
    """Extract reads from BAM file and output to FASTQ

    Args:
        aln_bam (str): path to the alignment bam file
        reference (str): desired reference to extract reads for
        contigs_df (pd.DataFrame): dataframe of contigs
        include_unmapped (bool): whether to include unmapped reads
        platform (str): sequencing platform
        cpus (int): number of cores to use
        output_root (str): Path to the output root
    """
    logging.info("Outputting FASTQs")
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
    if platform == "illumina":
        command += (
            f"-1 {output_root}reads_1.fastq.gz -2 {output_root}reads_2.fastq.gz"
            f" -0 /dev/null -s /dev/null {sorted_ref_bam}"
        )
    else:
        command += f"-0 {output_root}reads.fastq.gz {sorted_ref_bam}"
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
    manifest_contigs: pd.DataFrame,
    query: list[str],
    query_contig_stats: pd.DataFrame | None,
    ref_to_extract: str | None,
    platform: str,
    cpus: int,
    output_root: str,
):
    """Map reads to a manifest and produce a comparison of coverage and alignment stats

    Args:
        manifest (str): Path to the manifest file
        manifest_contigs (pd.DataFrame): Dataframe of manifest_contigs
        query (list[str]): List of paths to the read fastqs
        query_contig_stats (pd.DataFrame): Dataframe of contig stats if using assembled contigs
        ref_to_extract (str | None): Reference to extract reads for
        platform (str): Sequencing platform
        cpus (int): Number of cores to use
        output_root (str): Path to the output root
    """
    logging.info("Running competitive mapping")
    aln_bam = map_reads(manifest, query, platform, cpus, output_root)

    coverage_df = get_coverage_stats(aln_bam, manifest_contigs, cpus, output_root)

    overall_stats, aln_stats = get_aln_stats(
        aln_bam, manifest_contigs, query_contig_stats
    )
    df = pd.merge(coverage_df, aln_stats, on="genome_name")

    if query_contig_stats is not None:
        # ensure that contig_name is a string
        query_contig_stats["contig_name"] = query_contig_stats["contig_name"].astype(
            str
        )

        # Drop columns that have not been accurately calculated
        df.drop(
            columns=["meandepth", "meandepth_including_secondary", "numreads"],
            inplace=True,
        )

        df = pd.merge(
            df,
            estimate_mean_depth(aln_bam, manifest_contigs, query_contig_stats),
            on="genome_name",
        )

    logging.info("Writing output files")

    # Can update numreads to actually reflect reads
    # As samtools coverage actually counts alignments
    df["numreads"] = df["primary_reads"] + df["supplementary_reads"]

    df.sort_values(
        ["meandepth", "meandepth_including_secondary", "genome_name"],
        ascending=[False, False, True],
        inplace=True,
    )

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
    df.loc[-1] = new_row  # type: ignore

    # incorporate species information if available
    if "species" in manifest_contigs.columns:
        species_lookup = manifest_contigs.set_index("reference")["species"].to_dict()
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

    if ref_to_extract:
        output_fastqs(
            aln_bam,
            ref_to_extract,
            manifest_contigs,
            include_unmapped=True,
            platform=platform,
            cpus=cpus,
            output_root=output_root,
        )

    logging.info("Finished competitive mapping")


def cli_entry_point():
    """Entry point for the CLI"""
    parser = argparse.ArgumentParser(description="Process sylph report and genomes.")

    parser.add_argument("--manifest", required=True, help="Path to the manifest file")
    parser.add_argument(
        "--manifest_contigs",
        required=True,
        help="Path to the contigs file for the manifest",
    )
    parser.add_argument(
        "--platform",
        help="Sequencing platform (illumina/ont/fasta). Use fasta for assembled contigs",
        default="illumina",
    )
    parser.add_argument(
        "--query",
        required=True,
        help="Path to the fastq/fasta files to map against manifest",
        nargs="+",
    )
    parser.add_argument(
        "--query_contig_stats",
        required=False,
        help="If using assembled contigs, provide stats file for numreads and meandepth of each contig",
    )
    parser.add_argument(
        "--ref_to_extract",
        help="Reference to extract matching reads to output fastq",
        default=None,
    )
    parser.add_argument("--cpus", help="Number of CPUs to use", default=4, type=int)
    parser.add_argument("--output_root", required=True, help="Path to the output files")

    args = parser.parse_args()

    manifest_contigs = pd.read_csv(args.manifest_contigs)

    if manifest_contigs.empty:
        logging.warning("No manifest contigs found in input. Producing empty results.")
        produce_empty_outputs(args.output_root)
        return

    query_contig_stats = (
        pd.read_csv(args.query_contig_stats) if args.query_contig_stats else None
    )

    run_competitive_mapping(
        args.manifest,
        manifest_contigs,
        args.query,
        query_contig_stats,
        args.ref_to_extract,
        args.platform,
        args.cpus,
        args.output_root,
    )


if __name__ == "__main__":
    cli_entry_point()
