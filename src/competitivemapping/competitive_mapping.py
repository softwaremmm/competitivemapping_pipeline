"""Run Competitive Mapping against manifest and aggregate the results"""

import argparse
import dataclasses
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


@dataclasses.dataclass
class Config:
    """
    Config class for the module

    Args:
        seq_platform (str): Sequencing platform
        cpus (int): Number of cores to use
        output_root (str): Path to the output root
    """

    seq_platform: str
    cpus: int
    output_root: str


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
    aln_bam: str, contigs_df: pd.DataFrame, config: Config
) -> pd.DataFrame:
    """Get coverage stats using samtools coverage

    Args:
        aln_bam (str): path to the alignment bam file
        contigs_df (pd.DataFrame): dataframe of contigs
        config (Config): Config object

    Returns:
        pd.DataFrame: summary of coverage metrics for each reference
    """
    logging.info("Getting coverage stats")

    primary_coverage = f"{config.output_root}coverage_primary.tsv"
    full_coverage = f"{config.output_root}coverage_full.tsv"

    args = [
        (aln_bam, "", primary_coverage),
        (
            aln_bam,
            "--excl-flags 1540",
            full_coverage,
        ),  # 1540 just excludes unmapped, and failed
    ]

    with ProcessPoolExecutor(max_workers=config.cpus) as executor:
        _ = list(executor.map(run_samtools_coverage, args))

    coverage_df = process_coverage(primary_coverage, full_coverage, contigs_df)
    return coverage_df


def output_fastqs(
    aln_bam: str,
    references: list[str],
    contigs_df: pd.DataFrame,
    include_unmapped: bool,
    config: Config,
):
    """Extract reads from BAM file and output to FASTQ

    Args:
        aln_bam (str): path to the alignment bam file
        references (list[str]): list of desired references to extract reads for
        contigs_df (pd.DataFrame): dataframe of contigs
        include_unmapped (bool): whether to include unmapped reads
        config (Config): Config object
    """
    logging.info("Outputting FASTQs")
    rnames = contigs_df[contigs_df["reference"].isin(references)]["rname"].tolist()
    if include_unmapped:
        rnames.append('"*"')
    logging.info(f"Extracting reads for reference {references} using rnames: {rnames}")

    sorted_ref_bam = f"{config.output_root}output_aln.bam"

    subprocess.run(
        f"samtools index {aln_bam}",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    # Note: At some point can change to much simpler, but will slightly alter ordering of bam file
    rnames_string = " ".join(rnames)
    command = f"samtools view -h {aln_bam} -u {rnames_string} | samtools sort -n -@ {config.cpus} -o {sorted_ref_bam}"
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
    command = f"samtools fastq --excl-flags 0x100 -@ {config.cpus} "
    if config.seq_platform == "ont":
        command += f"-0 {config.output_root}reads.fastq.gz {sorted_ref_bam}"
    else:
        command += (
            f"-1 {config.output_root}reads_1.fastq.gz -2 {config.output_root}reads_2.fastq.gz"
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


def competitive_mapping_analysis(
    aln_bam: str,
    contigs: pd.DataFrame,
    ref_for_fastq: str | None,
    config: Config,
) -> tuple[pd.DataFrame, dict]:
    """Map reads to a manifest and produce a comparison of coverage and alignment stats

    Args:
        aln_bam (str): Path to the alignment bam file
        contigs (pd.DataFrame): Dataframe of contigs
        ref_for_fastq (str | None): Reference to extract reads for
        config (Config): Config object

    Returns:
        tuple[pd.DataFrame, dict]: Dataframe of species comparison and JSON report
    """
    logging.info("Running competitive mapping analysis")
    coverage_df = get_coverage_stats(aln_bam, contigs, config)

    overall_stats, aln_stats = get_aln_stats(aln_bam, contigs)
    df = pd.merge(coverage_df, aln_stats, on="genome_name")

    logging.info("Writing output files")

    # Can update numreads to actually reflect reads
    # As samtools coverage actually counts alignments
    df["numreads"] = df["primary_reads"] + df["supplementary_reads"]

    final_columns = [c for c in FINAL_COLUMNS if c in df.columns]
    df = df[final_columns].copy()

    # add final row with unmapped read count
    new_row: dict[str, typing.Any] = {col: 0 for col in final_columns}
    new_row["genome_name"] = "unmapped"
    for col in "numreads", "primary_reads", "total_reads":
        new_row[col] = overall_stats["unmapped_reads"]
    df.loc[-1] = pd.Series(new_row)

    # incorporate species information if available
    if "species" in contigs.columns:
        species_lookup = contigs.set_index("reference")["species"].to_dict()
        species_lookup["unmapped"] = "unmapped"
        df["species"] = df["genome_name"].map(species_lookup)
        df = df[["species"] + final_columns]

    df.to_csv(f"{config.output_root}species_comparison.csv", index=False)
    # remove the last row with unmapped read count

    output_json = {
        "total_read_counts": overall_stats,
        "references": (df.iloc[:-1]).to_dict(orient="records"),
        "definitions": COLUMN_EXPLANATIONS,
    }
    with open(
        f"{config.output_root}species_comparison.json", "w", encoding="utf-8"
    ) as file:
        json.dump(output_json, file, indent=4)

    if ref_for_fastq:
        refs = [ref.strip() for ref in ref_for_fastq.split(",")]
        output_fastqs(
            aln_bam,
            refs,
            contigs,
            include_unmapped=True,
            config=config,
        )

    logging.info("Finished competitive mapping")
    return df, output_json


def cli_entry_point():
    """Entry point for the CLI"""
    parser = argparse.ArgumentParser(description="Process sylph report and genomes.")

    parser.add_argument(
        "-i", "--input_bam", required=True, help="Path to the manifest file"
    )
    parser.add_argument("--contigs", required=True, help="Path to the contigs file")
    parser.add_argument(
        "--seq_platform", help="Sequencing platform", default="illumina"
    )
    parser.add_argument("--cpus", help="Number of CPUs to use", default=4, type=int)
    parser.add_argument(
        "--ref_for_fastq", help="Reference to extract reads for", default=None
    )
    parser.add_argument(
        "-o", "--output_root", required=True, help="Path to the output files"
    )

    args = parser.parse_args()

    contigs = pd.read_csv(args.contigs)

    if contigs.empty:
        logging.warning("No contigs found in input file. Producing empty results.")
        produce_empty_outputs(args.output_root)
        return

    config = Config(
        args.seq_platform,
        args.cpus,
        args.output_root,
    )

    competitive_mapping_analysis(
        args.input_bam,
        contigs,
        args.ref_for_fastq,
        config,
    )


if __name__ == "__main__":
    cli_entry_point()
