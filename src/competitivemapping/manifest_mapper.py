"""Maps reads against a manifest to produce bam file"""
import argparse
import dataclasses
import gzip
import logging
import os
import subprocess

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
        sort_by_name (bool): Whether to sort the output by name
        filter_secondary (bool): Whether to filter secondary alignments
        filter_supplementary (bool): Whether to filter supplementary alignments
        equality_in_cigar (bool): Whether to add equality in cigar string
    """

    seq_platform: str
    cpus: int
    sort_by_name: bool
    filter_secondary: bool
    filter_supplementary: bool
    equality_in_cigar: bool


def map_reads(manifest: str, reads: list[str], config: Config, output_bam: str) -> str:
    """Run minimap2 to map reads against manifest, returning sorted bam.

    Args:
        manifest (str): Path to the manifest file
        reads (list[str]): List of paths to the read fastqs
        config (Config): Config object

    Returns:
        str: path to the alignment bam file
    """
    logging.info("Mapping reads using config: %s", config)

    # -I100G means that more of manifest is loaded into memory
    command = f"minimap2 -I100G -t {config.cpus}"

    if config.seq_platform == "ont":
        command += " -ax map-ont"
    else:
        command += " -ax sr"

    if config.filter_secondary:
        command += " --secondary no"
    else:
        command += " --secondary yes -N 1000"

    # eqx adds more info to the cigar string
    if config.equality_in_cigar:
        command += " --eqx"
    # -p 0.4 lowers the bar for secondary alignments. May want in the future

    command += f" {manifest} {' '.join(reads)}"

    command += f"| samtools sort -@ {config.cpus}"

    if config.sort_by_name:
        command += " -n"

    # Note that if a read is secondary to a supplementary read, it will only be marked as secondary
    if config.filter_supplementary:
        command += " | samtools view -F 2048"

    command += f" -o {output_bam}"

    subprocess.run(
        command,
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    logging.info("Mapping complete")
    return output_bam


def is_manifest_empty(manifest: str) -> bool:
    """Check if the manifest file is empty"""
    if manifest.endswith(".gz"):
        with gzip.open(manifest, "rb") as f:
            # Check if the file is empty
            return f.read(1) == b""

    return os.path.getsize(manifest) == 0


def cli_entry_point():
    """Entry point for the CLI"""
    parser = argparse.ArgumentParser(description="Map reads against manifest.")

    parser.add_argument("--manifest", required=True, help="Path to the manifest file")
    parser.add_argument(
        "-r", "--reads", required=True, help="Path to the reads", nargs="+"
    )
    parser.add_argument(
        "--seq_platform", help="Sequencing platform", default="illumina"
    )
    parser.add_argument(
        "-c", "--cpus", help="Number of CPUs to use", default=4, type=int
    )
    parser.add_argument(
        "-n", "--sort_by_name", help="Sort by name", action="store_true", default=False
    )
    parser.add_argument(
        "--equality_in_cigar",
        help="Add equality in cigar string",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--filter_secondary",
        help="Remove secondary alignments",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--filter_supplementary",
        help="Remove supplementary alignments",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path for output bam file"
    )

    args = parser.parse_args()

    # Check if manifest is an empty file
    if is_manifest_empty(args.manifest):
        logging.error(f"Manifest file {args.manifest} is empty. No outputs produced.")
        return

    config = Config(
        args.seq_platform,
        args.cpus,
        args.sort_by_name,
        args.filter_secondary,
        args.filter_supplementary,
        args.equality_in_cigar,
    )

    map_reads(args.manifest, args.reads, config, args.output)


if __name__ == "__main__":
    cli_entry_point()
