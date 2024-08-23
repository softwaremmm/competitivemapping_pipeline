# pylint: disable=W1203

"""Check if the number of reads for a genome is above a threshold."""

import argparse
import json
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def check_threshold(json_file_path: str, threshold: int, genome_name: str = "M.tuberculosis") -> str:
    """Check if the number of reads is above a certain threshold.

    Args:
        json_file_path (str): Path to the JSON file.
        threshold (int): Threshold value.
        genome_name (str): Name of the genome. Default is "M.tuberculosis".

    Returns:
        str: "true" if the number of reads is above the threshold, "false" otherwise.
    """

    # Read and parse the JSON file
    with open(json_file_path, "r", encoding="utf-8") as file:
        data = json.load(file)

    # Extract the numreads value for the genome_name
    num_reads = None
    for reference in data.get("references", []):
        if reference.get("genome_name") == genome_name:
            num_reads = int(reference.get("numreads", None))
            break

    if not isinstance(num_reads, int):
        logger.info(f"Number of reads for {genome_name} not found in JSON file")
        return "false"

    logger.info(f"Number of reads for {genome_name}: {num_reads}")

    if num_reads > threshold:
        logger.info(f"Number of reads above the threshold ({num_reads} > {threshold})")
        return "true"

    return "false"


def cli_entry_point():
    """Main function for the script"""
    parser = argparse.ArgumentParser(description="Check if number of reads are above threshold")
    parser.add_argument(
        "--json_file_path",
        help="JSON file output by Competitive Mapping",
        required=True,
    )
    parser.add_argument("--read_threshold", help="Threshold number of reads", required=True)
    args = parser.parse_args()

    logger.info(f"Checking if number of reads are above threshold of {args.read_threshold} in {args.json_file_path}")

    # Using print statements as Nextflow expects output on stdout, without a newline
    # Nextflow expects the output to be in lowercase
    print(
        check_threshold(
            args.json_file_path,
            int(args.read_threshold),
        ),
        end="",
    )


if __name__ == "__main__":
    cli_entry_point()
