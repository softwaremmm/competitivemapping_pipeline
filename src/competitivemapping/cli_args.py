"""CLI"""
import argparse
from pathlib import Path


class Arguments:  # pylint: disable=too-few-public-methods
    """Class for holding command line arguments."""

    def __init__(self, argv: list):
        """Initialise a command line argument object.

        Args:
            argv (list): A list of command line arguments, usually `sys.argv[1:]`.
        """
        parser = argparse.ArgumentParser(description="Process competitive mapping output to create a JSON file")
        parser.add_argument(
            "--coverage",
            dest="coverage",
            help="Path to cov_<samplename>.tsv file created using the samtools coverage command",
        )
        parser.add_argument(
            "--species_list",
            dest="species_list",
            help="Path to species_list_<isodate>.csv file (reference data)",
        )
        parser.add_argument(
            "--output",
            default="output.json",
            dest="output",
            help="Path for output .json file",
        )

        args = parser.parse_args(argv)

        self.coverage = Path(args.coverage)
        self.species_list = Path(args.species_list)
        self.output = Path(args.output)
