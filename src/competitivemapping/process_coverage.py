# pylint: disable=logging-fstring-interpolation
"""Process output from Competitive Mapping"""

import argparse
import json
import logging
import sys
from importlib.resources import files
from pathlib import Path

import pandas as pd
from jsonschema import validate

logging.basicConfig(
    format="%(asctime)s — %(name)s — %(levelname)s — %(funcName)s:%(lineno)d — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S%z",
    level=logging.DEBUG,
)


def unmatched_rnames(coverage_table: pd.DataFrame, species_table: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Check if any rnames do not have a reference (species name)

    Args:
        coverage_table (pd.DataFrame): Output of `samtools coverage`
        species_table (pd.DataFrame): Table for rnames and references (reference data)

    Returns:
        [pd.Series, pd.Series]: A series listing those rnames in the species output, but
        not in the coverage table, and another listing those in the coverage output but not
        in the species table
    """

    all_rnames = coverage_table.merge(species_table, left_on="#rname", right_on="rname", how="outer")

    not_in_coverage = all_rnames[all_rnames["#rname"].isna()]["rname"]
    not_in_species = all_rnames[all_rnames["rname"].isna()]["#rname"]

    return not_in_coverage, not_in_species


def join_references(coverage_table: pd.DataFrame, species_table: pd.DataFrame) -> pd.DataFrame:
    """Look up references from rnames.

    Args:
        coverage_table (pd.DataFrame): Output of `samtools coverage`
        species_table (pd.DataFrame): Table for rnames and references (reference data)

    Returns:
        pd.DataFrame: As `coverage_table`, but with an extra column for reference
    """
    return coverage_table.merge(species_table, left_on="#rname", right_on="rname", how="left")


def aggregate_contigs(referenced_table: pd.DataFrame) -> pd.DataFrame:
    """Aggregate data over multiple contigs.

    Makes the assumption that all of the `totallength` values in a group
    are identical, and uses the first.

    Args:
        referenced_table (pd.DataFrame): Output of `samtools coverage` with
        references joined on

    Returns:
        pd.DataFrame: genome_name,length,coverage,numreads,meandepth
        Sorted in descending order of meandepth. Without genomes that
        have no reads.
    """
    aggregated = (
        referenced_table.groupby("reference")
        .apply(
            lambda contig: pd.Series(
                {
                    "length": contig.totallength.iloc[0],
                    "coverage": (contig.coverage * contig.endpos).sum() / contig.totallength.iloc[0],
                    "numreads": contig.numreads.sum(),
                    "meandepth": (contig.meandepth * contig.endpos / contig.totallength.iloc[0]).sum(),
                }
            )
        )
        .sort_values(by=["meandepth"], ascending=False)
        .reset_index()
        .rename(columns={"reference": "genome_name"})
    )

    aggregated["length"] = aggregated["length"].astype(int)
    aggregated["numreads"] = aggregated["numreads"].astype(int)

    return aggregated[aggregated["coverage"] > 0]


def lookup_and_aggregate(coverage_table: pd.DataFrame, species_table: pd.DataFrame) -> pd.DataFrame:
    """Determines the aggregated values from `samtools coverage` output where the manifest contains
    contigs.

    Args:
        coverage_table (pd.DataFrame): Output of `samtools coverage`
        species_table (pd.DataFrame): Table for rnames and references (reference data)

    Raises:
        ValueError: If there are contigs in the `samtools coverage` output that there
        aren't references for in the species table, calculations will be incorrect
        so an error is raised.

    Returns:
        pd.DataFrame: reference, overall_coverage, total_reads, totallength, mean_depth
    """
    not_in_coverage, not_in_species = unmatched_rnames(coverage_table, species_table)

    if not_in_species.size > 0:
        raise ValueError(f"Some #rnames could not be found in species list: {not_in_species.to_string()}")
    if not_in_coverage.size > 0:
        logging.info(f"Species contains rnames not in manifest: {not_in_coverage.to_string()}")

    joined = join_references(coverage_table, species_table)

    return aggregate_contigs(joined)


def validate_output(file_to_validate: Path):
    """Check if output is schema compliant (raises an error if not)

    Args:
        file_to_validate (Path): Output JSON file
    """
    with open(file_to_validate, "r", encoding="utf-8") as file:
        output = json.load(file)

    schema_path = files("competitivemapping").joinpath("competitivemapping.schema.json").read_text()

    schema = json.loads(schema_path)

    validate(instance=output, schema=schema)


def process_coverage(
    coverage_file: str | Path,
    secondary_coverage: str | Path | None,
    contigs: pd.DataFrame,
) -> pd.DataFrame:
    """Main function for producing summary Dataframe of coverage stats"""
    coverage_table = pd.read_table(coverage_file)

    aggregated = lookup_and_aggregate(coverage_table, contigs)

    if secondary_coverage:
        try:
            secondary_coverage_table = pd.read_table(secondary_coverage)
            secondary_aggregated = lookup_and_aggregate(secondary_coverage_table, contigs)
            secondary_aggregated = secondary_aggregated.rename(
                columns={
                    "coverage": "coverage_including_secondary",
                    "meandepth": "meandepth_including_secondary",
                }
            )[
                [
                    "genome_name",
                    "coverage_including_secondary",
                    "meandepth_including_secondary",
                ]
            ]

            aggregated = aggregated.merge(secondary_aggregated, on="genome_name", how="left")
        except pd.errors.EmptyDataError:
            # This should be fine, as it just means there is no secondary coverage file
            logging.info("No data found within secondary coverage file")

    return aggregated


def cli_entry_point() -> None:
    """CLI entry point."""
    args = Arguments(sys.argv[1:])

    contigs_df = pd.read_csv(args.species_list)
    aggregated = process_coverage(args.coverage, args.secondary_coverage, contigs_df)

    output = {}
    output["references"] = aggregated.to_dict(orient="records")
    with open(args.output, "w", encoding="utf-8") as file:
        json.dump(output, file, indent=4)

    validate_output(Path(args.output))


if __name__ == "__main__":
    cli_entry_point()


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
            help="Path to file created using the samtools coverage command",
        )
        parser.add_argument(
            "--secondary_coverage",
            dest="secondary_coverage",
            help="Path to file created using the samtools coverage command with secondary reads included",
            required=False,
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
        self.secondary_coverage = Path(args.secondary_coverage) if args.secondary_coverage else None
        self.species_list = Path(args.species_list)
        self.output = Path(args.output)
