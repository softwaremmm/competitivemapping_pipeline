# pylint: disable=logging-fstring-interpolation
"""Competitive Mapping"""

import logging
import sys
import pandas as pd

from competitivemapping.cli_args import Arguments

logging.basicConfig(
    format="%(asctime)s — %(name)s — %(levelname)s — %(funcName)s:%(lineno)d — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S%z",
    level=logging.DEBUG,
)


def unmatched_rnames(coverage_table: pd.DataFrame, species_table: pd.DataFrame) -> [pd.Series, pd.Series]:
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
    """Aggregate coverage over multiple contigs.

    Makes the assumption that all of the `totallength` values in a group
    are identical, and uses the first.

    Args:
        referenced_table (pd.DataFrame): Output of `samtools coverage` with
        references joined on

    Returns:
        pd.DataFrame: Aggregated coverage for each reference
    """
    return (
        referenced_table.groupby("reference")
        .apply(
            lambda x: pd.Series({"overall_coverage": (x["coverage"] * x["endpos"]).sum() / x["totallength"].iloc[0]})
        )
        .reset_index()
    )


def determine_overall_coverage(coverage_table: pd.DataFrame, species_table: pd.DataFrame) -> pd.DataFrame:
    """Determines the "overall coverage" from `samtools coverage` output where the manifest contains
    contigs over which coverage must be aggregated.

    Args:
        coverage_table (pd.DataFrame): Output of `samtools coverage`
        species_table (pd.DataFrame): Table for rnames and references (reference data)

    Raises:
        ValueError: If there are contigs in the `samtools coverage` output that there
        aren't references for in the species table, calculations will be incorrect
        so an error is raised.

    Returns:
        pd.DataFrame: Aggregated coverage for each reference
    """
    not_in_coverage, not_in_species = unmatched_rnames(coverage_table, species_table)

    if not_in_species.size > 0:
        raise ValueError(f"Some #rnames could not be found in species list: {not_in_species.to_string()}")
    if not_in_coverage.size > 0:
        logging.info(f"Species contains rnames not in manifest: {not_in_coverage.to_string()}")

    joined = join_references(coverage_table, species_table)

    return aggregate_contigs(joined)


def cli_entry_point() -> None:
    """CLI entry point."""
    Arguments(sys.argv[1:])
