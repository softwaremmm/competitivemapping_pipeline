"""Competitive Mapping"""

import sys
import pandas as pd

from competitivemapping.cli_args import Arguments


def join_references(coverage_table: pd.DataFrame, species_table: pd.DataFrame) -> pd.DataFrame:
    """Look up references from rnames.

    Args:
        coverage_table (pd.DataFrame): Output of `samtools coverage`
        species_table (pd.DataFrame): Table for rnames and references (reference data)

    Returns:
        pd.DataFrame: As `coverage_table`, but with an extra column for reference
    """
    return coverage_table.merge(species_table, left_on="#rname", right_on="rname", how="left")


def cli_entry_point() -> None:
    """CLI entry point."""
    Arguments(sys.argv[1:])
