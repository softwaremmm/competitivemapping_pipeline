# pylint: disable=logging-fstring-interpolation
"""Process output from Competitive Mapping"""

import json
import logging
from pathlib import Path
import sys
from importlib.resources import files

from jsonschema import validate
import pandas as pd

import competitivemapping.cli_args

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
        Sorted in descending order of coverage. Without genomes that
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
        .sort_values(by=["coverage"], ascending=False)
        .reset_index()
        .rename(columns={"reference": "genome_name"})
    )

    filtered = aggregated[aggregated["coverage"] > 0]

    return filtered


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


def cli_entry_point() -> None:
    """CLI entry point."""
    args = competitivemapping.cli_args.Arguments(sys.argv[1:])

    coverage_table = pd.read_table(args.coverage)
    species_table = pd.read_csv(args.species_list)

    aggregated = lookup_and_aggregate(coverage_table, species_table)

    if args.aln_summary:
        # Align summary has more detailed break down of number of reads/alns
        aln_summary = pd.read_csv(args.aln_summary)
        aggregated = aggregated.merge(aln_summary, on="genome_name", how="left")
        # numreads superseeded by total_reads. Could remove in future
        # aggregated["numreads"] = aggregated["total_reads"]

    aggregated.to_json(args.output, orient="records", indent=4)

    validate_output(Path(args.output))
