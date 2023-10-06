import json
import pytest

import pandas as pd

import competitivemapping.process_mapping as process_mapping


def test_not_in_species(coverage_long, species_short):
    not_in_coverage, not_in_species = process_mapping.unmatched_rnames(coverage_long, species_short)
    assert not_in_species.size > 0


def test_not_in_coverage(coverage_short, species_long):
    not_in_coverage, not_in_species = process_mapping.unmatched_rnames(coverage_short, species_long)
    assert not_in_coverage.size > 0


def test_full_match(coverage_long, species_long):
    not_in_coverage, not_in_species = process_mapping.unmatched_rnames(coverage_long, species_long)
    assert not_in_coverage.size == 0
    assert not_in_species.size == 0


def test_join_references(coverage_table, species_table, expected_joined):
    actual_joined = process_mapping.join_references(coverage_table, species_table)
    pd.testing.assert_frame_equal(actual_joined, expected_joined)


def test_aggregate_contigs(expected_joined, expected_aggregated):
    actual_aggregated = process_mapping.aggregate_contigs(expected_joined)
    pd.testing.assert_frame_equal(actual_aggregated, expected_aggregated)


def test_lookup_and_aggregate(coverage_table, species_table, expected_aggregated):
    actual_aggregated = process_mapping.determine_overall_coverage(coverage_table, species_table)
    pd.testing.assert_frame_equal(actual_aggregated, expected_aggregated)


def test_lookup_and_aggregate(coverage_table, species_short):
    with pytest.raises(ValueError):
        process_mapping.lookup_and_aggregate(coverage_table, species_short)


def test_cli_entry_point(coverage_table_path, species_table_path, tmp_path, mocker, expected_output):
    tmp_file = str(tmp_path / "output.json")
    args: list = [
        "process_mapping",
        "--coverage",
        coverage_table_path,
        "--species_list",
        species_table_path,
        "--output",
        tmp_file,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_mapping.cli_entry_point()

    with open(tmp_file, "r") as file:
        actual_output = json.load(file)

    assert actual_output == expected_output
