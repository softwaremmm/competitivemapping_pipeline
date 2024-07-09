import json
import pytest
import filecmp
from jsonschema import exceptions
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


def test_lookup_and_aggregate(coverage_table, species_short):
    with pytest.raises(ValueError):
        process_mapping.lookup_and_aggregate(coverage_table, species_short)


def test_validate_output(path_invalid_output):
    with pytest.raises(exceptions.ValidationError):
        process_mapping.validate_output(path_invalid_output)


def test_cli_entry_point(samples, species_table_path, tmp_path, mocker):
    tmp_file = str(tmp_path / "output.json")
    args = [
        "process_mapping",
        "--coverage",
        samples["coverage"],
        "--secondary_coverage",
        samples["secondary_coverage"],
        "--species_list",
        species_table_path,
        "--aln_summary",
        samples["aln_stats"],
        "--counts_summary",
        samples["summary"],
        "--output",
        tmp_file,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_mapping.cli_entry_point()

    assert filecmp.cmp(tmp_file, samples["report"])


def test_cli_entry_point_empty_seconday_coverage(samples, species_table_path, tmp_path, mocker):
    tmp_file = str(tmp_path / "output.json")
    args = [
        "process_mapping",
        "--coverage",
        samples["coverage"],
        "--secondary_coverage",
        "test_data/empty-secondary-coverage.tsv",
        "--species_list",
        species_table_path,
        "--aln_summary",
        samples["aln_stats"],
        "--counts_summary",
        samples["summary"],
        "--output",
        tmp_file,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_mapping.cli_entry_point()

    # Not passing secondary reads should give a different (valid) output
    assert not filecmp.cmp(tmp_file, samples["report"])

    # Regression that it should match the output without secondary reads
    assert filecmp.cmp(tmp_file, samples["report_no_secondary"])

def test_cli_entry_point_no_seconday_coverage(samples, species_table_path, tmp_path, mocker):
    tmp_file = str(tmp_path / "output.json")
    args = [
        "process_mapping",
        "--coverage",
        samples["coverage"],
        "--species_list",
        species_table_path,
        "--aln_summary",
        samples["aln_stats"],
        "--counts_summary",
        samples["summary"],
        "--output",
        tmp_file,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_mapping.cli_entry_point()

    # Not passing secondary reads should give a different (valid) output
    assert not filecmp.cmp(tmp_file, samples["report"])

    # Regression that it should match the output without secondary reads
    assert filecmp.cmp(tmp_file, samples["report_no_secondary"])
