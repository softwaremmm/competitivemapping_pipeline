import filecmp

import pandas as pd
import pytest
from jsonschema import exceptions
from test_utils import check_file

import competitivemapping.process_coverage as process_coverage


def test_not_in_species(coverage_long, species_short):
    not_in_coverage, not_in_species = process_coverage.unmatched_rnames(
        coverage_long, species_short
    )
    assert not_in_species.size > 0


def test_not_in_coverage(coverage_short, species_long):
    not_in_coverage, not_in_species = process_coverage.unmatched_rnames(
        coverage_short, species_long
    )
    assert not_in_coverage.size > 0


def test_full_match(coverage_long, species_long):
    not_in_coverage, not_in_species = process_coverage.unmatched_rnames(
        coverage_long, species_long
    )
    assert not_in_coverage.size == 0
    assert not_in_species.size == 0


def test_join_references(coverage_table, species_table, expected_joined):
    actual_joined = process_coverage.join_references(coverage_table, species_table)
    if not actual_joined.equals(expected_joined):
        # Since different to expectation, save to file for inspection
        actual_joined.to_csv("tests/test_outputs/actual_joined.csv", index=True)
    pd.testing.assert_frame_equal(actual_joined, expected_joined)


def test_aggregate_contigs(expected_joined, expected_aggregated):
    actual_aggregated = process_coverage.aggregate_contigs(expected_joined)
    if not actual_aggregated.equals(expected_aggregated):
        # Since different to expectation, save to file for inspection
        actual_aggregated.to_csv("tests/test_outputs/actual_aggregated.csv", index=True)
    pd.testing.assert_frame_equal(actual_aggregated, expected_aggregated)


def test_lookup_and_aggregate(coverage_table, species_short):
    with pytest.raises(ValueError):
        process_coverage.lookup_and_aggregate(coverage_table, species_short)


def test_validate_output(path_invalid_output):
    with pytest.raises(exceptions.ValidationError):
        process_coverage.validate_output(path_invalid_output)


def test_cli_entry_point(samples, species_table_path, tmp_path, mocker):
    tmp_file = str(tmp_path / "output.json")
    args = [
        "process_coverage",
        "--coverage",
        samples["coverage"],
        "--secondary_coverage",
        samples["secondary_coverage"],
        "--species_list",
        species_table_path,
        "--output",
        tmp_file,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_coverage.cli_entry_point()

    check_file(samples["coverage_summary"], tmp_file)
    # assert filecmp.cmp(tmp_file, samples["report"])


def test_cli_entry_point_empty_seconday_coverage(
    samples, species_table_path, tmp_path, mocker
):
    tmp_file = str(tmp_path / "output.json")
    args = [
        "process_coverage",
        "--coverage",
        samples["coverage"],
        "--secondary_coverage",
        "test_data/empty-secondary-coverage.tsv",
        "--species_list",
        species_table_path,
        "--output",
        tmp_file,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_coverage.cli_entry_point()

    # Not passing secondary reads should give a different (valid) output
    assert not filecmp.cmp(tmp_file, samples["coverage_summary"])

    check_file(samples["coverage_summary_no_secondary"], tmp_file)
    # assert filecmp.cmp(tmp_file, samples["report_no_secondary"])


def test_cli_entry_point_no_seconday_coverage(
    samples, species_table_path, tmp_path, mocker
):
    tmp_file = str(tmp_path / "output.json")
    args = [
        "process_coverage",
        "--coverage",
        samples["coverage"],
        "--species_list",
        species_table_path,
        "--output",
        tmp_file,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_coverage.cli_entry_point()

    check_file(samples["coverage_summary_no_secondary"], tmp_file)
    # assert filecmp.cmp(tmp_file, samples["report_no_secondary"])
