import pandas as pd
import pytest
import filecmp

import competitivemapping.process_aln_stats as process_aln_stats


def test_illumina(
    illumina_bam, illumina_aln_stats, species_table_path, tmp_path, mocker
):
    tmp_file = str(tmp_path / "aln.csv")
    args = [
        "process_aln_stats",
        "--bam",
        illumina_bam,
        "--species_list",
        species_table_path,
        "--output",
        tmp_file,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_aln_stats.cli_entry_point()

    assert filecmp.cmp(tmp_file, illumina_aln_stats)


def test_ont(ont_bam, ont_aln_stats, species_table_path, tmp_path, mocker):
    tmp_file = str(tmp_path / "aln.csv")
    args = [
        "process_aln_stats",
        "--bam",
        ont_bam,
        "--species_list",
        species_table_path,
        "--output",
        tmp_file,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_aln_stats.cli_entry_point()

    assert filecmp.cmp(tmp_file, ont_aln_stats)


def test_no_secondary(
    illumina_bam, species_table_path, tmp_path, mocker
):
    tmp_file = str(tmp_path / "aln.csv")
    args = [
        "process_aln_stats",
        "--bam",
        illumina_bam,
        "--species_list",
        species_table_path,
        "--output",
        tmp_file,
        "--exclude_secondary",
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_aln_stats.cli_entry_point()
    df = pd.read_csv(tmp_file)
    assert (df["secondary_reads"] == 0).all()
