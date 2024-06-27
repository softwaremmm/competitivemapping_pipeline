import pandas as pd
import pytest
import filecmp

import competitivemapping.process_aln_stats as process_aln_stats


def test_summarise_by_chrom():
    chroms = {"a", "b", "c"}
    read_info = {
        "read1": {
            "primary": "a",
            "secondary": [],
            "supplementary": [],
        },
        "read2": {
            "primary": "a",
            "secondary": [],
            "supplementary": ["a", "a"],
        },
        "read3": {
            "primary": "a",
            "secondary": ["a"],
            "supplementary": [],
        },
        "read4": {
            "primary": "b",
            "secondary": ["c", "c"],
            "supplementary": ["a"],
        },
        "read5": {
            "primary": "",
            "secondary": ["b"],
            "supplementary": ["c"],
        },
    }

    expected_dict = {
        "a": {
            "total_reads": 4,  # Total number of reads which map (in any way) to chrom
            "total_alns": 7,  # Total number of alignments, so will count supplementary separately
            "exclusive_reads": 3,  # reads which only maps to this chrom
            "primary_reads": 3,  # reads which map best to this chrom
            "secondary_reads": 0,  # reads which map to this chrom but not with primary
            "supplementary_reads": 1,  # reads which have supplementary alignments but not primary
            "supplementary_alns": 3,  # number of supplementary alignments
        },
        "b": {
            "total_reads": 2,
            "total_alns": 2,
            "exclusive_reads": 0,
            "primary_reads": 1,
            "secondary_reads": 1,
            "supplementary_reads": 0,
            "supplementary_alns": 0,
        },
        "c": {
            "total_reads": 2,
            "total_alns": 3,
            "exclusive_reads": 0,
            "primary_reads": 0,
            "secondary_reads": 1,
            "supplementary_reads": 1,
            "supplementary_alns": 1,
        },
    }
    df = pd.DataFrame.from_dict(expected_dict, orient="index")
    df.reset_index(inplace=True, names="genome_name")

    result = process_aln_stats.summarise_by_chrom(chroms, read_info)
    assert result.reset_index(drop=True).equals(df.reset_index(drop=True))


def test_aln_stats(samples, species_table_path, tmp_path, mocker):
    tmp_stats = str(tmp_path / "aln.csv")
    tmp_summary = str(tmp_path / "summary.csv")
    args = [
        "process_aln_stats",
        "--bam",
        samples["bam"],
        "--species_list",
        species_table_path,
        "--output",
        tmp_stats,
        "--output_summary",
        tmp_summary,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_aln_stats.cli_entry_point()

    assert filecmp.cmp(tmp_stats, samples["aln_stats"])
    assert filecmp.cmp(tmp_summary, samples["summary"])


def test_no_secondary(samples, species_table_path, tmp_path, mocker):
    tmp_stats = str(tmp_path / "aln.csv")
    tmp_summary = str(tmp_path / "summary.csv")
    args = [
        "process_aln_stats",
        "--bam",
        samples["bam"],
        "--species_list",
        species_table_path,
        "--output",
        tmp_stats,
        "--output_summary",
        tmp_summary,
        "--exclude_secondary",
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    process_aln_stats.cli_entry_point()

    df = pd.read_csv(tmp_stats)
    assert (df["secondary_reads"] == 0).all()
