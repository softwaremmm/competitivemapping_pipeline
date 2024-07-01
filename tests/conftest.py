import json
from pathlib import Path
import pytest

import pandas as pd


@pytest.fixture
def coverage_long() -> pd.DataFrame:
    return pd.read_table("test_data/matching/coverage_long.tsv")


@pytest.fixture
def coverage_short() -> pd.DataFrame:
    return pd.read_table("test_data/matching/coverage_short.tsv")


@pytest.fixture
def species_long() -> pd.DataFrame:
    return pd.read_csv("test_data/matching/species_long.csv")


@pytest.fixture
def species_short() -> pd.DataFrame:
    return pd.read_csv("test_data/matching/species_short.csv")


@pytest.fixture
def coverage_table_path() -> str:
    return "test_data/cov_WTCHG_885333_73205296.tsv"


@pytest.fixture
def coverage_table(coverage_table_path) -> pd.DataFrame:
    return pd.read_table(coverage_table_path)


@pytest.fixture
def species_table_path() -> str:
    return "test_data/species_list_manifest_20231001.csv"


@pytest.fixture
def species_table(species_table_path) -> pd.DataFrame:
    return pd.read_csv(species_table_path)


@pytest.fixture
def expected_joined() -> pd.DataFrame:
    return pd.read_csv("test_data/expected_joined.csv")


@pytest.fixture
def expected_aggregated() -> pd.DataFrame:
    return pd.read_csv("test_data/expected_aggregated.csv")


@pytest.fixture
def path_invalid_output() -> pd.DataFrame:
    return Path("test_data/species_comparison_report_invalid.json")


@pytest.fixture
def expected_output() -> pd.DataFrame:
    with open("test_data/species_comparison_report.json", "r", encoding="utf-8") as file:
        output = json.load(file)
    return output


TB_10K = {
    "bam": "test_data/TB_10k/sorted_aln.bam",
    "aln_stats": "test_data/TB_10k/aln_stats.csv",
    "summary": "test_data/TB_10k/summary.json",
    "coverage": "test_data/TB_10k/coverage.tsv",
    "secondary_coverage": "test_data/TB_10k/coverage.secondary.tsv",
    "report": "test_data/TB_10k/species_comparison_report.json",
    "report_no_secondary": "test_data/TB_10k/species_comparison_report_no_secondary.json",
}

TB_ONT = {
    "bam": "test_data/TB_ont/sorted_aln.bam",
    "aln_stats": "test_data/TB_ont/aln_stats.csv",
    "summary": "test_data/TB_ont/summary.json",
    "coverage": "test_data/TB_ont/coverage.tsv",
    "secondary_coverage": "test_data/TB_ont/coverage.secondary.tsv",
    "report": "test_data/TB_ont/species_comparison_report.json",
    "report_no_secondary": "test_data/TB_ont/species_comparison_report_no_secondary.json",
}

CHLORO_10k = {
    "bam": "test_data/chloro_10k/sorted_aln.bam",
    "aln_stats": "test_data/chloro_10k/aln_stats.csv",
    "summary": "test_data/chloro_10k/summary.json",
    "coverage": "test_data/chloro_10k/coverage.tsv",
    "secondary_coverage": "test_data/chloro_10k/coverage.secondary.tsv",
    "report": "test_data/chloro_10k/species_comparison_report.json",
    "report_no_secondary": "test_data/chloro_10k/species_comparison_report_no_secondary.json",
}


@pytest.fixture(
    params=[
        TB_10K,
        TB_ONT,
        CHLORO_10k,
    ],
    ids=[
        "illumina",
        "ont",
        "chloro",
    ],
)
def samples(request) -> dict:
    return request.param
