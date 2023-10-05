import pandas as pd
import pytest


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
def coverage_table() -> pd.DataFrame:
    return pd.read_table("test_data/cov_WTCHG_885333_73205296.tsv")


@pytest.fixture
def species_table() -> pd.DataFrame:
    return pd.read_csv("test_data/species_list_manifest_20231001.csv")


@pytest.fixture
def expected_joined() -> pd.DataFrame:
    return pd.read_csv("test_data/expected_joined.csv")


@pytest.fixture
def expected_aggregated() -> pd.DataFrame:
    return pd.read_csv("test_data/expected_aggregated.csv")
