import json
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
def expected_output() -> pd.DataFrame:
    with open("test_data/expected_output.json", "r") as file:
        output = json.load(file)
    return output
