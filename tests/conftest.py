import json
from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def coverage_long() -> pd.DataFrame:
    return pd.read_table("test_data/process_coverage_examples/coverage_long.tsv")


@pytest.fixture
def coverage_short() -> pd.DataFrame:
    return pd.read_table("test_data/process_coverage_examples/coverage_short.tsv")


@pytest.fixture
def species_long() -> pd.DataFrame:
    return pd.read_csv("test_data/process_coverage_examples/species_long.csv")


@pytest.fixture
def species_short() -> pd.DataFrame:
    return pd.read_csv("test_data/process_coverage_examples/species_short.csv")


@pytest.fixture
def coverage_table() -> pd.DataFrame:
    return pd.read_table(
        "test_data/process_coverage_examples/cov_WTCHG_885333_73205296.tsv"
    )


@pytest.fixture
def species_table_path() -> str:
    return "test_data/species_list_manifest_20250324.csv"


@pytest.fixture
def species_table(species_table_path) -> pd.DataFrame:
    return pd.read_csv(species_table_path)


@pytest.fixture
def manifest() -> str:
    return "data/manifest/manifest_20231001"


@pytest.fixture
def expected_joined() -> pd.DataFrame:
    return pd.read_csv("test_data/process_coverage_examples/expected_joined.csv")


@pytest.fixture
def expected_aggregated() -> pd.DataFrame:
    return pd.read_csv("test_data/process_coverage_examples/expected_aggregated.csv")


@pytest.fixture
def path_invalid_output() -> Path:
    return Path(
        "test_data/process_coverage_examples/species_comparison_report_invalid.json"
    )


@pytest.fixture
def empty_secondary_cov() -> str:
    return "test_data/process_coverage_examples/empty-secondary-coverage.tsv"


@pytest.fixture
def expected_output() -> pd.DataFrame:
    with open(
        "test_data/species_comparison_report.json", "r", encoding="utf-8"
    ) as file:
        output = json.load(file)
    return output


TB_10K = {
    "sample": "tb_10k",
    "reads": [
        "test_data/samples/illumina/TB/TB_1.fastq.gz",
        "test_data/samples/illumina/TB/TB_2.fastq.gz",
    ],
    "bam": "test_data/samples/illumina/TB/sorted_aln.bam",
    "aln_stats": "test_data/samples/illumina/TB/aln_stats.csv",
    "summary": "test_data/samples/illumina/TB/summary.json",
    "coverage": "test_data/samples/illumina/TB/coverage.tsv",
    "secondary_coverage": "test_data/samples/illumina/TB/coverage.secondary.tsv",
    "coverage_summary": "test_data/samples/illumina/TB/coverage_summary.json",
    "coverage_summary_no_secondary": "test_data/samples/illumina/TB/coverage_summary.no_secondary.json",
    "species_comparison": "test_data/samples/illumina/TB/species_comparison_report.json",
    "species_comparison_csv": "test_data/samples/illumina/TB/species_comparison.csv",
    "sylph_report": "test_data/samples/illumina/TB/sylph.tsv",
    "sylph_species_comparison": "test_data/samples/illumina/TB/sylph_species_comparison.json",
    "sylph_csv_comparison": "test_data/samples/illumina/TB/sylph_species_comparison.csv",
    "sylph_species_comparison_all_genera": "test_data/samples/illumina/TB/sylph_species_comparison_all_genera.json",
}

TB_ONT = {
    "sample": "tb_ont",
    "reads": ["test_data/samples/ont/TB/ont_lineage4_10.fastq.gz"],
    "bam": "test_data/samples/ont/TB/sorted_aln.bam",
    "aln_stats": "test_data/samples/ont/TB/aln_stats.csv",
    "summary": "test_data/samples/ont/TB/summary.json",
    "coverage": "test_data/samples/ont/TB/coverage.tsv",
    "secondary_coverage": "test_data/samples/ont/TB/coverage.secondary.tsv",
    "coverage_summary": "test_data/samples/ont/TB/coverage_summary.json",
    "coverage_summary_no_secondary": "test_data/samples/ont/TB/coverage_summary.no_secondary.json",
    "species_comparison": "test_data/samples/ont/TB/species_comparison_report.json",
    "species_comparison_csv": "test_data/samples/ont/TB/species_comparison.csv",
    "sylph_report": "test_data/samples/ont/TB/sylph.tsv",
    "sylph_species_comparison": "test_data/samples/ont/TB/sylph_species_comparison.json",
    "sylph_csv_comparison": "test_data/samples/ont/TB/sylph_species_comparison.csv",
    "sylph_species_comparison_all_genera": "test_data/samples/ont/TB/sylph_species_comparison_all_genera.json",
}

CHLORO_10k = {
    "sample": "chloro_10k",
    "reads": [
        "test_data/samples/illumina/chloro/chloro_1.fastq.gz",
        "test_data/samples/illumina/chloro/chloro_2.fastq.gz",
    ],
    "bam": "test_data/samples/illumina/chloro/sorted_aln.bam",
    "aln_stats": "test_data/samples/illumina/chloro/aln_stats.csv",
    "summary": "test_data/samples/illumina/chloro/summary.json",
    "coverage": "test_data/samples/illumina/chloro/coverage.tsv",
    "secondary_coverage": "test_data/samples/illumina/chloro/coverage.secondary.tsv",
    "coverage_summary": "test_data/samples/illumina/chloro/coverage_summary.json",
    "coverage_summary_no_secondary": "test_data/samples/illumina/chloro/coverage_summary.no_secondary.json",
    "species_comparison": "test_data/samples/illumina/chloro/species_comparison_report.json",
    "species_comparison_csv": "test_data/samples/illumina/chloro/species_comparison.csv",
    "sylph_report": "test_data/samples/illumina/chloro/sylph.tsv",
    "sylph_species_comparison": "test_data/samples/illumina/chloro/sylph_species_comparison.json",
    "sylph_csv_comparison": "test_data/samples/illumina/chloro/sylph_species_comparison.csv",
    "sylph_species_comparison_all_genera": "test_data/samples/illumina/chloro/sylph_species_comparison_all_genera.json",
}

EMPTY_SYLPH = {
    "sample": "empty",
    "reads": [
        "test_data/samples/illumina/chloro/chloro_1.fastq.gz",
        "test_data/samples/illumina/chloro/chloro_2.fastq.gz",
    ],
    "bam": "path_to_nothing.bam",
    "sylph_report": "test_data/empty_sylph/sylph.tsv",
    "manifest": "test_data/empty_sylph/empty_manifest.fasta.gz",
    "contigs": "test_data/empty_sylph/empty_contigs.csv",
    "sylph_species_comparison": "test_data/empty_sylph/sylph_species_comparison.json",
    "sylph_csv_comparison": "test_data/empty_sylph/sylph_species_comparison.csv",
}


@pytest.fixture()
def empty_sylph():
    return EMPTY_SYLPH


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


# Fixtures for the test cases in test_check_read_count.py


@pytest.fixture
def check_read_count_files() -> dict[str, Path]:
    return {
        "enough": Path("test_data/check_read_count/enough_tb.json"),
        "not_enough": Path("test_data/check_read_count/not_enough_tb.json"),
        "no_mb_value": Path("test_data/check_read_count/no_mb_value.json"),
        "no_genome_name_key": Path(
            "test_data/check_read_count/no_genome_name_key.json"
        ),
    }


@pytest.fixture
def test_outputs_dir() -> Path:
    return Path("tests/test_outputs/")


@pytest.fixture
def sylph_rep_paths() -> str:
    return "data/sylph/gtdb_genomes_reps_r220"


@pytest.fixture
def sylph_metadata() -> str:
    return "data/sylph/gtdb_r220_metadata.tsv"


@pytest.fixture
def sylph_db_A() -> dict:
    return {
        "db": "test_data/sylph_databases/A/db.syldb",
        "taxonomy": "test_data/sylph_databases/A/taxonomy.tsv",
        "metadata": "test_data/sylph_databases/A/metadata.csv",
        "genomes_dir": "test_data/sylph_databases/A",
        "sylph_report": "test_data/sylph_databases/A/sylph.tsv",
        "manifest": "test_data/sylph_databases/A/manifest.fasta.gz",
        "contigs": "test_data/sylph_databases/A/contigs.csv",
        "manifest_with_whole_genus": "test_data/sylph_databases/A/manifest_with_whole_genus.fasta.gz",
        "contigs_with_whole_genus": "test_data/sylph_databases/A/contigs_with_whole_genus.csv",
    }


@pytest.fixture
def sylph_db_B() -> dict:
    return {
        "db": "test_data/sylph_databases/B/db.syldb",
        "taxonomy": "test_data/sylph_databases/B/taxonomy.tsv",
        "genomes_dir": "test_data/sylph_databases/B",
        "sylph_report": "test_data/sylph_databases/B/sylph.tsv",
        "manifest": "test_data/sylph_databases/B/manifest.fasta.gz",
        "contigs": "test_data/sylph_databases/B/contigs.csv",
    }


@pytest.fixture
def extract_read_files() -> dict:
    prefix = "test_data/extract_reads/"
    return {
        "input": [f"{prefix}reads_1.fq", f"{prefix}reads_2.fq"],
        "alns": f"{prefix}alns.csv",
        "contigs": f"{prefix}contigs.csv",
        "reference": "ref1",
        "expectation": [f"{prefix}output_1.fq", f"{prefix}output_2.fq"],
    }
