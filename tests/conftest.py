from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def myco_manifest() -> dict:
    return {
        "manifest": "test_data/myco_manifest/manifest.fasta.gz",
        "contigs": "test_data/myco_manifest/contigs.csv",
    }


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
    return "test_data/process_coverage_examples/species_list_manifest_20250324.csv"


@pytest.fixture
def species_table() -> pd.DataFrame:
    return pd.read_csv(
        "test_data/process_coverage_examples/species_list_manifest_20250324.csv"
    )


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
def sylph_db_A() -> dict:
    return {
        "db": "test_data/sylph_databases/A/db.syldb",
        "taxonomy": "test_data/sylph_databases/A/taxonomy.tsv",
        "metadata": "test_data/sylph_databases/A/metadata.csv",
        "genomes_dir": "test_data/sylph_databases/A",
        "sylph_report": "test_data/sylph_databases/A/sylph.tsv",
        "manifest": "test_data/sylph_databases/A/manifest.fasta.gz",
        "contigs": "test_data/sylph_databases/A/contigs.csv",
        "contigs_with_ani_group": "test_data/sylph_databases/A/contigs_with_ani_group.csv",
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
def sylph_db_C() -> dict:
    return {
        "db": "test_data/sylph_databases/C/db.syldb",
        "taxonomy": "test_data/sylph_databases/C/taxonomy.tsv",
        "genomes_dir": "test_data/sylph_databases/C",
        "sylph_report": "test_data/sylph_databases/C/sylph.tsv",
        "manifest": "test_data/sylph_databases/C/manifest.fasta.gz",
        "contigs": "test_data/sylph_databases/C/contigs.csv",
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
}

TB_ONT = {
    "sample": "tb_ont",
    "reads": ["test_data/samples/ont/TB/lin4_10.fastq.gz"],
    "bam": "test_data/samples/ont/TB/sorted_aln.bam",
    "aln_stats": "test_data/samples/ont/TB/aln_stats.csv",
    "summary": "test_data/samples/ont/TB/summary.json",
    "coverage": "test_data/samples/ont/TB/coverage.tsv",
    "secondary_coverage": "test_data/samples/ont/TB/coverage.secondary.tsv",
    "coverage_summary": "test_data/samples/ont/TB/coverage_summary.json",
    "coverage_summary_no_secondary": "test_data/samples/ont/TB/coverage_summary.no_secondary.json",
    "species_comparison": "test_data/samples/ont/TB/species_comparison_report.json",
    "species_comparison_csv": "test_data/samples/ont/TB/species_comparison.csv",
}

CHLORO_10K = {
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
}

EMPTY_FILES = {
    "sample": "empty",
    "reads": [
        "test_data/samples/illumina/chloro/chloro_1.fastq.gz",
        "test_data/samples/illumina/chloro/chloro_2.fastq.gz",
    ],
    "sylph_query": "test_data/empty/sylph_query.tsv",
    "sylph_profile": "test_data/empty/sylph_profile.tsv",
    "manifest": "test_data/empty/manifest.fasta.gz",
    "contigs": "test_data/empty/contigs.csv",
    "species_comparison": "test_data/empty/species_comparison.json",
    "csv_comparison": "test_data/empty/species_comparison.csv",
}


# Reads are genuine
@pytest.fixture()
def empty_files():
    return EMPTY_FILES


@pytest.fixture(
    params=[
        TB_10K,
        TB_ONT,
        CHLORO_10K,
    ],
    ids=[
        "illumina",
        "ont",
        "chloro",
    ],
)
def samples(request) -> dict:
    return request.param
