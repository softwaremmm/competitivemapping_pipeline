"""
Testing for manifest builder.
Note that majority of tests are in test_sylph.py
as manifest is a large intermediate file for competitive mapping.
"""

import os
from gzip import open as gzopen

import pandas as pd
from test_utils import check_file

from competitivemapping import manifest_builder


def test_get_base_species_name():
    assert (
        manifest_builder.get_base_species_name("Streptococcus mitis")
        == "Streptococcus mitis"
    )
    assert (
        manifest_builder.get_base_species_name("Haemophilus_D parainfluenzae_Y")
        == "Haemophilus_D parainfluenzae"
    )
    assert (
        manifest_builder.get_base_species_name("Aerococcus urinae_E")
        == "Aerococcus urinae"
    )
    assert (
        manifest_builder.get_base_species_name("Haemophilus_D sp015255025")
        == "Haemophilus_D sp015255025"
    )


def test_select_extra_species():
    df = pd.DataFrame(
        {
            "species": [
                "sp12345678",
                "Streptococcus mitis_CU",
                "Streptococcus mitis",
                "Streptococcus oralis_V",
                "Streptococcus oralis_E",
                "Streptococcus halitosis",
                "Streptococcus sp943736975",
            ]
        }
    )

    assert set(
        manifest_builder.select_extra_species(["other"], df)["species"].tolist()
    ) == {
        "Streptococcus mitis",
        "Streptococcus oralis_E",
        "Streptococcus halitosis",
    }

    assert set(
        manifest_builder.select_extra_species(["Streptococcus mitis_CU", "other"], df)[
            "species"
        ].tolist()
    ) == {
        "Streptococcus oralis_E",
        "Streptococcus halitosis",
    }


def test_get_genome_paths(test_outputs_dir):
    test_dir = os.path.join(test_outputs_dir, "test_get_genome_paths")
    os.makedirs(test_dir, exist_ok=True)
    df = pd.DataFrame(
        {
            "filename": [
                "species_A.fna.gz",
                "species_B.fna.gz",
            ],
            "path": [
                "genomes/",
                "genomes/",
            ],
        }
    )
    expectation = pd.DataFrame(
        {
            "filename": [
                "species_A.fna.gz",
                "species_B.fna.gz",
            ],
            "path": [
                os.path.join(test_dir, "genomes/", "species_A.fna.gz"),
                os.path.join(test_dir, "genomes/", "species_B.fna.gz"),
            ],
        }
    )

    paths_file = os.path.join(test_dir, "genome_paths.csv")
    df.to_csv(
        paths_file,
        sep="\t",
        index=False,
        header=False,
    )

    result_df = manifest_builder.get_genome_paths(paths_file)

    if not expectation.equals(result_df):
        # If the dataframes are not equal, we can print the differences
        expectation.to_csv(
            os.path.join(test_dir, "expected_genome_paths.csv"),
            sep="\t",
            index=False,
        )
        result_df.to_csv(
            os.path.join(test_dir, "result_genome_paths.csv"),
            sep="\t",
            index=False,
        )

    assert expectation.equals(result_df)


def test_empty_sylph(
    empty_sylph, sylph_rep_paths, sylph_metadata, test_outputs_dir, mocker
):
    output_root = os.path.join(
        test_outputs_dir, "build_manifest", empty_sylph["sample"] + "."
    )
    os.makedirs(os.path.join(test_outputs_dir, "build_manifest"), exist_ok=True)

    mocker.patch(
        "sys.argv",
        [
            "manifest_builder",
            "--sylph_report",
            empty_sylph["sylph_report"],
            "--genome_path_files",
            sylph_rep_paths,
            "--taxonomy_files",
            sylph_metadata,
            "--output_root",
            output_root,
            "--cpus",
            "4",
        ],
    )

    manifest_builder.cli_entry_point()

    # compare contents of the manifest file which is gzipped
    with (
        gzopen(empty_sylph["manifest"], "rt") as f_expectation,
        gzopen(output_root + "manifest.fasta.gz", "rt") as f_result,
    ):
        expectation = f_expectation.read()
        result = f_result.read()
        assert expectation == result

    check_file(empty_sylph["contigs"], output_root + "contigs.csv")
