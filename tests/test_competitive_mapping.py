import os
from gzip import open as gzopen

from Bio import SeqIO
from test_utils import check_file

from competitivemapping import competitive_mapping


def test_cli_entry_point_manifest(
    samples, species_table_path, manifest, test_outputs_dir, mocker
):
    output_root = os.path.join(test_outputs_dir, samples["sample"] + ".")

    seq_platform = "ont" if len(samples["reads"]) == 1 else "illumina"

    args = [
        "competitive_mapping",
        "manifest",
        "--reads",
        " ".join(samples["reads"]),
        "--seq_platform",
        seq_platform,
        "--ref_for_fastq",
        "M.tuberculosis",
        "--manifest",
        manifest,
        "--contigs",
        species_table_path,
        "--output_root",
        output_root,
        "--cpus",
        "4",
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    competitive_mapping.cli_entry_point()

    check_file(samples["species_comparison"], output_root + "species_comparison.json")
    check_file(samples["sylph_csv_comparison"], output_root + "species_comparison.csv")

    def check_length(path, length):
        assert len(list(SeqIO.parse(gzopen(path, "rt"), format="fastq"))) == length

    match samples["sample"]:
        case "chloro_10k":
            check_length(output_root + "reads_1.fastq.gz", 2736)
            check_length(output_root + "reads_2.fastq.gz", 2736)
        case "tb_10k":
            check_length(output_root + "reads_1.fastq.gz", 9648)
            check_length(output_root + "reads_2.fastq.gz", 9648)
        case "tb_ont":
            check_length(output_root + "reads.fastq.gz", 1042)


def test_cli_entry_point_sylph(
    samples, sylph_db, sylph_metadata, test_outputs_dir, mocker
):
    output_root = os.path.join(test_outputs_dir, "cm_sylph", samples["sample"] + ".")
    os.makedirs(os.path.join(test_outputs_dir, "cm_sylph"), exist_ok=True)

    seq_platform = "ont" if len(samples["reads"]) == 1 else "illumina"

    args = [
        "competitive_mapping",
        "sylph",
        "--reads",
        " ".join(samples["reads"]),
        "--seq_platform",
        seq_platform,
        "--ref_for_fastq",
        "GCF_000195955.2",  # M.tuberculosis
        "--sylph_report",
        samples["sylph_report"],
        "--genomes",
        sylph_db,
        "--db_metadata",
        sylph_metadata,
        "--output_root",
        output_root,
        "--cpus",
        "4",
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    competitive_mapping.cli_entry_point()

    check_file(
        samples["sylph_species_comparison"], output_root + "species_comparison.json"
    )
    check_file(samples["sylph_csv_comparison"], output_root + "species_comparison.csv")

    def check_length(path, length):
        assert len(list(SeqIO.parse(gzopen(path, "rt"), format="fastq"))) == length

    match samples["sample"]:
        case "chloro_10k":
            check_length(output_root + "reads_1.fastq.gz", 14)
            check_length(output_root + "reads_2.fastq.gz", 14)
        case "tb_10k":
            check_length(output_root + "reads_1.fastq.gz", 9940)
            check_length(output_root + "reads_2.fastq.gz", 9940)
        case "tb_ont":
            check_length(output_root + "reads.fastq.gz", 1044)


def test_empty_sylph(empty_sylph, sylph_db, sylph_metadata, test_outputs_dir, mocker):
    output_root = os.path.join(
        test_outputs_dir, "cm_sylph", empty_sylph["sample"] + "."
    )
    os.makedirs(os.path.join(test_outputs_dir, "cm_sylph"), exist_ok=True)

    seq_platform = "ont" if len(empty_sylph["reads"]) == 1 else "illumina"

    args = [
        "competitive_mapping",
        "sylph",
        "--reads",
        " ".join(empty_sylph["reads"]),
        "--seq_platform",
        seq_platform,
        "--ref_for_fastq",
        "GCF_000195955.2",  # M.tuberculosis
        "--sylph_report",
        empty_sylph["sylph_report"],
        "--genomes",
        sylph_db,
        "--db_metadata",
        sylph_metadata,
        "--output_root",
        output_root,
        "--cpus",
        "4",
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    competitive_mapping.cli_entry_point()

    check_file(
        empty_sylph["sylph_species_comparison"], output_root + "species_comparison.json"
    )
    check_file(
        empty_sylph["sylph_csv_comparison"], output_root + "species_comparison.csv"
    )
