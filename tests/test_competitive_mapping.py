import os
from gzip import open as gzopen

from Bio import SeqIO
from test_utils import check_file

from competitivemapping import competitive_mapping, manifest_builder


def test_cli_entry_point(
    samples, species_table_path, manifest, test_outputs_dir, mocker
):
    output_root = os.path.join(test_outputs_dir, samples["sample"] + ".")

    seq_platform = "ont" if len(samples["reads"]) == 1 else "illumina"

    args = [
        "competitive_mapping",
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
    check_file(
        samples["species_comparison_csv"], output_root + "species_comparison.csv"
    )

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


def test_empty_manifest(empty_sylph, test_outputs_dir, mocker):
    output_root = os.path.join(
        test_outputs_dir, "cm_sylph", empty_sylph["sample"] + "."
    )
    os.makedirs(os.path.join(test_outputs_dir, "cm_sylph"), exist_ok=True)

    seq_platform = "ont" if len(empty_sylph["reads"]) == 1 else "illumina"

    args = [
        "competitive_mapping",
        "--manifest",
        empty_sylph["manifest"],
        "--contigs",
        empty_sylph["contigs"],
        "--reads",
        " ".join(empty_sylph["reads"]),
        "--seq_platform",
        seq_platform,
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


def test_sylph_manifest(
    samples, sylph_rep_paths, sylph_metadata, test_outputs_dir, mocker
):
    output_root = os.path.join(test_outputs_dir, "cm_sylph", samples["sample"] + ".")
    os.makedirs(os.path.join(test_outputs_dir, "cm_sylph"), exist_ok=True)

    seq_platform = "ont" if len(samples["reads"]) == 1 else "illumina"

    # make manifest first
    mocker.patch(
        "sys.argv",
        [
            "manifest_builder",
            "--sylph_report",
            samples["sylph_report"],
            "--genome_dirs",
            sylph_rep_paths,
            "--metadata_files",
            sylph_metadata,
            "--output_root",
            output_root,
            "--cpus",
            "4",
        ],
    )
    manifest_builder.cli_entry_point()

    mocker.patch(
        "sys.argv",
        [
            "competitive_mapping",
            "--manifest",
            output_root + "manifest.fasta.gz",
            "--contigs",
            output_root + "contigs.csv",
            "--reads",
            " ".join(samples["reads"]),
            "--seq_platform",
            seq_platform,
            "--ref_for_fastq",
            "GCF_000195955.2",  # M.tuberculosis
            "--output_root",
            output_root,
            "--cpus",
            "10",
        ],
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


def test_sylph_manifest_all_genera(
    samples, sylph_rep_paths, sylph_metadata, test_outputs_dir, mocker
):
    output_root = os.path.join(
        test_outputs_dir, "cm_sylph_all_genera", samples["sample"] + "."
    )
    os.makedirs(os.path.join(test_outputs_dir, "cm_sylph_all_genera"), exist_ok=True)

    seq_platform = "ont" if len(samples["reads"]) == 1 else "illumina"

    # make manifest first
    mocker.patch(
        "sys.argv",
        [
            "manifest_builder",
            "--include_whole_genus",
            "--sylph_report",
            samples["sylph_report"],
            "--genome_dirs",
            sylph_rep_paths,
            "--metadata_files",
            sylph_metadata,
            "--output_root",
            output_root,
            "--cpus",
            "4",
        ],
    )
    manifest_builder.cli_entry_point()

    mocker.patch(
        "sys.argv",
        [
            "competitive_mapping",
            "--manifest",
            output_root + "manifest.fasta.gz",
            "--contigs",
            output_root + "contigs.csv",
            "--reads",
            " ".join(samples["reads"]),
            "--seq_platform",
            seq_platform,
            "--ref_for_fastq",
            "GCF_000195955.2",  # M.tuberculosis
            "--output_root",
            output_root,
            "--cpus",
            "10",
        ],
    )
    competitive_mapping.cli_entry_point()

    check_file(
        samples["sylph_species_comparison_all_genera"],
        output_root + "species_comparison.json",
    )
