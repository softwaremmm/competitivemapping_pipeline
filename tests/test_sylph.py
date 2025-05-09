"""Test the sylph process which requires building, mapping and analysing.
Only done end to end due to large intermediate files."""
import os
from gzip import open as gzopen

from Bio import SeqIO
from test_utils import check_file, get_cpus

from competitivemapping import competitive_mapping, manifest_builder, manifest_mapper


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
            str(get_cpus()),
        ],
    )
    manifest_builder.cli_entry_point()

    # Then need to map reads
    mocker.patch(
        "sys.argv",
        [
            "manifest_mapper",
            "--manifest",
            output_root + "manifest.fasta.gz",
            "--reads",
            " ".join(samples["reads"]),
            "--seq_platform",
            seq_platform,
            "--cpus",
            str(get_cpus()),
            "-o",
            output_root + "aln.bam",
        ],
    )
    manifest_mapper.cli_entry_point()

    mocker.patch(
        "sys.argv",
        [
            "competitive_mapping",
            "--input_bam",
            output_root + "aln.bam",
            "--contigs",
            output_root + "contigs.csv",
            "--seq_platform",
            seq_platform,
            "--ref_for_fastq",
            "GCF_000195955.2",  # M.tuberculosis
            "--output_root",
            output_root,
            "--cpus",
            str(get_cpus()),
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
            check_length(output_root + "reads_1.fastq.gz", 10000)
            check_length(output_root + "reads_2.fastq.gz", 10000)
        case "tb_ont":
            check_length(output_root + "reads.fastq.gz", 1000)


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
            str(get_cpus()),
        ],
    )
    manifest_builder.cli_entry_point()

    # Then need to map reads
    mocker.patch(
        "sys.argv",
        [
            "manifest_mapper",
            "--manifest",
            output_root + "manifest.fasta.gz",
            "--reads",
            " ".join(samples["reads"]),
            "--seq_platform",
            seq_platform,
            "--cpus",
            str(get_cpus()),
            "-o",
            output_root + "aln.bam",
        ],
    )
    manifest_mapper.cli_entry_point()

    mocker.patch(
        "sys.argv",
        [
            "competitive_mapping",
            "--input_bam",
            output_root + "aln.bam",
            "--contigs",
            output_root + "contigs.csv",
            "--seq_platform",
            seq_platform,
            "--ref_for_fastq",
            "GCF_000195955.2",  # M.tuberculosis
            "--output_root",
            output_root,
            "--cpus",
            str(get_cpus()),
        ],
    )
    competitive_mapping.cli_entry_point()

    check_file(
        samples["sylph_species_comparison_all_genera"],
        output_root + "species_comparison.json",
    )
