import hashlib
import os
from gzip import open as gzopen

from Bio import SeqIO
from test_utils import check_file, get_cpus

from competitivemapping import competitive_mapping


def test_cli_entry_point(samples, species_table_path, test_outputs_dir, mocker):
    output_root = os.path.join(test_outputs_dir, samples["sample"] + ".")

    seq_platform = "ont" if len(samples["reads"]) == 1 else "illumina"

    args = [
        "competitive_mapping",
        "--input_bam",
        samples["bam"],
        "--seq_platform",
        seq_platform,
        "--ref_for_fastq",
        "M.tuberculosis",
        "--contigs",
        species_table_path,
        "--output_root",
        output_root,
        "--cpus",
        str(get_cpus()),
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

    def check_order(path, hash_expected):
        read_names = [
            record.id for record in SeqIO.parse(gzopen(path, "rt"), format="fastq")
        ]
        print(",".join(read_names))
        read_hash = hashlib.sha256(",".join(read_names).encode()).hexdigest()
        assert read_hash == hash_expected

    match samples["sample"]:
        case "chloro_10k":
            check_length(output_root + "reads_1.fastq.gz", 2736)
            check_length(output_root + "reads_2.fastq.gz", 2736)
            check_order(
                output_root + "reads_1.fastq.gz",
                "729a5634addf363b4848f38ffc1bc9a9edbd10538b45669a7ac2719a6e2d61c2",
            )
            check_order(
                output_root + "reads_2.fastq.gz",
                "729a5634addf363b4848f38ffc1bc9a9edbd10538b45669a7ac2719a6e2d61c2",
            )
        case "tb_10k":
            check_length(output_root + "reads_1.fastq.gz", 9711)
            check_length(output_root + "reads_2.fastq.gz", 9711)
        case "tb_ont":
            check_length(output_root + "reads.fastq.gz", 998)


def test_empty_contigs(empty_sylph, test_outputs_dir, mocker):
    output_root = os.path.join(
        test_outputs_dir, "cm_sylph", empty_sylph["sample"] + "."
    )
    os.makedirs(os.path.join(test_outputs_dir, "cm_sylph"), exist_ok=True)

    seq_platform = "ont" if len(empty_sylph["reads"]) == 1 else "illumina"

    args = [
        "competitive_mapping",
        "--input_bam",
        empty_sylph["bam"],
        "--contigs",
        empty_sylph["contigs"],
        "--seq_platform",
        seq_platform,
        "--output_root",
        output_root,
        "--cpus",
        str(get_cpus()),
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
