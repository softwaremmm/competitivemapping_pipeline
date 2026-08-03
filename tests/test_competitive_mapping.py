import hashlib
import os
import subprocess
import tempfile
from gzip import open as gzopen
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
from Bio import SeqIO
from test_utils import check_file, get_cpus

from competitivemapping import competitive_mapping
from competitivemapping.competitive_mapping import Config, output_fastqs


def test_cli_entry_point(samples, myco_manifest, test_outputs_dir, mocker):
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
        myco_manifest["contigs"],
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
            check_length(output_root + "reads_1.fastq.gz", 13)
            check_length(output_root + "reads_2.fastq.gz", 13)
            check_order(
                output_root + "reads_1.fastq.gz",
                "5a36937177ce534d20988350a583c445b3dca532f8c359c285da324f95f75346",
            )
            check_order(
                output_root + "reads_2.fastq.gz",
                "5a36937177ce534d20988350a583c445b3dca532f8c359c285da324f95f75346",
            )
        case "tb_10k":
            check_length(output_root + "reads_1.fastq.gz", 5235)
            check_length(output_root + "reads_2.fastq.gz", 5235)
        case "tb_ont":
            check_length(output_root + "reads.fastq.gz", 438)


def test_empty_contigs(empty_files, test_outputs_dir, mocker):
    output_root = os.path.join(test_outputs_dir, "cm", empty_files["sample"] + ".")
    os.makedirs(os.path.join(test_outputs_dir, "cm"), exist_ok=True)

    seq_platform = "ont" if len(empty_files["reads"]) == 1 else "illumina"

    args = [
        "competitive_mapping",
        "--input_bam",
        "path_to_nothing.bam",
        "--contigs",
        empty_files["contigs"],
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
        empty_files["species_comparison"], output_root + "species_comparison.json"
    )
    check_file(empty_files["csv_comparison"], output_root + "species_comparison.csv")


## Test cases for output_fastqs function


@pytest.fixture
def mock_contigs_df():
    """Create a mock contigs dataframe for testing"""
    return pd.DataFrame(
        {"reference": ["ref1", "ref2", "ref3"], "rname": ["r1", "r2", "r3"]}
    )


@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for test outputs"""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@patch("subprocess.run")
def test_output_fastqs_single_reference(mock_run, mock_contigs_df, temp_output_dir):
    """Test output_fastqs with a single reference"""
    # Setup
    mock_run.return_value = MagicMock()
    aln_bam = "test.bam"
    references = ["ref1"]

    config = Config(
        cpus=1,
        seq_platform="illumina",
        output_root=temp_output_dir,
    )

    # Execute
    output_fastqs(
        aln_bam=aln_bam,
        references=references,
        contigs_df=mock_contigs_df,
        include_unmapped=False,
        config=config,
    )

    # Verify
    # Check that samtools index was called
    mock_run.assert_any_call(
        f"samtools index {aln_bam}", shell=True, check=True, stdout=subprocess.PIPE
    )
    # Check that the correct rname was used
    mock_run.assert_any_call(
        f"samtools view -h {aln_bam} -u r1 | samtools sort -n -@ 1 -o {temp_output_dir}output_aln.bam",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )


@patch("subprocess.run")
def test_output_fastqs_multiple_references(mock_run, mock_contigs_df, temp_output_dir):
    """Test output_fastqs with multiple references"""
    # Setup
    mock_run.return_value = MagicMock()
    aln_bam = "test.bam"
    references = ["ref1", "ref2"]

    config = Config(
        cpus=1,
        seq_platform="illumina",
        output_root=temp_output_dir,
    )

    # Execute
    output_fastqs(
        aln_bam=aln_bam,
        references=references,
        contigs_df=mock_contigs_df,
        include_unmapped=False,
        config=config,
    )

    # Verify
    # Check that both rnames were used
    mock_run.assert_any_call(
        f"samtools view -h {aln_bam} -u r1 r2 | samtools sort -n -@ 1 -o {temp_output_dir}output_aln.bam",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )


@patch("subprocess.run")
def test_output_fastqs_with_unmapped(mock_run, mock_contigs_df, temp_output_dir):
    """Test output_fastqs with unmapped reads included"""
    # Setup
    mock_run.return_value = MagicMock()
    aln_bam = "test.bam"
    references = ["ref1"]

    config = Config(
        cpus=1,
        seq_platform="illumina",
        output_root=temp_output_dir,
    )

    # Execute
    output_fastqs(
        aln_bam=aln_bam,
        references=references,
        contigs_df=mock_contigs_df,
        include_unmapped=True,
        config=config,
    )

    # Verify
    # Check that unmapped reads were included
    mock_run.assert_any_call(
        f'samtools view -h {aln_bam} -u r1 "*" | samtools sort -n -@ 1 -o {temp_output_dir}output_aln.bam',
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )


@patch("subprocess.run")
def test_output_fastqs_nonexistent_reference(
    mock_run, mock_contigs_df, temp_output_dir
):
    """Test output_fastqs with references that don't exist in contigs"""
    # Setup
    mock_run.return_value = MagicMock()
    aln_bam = "test.bam"
    references = ["nonexistent_ref"]

    config = Config(
        cpus=1,
        seq_platform="illumina",
        output_root=temp_output_dir,
    )

    # Execute
    output_fastqs(
        aln_bam=aln_bam,
        references=references,
        contigs_df=mock_contigs_df,
        include_unmapped=False,
        config=config,
    )

    # Verify
    # Check that no rnames were used (empty list)
    assert not any("samtools view" in str(call) for call in mock_run.call_args_list)


@patch("subprocess.run")
def test_output_fastqs_ont_platform(mock_run, mock_contigs_df, temp_output_dir):
    """Test output_fastqs with ONT platform"""
    # Setup
    mock_run.return_value = MagicMock()
    aln_bam = "test.bam"
    references = ["ref1"]

    config = Config(
        cpus=1,
        seq_platform="ont",
        output_root=temp_output_dir,
    )

    # Execute
    output_fastqs(
        aln_bam=aln_bam,
        references=references,
        contigs_df=mock_contigs_df,
        include_unmapped=False,
        config=config,
    )

    # Verify
    # Check that ONT-specific command was used
    mock_run.assert_any_call(
        f"samtools fastq --excl-flags 0x100 -@ 1 -0 {temp_output_dir}reads.fastq.gz {temp_output_dir}output_aln.bam",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )
