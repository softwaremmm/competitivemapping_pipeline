import os

from test_utils import check_file

from competitivemapping.extract_reads import extract_reads


def test_extract_reads(extract_read_files, test_outputs_dir):
    os.makedirs(test_outputs_dir / "extract_reads", exist_ok=True)

    extract_reads(
        input_fastqs=extract_read_files["input"],
        alns_csv=extract_read_files["alns"],
        contigs_csv=extract_read_files["contigs"],
        reference=extract_read_files["reference"],
        require_both_in_pair=False,
        output_root=f"{test_outputs_dir}/extract_reads/output",
    )

    output_fastqs = [
        f"{test_outputs_dir}/extract_reads/output_{index}.fastq"
        for index in range(1, 3)
    ]

    expectation = extract_read_files["expectation"]
    for output_fq, expected_fq in zip(output_fastqs, expectation):
        check_file(output_fq, expected_fq)
