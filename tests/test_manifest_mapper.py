import os

from test_utils import check_bam_files

from competitivemapping import manifest_mapper


def test_if_manifest_empty(manifest, empty_sylph):
    """Test if the manifest is empty"""
    assert not manifest_mapper.is_manifest_empty(manifest)
    assert manifest_mapper.is_manifest_empty(empty_sylph["manifest"])


def test_manifest_mapper(samples, manifest, test_outputs_dir, mocker):
    output_dir = os.path.join(test_outputs_dir, "test_mapper")
    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, samples["sample"] + ".mapped.bam")

    seq_platform = "ont" if len(samples["reads"]) == 1 else "illumina"

    args = [
        "manifest_mapper",
        "--manifest",
        manifest,
        "--reads",
        " ".join(samples["reads"]),
        "--seq_platform",
        seq_platform,
        "--cpus",
        "40",
        "-o",
        outfile,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    manifest_mapper.cli_entry_point()

    check_bam_files(samples["bam"], outfile)


def test_empty_manifest(empty_sylph, test_outputs_dir, mocker):
    output_dir = os.path.join(test_outputs_dir, "test_mapper")
    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, "empty" + ".mapped.bam")
    # delete outfile if it exists
    if os.path.exists(outfile):
        os.remove(outfile)

    seq_platform = "ont" if len(empty_sylph["reads"]) == 1 else "illumina"

    args = [
        "manifest_mapper",
        "--manifest",
        empty_sylph["manifest"],
        "--reads",
        " ".join(empty_sylph["reads"]),
        "--seq_platform",
        seq_platform,
        "--cpus",
        "40",
        "-o",
        outfile,
    ]

    mocker.patch(
        "sys.argv",
        args,
    )

    manifest_mapper.cli_entry_point()
    # No output should be produced
    assert not os.path.exists(outfile)
