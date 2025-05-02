import filecmp
import gzip
import hashlib
import os
import shutil
import subprocess

REPLACE_EXPECTATION = False


def replace(expectation, result):
    """Replace the expectation file with the result file."""
    if REPLACE_EXPECTATION:
        shutil.copy(result, expectation)
    else:
        new_file_path = f"{expectation}.new"
        shutil.copy(result, new_file_path)


def check_file(expectation, result):
    """Check if files match, if not replace expectation or save new result for inspection."""

    # check if expectation missing
    if not os.path.exists(expectation):
        # Then just create an empty file
        with open(expectation, "w", encoding="utf-8") as file:
            file.write("")

    if not filecmp.cmp(expectation, result):
        replace(expectation, result)

    assert filecmp.cmp(expectation, result)


def md5(file_path):
    """Calculate the md5 hash of a file."""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def check_file_by_md5(expectation, result):
    """Check if files match by md5, if not replace expectation or save new result for inspection."""

    # check if expectation missing
    if not os.path.exists(expectation):
        with open(expectation, "w", encoding="utf-8") as file:
            file.write("")

    matches = md5(expectation) == md5(result)

    if not matches:
        replace(expectation, result)
    assert matches


def check_gzipped_file(expectation, result):
    """Check if gzipped files match, if not replace expectation or save new result for inspection."""

    # check if expectation missing
    if not os.path.exists(expectation):
        replace(expectation, result)
        assert False

    equal = True
    with gzip.open(expectation, "rb") as f1, gzip.open(result, "rb") as f2:
        while True:
            b1 = f1.read(4096)
            b2 = f2.read(4096)
            if b1 != b2:
                equal = False
                break
            if not b1:
                break

    if not equal:
        replace(expectation, result)
    assert equal


def get_bam_md5(file_path):
    """Get the md5 hash of a bam file excluding header."""
    # convert filepath to absolute path
    file_path = os.path.abspath(file_path)
    command = f"samtools view {file_path} | md5sum | cut -f1 -d' '"
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        check=True,
    )
    if result.returncode != 0:
        raise ValueError(f"Error getting md5 for {file_path}: {result.stderr}")
    return result.stdout.strip()


def check_bam_files(expectation, result):
    """Check if bam files match, if not replace expectation or save new result for inspection."""

    # check if expectation missing
    if not os.path.exists(expectation):
        replace(expectation, result)
        assert False

    # Check md5 of the bam files
    md5_expectation = get_bam_md5(expectation)
    md5_result = get_bam_md5(result)

    if md5_expectation != md5_result:
        replace(expectation, result)
        assert False
    assert md5_expectation == md5_result


def get_cpus():
    """Get the number of CPUs available."""
    if os.cpu_count():
        return os.cpu_count()
    return 4
