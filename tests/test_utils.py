import filecmp
import os

REPLACE_EXPECTATION = False


def check_file(expectation, result):
    """Check if files match, if not replace/save replacement."""

    # check if expectation missing
    if not os.path.exists(expectation):
        with open(expectation, "w", encoding="utf-8") as file:
            file.write("")

    if not filecmp.cmp(expectation, result):
        new_file_path = f"{expectation}.new"
        if REPLACE_EXPECTATION:
            new_file_path = expectation
        with (
            open(new_file_path, "w", encoding="utf-8") as new_file,
            open(result, "r", encoding="utf-8") as result_file,
        ):
            new_file.write(result_file.read())

    assert filecmp.cmp(expectation, result)
