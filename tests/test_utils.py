import filecmp

REPLACE_EXPECTATION = False


def check_file(expectation, result):
    if not filecmp.cmp(expectation, result):
        new_file_path = f"{expectation}.new"
        if REPLACE_EXPECTATION:
            new_file_path = expectation
        with open(new_file_path, "w", encoding="utf-8") as new_file, open(
            result, "r", encoding="utf-8"
        ) as result_file:
            new_file.write(result_file.read())

    assert filecmp.cmp(expectation, result)
