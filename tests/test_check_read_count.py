from competitivemapping.check_read_count import check_threshold


def test_check_enough(check_read_count_files):
    assert check_threshold(check_read_count_files["enough"], 100000) == "true"
    assert check_threshold(check_read_count_files["not_enough"], 100000) == "false"
    assert check_threshold(check_read_count_files["no_mb_value"], 100000) == "false"
    assert (
        check_threshold(check_read_count_files["no_genome_name_key"], 100000) == "false"
    )
