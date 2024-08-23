from competitivemapping.check_read_count import check_threshold


def test_check_enough_tb(enough_tb):
    assert check_threshold(enough_tb, 100000) == "true"


def test_check_not_enough_tb(not_enough_tb):
    assert check_threshold(not_enough_tb, 100000) == "false"


def test_check_no_mb_value(no_mb_value):
    assert check_threshold(no_mb_value, 100000) == "false"


def test_check_no_genome_name_key(no_genome_name_key):
    assert check_threshold(no_genome_name_key, 100000) == "false"
