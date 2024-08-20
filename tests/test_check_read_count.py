import pytest
from competitivemapping import check_read_count


def test_check_enough_tb(enough_tb):
    num_reads = check_read_count.count_reads(enough_tb)
    assert check_read_count.check_threshold(num_reads, 100000)


def test_check_not_enough_tb(not_enough_tb):
    num_reads = check_read_count.count_reads(not_enough_tb)
    assert not check_read_count.check_threshold(num_reads, 100000)


def test_check_no_mb_value(no_mb_value):
    with pytest.raises(ValueError):
        check_read_count.count_reads(no_mb_value)


def test_check_no_genome_name_key(no_genome_name_key):
    with pytest.raises(ValueError):
        check_read_count.count_reads(no_genome_name_key)
