import pandas as pd
import competitivemapping.process_mapping as process_mapping


def test_join_references(coverage_table, species_table, expected_joined):
    actual_joined = process_mapping.join_references(coverage_table, species_table)
    pd.testing.assert_frame_equal(actual_joined, expected_joined)
