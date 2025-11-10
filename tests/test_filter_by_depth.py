import os

from test_utils import check_file

from competitivemapping.filter_by_depth import filter_by_depth


def test_filter_by_depth(filter_by_depth_files, test_outputs_dir):
    os.makedirs(test_outputs_dir / "filter_by_depth", exist_ok=True)

    filter_by_depth(
        tie_break_report=filter_by_depth_files["tie_break_report"],
        species_list=filter_by_depth_files["species_list"],
        min_depth=filter_by_depth_files["min_depth"],
        output_root=f"{test_outputs_dir}/filter_by_depth/high_depth",
    )

    output_file = f"{test_outputs_dir}/filter_by_depth/high_depth_refs.txt"

    expectation = filter_by_depth_files["expectation"]
    check_file(output_file, expectation)
