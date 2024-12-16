import os

import pandas as pd

from competitivemapping import contig_mapping


def test_get_megahit_contig_stats(chloro_10k):
    df = contig_mapping.get_megahit_contig_stats([chloro_10k["contigs"]])
    pd.testing.assert_frame_equal(df, pd.read_csv(chloro_10k["contig_stats"]))


def test_competitive_map_contigs(
    chloro_10k, species_table_path, manifest, test_outputs_dir
):
    output_root = os.path.join(
        test_outputs_dir, "cm_contigs", chloro_10k["sample"] + "."
    )
    os.makedirs(os.path.join(test_outputs_dir, "cm_contigs"), exist_ok=True)

    df = contig_mapping.competitive_map_contigs(
        manifest,
        pd.read_csv(species_table_path),
        [chloro_10k["contigs"]],
        4,
        output_root,
    )
    df.to_csv(output_root + "contig_species_comparison.csv", index=False)

    pd.testing.assert_frame_equal(
        df.reset_index(drop=True), pd.read_csv(chloro_10k["contig_species_comparison"])
    )
