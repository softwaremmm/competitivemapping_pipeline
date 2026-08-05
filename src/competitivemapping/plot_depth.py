"""Produce depth plot for most common references in competitive mapping output."""

import warnings
import argparse
import pandas as pd
from plotnine import (
    ggplot,
    aes,
    geom_line,
    facet_wrap,
    theme_bw,
    labs,
    geom_vline,
    guides,
)


def bin_depths(df: pd.DataFrame, bin_size: int) -> pd.DataFrame:
    """Bin positions together and calculate mean depth for each bin."""
    df = df.copy()
    df["bin"] = (df["position"] // bin_size) * bin_size
    binned = (
        df.groupby(["contig", "bin"], as_index=False)
        .agg({"depth": "mean"})
        .rename(columns={"bin": "position"})
    )
    return binned


def smart_depth_sampling(
    df: pd.DataFrame, rate: int, fold_change: float, abs_change: int
) -> pd.DataFrame:
    """Sample depth values at a specified rate, but keep all positions where depth changes significantly."""
    sampled_points = []
    for contig, group in df.groupby("contig"):
        prev_depth = None
        prev_pos = None
        max_pos = max(group["position"])  # always include the last position
        for pos, depth in zip(group["position"], group["depth"]):
            if prev_depth is None:
                sampled_points.append((contig, pos, depth))
                prev_depth = depth
                prev_pos = pos
                continue

            if pos - prev_pos >= rate or pos == max_pos:
                sampled_points.append((contig, pos, depth))
                prev_depth = depth
                prev_pos = pos
                continue

            if depth == prev_depth:
                continue

            # always sample a change to or from zero depth
            if depth == 0 or prev_depth == 0:
                # also add a point before the change if it is not already sampled
                if pos - prev_pos > 1:
                    sampled_points.append((contig, pos - 1, prev_depth))
                sampled_points.append((contig, pos, depth))
                prev_depth = depth
                prev_pos = pos
                continue

            abs_diff = abs(depth - prev_depth)
            fold_diff = abs_diff / prev_depth
            if abs_diff >= abs_change and fold_diff >= fold_change:
                sampled_points.append((contig, pos, depth))
                prev_depth = depth
                prev_pos = pos

    return pd.DataFrame(sampled_points, columns=["contig", "position", "depth"])


def plot_depth(df: pd.DataFrame, contigs: pd.DataFrame, bin_size: int, outfile: str):
    """Produce multi-panel plot of depth per base for each reference."""

    # Calculate mean depth per contig and species to determine which species to plot
    contig_means = df.groupby("contig")["depth"].mean().reset_index()
    contig_means = contig_means.merge(
        contigs[["contig", "length", "species"]], on="contig", how="right"
    ).fillna(0)
    contig_means.to_csv("contig_means.csv", index=False)
    contig_means["total_bp"] = contig_means["length"] * contig_means["depth"]
    species_means = (
        contig_means.groupby("species")
        .agg({"total_bp": "sum", "length": "sum"})
        .reset_index()
    )
    species_means["mean_depth"] = species_means["total_bp"] / species_means["length"]
    species_means.to_csv("species_means.csv", index=False)
    species_means = species_means.sort_values(by="mean_depth", ascending=False)

    # only keep top species for plotting
    top_species = list(species_means.head(8)["species"])
    print(f"Top species for plotting: {top_species}")
    top_contigs = contigs[contigs["species"].isin(top_species)]["contig"].tolist()
    df = df[df["contig"].isin(top_contigs)].copy()

    # bin or smart sample the depth values for plotting
    if bin_size > 0:
        df = bin_depths(df, bin_size)
    else:
        df = smart_depth_sampling(df, rate=500, fold_change=0.1, abs_change=3)

    contigs = contigs[contigs["contig"].isin(top_contigs)].copy()
    contigs["species"] = pd.Categorical(
        contigs["species"], categories=top_species, ordered=True
    )

    # within a reference need to calculate cumulative length of contigs to plot them in order
    contigs["cumulative_length"] = contigs.groupby("species", observed=True)[
        "length"
    ].cumsum()
    contigs["contig_start"] = contigs["cumulative_length"] - contigs["length"]

    df = df.merge(
        contigs[["contig", "species", "contig_start"]],
        on="contig",
        how="left",
    )
    df["position"] = df["position"] + df["contig_start"]

    transitions = contigs[contigs["contig_start"] > 0][["species", "contig_start"]]

    p = (
        ggplot(df, aes(x="position", y="depth"))
        + geom_line(aes(group="contig"), size=0.2, color="#0015D1")
        + geom_vline(
            data=transitions,
            mapping=aes(xintercept="contig_start"),
            linetype="dashed",
            color="grey",
            alpha=0.5,
        )
        + facet_wrap(
            "species",
            scales="free",
            ncol=1,
        )
        + theme_bw()
        + labs(title="Read depth per contig", x="Position", y="Depth")
        + guides(color=False)
    )

    n_contigs = df["species"].nunique()
    height = max(6, n_contigs * 2)
    p.save(outfile, width=14, height=height, limitsize=False)


def cli_entry_point():
    """CLI entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Plot depth per base for each reference"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Depth TSV file with columns: contig, position, depth (like samtools depth output)",
    )
    parser.add_argument(
        "--contigs",
        type=str,
        required=True,
        help="Contigs csv relating contigs to reference with lengths",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output plot file",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=0,
        help="Bin size for plotting depth. Or use 0 for smart sampling (default 0)",
    )
    args = parser.parse_args()

    warnings.filterwarnings("ignore", message=".*ChainedAssignmentError.*")

    contigs = pd.read_csv(
        args.contigs, usecols=["reference", "rname", "length", "species"]
    ).rename(columns={"rname": "contig"})
    df = pd.read_csv(args.input, sep="\t", names=["contig", "position", "depth"])

    plot_depth(df, contigs, args.bin_size, args.output)
