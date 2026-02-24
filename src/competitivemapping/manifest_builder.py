# pylint: disable=too-many-locals
"""Build manifest and contigs dataframe from sylph report"""

import argparse
import dataclasses
import gzip
import logging
import multiprocessing
import os
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor

import networkx as nx
import pandas as pd
from Bio import SeqIO

multiprocessing.set_start_method("fork", force=True)

logging.basicConfig(
    format="%(asctime)s — %(relativeCreated)d — %(levelname)s — %(funcName)s:%(lineno)d — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S%z",
    level=logging.DEBUG,
)

REF_FILE_SUFFIXES = ["_genomic.fna.gz", ".fasta.gz"]
REF_FILE_SUFFIX_PATTERN = (
    "(" + "|".join(re.escape(suffix) for suffix in REF_FILE_SUFFIXES) + ")$"
)


@dataclasses.dataclass
class Config:
    """
    Config class for the module

    Args:
        cpus (int): Number of cores to use
        include_whole_genus (bool): Whether to include all genomes from genera found
        ani_threshold (float): ANI threshold for grouping genomes
        af_threshold (float): Alignment fraction threshold for grouping genomes
        force_ani_groups (bool): Whether to force use of ANI groups even if present in taxonomy file
        output_root (str): Path to the output root
        just_contigs (bool): Whether to skip manifest generation and only produce contigs.csv
    """

    cpus: int
    include_whole_genus: bool
    ani_threshold: float
    af_threshold: float
    force_ani_groups: bool
    output_root: str
    just_contigs: bool


def read_taxonomy_file(filepath: str) -> pd.DataFrame:
    """Read a taxonomy file and return a DataFrame"""
    # Note that pandas will automatically handle gzipped files
    sep = "\t" if filepath.endswith(".tsv") or filepath.endswith(".tsv.gz") else ","

    # Need to check if this is a headered file or a two column file such as GTDB/sylph provide by default

    first_row = pd.read_csv(filepath, sep=sep, nrows=1, header=None)

    # is "accession" in first row?
    if "accession" in first_row.values:
        usecols = ["accession", "taxonomy"]
        if "ani_group" in first_row.values:
            usecols.append("ani_group")
        return pd.read_csv(filepath, sep=sep, usecols=usecols)

    # Assume a headerless file
    return pd.read_csv(filepath, sep=sep, header=None, names=["accession", "taxonomy"])


def read_contigs(args: tuple[str, str]) -> list[dict[str, str]]:
    """Read contigs from a gzipped fasta file"""
    accession, filepath = args
    contigs = []
    with gzip.open(filepath, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contigs.append(
                {
                    "reference": accession,
                    "rname": record.id,
                    "length": len(record.seq),
                }
            )
    return contigs


def get_base_species_name(species: str) -> str:
    """removes _AB etc from species names if present.
    Should not effect genus names"""
    return re.sub(r"_[A-Z]+$", "", species)


def select_extra_species(
    sylph_species: list[str], potential_species_df: pd.DataFrame
) -> pd.DataFrame:
    """Filter extra species to include in manifest.
    Want to include one reference for each named species.
    So exlucde sp12345678 and only include one of <species>_A and <species>_B

    Args:
        sylph_species (list[str]): list of species found by sylph
        potential_species_df (pd.DataFrame): df with genomes from rest of genera.

    Returns:
        pd.DataFrame: Filtered dataframe
    """
    sylph_base_species = [get_base_species_name(species) for species in sylph_species]
    potential_species_df["base_species"] = potential_species_df["species"].apply(
        get_base_species_name
    )
    potential_species_df = potential_species_df[
        ~potential_species_df["base_species"].isin(sylph_base_species)
    ]

    # Remove species which have sp followed by 8 digits
    potential_species_df = potential_species_df[
        ~potential_species_df["species"].str.contains(r"sp\d{8}", na=False)
    ]

    # Now group by base_species and select the first alphabetically
    potential_species_df = (
        potential_species_df.sort_values("species")
        .groupby("base_species")
        .first()
        .reset_index()
    )

    return potential_species_df.copy()


def get_genome_paths(genome_dir: str) -> pd.DataFrame:
    """Produce df of genome paths using directory provided.
    Directory should contain a genome_paths.tsv file
    The paths in the table must be relative to the directory"""

    # File when downloaded actually seems to be space separated
    # using regex needs python engine, but file generally small so not a problem
    df = pd.read_csv(
        genome_dir + "/genome_paths.tsv",
        sep=r"\s",
        engine="python",
        header=None,
        names=["filename", "path"],
    )
    # use os.path.join to ensure correct path separator
    df["path"] = df.apply(
        lambda x: os.path.join(genome_dir, x["path"], x["filename"]), axis=1
    )
    return df


def get_ani_distances(
    selected_df: pd.DataFrame,
    config: Config,
) -> pd.DataFrame:
    """Get ANI distances between genomes using skani"""
    all_genomes = " ".join(selected_df["path"].tolist())
    command = f"skani triangle -E -t {config.cpus} --medium {all_genomes} > {config.output_root}ani_edge_list.tsv"
    subprocess.run(command, shell=True, check=True)
    ani_df = pd.read_csv(f"{config.output_root}ani_edge_list.tsv", sep="\t")

    # Want to replace path with accession
    path_to_accession = selected_df.set_index("path")["accession"].to_dict()
    ani_df["Ref"] = ani_df["Ref_file"].map(path_to_accession)
    ani_df["Query"] = ani_df["Query_file"].map(path_to_accession)

    return ani_df


def assign_ani_groups(contigs_df, ani_df, ani_threshold, af_threshold) -> pd.DataFrame:
    """Find groups of similar genomes based on ANI threshold.
    Returns copy of contigs_df with ani_group column added"""
    ani_df = ani_df[
        (ani_df["ANI"] >= ani_threshold)
        & (
            (ani_df["Align_fraction_ref"] >= af_threshold)
            | (ani_df["Align_fraction_query"] >= af_threshold)
        )
    ].copy()

    # build graph from filtered data
    graph = nx.from_pandas_edgelist(ani_df, "Ref", "Query")

    groups = {}
    group_index = 1
    for c in nx.connected_components(graph):
        for accession in c:
            groups[accession] = group_index
        group_index += 1

    new_df = contigs_df.copy()
    # for singletons ani_group should be none
    new_df["ani_group"] = new_df["reference"].map(groups, na_action="ignore")
    # set to Int64 to allow NaNs
    new_df["ani_group"] = new_df["ani_group"].astype("Int64")

    return new_df


def make_manifest(
    report_path: str,
    taxonomy_files: list[str],
    genome_dirs: list[str],
    config: Config,
) -> tuple[str, pd.DataFrame]:
    """Produce a multifasta manifest and contig df from a sylph report and genomes folder

    Args:
        report_path (str): Path to the sylph report
        taxonomy_files (list[str]): path to the db taxonomy files
        genome_dirs (list[str]): path to the directories with the genome fastas
        config (Config): Config object

    Returns:
        tuple[str, pd.DataFrame]: Path to the manifest file and a dataframe of contigs
    """
    logging.info("Creating manifest and reading contigs")
    manifest_file = f"{config.output_root}manifest.fasta.gz"

    # Check if the sylph report is empty
    if os.stat(report_path).st_size == 0 or pd.read_csv(report_path, sep="\t").empty:
        logging.warning("Sylph report is empty, producing empty outputs")
        with gzip.open(manifest_file, "wb") as _outfile:
            pass
        contigs_df = pd.DataFrame(
            columns=[
                "reference",
                "rname",
                "length",
                "totallength",
                "ani_group",
                "species",
            ]
        )
        return manifest_file, contigs_df

    sylph_df = pd.read_csv(report_path, sep="\t")
    sylph_df["accession"] = (
        sylph_df["Genome_file"]
        .str.split("/")
        .str[-1]
        .str.replace(REF_FILE_SUFFIX_PATTERN, "", regex=True)
    )
    sylph_accessions = sylph_df["accession"].tolist()

    genome_paths = pd.concat(get_genome_paths(genome_dir) for genome_dir in genome_dirs)
    genome_paths["accession"] = genome_paths["filename"].str.replace(
        REF_FILE_SUFFIX_PATTERN, "", regex=True
    )

    taxonomy_df = pd.concat([read_taxonomy_file(f) for f in taxonomy_files])

    # Can restrict to only representative genomes (those with a genome path)
    taxonomy_df = taxonomy_df[
        taxonomy_df["accession"].isin(genome_paths["accession"])
    ].copy()

    def select_taxa_level(taxonomy: str, key: str) -> str:
        parts = taxonomy.split(";")
        for taxa in parts:
            if taxa.startswith(key):
                return taxa
        return ""

    taxonomy_df["genus"] = taxonomy_df["taxonomy"].apply(
        lambda x: select_taxa_level(x, "g__").replace("g__", "")
    )
    taxonomy_df["species"] = taxonomy_df["taxonomy"].apply(
        lambda x: select_taxa_level(x, "s__").replace("s__", "")
    )

    if config.include_whole_genus:
        # Extend accessions to include genomes from rest of the genus(/genera)
        sylph_taxonomy_df = taxonomy_df[
            taxonomy_df["accession"].isin(sylph_accessions)
        ].copy()

        sylph_species = (
            taxonomy_df[taxonomy_df["accession"].isin(sylph_accessions)]["species"]
            .unique()
            .tolist()
        )

        found_genera = sylph_taxonomy_df["genus"].unique()

        potential_genomes = taxonomy_df[taxonomy_df["genus"].isin(found_genera)].copy()

        potential_genomes = select_extra_species(sylph_species, potential_genomes)

        # Now add these to the accessions
        accessions = sylph_accessions + potential_genomes["accession"].tolist()
    else:
        accessions = sylph_accessions

    # Now look up the genome paths
    selected_df = genome_paths[genome_paths["accession"].isin(accessions)]

    # Cat all genomes into a single file
    if not config.just_contigs:
        with open(manifest_file, "wb") as outfile:
            for filepath in selected_df["path"]:
                with open(filepath, "rb") as infile:
                    outfile.write(infile.read())

    # Read contigs in parallel
    with ProcessPoolExecutor(max_workers=config.cpus) as executor:
        results = list(
            executor.map(
                read_contigs, zip(selected_df["accession"], selected_df["path"])
            )
        )

    contigs_df = pd.DataFrame([contig for result in results for contig in result])
    contigs_df["totallength"] = contigs_df.groupby("reference")["length"].transform(
        "sum"
    )

    # add ani information if missing
    if "ani_group" in taxonomy_df.columns and not config.force_ani_groups:
        # If ani_group is already present, use it
        contigs_df = contigs_df.merge(
            taxonomy_df[["accession", "ani_group"]],
            left_on="reference",
            right_on="accession",
            how="left",
        )
        contigs_df.drop(columns=["accession"], inplace=True)
    elif config.ani_threshold > 0:
        ani_df = get_ani_distances(selected_df, config)
        contigs_df = assign_ani_groups(
            contigs_df, ani_df, config.ani_threshold, config.af_threshold
        )

    # add species information
    species_lookup = taxonomy_df.set_index("accession")["species"].to_dict()
    contigs_df["species"] = contigs_df["reference"].map(species_lookup)
    return manifest_file, contigs_df


def cli_entry_point():
    """Entry point for the CLI"""
    parser = argparse.ArgumentParser(description="Build manifest from sylph report.")
    parser.add_argument(
        "--sylph_report", required=True, help="Path to the sylph report TSV file"
    )
    parser.add_argument(
        "--taxonomy_files",
        required=True,
        help="Path to tsv files with taxonomy (and optionally ani) mapping. Can be gzipped",
        nargs="+",
    )
    # Need to provide parent directory to work with nextflow symlinks
    parser.add_argument(
        "--genome_dirs",
        required=True,
        help="Path to the directories containing the genome files."
        + " Must contain a genome_paths.tsv file like in gtdb_genomes_reps",
        nargs="+",
    )
    parser.add_argument(
        "--include_whole_genus",
        help="Include all genomes from genera found in the sylph report",
        action="store_true",
    )
    parser.add_argument(
        "--ani_threshold",
        help=(
            "ANI threshold for grouping genomes. 0 means do not group (default). "
            "Ignored if ani_group present in taxonomy file."
        ),
        default=0,
        type=float,
    )
    parser.add_argument(
        "--af_threshold",
        help=(
            "Alignment fraction threshold for grouping genomes."
            "Requires at least this in EITHER direction for a genome pair to be grouped."
        ),
        default=80,
        type=float,
    )
    parser.add_argument(
        "--force_ani_groups",
        help=(
            "Force use of ANI groups even if present in taxonomy file. "
            "Will recalculate ANI groups using provided thresholds."
        ),
        action="store_true",
    )
    parser.add_argument("--cpus", help="Number of CPUs to use", default=4, type=int)
    parser.add_argument("--output_root", required=True, help="Path to the output files")
    parser.add_argument(
        "--just_contigs",
        help="Only output contigs.csv and skip manifest generation",
        action="store_true",
    )

    args = parser.parse_args()

    sylph_report = args.sylph_report

    config = Config(
        cpus=int(args.cpus),
        include_whole_genus=args.include_whole_genus,
        ani_threshold=float(args.ani_threshold),
        af_threshold=float(args.af_threshold),
        force_ani_groups=args.force_ani_groups,
        output_root=args.output_root,
        just_contigs=args.just_contigs,
    )

    _manifest, contigs = make_manifest(
        sylph_report,
        args.taxonomy_files,
        args.genome_dirs,
        config,
    )

    contigs.to_csv(f"{args.output_root}contigs.csv", index=False)


if __name__ == "__main__":
    cli_entry_point()
