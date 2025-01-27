# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments
# pylint: disable=too-many-locals
"""Build manifest and contigs dataframe from sylph report"""

import argparse
import gzip
import logging
import multiprocessing
import os
import re
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from Bio import SeqIO

multiprocessing.set_start_method("fork", force=True)

logging.basicConfig(
    format="%(asctime)s — %(relativeCreated)d — %(levelname)s — %(funcName)s:%(lineno)d — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S%z",
    level=logging.DEBUG,
)


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


def select_close_genera(genera_list: list[str]) -> set[str]:
    """Add extra genera to the list if certain genera are present. Like Escherichia and Shigella

    Args:
        genera_list (list[str]): List of genera found in the sylph report

    Returns:
        set[str]: Set of genera to include
    """
    new_list = genera_list.copy()
    for group in [["Escherichia", "Shigella"]]:
        if any(genus in genera_list for genus in group):
            new_list.extend(group)
    # remove duplicates
    return set(new_list)


def select_extra_species(
    sylph_species: list[str], potential_species_df: pd.DataFrame
) -> pd.DataFrame:
    """Filter extra species to include in manifest.
    Want to include one reference for each named species.
    So exlucde sp12345678 and only include one of <species>_A and <species>_B

    Args:
        sylph_species (list[str]): list of species found by sylph
        potential_species_df (pd.DataFrame): metadata df with genomes from rest of genera.

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


def make_manifest(
    report_path: str,
    metadata_files: list[str],
    genome_dirs: list[str],
    include_whole_genus: bool,
    output_root: str,
    cpus: int,
) -> tuple[str, pd.DataFrame]:
    """Produce a multifasta manifest and contig df from a sylph report and genomes folder

    Args:
        report_path (str): Path to the sylph report
        metadata_files (list[str]): path to the db metadata files, with taxonomy info
        genome_dirs (list[str]): path to the directories with the genome fastas
        include_whole_genus (bool): Whether to include all genomes from genera found
        output_root (str): Path to the output root
        cpus (int): number of cores to use

    Returns:
        tuple[str, pd.DataFrame]: Path to the manifest file and a dataframe of contigs
    """
    logging.info("Creating manifest and reading contigs")
    manifest_file = f"{output_root}manifest.fasta.gz"

    # Check if the sylph report is empty
    if os.stat(report_path).st_size == 0 or pd.read_csv(report_path, sep="\t").empty:
        logging.warning("Sylph report is empty, producing empty outputs")
        with gzip.open(manifest_file, "wb") as _outfile:
            pass
        contigs_df = pd.DataFrame(
            columns=["reference", "rname", "length", "totallength", "species"]
        )
        return manifest_file, contigs_df

    sylph_df = pd.read_csv(report_path, sep="\t")
    sylph_df["accession"] = (
        sylph_df["Genome_file"]
        .str.split("/")
        .str[-1]
        .str.replace("_genomic.fna.gz", "")
    )
    sylph_accessions = sylph_df["accession"].tolist()

    def get_genome_paths(dir_path):
        """Produce df of genome paths for given directory.
        Assumes a genomes_paths.tsv file which species relative paths to genomes"""
        # File when downloaded actually seems to be space separated
        # using regex needs python engine, but file generally small so not a problem
        df = pd.read_csv(
            dir_path + "/genome_paths.tsv",
            sep=r"\s",
            engine="python",
            header=None,
            names=["filename", "path"],
        )
        df["path"] = df["path"].apply(lambda x: os.path.join(dir_path, x))
        df["path"] = df["path"] + "/" + df["filename"]
        return df

    genome_paths = pd.concat(get_genome_paths(g_dir) for g_dir in genome_dirs)
    genome_paths["accession"] = genome_paths["filename"].str.replace(
        "_genomic.fna.gz", ""
    )

    metadata_df = pd.concat(
        [
            pd.read_csv(f, sep="\t", header=None, names=["accession", "taxonomy"])
            for f in metadata_files
        ]
    )

    # Can restrict to only representative genomes (those with a genome path)
    metadata_df = metadata_df[
        metadata_df["accession"].isin(genome_paths["accession"])
    ].copy()

    def select_taxa_level(taxonomy: str, key: str) -> str:
        parts = taxonomy.split(";")
        for taxa in parts:
            if taxa.startswith(key):
                return taxa
        return ""

    metadata_df["genus"] = metadata_df["taxonomy"].apply(
        lambda x: select_taxa_level(x, "g__").replace("g__", "")
    )
    metadata_df["species"] = metadata_df["taxonomy"].apply(
        lambda x: select_taxa_level(x, "s__").replace("s__", "")
    )

    if include_whole_genus:
        # Extend accessions to include genomes from rest of the genus(/genera)
        sylph_metadata_df = metadata_df[
            metadata_df["accession"].isin(sylph_accessions)
        ].copy()

        sylph_species = sylph_metadata_df["species"].unique().tolist()

        found_genera = select_close_genera(sylph_metadata_df["genus"].unique().tolist())

        potential_genomes = metadata_df[metadata_df["genus"].isin(found_genera)].copy()

        potential_genomes = select_extra_species(sylph_species, potential_genomes)

        # Now add these to the accessions
        accessions = sylph_accessions + potential_genomes["accession"].tolist()
    else:
        accessions = sylph_accessions

    # Now look up the genome paths
    selected_df = genome_paths[genome_paths["accession"].isin(accessions)]

    # Cat all genomes into a single file
    with open(manifest_file, "wb") as outfile:
        for filepath in selected_df["path"]:
            with open(filepath, "rb") as infile:
                outfile.write(infile.read())

    # Read contigs in parallel
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        results = list(
            executor.map(
                read_contigs, zip(selected_df["accession"], selected_df["path"])
            )
        )

    contigs_df = pd.DataFrame([contig for result in results for contig in result])
    contigs_df["totallength"] = contigs_df.groupby("reference")["length"].transform(
        "sum"
    )

    # add species information from metadata
    species_lookup = metadata_df.set_index("accession")["species"].to_dict()
    contigs_df["species"] = contigs_df["reference"].map(species_lookup)
    return manifest_file, contigs_df


def cli_entry_point():
    """Entry point for the CLI"""
    parser = argparse.ArgumentParser(description="Build manifest from sylph report.")
    parser.add_argument(
        "--sylph_report", required=True, help="Path to the sylph report TSV file"
    )
    parser.add_argument(
        "--metadata_files",
        required=True,
        help="Path to the metadata files with assembly to taxonomy mapping",
        nargs="+",
    )
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
    parser.add_argument("--cpus", help="Number of CPUs to use", default=4, type=int)
    parser.add_argument("--output_root", required=True, help="Path to the output files")

    args = parser.parse_args()

    sylph_report = args.sylph_report

    _manifest, contigs = make_manifest(
        sylph_report,
        args.metadata_files,
        args.genome_dirs,
        args.include_whole_genus,
        args.output_root,
        args.cpus,
    )

    contigs.to_csv(f"{args.output_root}contigs.csv", index=False)


if __name__ == "__main__":
    cli_entry_point()
