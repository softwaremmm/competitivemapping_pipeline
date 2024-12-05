# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments
# pylint: disable=too-many-locals
"""Run Competitive Mapping and aggregate the results"""

import argparse
import gzip
import json
import logging
import multiprocessing
import os
import subprocess
import typing
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from Bio import SeqIO

from competitivemapping.process_aln_stats import get_alignment_stats
from competitivemapping.process_coverage import process_coverage

multiprocessing.set_start_method("fork", force=True)

logging.basicConfig(
    format="%(asctime)s — %(relativeCreated)d — %(levelname)s — %(funcName)s:%(lineno)d — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S%z",
    level=logging.DEBUG,
)


FINAL_COLUMNS = [
    "genome_name",
    "length",
    "numreads",
    "coverage",
    "meandepth",
    "coverage_including_secondary",
    "meandepth_including_secondary",
    "total_reads",
    "total_alns",
    "exclusive_reads",
    "primary_reads",
    "secondary_reads",
    "supplementary_reads",
    "supplementary_alns",
]

COLUMN_EXPLANATIONS = {
    "length": "Length of all contigs in reference",
    "numreads": "Number of reads which map well to this reference (primary + supplementary)",
    "coverage": "Percentage of reference covered by reads (not including secondary alignments)",
    "meandepth": "Mean depth of coverage of reference (not including secondary alignments)",
    "coverage_including_secondary": "Percentage of reference covered by reads (including secondary alignments)",
    "meandepth_including_secondary": "Mean depth of coverage of reference (including secondary alignments)",
    "total_reads": "Total number of reads which map (in any way) to chrom",
    "total_alns": "Total number of alignments, so will count supplementary separately",
    "exclusive_reads": "reads which only maps to this chrom",
    "primary_reads": "reads which map best to this chrom",
    "secondary_reads": "reads which map to this chrom but not with primary",
    "supplementary_reads": "reads which have supplementary alignments but not primary",
    "supplementary_alns": "number of supplementary alignments",
}


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
    """removes _AB etc from species names if present"""
    if "_" in species:
        return species.split("_")[0]
    return species


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

        sylph_species = (
            metadata_df[metadata_df["accession"].isin(sylph_accessions)]["species"]
            .unique()
            .tolist()
        )

        found_genera = sylph_metadata_df["genus"].unique()

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


def map_reads(
    manifest: str, reads: list[str], seq_platform: str, cpus: int, output_root: str
) -> str:
    """Run minimap2 to map reads against manifest, returning sorted bam.

    Args:
        manifest (str): Path to the manifest file
        reads (list[str]): List of paths to the read fastqs
        seq_platform (str): Sequencing platform
        cpus (int): number of cores to use
        output_root (str): Path to the output root

    Returns:
        str: path to the alignment bam file
    """
    logging.info("Mapping reads")
    aln_bam = f"{output_root}alignment.bam"

    command = f"minimap2 -t {cpus} --secondary yes -N 1000"
    if seq_platform == "ont":
        command += " -ax map-ont"
    else:
        command += " -ax sr"
    command += f" {manifest} {' '.join(reads)}"

    command += f"| samtools sort -@ {cpus} -o {aln_bam}"

    subprocess.run(
        command,
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    return aln_bam


def get_aln_stats(aln_bam: str, contigs_df: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    """Get stats from bam files using pysam.

    Args:
        aln_bam (str): path to the alignment bam file
        contigs_df (pd.DataFrame): dataframe of contigs

    Returns:
        tuple[dict, pd.DataFrame]: dict with overall stats and
            dataframe with alignment stats per reference
    """
    logging.info("Getting alignment stats")
    name_mapping = contigs_df.set_index("rname")["reference"].to_dict()
    overall_stats, aln_stats_df = get_alignment_stats(aln_bam, name_mapping)

    return overall_stats, aln_stats_df


def run_samtools_coverage(args: tuple[str, str, str]) -> str:
    """Run samtools coverage on a bam file

    Args:
        args (tuple[str, str, str]): (bam file, flags, output file)

    Returns:
        str: output file
    """
    aln_bam, flags, output_file = args
    subprocess.run(
        f"samtools coverage {aln_bam} {flags} -o {output_file}",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )
    return output_file


def get_coverage_stats(
    aln_bam: str, contigs_df: pd.DataFrame, cpus: int, output_root: str
) -> pd.DataFrame:
    """Get coverage stats using samtools coverage

    Args:
        aln_bam (str): path to the alignment bam file
        contigs_df (pd.DataFrame): dataframe of contigs
        cpus (int): number of cores to use (uses 2 max)
        output_root (str): Path to the output root

    Returns:
        pd.DataFrame: summary of coverage metrics for each reference
    """
    logging.info("Getting coverage stats")

    primary_coverage = f"{output_root}coverage_primary.tsv"
    full_coverage = f"{output_root}coverage_full.tsv"
    args = [
        (aln_bam, "", primary_coverage),
        (aln_bam, "--excl-flags 1540", full_coverage),
    ]

    with ProcessPoolExecutor(max_workers=cpus) as executor:
        _results = list(executor.map(run_samtools_coverage, args))

    coverage_df = process_coverage(primary_coverage, full_coverage, contigs_df)

    return coverage_df


def output_fastqs(
    aln_bam: str,
    reference: str,
    contigs_df: pd.DataFrame,
    include_unmapped: bool,
    seq_platform: str,
    cpus: int,
    output_root: str,
):
    """Extract reads from BAM file and output to FASTQ

    Args:
        aln_bam (str): path to the alignment bam file
        reference (str): desired reference to extract reads for
        contigs_df (pd.DataFrame): dataframe of contigs
        include_unmapped (bool): whether to include unmapped reads
        seq_platform (str): sequencing platform
        cpus (int): number of cores to use
        output_root (str): Path to the output root
    """
    logging.info("Outputting FASTQs")
    rnames = contigs_df[contigs_df["reference"] == reference]["rname"].tolist()
    if include_unmapped:
        rnames.append('"*"')

    sorted_ref_bam = f"{output_root}output_aln.bam"

    subprocess.run(
        f"samtools index {aln_bam}",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )
    # In future could add flag -P to always include read pairs
    # But note that this fails to fetch the pair for supplementary alignments
    subprocess.run(
        f"samtools view {aln_bam} -u {' '.join(rnames)} | samtools sort -@ {cpus} -o {sorted_ref_bam}",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )

    # Convert BAM output to FASTQ
    # Default excl-flag is 0x900 which is secondary (0x100) and supplementary (0x800)
    # So we need to exclude secondary alignments (0x100) only
    # Note that singletons will be discarded;
    # this is where only one of the read pair is in the bam or passes the flags

    # So for paired reads, only output if both are unmapped, or both map
    # (maybe with supplementary) to the reference
    command = f"samtools fastq --excl-flags 0x100 -@ {cpus} "
    if seq_platform == "ont":
        command += f"-0 {output_root}reads.fastq.gz {sorted_ref_bam}"
    else:
        command += (
            f"-1 {output_root}reads_1.fastq.gz -2 {output_root}reads_2.fastq.gz"
            f" -0 /dev/null -s /dev/null {sorted_ref_bam}"
        )
    subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE)


def produce_empty_outputs(output_root: str):
    """Create empty output files for species comparison json/csv"""
    report = {
        "total_read_counts": {
            "mapped_reads": 0,
            "unmapped_reads": 0,
        },
        "references": [],
        "definitions": COLUMN_EXPLANATIONS,
    }
    with open(f"{output_root}species_comparison.json", "w", encoding="utf-8") as file:
        json.dump(report, file, indent=4)

    with open(f"{output_root}species_comparison.csv", "w", encoding="utf-8") as file:
        file.write(",".join(FINAL_COLUMNS) + "\n")


def run_competitive_mapping(
    manifest: str,
    contigs: pd.DataFrame,
    reads: list[str],
    ref_for_fastq: str | None,
    seq_platform: str,
    cpus: int,
    output_root: str,
):
    """Map reads to a manifest and produce a comparison of coverage and alignment stats

    Args:
        manifest (str): Path to the manifest file
        contigs (pd.DataFrame): Dataframe of contigs
        reads (list[str]): List of paths to the read fastqs
        ref_for_fastq (str | None): Reference to extract reads for
        seq_platform (str): Sequencing platform
        cpus (int): Number of cores to use
        output_root (str): Path to the output root
    """
    logging.info("Running competitive mapping")
    aln_bam = map_reads(manifest, reads, seq_platform, cpus, output_root)

    coverage_df = get_coverage_stats(aln_bam, contigs, cpus, output_root)

    overall_stats, aln_stats = get_aln_stats(aln_bam, contigs)
    df = pd.merge(coverage_df, aln_stats, on="genome_name")

    logging.info("Writing output files")

    # Can update numreads to actually reflect reads
    # As samtools coverage actually counts alignments
    df["numreads"] = df["primary_reads"] + df["supplementary_reads"]

    for col in (
        "coverage",
        "meandepth",
        "coverage_including_secondary",
        "meandepth_including_secondary",
    ):
        if col in df.columns:
            df[col] = df[col].apply(lambda x: float(f"{x:.4g}"))

    df = df[FINAL_COLUMNS].copy()

    # add final row with unmapped read count
    new_row: dict[str, typing.Any] = {col: 0 for col in FINAL_COLUMNS}
    new_row["genome_name"] = "unmapped"
    for col in "numreads", "primary_reads", "total_reads":
        new_row[col] = overall_stats["unmapped_reads"]
    df.loc[-1] = new_row

    # incorporate species information if available
    if "species" in contigs.columns:
        species_lookup = contigs.set_index("reference")["species"].to_dict()
        species_lookup["unmapped"] = "unmapped"
        df["species"] = df["genome_name"].map(species_lookup)
        df = df[["species"] + FINAL_COLUMNS]

    df.to_csv(f"{output_root}species_comparison.csv", index=False)
    # remove the last row with unmapped read count
    df = df.iloc[:-1]

    output = {
        "total_read_counts": overall_stats,
        "references": df.to_dict(orient="records"),
        "definitions": COLUMN_EXPLANATIONS,
    }
    with open(f"{output_root}species_comparison.json", "w", encoding="utf-8") as file:
        json.dump(output, file, indent=4)

    if ref_for_fastq:
        output_fastqs(
            aln_bam,
            ref_for_fastq,
            contigs,
            include_unmapped=True,
            seq_platform=seq_platform,
            cpus=cpus,
            output_root=output_root,
        )

    logging.info("Finished competitive mapping")


def run_dynamic_competitive_mapping(
    sylph_report: str,
    metadata_files: list[str],
    genome_dirs: list[str],
    include_whole_genus: bool,
    reads: list[str],
    ref_for_fastq: str | None,
    seq_platform: str,
    cpus: int,
    output_root: str,
):
    """Create manifest from sylph report then competitive mapping

    Args:
        sylph_report (str): path to the sylph report
        metadata_files (list[str]): path to the db metadata files, with taxonomy info
        genome_dirs (list[str]): path to the directories with the genome fastas
        include_whole_genus (bool): whether to include all genomes from genera found
        reads (list[str]): list of paths to the read fastqs
        ref_for_fastq (str | None): reference to extract reads for
        seq_platform (str): sequencing platform
        cpus (int): number of cores to use
        output_root (str): path to the output root
    """

    # check if sylph_report is empty
    if os.stat(sylph_report).st_size == 0 or pd.read_csv(sylph_report, sep="\t").empty:
        logging.warning("Sylph report is empty")
        produce_empty_outputs(output_root)
        return

    manifest, contigs = make_manifest(
        sylph_report,
        metadata_files,
        genome_dirs,
        include_whole_genus,
        output_root,
        cpus,
    )

    run_competitive_mapping(
        manifest, contigs, reads, ref_for_fastq, seq_platform, cpus, output_root
    )


def cli_entry_point():
    """Entry point for the CLI"""
    parser = argparse.ArgumentParser(description="Process sylph report and genomes.")

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "--reads", required=True, help="Path to the reads", nargs="+"
    )
    common_parser.add_argument(
        "--seq_platform", help="Sequencing platform", default="illumina"
    )
    common_parser.add_argument(
        "--cpus", help="Number of CPUs to use", default=4, type=int
    )
    common_parser.add_argument(
        "--ref_for_fastq", help="Reference to extract reads for", default=None
    )
    common_parser.add_argument(
        "--output_root", required=True, help="Path to the output files"
    )

    # add subcommand "with_manifest" to run_dynamic_competitive_mapping
    subparsers = parser.add_subparsers(dest="subcommand", required=True)
    manifest_parser = subparsers.add_parser(
        "manifest",
        parents=[common_parser],
        help="Run competitive mapping with set manifest and contigs list",
    )
    manifest_parser.add_argument(
        "--manifest", required=True, help="Path to the manifest file"
    )
    manifest_parser.add_argument(
        "--contigs", required=True, help="Path to the contigs file"
    )

    sylph_parser = subparsers.add_parser(
        "sylph",
        parents=[common_parser],
        help="Run competitive mapping dynamically based on sylph report",
    )
    sylph_parser.add_argument(
        "--sylph_report", required=True, help="Path to the sylph report TSV file"
    )
    sylph_parser.add_argument(
        "--metadata_files",
        required=True,
        help="Path to the metadata files with assembly to taxonomy mapping",
        nargs="+",
    )
    sylph_parser.add_argument(
        "--genome_dirs",
        required=True,
        help="Path to the directories containing the genome files."
        + " Must contain a genome_paths.tsv file like in gtdb_genomes_reps",
        nargs="+",
    )
    sylph_parser.add_argument(
        "--include_whole_genus",
        help="Include all genomes from genera found in the sylph report",
        action="store_true",
    )

    args = parser.parse_args()

    if args.subcommand == "manifest":
        run_competitive_mapping(
            args.manifest,
            pd.read_csv(args.contigs),
            args.reads,
            args.ref_for_fastq,
            args.seq_platform,
            args.cpus,
            args.output_root,
        )
    elif args.subcommand == "sylph":
        run_dynamic_competitive_mapping(
            args.sylph_report,
            args.metadata_files,
            args.genome_dirs,
            args.include_whole_genus,
            args.reads,
            args.ref_for_fastq,
            args.seq_platform,
            args.cpus,
            args.output_root,
        )


if __name__ == "__main__":
    cli_entry_point()
