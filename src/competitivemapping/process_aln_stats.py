"""Script to process a bam file and give summary of
number of reads and alignmments for each reference"""

import typing
import argparse
from pysam import AlignmentFile  # pylint: disable = no-name-in-module
import pandas as pd


def get_name_mapping(names_file: str) -> dict:
    """Produce dictionary of reference code to their human names

    Args:
        names_file (str): path to species list file

    Returns:
        dict: dict with mapping for reference code to human names
    """

    df = pd.read_csv(names_file).set_index("rname")
    name_mapping = df["reference"].to_dict()
    return name_mapping


def get_alignment_stats(
    bam_file,
    name_mapping,
    exclude_secondary=False,
    exclude_supplementary=False,
    exclude_mq0=False,
) -> pd.DataFrame:
    """Iterate through all alignments in the bam file and produce a summary by reference

    Args:
        bam_file (_type_): bam file to process
        name_mapping (_type_): dict with mapping for reference code to human names
        exclude_secondary (bool, optional): exclude secondary alignments. Defaults to False.
        exclude_supplementary (bool, optional): exclude supplementary alignments. Defaults to False.
        exclude_mq0 (bool, optional): exclude alignments with map quality 0. Defaults to False.

    Raises:
        ValueError: if a read is both secondary and supplementary
        ValueError: if a read has multiple primary alignments

    Returns:
        pd.DataFrame: table with summary of alignments by reference
    """
    with AlignmentFile(bam_file, "rb") as bam:  # ignore: no-member
        # This is a mapping from the chomosome id to the human readable name
        # This will also have the effect of combining contigs of multi-chromosome references
        chrom_id_to_name = {}
        # read_info = {query: {"primary": "", "secondary": [], "supplementary": []}}
        read_info: dict[str, dict[str, typing.Any]] = {}
        chroms = set()

        for read in bam:
            if read.is_unmapped or read.is_qcfail or read.is_duplicate:
                continue

            if exclude_secondary and read.is_secondary:
                continue

            if exclude_supplementary and read.is_supplementary:
                continue

            if exclude_mq0 and read.mapping_quality <= 0:
                continue

            chrom_id = read.reference_id
            if chrom_id not in chrom_id_to_name:
                chrom_id_to_name[chrom_id] = name_mapping[bam.get_reference_name(chrom_id)]
            chrom_name = chrom_id_to_name[chrom_id]
            chroms.add(chrom_name)

            query = f"{read.query_name}_{2 if read.is_read2 else 1}"
            if query not in read_info:
                read_info[query] = {
                    "primary": "",
                    "secondary": [],
                    "supplementary": [],
                }

            if read.is_secondary:
                read_info[query]["secondary"].append(chrom_name)
            elif read.is_supplementary:
                read_info[query]["supplementary"].append(chrom_name)
            else:
                assert read_info[query]["primary"] == "", f"Multiple Primary Reads! \n{query=}\n{read=}"
                read_info[query]["primary"] = chrom_name

        # _summarise_reads(read_info).to_csv("reads_summary.csv")
        return summarise_by_chrom(chroms, read_info)


def _summarise_reads(read_info: dict[str, dict[str, typing.Any]]) -> pd.DataFrame:
    """Summarise the alignment information by read"""
    read_counts = {
        query: {
            "primary": read_dict["primary"],
            "n_secondary": len(read_dict["secondary"]),
            "n_supplementary": len(read_dict["supplementary"]),
            "secondary": ",".join(read_dict["secondary"]),
            "supplementary": ",".join(read_dict["supplementary"]),
        }
        for query, read_dict in read_info.items()
    }
    df = pd.DataFrame.from_dict(read_counts, orient="index")
    df.reset_index(inplace=True, names="query")
    return df


def summarise_by_chrom(chroms: set[str], read_info: dict[str, dict[str, typing.Any]]) -> pd.DataFrame:
    """Summarise the alignment information by reference"""
    chrom_info = {
        chrom: {
            "unique_alns": 0,  # read maps only once, and is to this chrom
            "unique_reads": 0,  # read only maps to this chrom
            "primary_alns": 0,
            "total_reads": 0,
            "total_alns": 0,
        }
        for chrom in chroms
    }
    for _, read_dict in read_info.items():
        primary = read_dict["primary"]
        secondary = read_dict["secondary"]
        supplementary = read_dict["supplementary"]

        if primary == "":
            for c in set(secondary) | set(supplementary):
                chrom_info[c]["total_reads"] += 1
            for c in secondary + supplementary:
                chrom_info[c]["total_alns"] += 1
            continue

        rest = set(secondary) | set(supplementary)

        chrom_info[primary]["primary_alns"] += 1
        if len(rest) == 0:
            chrom_info[primary]["unique_alns"] += 1
            chrom_info[primary]["unique_reads"] += 1
        elif rest == {primary}:
            chrom_info[primary]["unique_reads"] += 1

        for c in rest | {primary}:
            chrom_info[c]["total_reads"] += 1
        for c in secondary + supplementary + [primary]:
            chrom_info[c]["total_alns"] += 1

    df = pd.DataFrame.from_dict(chrom_info, orient="index")
    df.reset_index(inplace=True, names="genome_name")
    df.sort_values(by=["total_reads"], ascending=False, inplace=True)
    return df


def cli_entry_point():
    """Main function for the script"""
    parser = argparse.ArgumentParser(description="Process bam file")
    parser.add_argument("--bam", help="BAM file to process", required=True)
    parser.add_argument("--species_list", help="Reference names file", required=True)
    parser.add_argument("--output", help="Output file", required=True)
    parser.add_argument("--exclude_secondary", action="store_true", default=False)
    parser.add_argument("--exclude_supplementary", action="store_true", default=False)
    parser.add_argument("--exclude_mq0", action="store_true", default=False)
    args = parser.parse_args()

    name_mapping = get_name_mapping(args.species_list)
    df = get_alignment_stats(
        args.bam,
        name_mapping,
        exclude_secondary=args.exclude_secondary,
        exclude_supplementary=args.exclude_supplementary,
        exclude_mq0=args.exclude_mq0,
    )
    df.to_csv(args.output, index=False)


if __name__ == "__main__":
    cli_entry_point()
