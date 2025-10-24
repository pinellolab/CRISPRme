"""
This module provides functionality for merging target data from input files based 
on specified criteria. It includes functions for parsing command line arguments, 
processing target data, and writing results to output files.

Key functions include:
- `parse_commandline`: Parses and validates command line arguments for target 
    merging configuration.
- `split_target_row`: Splits a target row string into its individual components.
- `update_target_fields`: Updates specified fields in the target list by appending 
    corresponding values from another list.
- `distribute_targets`: Distributes targets between reference and variant targets 
    from a given cluster.
- `target_only_var`: Updates a target list to indicate whether it is a variant-only 
    target.
- `remove_duplicate_targets`: Removes duplicate values from specified fields in a 
    target list.
- `unfold_variant_targets`: Recovers and processes all variant targets from a given 
    dictionary.
- `sorting_score`: Generates a sorting key function based on specified criteria 
    for sorting.
- `sorting_fewest`: Creates a sorting key function based on the fewest specified 
    criteria.
- `initialize_sorting_criteria`: Initializes a sorting criteria function based on 
    the provided parameters.
- `retrieve_best_target`: Identifies and retrieves the best target from a given 
    cluster of targets.
- `merge_targets`: Merges target data from an input file and writes the best 
    targets to an output file.
- `main`: The entry point of the module that orchestrates the merging process.

This module is designed to facilitate the analysis of genomic target data, allowing 
users to efficiently merge and sort targets based on various criteria.
"""

from typing import List, Tuple, Dict, Callable
from time import time
from io import TextIOWrapper

import sys
import os

SORTCRITERIA = {"mm": 2, "bulges": 1, "mm+bulges": 0}


def parse_commandline(
    args: List[str],
) -> Tuple[str, str, int, int, int, int, int, int, int, str, List[str]]:
    """
    Parses command line arguments for target merging configuration.
    This function validates the input arguments and returns the necessary
    parameters for processing target files.

    Args:
        args (List[str]): A list of command line arguments where:
            - args[0] is the targets file name.
            - args[1] is the output file name.
            - args[2] is the targets merge range in base pairs.
            - args[3] is the chromosome column index (1-based).
            - args[4] is the position column index (1-based).
            - args[5] is the mm+bulges column index (1-based).
            - args[6] is the guide column index (1-based).
            - args[7] is the SNP info column index (1-based).
            - args[8] is the score column index (1-based).
            - args[9] is the sorting pivot (score or mm+bulges).
            - args[10] is a comma-separated list of sorting criteria.

    Returns:
        Tuple[str, str, int, int, int, int, int, int, int, str, List[str]]:
        A tuple containing the parsed parameters:
            - targets_fname: The targets file name.
            - outfname: The output file name.
            - rangebp: The targets merge range in base pairs.
            - chromidx: The chromosome column index (0-based).
            - posidx: The position column index (0-based).
            - mmbidx: The mm+bulges column index (0-based).
            - guideidx: The guide column index (0-based).
            - snpidx: The SNP info column index (0-based).
            - scoreidx: The score column index (0-based).
            - pivot: The sorting pivot.
            - sortcrit: A list of sorting criteria.

    Raises:
        FileNotFoundError: If the targets file cannot be found.
        ValueError: If the merge range is less than 1 or if invalid sort criteria are provided.
    """

    targets_fname = args[0]  # targets file
    if not os.path.isfile(targets_fname):
        raise FileNotFoundError(f"Unable to locate {targets_fname}")
    outfname = args[1]  # output file
    rangebp = int(args[2])  # targets merge range (bp)
    if rangebp < 1:
        raise ValueError(f"Forbidden targets merge range ({rangebp})")
    chromidx = int(args[3]) - 1  # chromosome col idx
    posidx = int(args[4]) - 1  # position col idx
    mmbidx = int(args[5]) - 1  # mm+bulges col idx
    guideidx = int(args[6]) - 1  # guide col idx
    snpidx = int(args[7]) - 1  # snp info col idx
    scoreidx = int(args[8]) - 1  # score col idx
    pivot = args[9]  # targets sorting pivot (score or mm+bulges)
    # comma-separated list of criteria to use while sorting targets
    sortcrit = args[10].split(",")
    if len(sortcrit) > 3 or any(c not in SORTCRITERIA for c in sortcrit):
        offending_vals = ",".join([c for c in sortcrit if c not in SORTCRITERIA])
        raise ValueError(f"Forbidden sort criteria selected: {offending_vals}")
    return (
        targets_fname,
        outfname,
        rangebp,
        chromidx,
        posidx,
        mmbidx,
        guideidx,
        snpidx,
        scoreidx,
        pivot,
        sortcrit,
    )


def split_target_row(
    target_row: str, guideidx: int, chromidx: int, posidx: int
) -> Tuple[str, str, int]:
    """
    Splits a target row string into its individual components.
    This function retrieves the guide, chromosome, and position from a target
    row based on specified indices.

    Args:
        target_row (str): A string representing a single target row, with fields
            separated by whitespace.
        guideidx (int): The index of the guide field in the target row.
        chromidx (int): The index of the chromosome field in the target row.
        posidx (int): The index of the position field in the target row.

    Returns:
        Tuple[str, str, int]: A tuple containing the extracted guide, chromosome,
            and position:
            - guide (str): The guide extracted from the target row.
            - chromosome (str): The chromosome extracted from the target row.
            - position (int): The position extracted from the target row, converted
                to an integer.
    """

    fields = target_row.strip().split()
    return fields[guideidx], fields[chromidx], int(fields[posidx])


def update_target_fields(
    target: List[str], fields: List[str], samplesidx: int, snpid_idx: int, afidx: int
) -> List[str]:
    """Update specified fields in the target list by appending corresponding
    values from the fields list.

    This function modifies the target list by concatenating values from the fields
    list at given indices, which is useful for aggregating information related to
    samples, SNP IDs, and allele frequencies.

    Args:
        target (List[str]): The list of target values to be updated.
        fields (List[str]): The list of new values to append to the target.
        samplesidx (int): The index in the target list for sample information.
        snpid_idx (int): The index in the target list for SNP ID information.
        afidx (int): The index in the target list for allele frequency information.

    Returns:
        List[str]: The updated target list with concatenated values.
    """

    target[samplesidx] = f"{target[samplesidx]},{fields[samplesidx]}"
    target[snpid_idx] = f"{target[snpid_idx]},{fields[snpid_idx]}"
    target[afidx] = f"{target[afidx]},{fields[afidx]}"
    return target


def distribute_targets(
    cluster: List[str],
    snpidx: int,
    posidx: int,
    snpid_idx: int,
    samplesidx: int,
    afidx: int,
) -> Tuple[List[List[str]], Dict[str, List[List[str]]]]:
    """
    Distributes targets between reference and variant targets from a given cluster.
    It merges identical targets found in different datasets into a structured
    format.

    This function processes a list of target strings, categorizing them into
    reference targets and variant targets based on specific indices. Reference
    targets are collected in a list, while variant targets are stored in a dictionary,
    allowing for the merging of identical targets across datasets.

    Args:
        cluster (List[str]): A list of target strings to be processed.
        snpidx (int): The index indicating the SNP status of the target.
        posidx (int): The index indicating the position of the target.
        snpid_idx (int): The index for SNP IDs in the target fields.
        samplesidx (int): The index for sample information in the target fields.
        afidx (int): The index for allele frequencies in the target fields.

    Returns:
        Tuple[List[List[str]], Dict[str, List[List[str]]]]: A tuple containing a
            list of reference
        targets and a dictionary of variant targets, where each key corresponds
            to a unique target and its value is a list of associated target fields.

    Raises:
        ValueError: If the input data is malformed or indices are out of range.
    """

    # distribute targets between reference and variant targets
    # dict used to merge identical targets found in different datasets
    reftargets, vartargets = [], {}
    for target in cluster:
        fields = target.strip().split()  # retrieve target fields
        if fields[snpidx] == "n":  # target found in reference
            reftargets.append(fields)
        else:  # target found in variant genomes
            targetkey = f"{fields[posidx]}_{fields[snpidx]}"
            current_target = vartargets.get(targetkey)
            if current_target:
                # update current sample list, snp ids, and allele freqs
                vartargets[targetkey][0] = update_target_fields(
                    current_target[0], fields, samplesidx, snpid_idx, afidx
                )
            else:
                vartargets[targetkey] = [fields]  # first target at position
    return reftargets, vartargets


def target_only_var(target: List[str], varonly: bool) -> List[str]:
    """
    Updates a target list to indicate whether it is a variant-only target.
    This function modifies the target's status based on the provided flag.

    The function checks the `varonly` flag and updates the 13th element of the
    target list to "y" if the flag is set to True, indicating that the target is
    only found in variant genomes. It returns the modified target list.

    Args:
        target (List[str]): A list representing the target information.
        varonly (bool): A flag indicating whether the target is variant-only.

    Returns:
        List[str]: The updated target list with the appropriate status.
    """

    # set apprpriate flag if no target in reference in current cluster
    target[12] = "y" if varonly else target[12]
    return target


def remove_duplicate_targets(
    target: List[str], snpidx: int, snpid_idx: int, afidx: int, samplesidx: int
) -> List[str]:
    """
    Removes duplicate values from specified fields in a target list.
    This function ensures that SNP IDs, SNP information, allele frequencies, and
    sample data contain only unique entries.

    The function takes a target list and specified indices for SNPs, SNP IDs,
    allele frequencies, and samples. It processes each of these fields to eliminate
    duplicates by converting them into sets and then back into comma-separated
    strings, returning the modified target list.

    Args:
        target (List[str]): A list representing the target information.
        snpidx (int): The index for SNP information in the target list.
        snpid_idx (int): The index for SNP IDs in the target list.
        afidx (int): The index for allele frequencies in the target list.
        samplesidx (int): The index for sample information in the target list.

    Returns:
        List[str]: The updated target list with duplicates removed from specified
            fields.
    """

    # remove duplicate values from snp ids, snp info, allele freqs, and samples
    target[snpidx] = ",".join(set(target[snpidx].split(",")))
    target[snpid_idx] = ",".join(set(target[snpid_idx].split(",")))
    target[afidx] = ",".join(set(target[afidx].split(",")))
    target[samplesidx] = ",".join(set(target[samplesidx].split(",")))
    return target


def unfold_variant_targets(
    vartargets: Dict[str, List[List[str]]],
    varonly: bool,
    snpidx: int,
    snpid_idx: int,
    afidx: int,
    samplesidx: int,
) -> List[List[str]]:
    """
    Recovers and processes all variant targets from a given dictionary.
    This function compiles variant targets into a single list while applying
    necessary transformations and removing duplicates.

    The function iterates through the provided dictionary of variant targets,
    applying the `target_only_var` function to each target based on the `varonly`
    flag. It then removes duplicates from the resulting list of targets using the
    `remove_duplicate_targets` function, returning a cleaned list of variant targets.

    Args:
        vartargets (Dict[str, List[List[str]]]): A dictionary of variant targets,
            where each key corresponds to a unique target identifier and the value
            is a list of target fields.
        varonly (bool): A flag indicating whether to mark targets as variant-only.
        snpidx (int): The index for SNP information in the target list.
        snpid_idx (int): The index for SNP IDs in the target list.
        afidx (int): The index for allele frequencies in the target list.
        samplesidx (int): The index for sample information in the target list.

    Returns:
        List[List[str]]: A list of processed variant targets with duplicates removed.
    """

    # recover all variant targets and store in a list
    vartargets_list = []
    for targets in vartargets.values():
        vartargets_list.extend([target_only_var(t, varonly) for t in targets])
    # remove duplicate values in targets
    return [
        remove_duplicate_targets(t, snpidx, snpid_idx, afidx, samplesidx)
        for t in vartargets_list
    ]


def sorting_score(criteria: List[str], score_idx: int, mmbidx: int) -> Callable:
    """
    Generates a sorting key function based on specified criteria for sorting.
    This function allows for dynamic sorting of items based on one to three criteria,
    prioritizing scores and additional metrics.

    The function returns a callable that can be used as a key in sorting operations.
    Depending on the number of criteria provided, it constructs a tuple that includes
    the negative score (to sort in descending order) and additional metrics derived
    from the specified indices, ensuring that items are sorted according to the defined
    priorities.

    Args:
        criteria (List[str]): A list of criteria used for sorting.
        score_idx (int): The index of the score in the items to be sorted.
        mmbidx (int): The base index used to calculate additional metrics from
            the criteria.

    Returns:
        Callable: A function that takes an item and returns a tuple for sorting purposes.
    """

    if len(criteria) == 1:  # single criterion
        return lambda x: (
            -float(x[score_idx]),
            int(x[mmbidx - SORTCRITERIA[criteria[0]]]),
        )
    elif len(criteria) == 2:
        return lambda x: (
            -float(x[score_idx]),
            int(x[mmbidx - SORTCRITERIA[criteria[0]]]),
            int(x[mmbidx - SORTCRITERIA[criteria[1]]]),
        )
    # base case (all three )
    return lambda x: (
        -float(x[score_idx]),
        int(x[mmbidx - SORTCRITERIA[criteria[0]]]),
        int(x[mmbidx - SORTCRITERIA[criteria[1]]]),
        int(x[mmbidx - SORTCRITERIA[criteria[2]]]),
    )


def sorting_fewest(criteria: List[str], mmbidx: int) -> Callable:
    """
    Creates a sorting key function based on the fewest specified criteria.
    This function allows for sorting items by one to three criteria, focusing on
    the values derived from the specified indices.

    The function returns a callable that can be used as a key in sorting operations.
    Depending on the number of criteria provided, it constructs a tuple of integer
    values from the specified indices, enabling sorting based on the defined
    priorities.

    Args:
        criteria (List[str]): A list of criteria used for sorting.
        mmbidx (int): The base index used to calculate values from the criteria.

    Returns:
        Callable: A function that takes an item and returns a tuple for sorting
            purposes.
    """

    if len(criteria) == 1:  # one criterion
        return lambda x: (int(x[mmbidx - SORTCRITERIA[criteria[0]]]))
    elif len(criteria) == 2:
        return lambda x: (
            int(x[mmbidx - SORTCRITERIA[criteria[0]]]),
            int(x[mmbidx - SORTCRITERIA[criteria[1]]]),
            # int(x[mmbidx - 2]),
            # int(x[mmbidx - 1]),
        )
    # base case (all three )
    return lambda x: (
        int(x[mmbidx - SORTCRITERIA[criteria[0]]]),
        int(x[mmbidx - SORTCRITERIA[criteria[1]]]),
        int(x[mmbidx - SORTCRITERIA[criteria[2]]]),
    )


def initialize_sorting_criteria(
    criteria: List[str], scoreidx: int, mmbidx: int, score: bool
) -> Callable:
    """
    Initializes a sorting criteria function based on the provided parameters.
    This function determines whether to sort by score or by the fewest criteria
    based on the input flag.

    Depending on the value of the `score` flag, the function returns either a
    sorting function that prioritizes scores or one that focuses on the fewest
    specified criteria. This allows for flexible sorting behavior based on the
    user's needs.

    Args:
        criteria (List[str]): A list of criteria used for sorting.
        scoreidx (int): The index of the score in the items to be sorted.
        mmbidx (int): The base index used to calculate additional metrics from
            the criteria.
        score (bool): A flag indicating whether to sort by score or by the fewest
            criteria.

    Returns:
        Callable: A function that can be used as a key in sorting operations.
    """

    if score:
        return sorting_score(criteria, scoreidx, mmbidx)
    return sorting_fewest(criteria, mmbidx)


def retrieve_best_target(
    cluster: List[str],
    snpidx: int,
    posidx: int,
    guideidx: int,
    scoreidx: int,
    mmbidx: int,
    pivot: str,
    sorting_criteria: List[str],
    outfile: TextIOWrapper,
    outfile_disc: TextIOWrapper,
) -> None:
    """
    Identifies and retrieves the best target from a given cluster of targets.
    This function processes the targets based on specified criteria and outputs
    the best target along with alternative alignments to the provided output files.

    The function first distributes the targets into reference and variant categories,
    then unfolds the variant targets into a list. It sorts both reference and variant
    targets according to the specified sorting criteria, determines the best target
    based on the pivot condition, and writes the results to the specified output files,
    including the count of remaining targets.

    Args:
        cluster (List[str]): A list of target strings representing the cluster.
        snpidx (int): The index indicating the SNP status of the target.
        posidx (int): The index indicating the position of the target.
        guideidx (int): The index for guide information in the target fields.
        scoreidx (int): The index for scores in the target fields.
        mmbidx (int): The base index used to calculate additional metrics from
            the criteria.
        pivot (str): A string indicating the pivot for sorting (e.g., "score").
        sorting_criteria (List[str]): A list of criteria used for sorting the
            targets.
        outfile (TextIOWrapper): The output file for the best target.
        outfile_disc (TextIOWrapper): The output file for alternative alignments.

    Returns:
        None: This function does not return a value but writes output to the
            specified files.
    """

    if not cluster:  # opening the first cluster, it will be empty
        return  # do nothing
    reftargets, vartargets = distribute_targets(
        cluster, snpidx, posidx, snpidx - 2, guideidx - 2, snpidx - 1
    )
    varonly = not reftargets  # check if found only variant targets
    # retrieve variant targets in list
    vartargets = unfold_variant_targets(
        vartargets, varonly, snpidx, snpidx - 2, snpidx - 1, guideidx - 2
    )
    # sort targets using the criteria specified in input
    score = pivot == "score"
    if reftargets:
        reftargets = sorted(
            reftargets,
            key=initialize_sorting_criteria(sorting_criteria, scoreidx, mmbidx, score),
        )
    if vartargets:
        vartargets = sorted(
            vartargets,
            key=initialize_sorting_criteria(sorting_criteria, scoreidx, mmbidx, score),
        )
    if varonly:
        target = vartargets.pop(0)  # retrieve best target
        # count the targets remaining in the cluster
        target[scoreidx - 1] = str(len(vartargets))
    elif reftargets and vartargets:
        if score:  # check on score
            target = (
                vartargets.pop(0)
                if float(vartargets[0][scoreidx]) > float(reftargets[0][scoreidx])
                else reftargets.pop(0)
            )
        elif int(vartargets[0][mmbidx]) < int(reftargets[0][mmbidx]):
            target = vartargets.pop(0)
        else:
            target = reftargets.pop(0)
        # count the targets remaining in the cluster
        target[scoreidx - 1] = str(len(reftargets) + len(vartargets))
    else:
        target = reftargets.pop(0)
        target[scoreidx - 1] = str(len(reftargets))
    outfile.write("\t".join(target) + "\n")  # report the best target
    # write alternative alignments
    for target in reftargets + vartargets:
        target[scoreidx - 1] = str(len(reftargets) + len(vartargets))
        outfile_disc.write("\t".join(target) + "\n")  # report the alternative target


def merge_targets(
    inargs: Tuple[str, str, int, int, int, int, int, int, int, str, List[str]]
) -> None:
    """
    Merges target data from an input file and writes the best targets to an output
    file. This function processes clusters of targets based on specified criteria
    and handles discarded samples in a separate output file.

    The function reads target data from the input file, grouping targets into
    clusters based on guide, chromosome, and position. It retrieves the best target
    from each cluster and writes the results to the specified output files, ensuring
    that discarded samples are also recorded.

    Args:
        inargs (Tuple[str, str, int, int, int, int, int, int, int, str, List[str]]):
            A tuple containing input parameters, including input and output file
            names, indices for various target fields, and sorting criteria.

    Returns:
        None: This function does not return a value but writes output to the
            specified files.
    """

    outfname_disc = f"{inargs[1]}.discarded_samples"  # discarded targets file
    # initialize variables used during merge
    prevpos, prevguide, prevchrom, cluster = -(inargs[2] + 1), "", "", []
    with open(inargs[0], mode="r") as infile:
        with open(inargs[1], mode="w") as outfile:
            with open(outfname_disc, mode="w") as outfile_disc:
                # header placed in both outfiles
                header = infile.readline()
                outfile_disc.write(header)
                outfile.write(header)
                for line in infile:  # start reading targets
                    # retrieve guide chromosome and position
                    guide, chrom, pos = split_target_row(
                        line, inargs[6], inargs[3], inargs[4]
                    )
                    # open new targets cluster and retrieve the best target from previous cluster
                    if (
                        prevguide != guide
                        or prevchrom != chrom
                        or (pos - prevpos) > inargs[2]
                    ):
                        retrieve_best_target(
                            cluster,
                            inargs[7],
                            inargs[4],
                            inargs[6],
                            inargs[8],
                            inargs[5],
                            inargs[9],
                            inargs[10],
                            outfile,
                            outfile_disc,
                        )
                        cluster = [line]
                    else:  # append target data to current cluster
                        cluster.append(line)
                    # update lookup variables
                    prevpos, prevguide, prevchrom = pos, guide, chrom
                retrieve_best_target(
                    cluster,
                    inargs[7],
                    inargs[4],
                    inargs[6],
                    inargs[8],
                    inargs[5],
                    inargs[9],
                    inargs[10],
                    outfile,
                    outfile_disc,
                )  # process the last cluster


def main():
    # read input args
    inargs = parse_commandline(sys.argv[1:])
    start = time()
    merge_targets(inargs)
    sys.stdout.write(f"Targets merge completed in {(time() - start):.2f}s\n")


if __name__ == "__main__":
    main()
