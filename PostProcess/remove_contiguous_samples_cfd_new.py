"""
This script parses the input arguments and removes contiguous samples from a file based on specified criteria.

The script expects the following input arguments:
- infname (str): The input file name.
- outfname (str): The output file name.
- rangebp (int): The maximum range in base pairs for samples to be considered contiguous.
- chrom_idx (int): The index of the chromosome column in the input file.
- position_idx (int): The index of the position column in the input file.
- mm_bul_count_idx (int): The index of the mismatch and bulge count column in the input file.
- guide_idx (int): The index of the guide column in the input file.
- snp_info_idx (int): The index of the SNP info column in the input file.
- score_idx (int): The index of the score column in the input file.
- sort_criterion (str): The sorting criterion for the output file.

The script reads the input file, parses the specified columns, and removes contiguous samples based on the given range. The resulting data is written to the output file.

Examples:
    $ python remove_contiguous_samples_cfd.py input.txt output.txt 10 2 3 4 5 6 7 8 score
"""

from typing import List, Tuple, Dict, Callable, Union
from io import TextIOWrapper
from time import time

import sys
import os

INPUT_ARG_COUNT = 12
FLAG_ALT_ONLY = 12
SORTING_PIVOTS = ["score", "total"]
SORTING_CRITERIA = {"mm+bulges": 0, "mm": 2, "bulges": 1}


class MergeTargets:
    """
    Represents a class for merging targets with specified criteria for sorting
    and indexing.

    Args:
        args (List[str]): List of arguments containing input parameters for
            target merging.

    Returns:
        None

    Raises:
        FileNotFoundError: If the input targets file is not found.
        ValueError: If there are issues with the specified merge range, column
            indices, or sorting criteria.
    """

    def __init__(self, args: List[str]) -> None:
        self._infname = args[0]  # input targets
        if not os.path.exists(self._infname) or not os.path.isfile(self._infname):
            raise FileNotFoundError(f"{self._infname} not found")
        self._outfname = args[1]  # output merged targets
        self._rangebp = int(args[2])  # merge bp range
        if self._rangebp <= 0:
            raise ValueError(f"Invalid merge range ({self._rangebp})")
        self._chromidx = int(args[3]) - 1  # chromosome index
        if self._chromidx != 4:
            raise ValueError(
                f"Chromosome data is expected on column 5, got {self._chromidx}"
            )
        self._posidx = int(args[4]) - 1  # position index
        if self._posidx != 6:
            raise ValueError(
                f"Position data is expected on column 7, got {self._posidx}"
            )
        self._mmb_idx = int(args[5]) - 1  # mm+bulges index
        if self._mmb_idx != 10:
            raise ValueError(
                f"MM+Bulges data is expected on column 11, got {self._mmb_idx}"
            )
        self._guide_idx = int(args[6]) - 1  # guide index
        if self._guide_idx != 15:
            raise ValueError(
                f"Guide data is expected on column 16, got {self._guide_idx}"
            )
        self._snp_idx = int(args[7]) - 1  # snp info index
        if self._snp_idx != 18:
            raise ValueError(f"SNP data is expected on column 19, got {self._snp_idx}")
        self._score_idx = int(args[8]) - 1  # score index
        if self._score_idx != 20:
            raise ValueError(
                f"Score data is expected on column 21, got {self._score_idx}"
            )
        self._sort_criterion = args[9]  # sorting criterion (pivot)
        if self._sort_criterion not in SORTING_PIVOTS:
            raise ValueError(
                f"Allowed sort pivots: {SORTING_PIVOTS}, got {self._sort_criterion}"
            )
        self._sorting_criteria_scoring = args[10]  # sorting criteria (score is pivot)
        self._sorting_criteria = args[11]  # sorting criteria (mm+bulges is pivot)

    @property
    def targets_fname(self) -> str:
        return self._infname

    @property
    def targets_fname_merged(self) -> str:
        return self._outfname

    @property
    def targets_fname_discarded(self) -> str:
        return f"{self._outfname}.discarded_samples"

    @property
    def rangebp(self) -> int:
        return self._rangebp

    @property
    def chromidx(self) -> int:
        return self._chromidx

    @property
    def posidx(self) -> int:
        return self._posidx

    @property
    def mmbidx(self) -> int:
        return self._mmb_idx

    @property
    def guideidx(self) -> int:
        return self._guide_idx

    @property
    def snpidx(self) -> int:
        return self._snp_idx

    @property
    def scoreidx(self) -> int:
        return self._score_idx

    @property
    def sort_pivot(self) -> str:
        return self._sort_criterion

    @property
    def sorting_criteria_scoring(self) -> str:
        return self._sorting_criteria_scoring

    @property
    def sorting_criteria(self) -> str:
        return self._sorting_criteria


def parse_input_args(args: List[str]) -> MergeTargets:
    """
    Parses the input arguments to create an instance of MergeTargets for handling
    target merging based on the provided arguments.

    Args:
        args (List[str]): List of input arguments to be parsed for target merging.

    Returns:
        MergeTargets: An instance of MergeTargets class for further processing.
    """

    if not args:
        raise ValueError("No input argument provided")
    if len(args) != INPUT_ARG_COUNT:
        raise ValueError(f"Expected {INPUT_ARG_COUNT} arguments, got {len(args)}")
    return MergeTargets(args)  # parse input arguments


def initialize_targets_cluster(rangebp: int) -> Tuple[int, str, str, str, List[str]]:
    """
    Initializes the targets cluster by setting the previous cluster position,
    guide, chromosome, SNP, and an open starting cluster.

    Args:
        rangebp (int): The range of base pairs for the cluster.

    Returns:
        Tuple[int, str, str, str, List[str]]: Tuple containing the initialized
            cluster parameters.
    """

    pos_prev = -(rangebp + 1)  # previous cluster position
    guide_prev = ""  # previous cluster's guide
    chrom_prev = ""  # previous cluster's  chromosome
    snp_prev = ""  # previous cluster snp
    cluster = []  # open starting cluster
    return pos_prev, guide_prev, chrom_prev, snp_prev, cluster


def parse_target(target: str) -> List[str]:
    """
    Parses the fields of a target string and returns a list of the parsed fields.

    Args:
        target (str): The target string to be parsed.

    Returns:
        List[str]: List of parsed fields extracted from the target string.
    """

    # parse target's fields
    return target.strip().split()


def open_targets_cluster(
    guide: str,
    guide_prev: str,
    chrom: str,
    chrom_prev: str,
    pos: int,
    pos_prev: int,
    rangebp: int,
) -> bool:
    """
    Determines whether a new targets cluster should be opened based on changes
    in guide, chromosome, and position within a specified range of base pairs.

    Args:
        guide (str): Current guide.
        guide_prev (str): Previous guide.
        chrom (str): Current chromosome.
        chrom_prev (str): Previous chromosome.
        pos (int): Current position.
        pos_prev (int): Previous position.
        rangebp (int): Range of base pairs for cluster merging.

    Returns:
        bool: True if a new targets cluster should be opened, False otherwise.
    """

    return guide != guide_prev or chrom != chrom_prev or pos - pos_prev > rangebp


def cluster_targets_by_pos_snp(
    cluster: List[List[str]], guideidx: int, snpidx: int, posidx: int
) -> Tuple[List[List[str]], Dict[Tuple[int, str], List[List[str]]]]:
    """
    Clusters targets based on position and SNP information, separating them into
    reference and alternative targets.

    Args:
        cluster (List[List[str]]): List of target clusters.
        guideidx (int): Index of the guide.
        snpidx (int): Index of the SNP information.
        posidx (int): Index of the position.

    Returns:
        Tuple[List[List[str]], Dict[Tuple[int, str], List[List[str]]]: A tuple
            containing the reference targets list and a dictionary of alternative
            targets grouped by position and SNP information.
    """

    # initialize reference and alternative targets lists
    reference_targets, alternative_targets = [], {}
    samplesidx, snpididx, afidx = guideidx - 2, snpidx - 2, snpidx - 1
    for target in cluster:
        if target[snpidx] == "n":  # reference target
            reference_targets.append(target)
            continue
        # alternative target
        pos_t, snp_info_t = target[posidx], target[snpidx]
        if (pos_t, snp_info_t) in alternative_targets:
            alternative_targets[(pos_t, snp_info_t)][0][
                samplesidx
            ] += f",{target[samplesidx]}"  # add samples
            alternative_targets[(pos_t, snp_info_t)][0][
                snpididx
            ] += f",{target[snpididx]}"  # add snp ids
            alternative_targets[(pos_t, snp_info_t)][0][
                afidx
            ] += f",{target[afidx]}"  # add allele frequencies
        else:
            alternative_targets[(pos_t, snp_info_t)] = [target]
    return reference_targets, alternative_targets


def recover_alternative_targets_list(
    alternative_targets: Dict[Tuple[int, str], List[List[str]]], noref: bool
) -> List[List[str]]:
    """
    Recovers the list of alternative targets from a dictionary of alternative
    targets grouped by position and SNP information.

    Args:
        alternative_targets (Dict[Tuple[int, str], List[List[str]]): Dictionary
            of alternative targets grouped by position and SNP information.
        noref (bool): Flag indicating whether to include reference targets.

    Returns:
        List[List[str]]: List of alternative targets with updated flags.
    """

    alternative_targets_list = []  # initialize the list
    for v in alternative_targets.values():
        for target in v:
            # check whether the target is uniquely found on alternative sequences
            target[FLAG_ALT_ONLY] = "y" if noref else target[FLAG_ALT_ONLY]
            alternative_targets_list.append(target)
    return alternative_targets_list


def remove_duplicates(target: List[str], idx: int) -> str:
    """
    Removes duplicates from a target list at the specified index and returns the
    unique values as a comma-separated string.

    Args:
        target (List[str]): List of target values.
        idx (int): Index of the target list to remove duplicates from.

    Returns:
        str: Comma-separated string of unique values after removing duplicates.
    """

    return ",".join(set(target[idx].split(",")))


def unique_values(
    targets: List[List[str]], snpidx: int, guideidx: int
) -> List[List[str]]:
    """
    Returns a list of targets with unique values for SNP information, SNP ID,
    allele frequency, and samples.

    Args:
        targets (List[List[str]]): List of target clusters.
        snpidx (int): Index of the SNP information.
        guideidx (int): Index of the guide.

    Returns:
        List[List[str]]: List of targets with unique values for specified fields.
    """

    snpididx, afidx, samplesidx = snpidx - 2, snpidx - 1, guideidx - 2
    targets_filt = []
    # remove duplicate values
    for target in targets:
        target[snpidx] = remove_duplicates(target, snpidx)  # snp info
        target[snpididx] = remove_duplicates(target, snpididx)  # snp id
        target[afidx] = remove_duplicates(target, afidx)  # af
        target[samplesidx] = remove_duplicates(target, samplesidx)  # samples
        targets_filt.append(target)  # update target
    return targets_filt


def construct_cluster(
    cluster: List[List[str]],
    guideidx: int,
    snpidx: int,
    posidx: int,
    scoreidx: int,
    mmbidx: int,
    sort_pivot: str,
    sorting_criteria: str,
    sorting_criteria_scoring: str,
    outfile: TextIOWrapper,
    outfiledisc: TextIOWrapper,
) -> None:
    """
    Constructs target clusters by processing and sorting reference and
    alternative targets based on specified criteria, and writes the results to
    output files.

    Args:
        cluster (List[List[str]]): List of target clusters.
        guideidx (int): Index of the guide.
        snpidx (int): Index of the SNP information.
        posidx (int): Index of the position.
        scoreidx (int): Index of the score.
        mmbidx (int): Index of the MM+Bulges.
        sort_pivot (str): Sorting pivot.
        sorting_criteria (str): Sorting criteria.
        sorting_criteria_scoring (str): Sorting criteria for scoring.
        outfile (TextIOWrapper): Output file for reported alignments.
        outfiledisc (TextIOWrapper): Output file for alternative alignments.

    Returns:
        None
    """

    if not cluster:  # avoid crashes when cluster is empty
        return
    # recover reference and alternative targets
    reference_targets, alternative_targets = cluster_targets_by_pos_snp(
        cluster, guideidx, snpidx, posidx
    )
    noref = (
        not reference_targets
    )  # check whether only alternative targets have been found
    # recover alternative targets list
    alternative_targets = recover_alternative_targets_list(alternative_targets, noref)
    # remove duplicates values on snp info, snp id, af, and samples columns
    alternative_targets = unique_values(alternative_targets, snpidx, guideidx)
    # sort targets
    score = sort_pivot == SORTING_PIVOTS[0]
    sorting_criteria = sorting_criteria_scoring if score else sorting_criteria
    criteria = initialize_sorting_criteria(sorting_criteria, score, scoreidx, mmbidx)
    if reference_targets:
        reference_targets = sort_targets(reference_targets, criteria)
    if alternative_targets:
        alternative_targets = sort_targets(alternative_targets, criteria)
    # write targets to reported or alternative alignments files
    write_best_targets(
        reference_targets,
        alternative_targets,
        noref,
        score,
        scoreidx,
        mmbidx,
        outfile,
        outfiledisc,
    )


def write_best_targets(
    reference_targets: List[List[str]],
    alternative_targets: List[List[str]],
    noref: bool,
    score: bool,
    scoreidx: int,
    mmbidx: int,
    outfile: TextIOWrapper,
    outfiledisc: TextIOWrapper,
) -> None:
    """
    Writes the best targets to the reported alignments file and alternative
    alignments file based on specified criteria.

    Args:
        reference_targets (List[List[str]]): List of reference targets.
        alternative_targets (List[List[str]]): List of alternative targets.
        noref (bool): Flag indicating if no reference target is found.
        score (bool): Flag indicating if scoring is used for comparison.
        scoreidx (int): Index of the score.
        mmbidx (int): Index of the MM+Bulges.
        outfile (TextIOWrapper): Output file for reported alignments.
        outfiledisc (TextIOWrapper): Output file for alternative alignments.

    Returns:
        None
    """

    if noref:  # no reference target found
        target = alternative_targets.pop(0)  # pop best target
        target[scoreidx - 1] = str(len(alternative_targets))
    elif reference_targets and alternative_targets:  # targets found in both
        if score:
            target = (
                alternative_targets.pop(0)
                if float(reference_targets[0][scoreidx])
                < float(alternative_targets[0][scoreidx])
                else reference_targets.pop(0)
            )
        elif int(alternative_targets[0][mmbidx]) < int(reference_targets[0][mmbidx]):
            target = alternative_targets.pop(0)
        else:
            target = reference_targets.pop(0)
        target[scoreidx - 1] = str(len(reference_targets) + len(alternative_targets))
    else:  # no alternative target found
        target = reference_targets.pop(0)
        target[scoreidx - 1] = str(len(reference_targets))
    target = "\t".join(target)
    outfile.write(f"{target}\n")
    # write the alternative alignments targets
    for target in reference_targets + alternative_targets:
        target[scoreidx - 1] = str(len(reference_targets) + len(alternative_targets))
        target = "\t".join(target)
        outfiledisc.write(f"{target}\n")


def sort_targets(targets: List[List[str]], criteria: Callable) -> List[List[str]]:
    """
    Sorts the list of targets based on the specified criteria function.

    Args:
        targets (List[List[str]]): List of targets to be sorted.
        criteria (Callable): Sorting criteria function.

    Returns:
        List[List[str]]: Sorted list of targets based on the provided criteria.
    """

    return sorted(targets, key=criteria)


def sorting_score(criteria: List[str], score_idx: int, mmb_idx: int) -> Callable:
    """
    Returns a sorting function based on the specified criteria for scoring,
    MM+Bulges index, and multiple criteria.

    Args:
        criteria (List[str]): List of sorting criteria.
        score_idx (int): Index of the score.
        mmb_idx (int): Index of the MM+Bulges.

    Returns:
        Callable: Sorting function based on the provided criteria.
    """

    if len(criteria) == 1:  # single criterion
        return lambda x: (
            -float(x[score_idx]),
            int(x[mmb_idx - SORTING_CRITERIA[criteria[0]]]),
        )
    elif len(criteria) == 2:
        return lambda x: (
            -float(x[score_idx]),
            int(x[mmb_idx - SORTING_CRITERIA[criteria[0]]]),
            int(x[mmb_idx - SORTING_CRITERIA[criteria[1]]]),
        )
    # base case (all three )
    return lambda x: (
        -float(x[score_idx]),
        int(x[mmb_idx - SORTING_CRITERIA[criteria[0]]]),
        int(x[mmb_idx - SORTING_CRITERIA[criteria[1]]]),
        int(x[mmb_idx - SORTING_CRITERIA[criteria[2]]]),
    )


def sorting_fewest(criteria: List[str], mmb_idx: int) -> Callable:
    """
    Returns a sorting function based on the specified criteria for the fewest
    MM+Bulges index.

    Args:
        criteria (List[str]): List of sorting criteria.
        mmb_idx (int): Index of the MM+Bulges.

    Returns:
        Callable: Sorting function based on the provided criteria for the fewest
        MM+Bulges index.
    """

    if len(criteria) == 1:  # one criterion
        return lambda x: (int(x[mmb_idx - SORTING_CRITERIA[criteria[0]]]))
    elif len(criteria) == 2:
        return lambda x: (
            int(x[mmb_idx - SORTING_CRITERIA[criteria[0]]]),
            int(x[mmb_idx - SORTING_CRITERIA[criteria[1]]]),
        )
    # base case (all three )
    return lambda x: (
        int(x[mmb_idx - SORTING_CRITERIA[criteria[0]]]),
        int(x[mmb_idx - SORTING_CRITERIA[criteria[1]]]),
        int(x[mmb_idx - SORTING_CRITERIA[criteria[2]]]),
    )


def initialize_sorting_criteria(
    sorting_criteria: str, score: bool, score_idx: int, mmb_idx: int
) -> Callable:
    """
    Initializes the sorting criteria function based on the specified sorting
    criteria, score flag, score index, and MM+Bulges index.

    Args:
        sorting_criteria (str): Comma-separated string of sorting criteria.
        score (bool): Flag indicating if scoring is used for sorting.
        score_idx (int): Index of the score.
        mmb_idx (int): Index of the MM+Bulges.

    Returns:
        Callable: Sorting criteria function based on the provided parameters.
    """

    criteria = sorting_criteria.split(",")
    if len(criteria) > 3:
        raise ValueError("Mismatching sorting criteria selected")
    if any(c not in SORTING_CRITERIA for c in criteria):
        raise ValueError("Unknown sorting criteria")
    if score:
        return sorting_score(criteria, score_idx, mmb_idx)
    return sorting_fewest(criteria, mmb_idx)


def split_targets(input_args: MergeTargets) -> None:
    """
    Splits and processes the input targets file into best and alternative targets
    files based on specified criteria.

    Args:
        input_args (MergeTargets): Object containing input arguments for target
            splitting and processing.

    Returns:
        None

    Raises:
        OSError: If an error occurs during the process of splitting and processing
            the targets.
    """

    # open best, and alternative targets files
    try:
        with open(input_args.targets_fname, mode="r") as infile, open(
            input_args.targets_fname_merged, mode="w"
        ) as outfile, open(input_args.targets_fname_discarded, mode="w") as discfile:
            # write header to outfiles
            header = infile.readline().strip()
            outfile.write(f"{header}\n")
            discfile.write(f"{header}\n")
            # begin dividing targets
            (
                pos_prev,
                guide_prev,
                chrom_prev,
                snp_prev,
                cluster,
            ) = initialize_targets_cluster(input_args.rangebp)
            for line in infile:
                target = parse_target(line)
                guide, chrom, pos, snp = (
                    target[input_args.guideidx],
                    target[input_args.chromidx],
                    int(target[input_args.posidx]),
                    target[input_args.snpidx],
                )
                if open_targets_cluster(
                    guide,
                    guide_prev,
                    chrom,
                    chrom_prev,
                    pos,
                    pos_prev,
                    input_args.rangebp,
                ):
                    construct_cluster(
                        cluster,
                        input_args.guideidx,
                        input_args.snpidx,
                        input_args.posidx,
                        input_args.scoreidx,
                        input_args.mmbidx,
                        input_args.sort_pivot,
                        input_args.sorting_criteria,
                        input_args.sorting_criteria_scoring,
                        outfile,
                        discfile,
                    )
                    cluster = [target]
                else:
                    cluster.append(target)
                guide_prev = guide
                pos_prev = pos
                chrom_prev = chrom
                snp_prev = snp
            construct_cluster(
                cluster,
                input_args.guideidx,
                input_args.snpidx,
                input_args.posidx,
                input_args.scoreidx,
                input_args.mmbidx,
                input_args.sort_pivot,
                input_args.sorting_criteria,
                input_args.sorting_criteria_scoring,
                outfile,
                discfile,
            )
    except IOError as e:
        raise OSError("An error occurred while merging contiguous targets") from e


def merge_targets() -> None:
    """
    Merges targets by parsing input arguments, splitting the targets, and
    displaying the completion time.

    Returns:
        None
    """

    start = time()  # start time point
    input_args = parse_input_args(sys.argv[1:])
    split_targets(input_args)
    sys.stdout.write(f"Merge completed in {(time() - start)}s\n")


if __name__ == "__main__":
    merge_targets()
