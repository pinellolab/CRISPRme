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


from typing import List, Any, Tuple, Dict, Callable
from io import TextIOWrapper
from time import time

import sys

INPUT_ARG_COUNT = 12  # expected input args
FLAG_ALT_ONLY = 12  # target column for ALT target
SORTING_CRITERIA = {"mm+bulges": 0, "mm": 2, "bulges": 1}


def parse_input_args(args: List[str]) -> List[Any]:
    """
    Parses the input arguments and returns a list of parsed values.

    Args:
    args (List[str]): The input arguments to be parsed.

    Returns:
    List[Any]: A list of parsed values.

    Raises:
    ValueError: If no input argument is provided or if the number of input arguments is not equal to INPUT_ARG_COUNT.
    ValueError: If the input argument types are invalid.

    Examples:
    >>> parse_input_args(['input.txt', 'output.txt', '10', '2', '3', '4', '5', '6', '7', '8', 'score', 'mm', 'mm+bulges,mm'])
    ['input.txt', 'output.txt', 10, 1, 2, 3, 4, 5, 6, 'sort', 'mm', 'mm+bulges,mm']
    """

    if not args:
        raise ValueError("No input argument provided")
    if len(args) != INPUT_ARG_COUNT:
        raise ValueError(f"Expected {INPUT_ARG_COUNT} input arguments, got {len(args)}")
    try:
        # Use more descriptive variable names
        (
            infname,
            outfname,
            rangebp,
            chrom_idx,
            position_idx,
            mm_bul_count_idx,
            guide_idx,
            snp_info_idx,
            score_idx,
            sort_criterion,
            sorting_criteria_scoring,
            sorting_criteria,
        ) = (
            args[0],
            args[1],
            int(args[2]),
            int(args[3]) - 1,
            int(args[4]) - 1,
            int(args[5]) - 1,
            int(args[6]) - 1,
            int(args[7]) - 1,
            int(args[8]) - 1,
            args[9],
            args[10],
            args[11],
        )
    except (ValueError, TypeError) as e:
        raise ValueError("Invalid input argument types") from e
    return [
        infname,
        outfname,
        rangebp,
        chrom_idx,
        position_idx,
        mm_bul_count_idx,
        guide_idx,
        snp_info_idx,
        score_idx,
        sort_criterion,
        sorting_criteria_scoring,
        sorting_criteria,
    ]


def read_raw_targets(targets_file: str) -> Tuple[str, List[List[str]]]:
    """
    Reads and parses a targets file.

    Args:
        targets_file (str): The path to the targets file.

    Returns:
        Tuple[str, List[List[str]]]: A tuple containing the header of the targets file and the parsed data.

    Raises:
        OSError: If parsing the targets file fails.

    Examples:
        >>> read_raw_targets('targets.txt')
        ('Header line', [['data1', 'data2'], ['data3', 'data4']])
    """

    try:
        with open(targets_file, mode="r") as infile:
            header = infile.readline()  # recover targets file header
            targets_data = [line.split() for line in infile]  # parse line fields
    except OSError as e:
        raise OSError(f"Parsing targets file {targets_file} failed!") from e
    return header, targets_data


def open_targets_cluster(
    guide: str,
    prev_guide: str,
    chrom: str,
    prev_chrom: str,
    pos: int,
    prev_pos: int,
    rangebp: int,
) -> bool:
    """
    Determines whether to open a new targets cluster based on the given conditions.

    Args:
        guide (str): The current guide.
        prev_guide (str): The previous guide.
        chrom (str): The current chromosome.
        prev_chrom (str): The previous chromosome.
        pos (int): The current position.
        prev_pos (int): The previous position.
        rangebp (int): The maximum range in base pairs for samples to be considered contiguous.

    Returns:
        bool: True if a new targets cluster should be opened, False otherwise.
    """

    # if the condition is satisfied open a new tragets cluster
    return guide != prev_guide or chrom != prev_chrom or pos - prev_pos > rangebp


def merge_targets_by_snp(
    cluster: List[List[str]], snp_info_idx: int, pos_idx: int, guide_idx: int
) -> Tuple[List[List[str]], Dict[Tuple[str, str], List[List[str]]]]:
    """
    Merges targets in a cluster based on SNP information.

    Args:
        cluster (List[List[str]]): The cluster of targets to be merged.
        snp_info_idx (int): The index of the SNP info column in the target data.
        pos_idx (int): The index of the position column in the target data.
        guide_idx (int): The index of the guide column in the target data.

    Returns:
        Tuple[List[List[str]], Dict[Tuple[str, str], List[List[str]]]]: A tuple containing the merged reference targets and a dictionary of merged variant targets.

    Raises:
        ValueError: If the cluster is empty.

    Examples:
        >>> merge_targets_by_snp([['target1', 'n', 'pos1'], ['target2', 'A', 'pos2'], ['target3', 'A', 'pos1']], 1, 2, 3)
        ([['target1', 'n', 'pos1']], {('pos2', 'A'): [['target2', 'A', 'pos2']], ('pos1', 'A'): [['target3', 'A', 'pos1']]}])
    """

    if not cluster:
        raise ValueError("Empty cluster, nothing to merge")
    reference_targets = []  # REF targets list
    variants = {}  # dictionary storing ALT targets positions and info
    for target in cluster:
        if target[snp_info_idx] == "n":  # reference target
            reference_targets.append(target)
        else:  # ALT target
            pos_t, snp_info_t = target[pos_idx], target[snp_info_idx]
            if (pos_t, snp_info_t) in variants:  # merge samples with identical targets
                variants[(pos_t, snp_info_t)][0][guide_idx - 2] += (
                    "," + target[guide_idx - 2]
                )
                variants[(pos_t, snp_info_t)][0][snp_info_idx - 2] += (
                    "," + target[snp_info_idx - 2]
                )
                variants[(pos_t, snp_info_t)][0][snp_info_idx - 1] += (
                    "," + target[snp_info_idx - 1]
                )
            else:  # add new ALT target
                variants[(pos_t, snp_info_t)] = [target]
    return reference_targets, variants


def recover_alt_targets(
    variants: Dict[Tuple[str, str], List[List[str]]], noref: bool
) -> List[List[str]]:
    """
    Recovers alternative targets from a dictionary of variants.

    Args:
        variants (Dict[Tuple[str, str], List[List[str]]]): A dictionary containing merged variant targets.
        noref (bool): Flag indicating whether to mark non-reference targets.

    Returns:
        List[List[str]]: A list of recovered alternative targets.

    Examples:
        >>> recover_alt_targets({('pos2', 'A'): [['target2', 'A', 'pos2']], ('pos1', 'A'): [['target3', 'A', 'pos1']]}, True)
        [['target2', 'A', 'pos2'], ['target3', 'A', 'pos1']]
    """

    alternative_targets = []
    for value in variants.values():
        for target in value:
            target[FLAG_ALT_ONLY] = "y" if noref else target[FLAG_ALT_ONLY]
            alternative_targets.append(target)
    return alternative_targets


def remove_duplicate_data(target: List[str], idx: int) -> str:
    """
    Removes duplicate data from a target list.

    Args:
        target (List[str]): The target list.
        idx (int): The index of the target element to remove duplicates from.

    Returns:
        str: The target element with duplicates removed.

    Examples:
        >>> remove_duplicate_data(['target1', 'A,A,A,A,A,A,A,A,A'], 1)
        'A'
    """

    return ",".join(set(target[idx].split(",")))


def remove_duplicate_alt_targets(
    targets_alt: List[List[str]], redundant_idxs: List[int]
) -> List[List[str]]:
    """
    Removes duplicate data from alternative targets.

    Args:
        targets_alt (List[List[str]]): The list of alternative targets.
        redundant_idxs (List[int]): The list of indices indicating the elements to remove duplicates from.

    Returns:
        List[List[str]]: The list of alternative targets with duplicates removed.

    Examples:
        >>> remove_duplicate_alt_targets([['target1', 'A,A,A,A,A,A,A,A,A'], ['target2', 'C,C,C,C,C,C,C,C,C']], [1])
        [['target1', 'A'], ['target2', 'C']]
    """

    targets_alt_polished = []  # remove duplicate data from ALT target
    for target in targets_alt:
        for idx in redundant_idxs:
            target[idx] = remove_duplicate_data(target, idx)
        targets_alt_polished.append(target)
    assert len(targets_alt_polished) == len(targets_alt)
    return targets_alt_polished


def sorting_score(
    criteria: List[str], score_idx: int, mm_bul_count_idx: int
) -> Callable:
    """
    Creates a sorting key function based on the given criteria. Score has highest
    priority and is always included in the sorting function

    Args:
        criteria (List[str]): The sorting criteria.
        score_idx (int): The index of the score column in the data.
        mm_bul_count_idx (int): The index of the mismatch and bulge count column in the data.

    Returns:
        Callable: A sorting key function that can be used with the `sorted` function.

    Examples:
        >>> sorting_score(['mm+bulges', 'mm'], 8, 5)
        <function sorting_score.<locals>.<lambda> at 0x7f1234567890>
    """

    if len(criteria) == 1:  # one criterion
        return lambda x: (
            -float(x[score_idx]),
            int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[0]]]),
        )
    elif len(criteria) == 2:
        return lambda x: (
            -float(x[score_idx]),
            int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[0]]]),
            int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[1]]]),
        )
    # base case (all three )
    return lambda x: (
        -float(x[score_idx]),
        int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[0]]]),
        int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[1]]]),
        int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[2]]]),
    )


def sorting_fewest(criteria: List[str], mm_bul_count_idx: int) -> Callable:
    """
    Creates a sorting key function based on the input criteria. The sorting
    criteria include mismatches, bulges, and mismatches+bulges

    Args:
        criteria (List[str]): The sorting criteria.
        mm_bul_count_idx (int): The index of the mismatch and bulge count column in the data.

    Returns:
        Callable: A sorting key function that can be used with the `sorted` function.

    Examples:
        >>> sorting_fewest(['mm+bulges', 'mm'], 5)
        <function sorting_fewest.<locals>.<lambda> at 0x7f1234567890>
    """

    if len(criteria) == 1:  # one criterion
        return lambda x: (int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[0]]]))
    elif len(criteria) == 2:
        return lambda x: (
            int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[0]]]),
            int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[1]]]),
        )
    # base case (all three )
    return lambda x: (
        int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[0]]]),
        int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[1]]]),
        int(x[mm_bul_count_idx - SORTING_CRITERIA[criteria[2]]]),
    )


def define_sorting_criteria(
    sorting_criteria: str, score: bool, score_idx: int, mm_bul_count_idx: int
) -> Callable:
    """
    Defines the sorting criteria based on the given parameters.

    Args:
        sorting_criteria (str): The sorting criteria as a comma-separated string.
        score (bool): Flag indicating whether to prioritize sorting by score.
        score_idx (int): The index of the score column in the data.
        mm_bul_count_idx (int): The index of the mismatch and bulge count column in the data.

    Returns:
        Callable: A sorting key function based on the defined criteria.

    Raises:
        ValueError: If the number of sorting criteria exceeds 3.
        ValueError: If unknown sorting criteria are provided.

    Examples:
        >>> define_sorting_criteria('mm+bulges,mm', True, 8, 5)
        <function sorting_score.<locals>.<lambda> at 0x7f1234567890>
    """

    criteria = sorting_criteria.split(",")
    # if score:
    #     print("score")
    # else:
    #     print("fewest")
    # print(criteria)
    if len(criteria) > 3:
        raise ValueError("Mismatching sorting criteria selected")
    if any(c not in SORTING_CRITERIA for c in criteria):
        raise ValueError("Unknown sorting criteria")
    if score:  # sort by criteria (priority to score)
        return sorting_score(criteria, score_idx, mm_bul_count_idx)
    return sorting_fewest(criteria, mm_bul_count_idx)  # sort by criteria


def report_best_targets(
    targets: List[List[str]], score_idx: int, outfile: TextIOWrapper
) -> List[List[str]]:
    """
    Reports the best targets and writes them to an output file.

    Args:
        targets (List[List[str]]): The list of targets.
        score_idx (int): The index of the score column in the data.
        outfile (TextIOWrapper): The output file to write the best target.

    Returns:
        List[List[str]]: The remaining targets after removing the best target.

    Examples:
        >>> report_best_targets([['target1', 'A', 'score1'], ['target2', 'C', 'score2']], 2, outfile)
        [['target2', 'C', 'score2']]
    """

    # count remaining targets
    targets[0][score_idx - 1] = str(len(targets) - 1)
    reportline = "\t".join(targets[0])  # write best target to merge file
    outfile.write(f"{reportline}\n")
    best = targets.pop(0)  # pop best target from list
    return targets


def report_best_targets_ref_alt(
    targets_ref: List[List[str]],
    targets_alt: List[List[str]],
    score: bool,
    mm_bul_count_idx: int,
    score_idx: int,
    outfile: TextIOWrapper,
) -> Tuple[List[List[str]], List[List[str]]]:
    """
    Reports the best targets (reference and alternative) and writes them to an output file.

    Args:
        targets_ref (List[List[str]]): The list of reference targets.
        targets_alt (List[List[str]]): The list of alternative targets.
        score (bool): Flag indicating whether to prioritize sorting by score.
        mm_bul_count_idx (int): The index of the mismatch and bulge count column in the target data.
        score_idx (int): The index of the score column in the target data.
        outfile (TextIOWrapper): The output file to write the best targets.

    Returns:
        Tuple[List[List[str]], List[List[str]]]: A tuple containing the remaining reference targets and the remaining alternative targets.

    Examples:
        >>> report_best_targets_ref_alt([['target1', 'A', 'score1']], [['target2', 'C', 'score2']], True, 6, 5, outfile)
    """

    # compare values (score or mm+bulges) to determine the best target
    if score:  # compare scores
        if float(targets_ref[0][score_idx]) >= float(targets_alt[0][score_idx]):
            targets_ref[0][score_idx - 1] = str(
                len(targets_ref) + len(targets_alt) - 1
            )  # recover remaining targets
            reportline = "\t".join(targets_ref[0])
            outfile.write(f"{reportline}\n")
            best = targets_ref.pop(0)  # remove best target
        else:
            targets_alt[0][score_idx - 1] = str(
                len(targets_ref) + len(targets_alt) - 1
            )  # recover remaining targets
            reportline = "\t".join(targets_alt[0])
            outfile.write(f"{reportline}\n")
            best = targets_alt.pop(0)  # remove best target
    elif int(targets_ref[0][mm_bul_count_idx]) <= int(targets_alt[0][mm_bul_count_idx]):
        targets_ref[0][score_idx - 1] = str(
            len(targets_ref) + len(targets_alt) - 1
        )  # recover remaining targets
        reportline = "\t".join(targets_ref[0])
        outfile.write(f"{reportline}\n")
        best = targets_ref.pop(0)  # remove best target
    else:
        targets_alt[0][score_idx - 1] = str(
            len(targets_ref) + len(targets_alt) - 1
        )  # recover remaining targets
        reportline = "\t".join(targets_alt[0])
        outfile.write(f"{reportline}\n")
        best = targets_alt.pop(0)  # remove best target
    return targets_ref, targets_alt


def report_remaining_targets(
    targets_ref: List[List[str]],
    targets_alt: List[List[str]],
    score_idx: int,
    outfile: TextIOWrapper,
) -> None:
    """
    Reports the discarded targets (reference and alternative) and writes them to an output file.

    Args:
        targets_ref (List[List[str]]): The list of remaining reference targets.
        targets_alt (List[List[str]]): The list of remaining alternative targets.
        score_idx (int): The index of the score column in the data.
        outfile (TextIOWrapper): The output file to write the remaining targets.

    Returns:
        None

    Examples:
        >>> report_remaining_targets([['target1', 'A', 'score1']], [['target2', 'C', 'score2']], 3, outfile)
    """

    for i, target in enumerate(targets_ref):
        targets_ref[i][score_idx - 1] = str(len(targets_ref) + len(targets_alt) - 1)
        reportline = "\t".join(target)
        outfile.write(f"{reportline}\n")
    for i, target in enumerate(targets_alt):
        targets_alt[i][score_idx - 1] = str(len(targets_ref) + len(targets_alt) - 1)
        reportline = "\t".join(target)
        outfile.write(f"{reportline}\n")


def recover_best_targets(
    cluster: List[List[str]],
    snp_info_idx: int,
    pos_idx: int,
    guide_idx: int,
    score_idx: int,
    mm_bul_count_idx: int,
    criterion: str,
    sorting_criteria_scoring: str,
    sorting_criteria: str,
    outfile: TextIOWrapper,
    outfile_discarded: TextIOWrapper,
) -> None:
    """
    Recovers the best targets from a cluster and writes them to an output file.

    Args:
        cluster (List[List[str]]): The cluster of targets.
        snp_info_idx (int): The index of the SNP info column in the target data.
        pos_idx (int): The index of the position column in the target data.
        guide_idx (int): The index of the guide column in the target data.
        score_idx (int): The index of the score column in the target data.
        mm_bul_count_idx (int): The index of the mismatch and bulge count column in the target data.
        criterion (str): The criterion for selecting the best targets (score or fewest).
        sorting_criteria_scoring (str): The sorting criteria for the targets (scoring has highest priority).
        sorting_criteria (str): The sorting criteria for the targets.
        outfile (TextIOWrapper): The output file to write the best targets.
        outfile_discarded (TextIOWrapper): The output file to write the discarded targets.

    Returns:
        None

    Examples:
        >>> recover_best_targets(cluster, 2, 3, 4, 5, 6, 'score', 'mm+bulges,mm', outfile, outfile_discarded)
    """

    if not cluster:
        return  # avoids potential crash at first iteration when opened new cluster
    targets_ref, targets_alt_dict = merge_targets_by_snp(
        cluster, snp_info_idx, pos_idx, guide_idx
    )
    noref = not targets_ref  # only ALT targets
    targets_alt = remove_duplicate_alt_targets(
        recover_alt_targets(targets_alt_dict, noref),
        [snp_info_idx, snp_info_idx - 2, snp_info_idx - 1, guide_idx - 2],
    )  # recover ALT targets
    # sort by CFD or CRISTA score (no impact when using non Cas9 proteins)
    if criterion == "score":
        # recover sorting criteria order (score has priority)
        criteria = define_sorting_criteria(
            sorting_criteria_scoring, True, score_idx, mm_bul_count_idx
        )
    else:  # sort on mm and bulges only (fewest)
        # recover sorting criteria order (no score)
        criteria = define_sorting_criteria(
            sorting_criteria, False, score_idx, mm_bul_count_idx
        )
    if targets_ref:
        targets_ref = sorted(targets_ref, key=criteria)
    if targets_alt:
        targets_alt = sorted(targets_alt, key=criteria)
    # recover best targets
    if noref:  # only alt targets
        targets_alt = report_best_targets(targets_alt, score_idx, outfile)
    elif targets_ref and targets_alt:  # both ref and alt targets
        score = criterion == "score"
        targets_ref, targets_alt = report_best_targets_ref_alt(
            targets_ref, targets_alt, score, mm_bul_count_idx, score_idx, outfile
        )
    else:  # only ref targets
        targets_ref = report_best_targets(targets_ref, score_idx, outfile)
    report_remaining_targets(targets_ref, targets_alt, score_idx, outfile_discarded)


def merge_targets() -> None:
    """
    Merges targets based on specified criteria and writes the merged targets to output files.

    Returns:
        None

    Raises:
        OSError: If the targets merge fails.

    Examples:
        >>> merge_targets()
    """

    start = time()  # merging start time
    input_args = parse_input_args(sys.argv[1:])
    # recover raw targets file header and data
    header, targets_data_raw = read_raw_targets(input_args[0])
    # define fnames for merged target files (kept and discarded)
    targets_file_merged = input_args[1]
    targets_file_merged_discarded = f"{targets_file_merged}.discarded_samples"
    # initialize target clustering values
    prev_pos, prev_guide, prev_chrom = -(input_args[2] + 1), "", ""
    cluster = []  # targets clusters
    try:
        with open(targets_file_merged, mode="w") as outfile:
            with open(targets_file_merged_discarded, mode="w") as outfile_discarded:
                # write headers to target merged files
                outfile.write(header)
                outfile_discarded.write(header)
                for target_data in targets_data_raw:
                    if open_targets_cluster(
                        target_data[input_args[6]],
                        prev_guide,
                        target_data[input_args[3]],
                        prev_chrom,
                        int(target_data[input_args[4]]),
                        prev_pos,
                        input_args[2],
                    ):
                        recover_best_targets(
                            cluster,
                            input_args[7],
                            input_args[4],
                            input_args[6],
                            input_args[8],
                            input_args[5],
                            input_args[9],
                            input_args[10],
                            input_args[11],
                            outfile,
                            outfile_discarded,
                        )
                        cluster = [target_data]
                    else:
                        cluster.append(target_data)
                    prev_guide, prev_pos, prev_chrom = (
                        target_data[input_args[6]],
                        int(target_data[input_args[4]]),
                        target_data[input_args[3]],
                    )
                recover_best_targets(
                    cluster,
                    input_args[7],
                    input_args[4],
                    input_args[6],
                    input_args[8],
                    input_args[5],
                    input_args[9],
                    input_args[10],
                    input_args[11],
                    outfile,
                    outfile_discarded,
                )
    except OSError as e:
        raise OSError(f"Targets merge failed for {input_args[0]}") from e
    sys.stdout.write(f"Merge completed in: {(time() - start):.2f}s\n")


if __name__ == "__main__":
    merge_targets()
