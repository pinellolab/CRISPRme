"""
"""

from typing import List, Any, Tuple, Dict, Callable
from io import TextIOWrapper
from time import time

import sys

INPUT_ARG_COUNT = 10  # expected input args
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
    >>> parse_input_args(['input.txt', 'output.txt', '10', '2', '3', '4', '5', '6', '7', '8', 'sort'])
    ['input.txt', 'output.txt', 10, 1, 2, 3, 4, 5, 6, 'sort', 'mm+bulges,mm']
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
        )
        sorting_criteria = "mm+bulges,mm"
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
        sorting_criteria,
    ]


def read_raw_targets(targets_file: str) -> Tuple[str, List[List[str]]]:
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
    # if the condition is satisfied open a new tragets cluster
    return guide != prev_guide or chrom != prev_chrom or pos - prev_pos > rangebp


def merge_targets_by_snp(
    cluster: List[List[str]], snp_info_idx: int, pos_idx: int, guide_idx: int
) -> Tuple[List[List[str]], Dict[Tuple[str, str], List[List[str]]]]:
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
    alternative_targets = []
    for value in variants.values():
        for target in value:
            target[FLAG_ALT_ONLY] = "y" if noref else target[FLAG_ALT_ONLY]
            alternative_targets.append(target)
    return alternative_targets


def remove_duplicate_data(target: List[str], idx: int) -> List[str]:
    target[idx] = ",".join(set(target[idx].split(",")))
    return target


def remove_duplicate_alt_targets(
    targets_alt: List[List[str]], redundant_idxs: List[int]
):
    targets_alt_polished = []  # remove duplicate data from ALT target
    for target in targets_alt:
        for idx in redundant_idxs:
            target = remove_duplicate_data(target, idx)
        targets_alt_polished.append(target)
    assert len(targets_alt_polished) == len(targets_alt)
    return targets_alt_polished


def sorting_score(
    criteria: List[str], score_idx: int, mm_bul_count_idx: int
) -> Callable:
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
    criteria = sorting_criteria.split(",")
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
    # count remaining targets
    targets[0][score_idx - 1] = str(len(targets) - 1)
    outfile.write("\t".join(targets[0]))  # write best target to merge file
    best = targets.pop(0)  # pop best target from list
    return targets


def report_best_targets_ref_alt(
    targets_ref: List[List[str]],
    targets_alt: List[List[str]],
    cmp_idx: int,
    score_idx: int,
    outfile: TextIOWrapper,
) -> Tuple[List[List[str]], List[List[str]]]:
    # compare values (score or mm+bulges) to determine the best target
    if float(targets_ref[0][cmp_idx]) >= float(targets_alt[0][cmp_idx]):
        targets_ref[0][score_idx - 1] = str(
            len(targets_ref) + len(targets_alt) - 1
        )  # recover remaining targets
        outfile.write("\t".join(targets_ref[0]))
        best = targets_ref.pop(0)  # remove best target
    else:
        targets_alt[0][score_idx - 1] = str(
            len(targets_ref) + len(targets_alt) - 1
        )  # recover remaining targets
        outfile.write("\t".join(targets_alt[0]))
        best = targets_alt.pop(0)  # remove best target
    return targets_ref, targets_alt


def report_remaining_targets(
    targets_ref: List[List[str]],
    targets_alt: List[List[str]],
    score_idx: int,
    outfile: TextIOWrapper,
) -> None:
    for i, target in enumerate(targets_ref):
        targets_ref[i][score_idx - 1] = str(len(targets_ref) + len(targets_alt) - 1)
        outfile.write("\t".join(target))
    for i, target in enumerate(targets_alt):
        targets_alt[i][score_idx - 1] = str(len(targets_ref) + len(targets_alt) - 1)
        outfile.write("\t".join(target))


def recover_best_targets(
    cluster: List[List[str]],
    snp_info_idx: int,
    pos_idx: int,
    guide_idx: int,
    score_idx: int,
    mm_bul_count_idx: int,
    criterion: str,
    sorting_criteria: str,
    outfile: TextIOWrapper,
    outfile_discarded: TextIOWrapper,
):
    if not cluster:
        return  # avoids potential crash at first iteration when opened new cluster
    targets_ref, targets_alt_dict = merge_targets_by_snp(
        cluster, snp_info_idx, pos_idx, guide_idx
    )
    noref = not targets_ref  # only ALT targets
    targets_alt = remove_duplicate_alt_targets(
        recover_alt_targets(targets_alt_dict, noref),
        [snp_info_idx, snp_info_idx - 2, snp_info_idx - 1, guide_idx],
    )  # recover ALT targets
    # sort by CFD or CRISTA score (no impact when using non Cas9 proteins)
    if criterion == "score":
        # recover sorting criteria order (score has priority)
        criteria = define_sorting_criteria(
            sorting_criteria, True, score_idx, mm_bul_count_idx
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
        cmp_idx = score_idx if criterion == "score" else mm_bul_count_idx
        targets_ref, targets_alt = report_best_targets_ref_alt(
            targets_ref, targets_alt, cmp_idx, score_idx, outfile
        )
    else:  # only ref targets
        targets_ref = report_best_targets(targets_ref, score_idx, outfile)
    report_remaining_targets(targets_ref, targets_alt, score_idx, outfile_discarded)


def merge_targets():
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
                            outfile,
                            outfile_discarded,
                        )
                        cluster = [target_data]
                    else:
                        cluster.append(target_data)
                    prev_guide, prev_pos, prev_chrom = (
                        target_data[input_args[6]],
                        int(target_data[input_args[5]]),
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
                    outfile,
                    outfile_discarded,
                )
    except OSError as e:
        raise OSError(f"Targets merge failed for {input_args[0]}") from e
    sys.stdout.write(f"Merge completed in: {(time() - start):.2f}s\n")


if __name__ == "__main__":
    merge_targets()
