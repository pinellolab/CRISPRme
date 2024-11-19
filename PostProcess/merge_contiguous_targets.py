"""
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


def split_target_row(target_row: str, guideidx: int, chromidx: int, posidx: int):
    fields = target_row.strip().split()
    return fields[guideidx], fields[chromidx], int(fields[posidx])


def update_target_fields(
    target: List[str], fields: List[str], samplesidx: int, snpid_idx: int, afidx: int
) -> List[str]:
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
    # set apprpriate flag if no target in reference in current cluster
    target[12] = "y" if varonly else target[12]
    return target


def remove_duplicate_targets(
    target: List[str], snpidx: int, snpid_idx: int, afidx: int, samplesidx: int
):
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
):
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
):
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
):
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
