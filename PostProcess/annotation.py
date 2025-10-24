""" """

from pysam import TabixFile
from pysam.utils import SamtoolsError
from typing import List, Tuple, Optional
from time import time

import pysam
import sys
import os

TBI = "tbi"

def parse_commandline(args: List[str]) -> Tuple[str, str, str]:
    if len(args) != 3:  # check input arguments consistency
        raise ValueError("Too many/few input arguments for annotation script")
    offtargets_fname = args[0]  # offtargets table report
    if not os.path.isfile(offtargets_fname):
        raise FileNotFoundError(f"Cannot find off-targets file {offtargets_fname}")
    annotation_fname = args[1]
    if annotation_fname != "vuoto.txt" and not os.path.isfile(annotation_fname):
        raise FileNotFoundError(
            f"Cannot find annotation files {annotation_fname}"
        )
    offtargets_out_fname = args[2]  # annotated offtargets table report
    assert offtargets_out_fname
    return offtargets_fname, annotation_fname, offtargets_out_fname


def load_annotation_bed(annotation_fname: str) -> TabixFile:
    # check that tabix index is available for all annotation bed
    if not os.path.isfile(f"{annotation_fname}.{TBI}"):
        pysam.tabix_index(annotation_fname, force=True, preset="bed")
    try:  # return tabix indexes for each annotation bed
        return pysam.TabixFile(annotation_fname)
    except (SamtoolsError, Exception) as e:
        raise SamtoolsError(
            "An error occurred while loading Annotation BED files"
        ) from e


def compute_target_coords(fields: List[str]) -> Tuple[str, int, int]:
    start = int(fields[5])  # retrieve target start position
    stop = start + len(fields[1].replace("-", ""))  # compute target stop
    return fields[4], start, stop


def annotate_target(
    chrom: str,
    start: int,
    stop: int,
    annotation: TabixFile,
) -> str:
    target_anns = {
        feature.split()[3]
        for feature in annotation.fetch(chrom, start, stop)
    }  # retrieve annotations for current target
    return ",".join(sorted(target_anns))  # report as comma-separated list


def annotate_offtargets(
    offtargets_fname: str, annotation: Optional[TabixFile], empty: bool, offtargets_out_fname: str,
) -> None:
    try:
        with open(offtargets_fname, mode="r") as infile, open(
            offtargets_out_fname, mode="w"
        ) as outfile:
            header = infile.readline().strip()  # preserve header in annotated report
            outfile.write(f"{header}\n")
            for offtarget in infile:  # read offtargets
                fields = offtarget.split()  # split offtarget individual fields
                # annotate current off-target using input annotation datasets
                assert annotation
                fields[14] = (
                    "n" if empty else annotate_target(
                        *compute_target_coords(fields), annotation
                    )
                )
                offtarget_annotated = "\t".join(fields)
                outfile.write(f"{offtarget_annotated}\n")  # write annotated offtarget
    except (IOError, Exception) as e:
        raise OSError(f"Annotation failed on off-targets in {offtargets_fname}") from e


def main() -> None:
    # parse input command line arguments
    offtargets_fname, annotation, offtargets_out_fname = parse_commandline(sys.argv[1:])
    start = time()  # track annotation time
    empty = os.path.basename(annotation) == "vuoto.txt"
    # load annotation bed files
    annotation = None if empty else load_annotation_bed(annotation)
    # annotate offtargets using input bed files
    annotate_offtargets(offtargets_fname, annotation, empty, offtargets_out_fname)
    sys.stdout.write(f"Annotation completed in {time() - start:.2f}s\n")


if __name__ == "__main__":
    main()