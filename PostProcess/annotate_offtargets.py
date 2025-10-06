""" """

from pysam import TabixFile
from pysam.utils import SamtoolsError
from typing import List, Tuple
from time import time

import pysam
import sys
import os

TBI = "tbi"


def parse_commandline(args: List[str]) -> Tuple[str, List[str], List[str], str]:
    if len(args) != 4:  # check input arguments consistency
        raise ValueError("Too many/few input arguments for annotation script")
    offtargets_fname = args[0]  # offtargets table report
    if not os.path.isfile(offtargets_fname):
        raise FileNotFoundError(f"Cannot find off-targets file {offtargets_fname}")
    annotation_fnames = list(args[1].split(","))
    if annotation_fnames[0] != "empty.txt" and not all(os.path.isfile(fname) for fname in annotation_fnames):
        raise FileNotFoundError(
            f"Cannot find annotation files {' '.join(annotation_fnames)}"
        )
    annotation_colnames = list(args[2].split(","))
    if len(annotation_colnames) != len(annotation_fnames):
        raise ValueError(
            f"Mismatching number of annotation files ({len(annotation_fnames)}) and annotation column names ({annotation_colnames})"
        )
    offtargets_out_fname = args[3]  # annotated offtargets table report
    assert offtargets_out_fname
    return (
        offtargets_fname,
        annotation_fnames,
        annotation_colnames,
        offtargets_out_fname,
    )


def load_annotation_bed(annotation_fnames: List[str]) -> List[TabixFile]:
    # check that tabix index is available for all annotation bed
    for fname in annotation_fnames:  # index bed with samtools
        pysam.tabix_index(fname, force=True, preset="bed")
    try:  # return tabix indexes for each annotation bed
        return [pysam.TabixFile(fname) for fname in annotation_fnames]
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
    annotations: List[TabixFile],
    annotations_colnames: List[str],
) -> str:
    target_anns = {
        f"{feature.split()[3]}__{annotations_colnames[i]}"
        for i, annotation in enumerate(annotations)
        for feature in annotation.fetch(chrom, start, stop)
    }  # retrieve annotations for current target
    return ",".join(sorted(target_anns))  # report as comma-separated list


def annotate_offtargets(
    offtargets_fname: str,
    annotations: List[TabixFile],
    annotations_colnames: List[str],
    empty: bool,
    offtargets_out_fname: str,
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
                fields[14] = (
                    "n" if empty else annotate_target(
                        *compute_target_coords(fields),
                        annotations,
                        annotations_colnames,
                    )
                )
                offtarget_annotated = "\t".join(fields)
                outfile.write(f"{offtarget_annotated}\n")  # write annotated offtarget
    except (IOError, Exception) as e:
        raise OSError(f"Annotation failed on off-targets in {offtargets_fname}") from e


def main() -> None:
    # parse input command line arguments
    offtargets_fname, annotations, annotations_colnames, offtargets_out_fname = (
        parse_commandline(sys.argv[1:])
    )
    start = time()  # track annotation time
    empty = len(annotations) == 1 and annotations[0] == "empty.txt"
    # load annotation bed files
    annotations = [] if empty else load_annotation_bed(annotations)
    # annotate offtargets using input bed files
    annotate_offtargets(
        offtargets_fname, annotations, annotations_colnames, empty, offtargets_out_fname
    )
    sys.stdout.write(f"Annotation completed in {time() - start:.2f}s\n")


if __name__ == "__main__":
    main()
