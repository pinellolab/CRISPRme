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
    """Parses and validates command line arguments for the annotation script.

    Ensures the correct number of arguments are provided and that input files exist. 
    Returns the file paths for the offtargets, annotation, and output files.

    Args:
        args: List of command line arguments.

    Returns:
        Tuple containing the offtargets file name, annotation file name, and output 
            file name.

    Raises:
        ValueError: If the number of arguments is incorrect.
        FileNotFoundError: If any of the specified files do not exist.
    """
    if len(args) != 3:  # check input arguments consistency
        raise ValueError("Too many/few input arguments for annotation script")
    offtargets_fname = args[0]  # offtargets table report
    if not os.path.isfile(offtargets_fname):
        raise FileNotFoundError(f"Cannot find off-targets file {offtargets_fname}")
    annotation_fname = args[1]
    if os.path.basename(annotation_fname) != "vuoto.txt" and not os.path.isfile(annotation_fname):
        raise FileNotFoundError(
            f"Cannot find annotation files {annotation_fname}"
        )
    offtargets_out_fname = args[2]  # annotated offtargets table report
    assert offtargets_out_fname
    return offtargets_fname, annotation_fname, offtargets_out_fname


def load_annotation_bed(annotation_fname: str) -> TabixFile:
    """Loads a BED annotation file and ensures a Tabix index is present.

    Checks for the existence of a Tabix index for the given annotation BED file 
    and creates one if necessary. Returns a TabixFile object for querying the 
    annotation data.

    Args:
        annotation_fname: Path to the annotation BED file.

    Returns:
        TabixFile object for the annotation BED file.

    Raises:
        SamtoolsError: If the annotation BED file cannot be loaded or indexed.
    """
    # check that tabix index is available for all annotation bed
    pysam.tabix_index(annotation_fname, force=True, preset="bed")
    try:  # return tabix indexes for each annotation bed
        return pysam.TabixFile(annotation_fname)
    except (SamtoolsError, Exception) as e:
        raise SamtoolsError(
            "An error occurred while loading Annotation BED files"
        ) from e


def compute_target_coords(fields: List[str]) -> Tuple[str, int, int]:
    """Computes the genomic coordinates for a target from a list of fields.

    Returns the chromosome, start, and stop positions for the target based on 
    the provided fields.

    Args:
        fields: List of string fields containing target information.

    Returns:
        Tuple containing the chromosome (str), start (int), and stop (int) 
            positions.
    """
    start = int(fields[5])  # retrieve target start position
    stop = start + len(fields[1].replace("-", ""))  # compute target stop
    return fields[4], start, stop


def annotate_target(chrom: str, start: int, stop: int, annotation: TabixFile) -> str:
    """Annotates a genomic target using a Tabix-indexed annotation file.

    Retrieves and returns a comma-separated list of annotation features overlapping 
    the specified genomic region.

    Args:
        chrom: Chromosome name.
        start: Start position of the target.
        stop: Stop position of the target.
        annotation: TabixFile object for annotation data.

    Returns:
        Comma-separated string of annotation features for the target.
    """
    if chrom not in annotation.contigs:
        return "n"  # contig not in input annotation file 
    target_anns = {
        feature.split()[3]
        for feature in annotation.fetch(chrom, start, stop)
    }  # retrieve annotations for current target
    return ",".join(sorted(target_anns))  # report as comma-separated list


def annotate_offtargets(
    offtargets_fname: str, annotation: Optional[TabixFile], offtargets_out_fname: str,
) -> None:
    """Annotates a genomic target using a Tabix-indexed annotation file.

    Retrieves and returns a comma-separated list of annotation features overlapping 
    the specified genomic region.

    Args:
        chrom: Chromosome name.
        start: Start position of the target.
        stop: Stop position of the target.
        annotation: TabixFile object for annotation data.

    Returns:
        Comma-separated string of annotation features for the target.
    """
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
                    "n" if annotation is None else annotate_target(
                        *compute_target_coords(fields), annotation
                    )
                )
                offtarget_annotated = "\t".join(fields)
                outfile.write(f"{offtarget_annotated}\n")  # write annotated offtarget
    except (IOError, Exception) as e:
        raise OSError(f"Annotation failed on off-targets in {offtargets_fname}") from e


def main() -> None:
    """Annotates a list of off-targets using an optional annotation dataset.

    Reads off-targets from a file, annotates each entry with features from the 
    provided annotation dataset, and writes the results to an output file.

    Args:
        offtargets_fname: Path to the input file containing off-targets.
        annotation: TabixFile object for annotation data, or None if no annotation.
        offtargets_out_fname: Path to the output file for annotated off-targets.

    Raises:
        OSError: If annotation fails due to file I/O or processing errors.
    """
    # parse input command line arguments
    offtargets_fname, annotation, offtargets_out_fname = parse_commandline(sys.argv[1:])
    start = time()  # track annotation time
    empty = os.path.basename(annotation) == "vuoto.txt"
    # load annotation bed files
    annotation = None if empty else load_annotation_bed(annotation)
    # annotate offtargets using input bed files
    annotate_offtargets(offtargets_fname, annotation, offtargets_out_fname)
    sys.stdout.write(f"Annotation completed in {time() - start:.2f}s\n")


if __name__ == "__main__":
    main()