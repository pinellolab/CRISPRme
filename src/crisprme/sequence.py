"""
"""

from utils import exception_handler

from typing import Tuple

import pybedtools

def extract_sequence(seqname: str, coord: str, genome: str, debug: bool) -> Tuple[str, str]:
    """Retrieve the sequence corresponding to the provided genomic coordinates, 
    given in BED format

    :param seqname: sequence name
    :type seqname: str
    :param coord: BED like genomic coordinate
    :type coord: str
    :param genome: reference genome
    :type genome: str
    :param debug: debug mode
    :type debug: bool
    :return: extracted sequence
    :rtype: Tuple[str, str]
    """
    if not isinstance(seqname, str):
        exception_handler(
            TypeError, f"Expected {str.__name__}, got {type(seqname).__name__}", debug
        )
    if not isinstance(coord, str):
        exception_handler(
            TypeError, f"Expected {str.__name__}, got {type(coord).__name__}", debug
        )
    seqname = "_".join(seqname.split())  # remove blank or tabs if any
    coordinate = pybedtools.BedTool(f"{coord}", from_string=True)  # initialize the object
    # extract the sequence
    sequence = coordinate.sequence(fi=genome).strip()
    return seqname, sequence

def recover_guides(sequence: str, pam: str, guide_length: int, pam_at_beginning: bool, debug: bool) -> None:
    pass


