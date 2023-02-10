"""
"""

from sequence import extract_sequence
from utils import exception_handler

from typing import Tuple

import os

def parse_pam(pamfile: str, debug: bool) -> Tuple[str, int, bool]:
    """Recover the PAM sequence from the input PAM file

    :param pamfile: PAM file
    :type pamfile: str
    :param debug: debug mode
    :type debug: bool
    :return: PAM sequence and expected guide length
    :rtype: Tuple[str, int, bool]
    """
    try:
        with open(pamfile, mode="r") as infile:
            pam_split = infile.readline().strip().split()
    except OSError:
        exception_handler(OSError, "An error occurred while parsing the PAM file", debug)
    pam = pam_split[0]
    total_pam_length = len(pam)
    try:
        pam_index = int(pam_split[-1])
    except ValueError:
        exception_handler(
            ValueError, f"Unable to interpret PAM index ({pam_split[-1]})", debug
        )
    stop = abs(pam_index) if pam_index < 0 else -(pam_index)
    pam_at_beginning = pam_index < 0
    pam = pam[:stop] if pam_at_beginning else pam[stop:] 
    return pam, total_pam_length - len(pam), pam_at_beginning

def parse_sequence(sequencefile: str, genome: str, debug: bool) -> None:
    guides = []
    try:
        with open(sequencefile, mode="r") as infile:
            sequence_data = infile.read()  # read sequence file content
    except OSError:
        exception_handler(
            OSError, "AN error occurred while parsing the sequence file", debug
        )
    sequences = sequence_data.split(">")  # recover seqname and associated sequence
    for seqname_sequence in sequences:
        seqname, sequence = seqname_sequence.strip().split("\n")
        if "chr" in sequence:  # BED like coordinates
            coordinates = sequence.split("\n")  # search if several consecutive coords
            for coordinate in coordinates:
                coord = coordinate.split()  # expected chr start stop
                if len(coord) != 3:
                    exception_handler(
                        ValueError, "Expected BED like coordinates, but the requirement is not satisfied", debug
                    )
                coord = f"{coord[0]}\t{coord[1]}\t{coord[2]}"  # ensure they are tab-separated
                genomefile = os.path.join(genome, f"{coord[0]}.fa")
                assert os.path.isfile(genomefile)
                seqname, sequence = extract_sequence(seqname, coord, genomefile, debug)
