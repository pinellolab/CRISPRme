"""
"""

from sequence import extract_sequence, recover_guides
from utils import add_n, exception_handler

from typing import List, Tuple

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

def parse_sequence(sequencefile: str, genome: str, pam: str, guide_length: int, pam_at_beginning: bool, debug: bool) -> List[str]:
    """Extract the guide sequences from the input sequence file, which include 
    the specified PAM sequence

    :param sequencefile: sequences file
    :type sequencefile: str
    :param genome: genome directory
    :type genome: str
    :param pam: PAM
    :type pam: str
    :param guide_length: guide length
    :type guide_length: int
    :param pam_at_beginning: PAM occurs at guide beginning 
    :type pam_at_beginning: bool
    :param debug: debug mode
    :type debug: bool
    :return: guides 
    :rtype: List[str]
    """
    guides = []
    try:
        with open(sequencefile, mode="r") as infile:
            sequence_data = infile.read()  # read sequence file content
    except OSError:
        exception_handler(
            OSError, "An error occurred while parsing the sequence file", debug
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
                guides.extend(recover_guides(sequence, pam, guide_length, pam_at_beginning, debug))
        else:  # FASTA like sequence
            sequence = "".join(sequence.split())  # clean the sequence from speces
            sequence = sequence.strip()
            # recover the guides sequences
            guides.extend(recover_guides(sequence, pam, guide_length, pam_at_beginning, debug))
    # compute guides appending N as prefix/suffix
    try:
        temp_guides = [add_n(guide, len(pam), pam_at_beginning) for guide in guides]
    except RuntimeError:
        exception_handler(
            RuntimeError, "An error occurred while appending N to guide sequences", debug
        )
    # keep 1000000000 guides at most
    guides = temp_guides[:1000000000] if len(guides) > 1000000000 else temp_guides
    return guides
