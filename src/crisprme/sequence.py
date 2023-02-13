"""
"""

from utils import PAM_DICT, exception_handler

from itertools import product
from typing import List, Tuple
from Bio.Seq import Seq

import pybedtools
import re

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

def _generate_iupac_pam(pam: str) -> List[str]:
    """ (PRIVATE)
    Generate all possible IUPAC PAMs based on the input PAM characters

    :param pam: input PAM
    :type pam: str
    :return: IUPAC PAMs 
    :rtype: List[str]
    """
    product_lst = [PAM_DICT[c] for c in pam]
    iupac_pam = ["".join(e) for e in product(*product_lst)]
    return iupac_pam

def _extract_guides(sequence: str, iupac_pam: List[str], pam_len: int, guide_len: int, pam_at_beginning: bool, debug: bool) -> List[str]:
    """(PRIVATE)
    Extract guide sequences from the input sequence, matching the input PAM 
    sequences.

    :param sequence: input sequence
    :type sequence: str
    :param iupac_pam: list of PAM sequences
    :type iupac_pam: List[str]
    :param pam_len: pam length
    :type pam_len: int
    :param guide_len: guide length
    :type guide_len: int
    :param pam_at_beginning: PAM occurs at the beginning of the guide
    :type pam_at_beginning: bool
    :param debug: debug mode
    :type debug: bool
    :return: matching guide sequences
    :rtype: List[str]
    """
    guides = []  # guides list
    try:
        for pam in iupac_pam:
            # find PAM occurrences in the input sequence
            pam_indices = [m.start() for m in re.finditer(f'(?={pam})', sequence)]
            for idx in pam_indices:  # recover guide start and stop positions
                start, stop = (idx + pam_len, idx + pam_len + guide_len) if pam_at_beginning else (idx - guide_len, idx)
                if (pam_at_beginning and stop > len(sequence) - guide_len) or (not pam_at_beginning and start < 0):
                    continue  # out of bounds -> skip
                guides.append(sequence[start:stop])
    except RuntimeError:
        exception_handler(
            RuntimeError, "An error occurred while extracting guide sequences", debug
        )
    return guides

def recover_guides(sequence: str, pam: str, guide_length: int, pam_at_beginning: bool, debug: bool) -> List[str]:
    """Recover guide sequences from the input sequence. The guide sequences 
    search is lead using the input PAM sequence 

    :param sequence: input sequence
    :type sequence: str
    :param pam: input PAM 
    :type pam: str
    :param guide_length: guide length
    :type guide_length: int
    :param pam_at_beginning: PAM occurs at the beginning of the guide
    :type pam_at_beginning: bool
    :param debug: debug mode
    :type debug: bool
    :return: guide sequences list
    :rtype: List[str]
    """
    if not isinstance(sequence, str):
        exception_handler(
            TypeError, f"Expected {str.__name__}, got {type(sequence).__name__}", debug
        )
    if not isinstance(pam, str):
        exception_handler(
            TypeError, f"Expected {str.__name__}, got {type(pam).__name__}", debug
        )
    if not isinstance(guide_length, int):
        exception_handler(
            TypeError, f"Expected {int.__name__}, got {type(guide_length).__name__}", debug
        )
    # generate IUPAC PAMs (e.g. NNNNN NGG / CCN NNNNN)
    iupac_pam_fwd = _generate_iupac_pam(pam)  # forward
    iupac_pam_rev = _generate_iupac_pam(str(Seq(pam).reverse_complement()))  # reverse
    sequence = sequence.upper()  # force characters to upper case
    # extract guides sequences
    guides_fwd = _extract_guides(
        sequence, iupac_pam_fwd, len(pam), guide_length, pam_at_beginning, debug
    )  # forward strand guides
    guides_rev = _extract_guides(
        sequence, iupac_pam_rev, len(pam), guide_length, pam_at_beginning, debug
    )  # reverse strand guides
    guides = guides_fwd + guides_rev
    assert len(guides) == (len(guides_fwd) + len(guides_rev))
    return guides



