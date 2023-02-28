"""
"""

from verbosity_handler import write_verbosity
from utils import PAM_DICT, exception_handler

from pybedtools import BedTool
from itertools import product
from typing import List, Optional, Tuple
from Bio.Seq import Seq

import re
import os


def extract_sequence(
    seqname: str, coord: str, genome: str, debug: bool
) -> Tuple[str, str]:
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
    coordinate = BedTool(f"{coord}", from_string=True)  # initialize the object
    # extract the sequence
    try:
        coordinate = coordinate.sequence(fi=genome)
        with open(coordinate.seqfn, mode="r") as infile:
            sequences = [seq.strip() for seq in infile if not seq.startswith(">")]
    except RuntimeError:
        exception_handler(
            RuntimeError, "An error occurred while extracting guide sequence", debug
        )
    assert len(sequences) == 1  # each coordinate should correspond to one sequence
    sequence = sequences[0]
    return seqname, sequence


def _generate_iupac_pam(pam: str) -> List[str]:
    """(PRIVATE)
    Generate all possible IUPAC PAMs based on the input PAM characters

    :param pam: input PAM
    :type pam: str
    :return: IUPAC PAMs
    :rtype: List[str]
    """
    product_lst = [PAM_DICT[c] for c in pam]
    return ["".join(e) for e in product(*product_lst)]


def _extract_guides(
    sequence: str,
    iupac_pam: List[str],
    pam_len: int,
    guide_len: int,
    pam_at_beginning: bool,
    debug: bool,
    reverse: Optional[bool] = False,
) -> List[str]:
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
    :param reverse: reverse strand, defaults to False
    :type reverse: Optional[bool], optional
    :return: matching guide sequences
    :rtype: List[str]
    """
    guides = []  # guides list
    try:
        for pam in iupac_pam:
            # find PAM occurrences in the input sequence
            if pam_indices := [
                m.start() for m in re.finditer(f"(?={pam})", sequence)
            ]:
                for idx in pam_indices:  # recover guide start and stop positions
                    if reverse:  # reverse strand
                        start, stop = (
                            (idx - guide_len, idx)
                            if pam_at_beginning
                            else (idx + pam_len, idx + guide_len + pam_len)
                        )
                        if (pam_at_beginning and idx < guide_len) or (
                            not pam_at_beginning
                            and idx > len(sequence) - guide_len - pam_len
                        ):
                            continue  # out of bounds -> skip
                        guide = str(Seq(sequence[start:stop]).reverse_complement())
                    else:  # forward strand
                        start, stop = (
                            (idx + pam_len, idx + pam_len + guide_len)
                            if pam_at_beginning
                            else (idx - guide_len, idx)
                        )
                        if (pam_at_beginning and stop > len(sequence) - guide_len) or (
                            not pam_at_beginning and start < 0
                        ):
                            continue  # out of bounds -> skip
                        guide = sequence[start:stop]
                    guides.append(guide)
    except RuntimeError:
        exception_handler(
            RuntimeError, "An error occurred while extracting guide sequences", debug
        )
    return guides


def recover_guides(
    sequence: str, pam: str, guide_length: int, pam_at_beginning: bool, debug: bool
) -> List[str]:
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
            TypeError,
            f"Expected {int.__name__}, got {type(guide_length).__name__}",
            debug,
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
        sequence,
        iupac_pam_rev,
        len(pam),
        guide_length,
        pam_at_beginning,
        debug,
        reverse=True,
    )  # reverse strand guides
    guides = guides_fwd + guides_rev
    assert len(guides) == (len(guides_fwd) + len(guides_rev))
    return guides

def write_guidefile(guide: str, verbosity: int, debug: bool) -> str:
    """Store the guide to a TXT file

    :param guide: guide
    :type guide: str
    :param verbosity: verbosity level
    :type verbosity: int
    :param debug: debug mode
    :type debug: bool
    :raises OSError: raise on errors occurring while writing the guide file
    :return: guide filename
    :rtype: str
    """
    guidefile = f".{guide}.guide"
    write_verbosity(f"writing guide {guide} to {guidefile}", verbosity, 2, debug)
    try:
        with open(guidefile, mode="w") as handle:
            handle.write(f"{guide}\n")
    except OSError:
        exception_handler(OSError, f"An error occurred while writing guide file for {guide}", debug)
    assert os.path.isfile(guidefile)
    return guidefile