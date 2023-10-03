"""
"""

from .bitset import Bitset
from .utils import reverse_complement
from .encoder import encode_pam, encode_genome

from typing import List, Tuple


def extract_pam(
    guide: str, pam_length: int, guide_pam_len: int, pam_in_start: bool
) -> str:
    """
    Extracts the PAM (Protospacer Adjacent Motif) sequence from a given guide sequence.

    Args:
        guide (str): The guide sequence.
        pam_length (int): The length of the PAM sequence.
        guide_pam_len (int): The limit of the PAM sequence to extract.
        pam_in_start (bool): Flag indicating if the PAM is in front of the guide sequence.

    Returns:
        str: The extracted PAM sequence.

    Example:
        ```python
        guide = 'ATCG'
        pam_length = 2
        guide_pam_len = 4
        pam_in_start = True
        extracted_pam = extract_pam(guide, pam_length, guide_pam_len, pam_in_start)
        print(extracted_pam)  # Output: 'AT'
        ```
    """
    if pam_in_start:  # PAM in front of guide sequence
        return guide[:pam_length]
    # PAM at the end of guide sequence
    return guide[(guide_pam_len - pam_length) : pam_length]


def match_pam(
    genome: List[Bitset],
    pam: Tuple[List[Bitset], List[Bitset]],
    pam_length: int,
    mm_max: int,
) -> Tuple[bool, bool]:
    """
    Matches the PAM (Protospacer Adjacent Motif) sequence against a given genome.

    Args:
        genome (List[Bitset]): The genome sequence represented as a list of Bitset objects.
        pam (Tuple[List[Bitset], List[Bitset]]): The PAM sequence represented as a tuple of two lists of Bitset objects (forward and reverse).
        pam_length (int): The length of the PAM sequence.
        mm_max (int): The maximum number of allowed mismatches.

    Returns:
        Tuple[bool, bool]: A tuple of two booleans indicating if the forward and reverse PAM sequences match the genome within the allowed number of mismatches.

    Example:
        ```python
        genome = [Bitset(), Bitset()]
        pam = ([Bitset(), Bitset()], [Bitset(), Bitset()])
        pam_length = 2
        mm_max = 1
        match_fwd, match_rev = match_pam(genome, pam, pam_length, mm_max)
        print(match_fwd, match_rev)  # Output: True, False
        ```
    """
    mm_fwd, mm_rev = 0, 0
    for i in range(pam_length):
        if not (genome[i] & pam[0][i]).to_bool():  # forward
            mm_fwd += 1
        if not (genome[i] & pam[1][i]).to_bool():  # reverse
            mm_rev += 1
    return (mm_fwd < mm_max, mm_rev < mm_max)


def match_genome(
    genome: List[Bitset],
    pam: Tuple[List[Bitset], List[Bitset]],
    pam_length: int,
    guide_pam_len: int,
    mm_max: int,
    pam_in_start: bool,
) -> Tuple[List[int], List[int]]:
    """
    Matches the PAM (Protospacer Adjacent Motif) sequence against a given genome and returns the positions of the matches.

    Args:
        genome (List[Bitset]): The genome sequence represented as a list of Bitset objects.
        pam (Tuple[List[Bitset], List[Bitset]]): The PAM sequence represented as a tuple of two lists of Bitset objects (forward and reverse).
        pam_length (int): The length of the PAM sequence.
        guide_pam_len (int): The length of the guide sequence including the PAM.
        mm_max (int): The maximum number of allowed mismatches.
        pam_in_start (bool): Flag indicating if the PAM is at the start of the guide sequence.

    Returns:
        Tuple[List[int], List[int]]: A tuple of two lists representing the positions of the matches on the forward and reverse strands.

    Example:
        ```python
        genome = [Bitset(), Bitset()]
        pam = ([Bitset(), Bitset()], [Bitset(), Bitset()])
        pam_length = 2
        guide_pam_len = 4
        mm_max = 1
        pam_in_start = True
        positions_fwd, positions_rev = match_genome(genome, pam, pam_length, guide_pam_len, mm_max, pam_in_start)
        print(positions_fwd, positions_rev)  # Output: [0], []
        ```
    """
    positions_pam = ([], [])
    for i in range(len(genome) - pam_length + 1):
        matching_pam = match_pam(genome[i : (i + pam_length)], pam, pam_length, mm_max)
        if pam_in_start:  # pam at 3'
            if matching_pam[0] and i < (
                len(genome) - guide_pam_len + 1
            ):  # match on forward strand
                positions_pam[0].append(i)
            if matching_pam[1]:  # match on reverse strand
                pamidx = (i + pam_length - 1) - (guide_pam_len - 1)
                if pamidx >= 0:
                    positions_pam[1].append(pamidx)
        else:  # pam at 5'
            if matching_pam[0]:  # match on forward strand
                pamidx = (i + pam_length - 1) - (guide_pam_len - 1)
                if pamidx >= 0:
                    positions_pam[0].append(pamidx)
            if matching_pam[1] and i < (
                len(genome) - guide_pam_len + 1
            ):  # match on reverse strand
                positions_pam[1].append(i)
    return positions_pam


def search_pam(
    guide: str,
    genome: str,
    pam_length: int,
    guide_pam_len: int,
    mm: int,
    pam_in_start: bool,
) -> None:
    """
    Searches for a PAM (Protospacer Adjacent Motif) sequence in a given genome.

    Args:
        guide (str): The guide sequence.
        genome (str): The genome sequence.
        pam_length (int): The length of the PAM sequence.
        guide_pam_len (int): The length of the guide sequence including the PAM.
        mm (int): The maximum number of allowed mismatches.
        pam_in_start (bool): Flag indicating if the PAM is at the start of the guide sequence.

    Returns:
        None

    Example:
        ```python
        guide = "ATCG"
        genome = "GATCGATCG"
        pam_length = 2
        guide_pam_len = 4
        mm = 1
        pam_in_start = True
        search_pam(guide, genome, pam_length, guide_pam_len, mm, pam_in_start)
        ```
    """
    pam = extract_pam(guide, pam_length, guide_pam_len, pam_in_start)
    # encode PAM in bits for efficient search
    pam_bits = encode_pam(pam), encode_pam(reverse_complement(pam))
    # encode genome in bits
    genome_bits = encode_genome(genome)
    # match PAM against the genome
    pam_positions = match_genome(
        genome_bits, pam_bits, pam_length, guide_pam_len, mm, pam_in_start
    )
