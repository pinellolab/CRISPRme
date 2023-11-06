"""
"""

from utils import reverse_complement
from encoder import encode_pam, encode_genome
from search import match_genome
from pam import PAM
from reporter import recover_guides

from typing import Tuple, List

import pysam


def read_pam(pamfile: str) -> PAM:
    return PAM(pamfile)


def read_genome(genome: str) -> pysam.FastaFile:
    return pysam.FastaFile(genome)  # load FastaFile object


def read_coordinates(bedfile: str) -> List[Tuple[str, int, int]]:
    try:
        with open(bedfile, mode="r") as infile:
            coordinates = [
                (fields[0], int(fields[1]), int(fields[2]))
                for line in infile
                for fields in [line.strip().split()]
            ]
    except IOError as e:
        raise IOError("BED file parsing failed!") from e
    return coordinates


def extract_sequences(
    genome: pysam.FastaFile, coords: List[Tuple[str, int, int]]
) -> str:
    try:
        sequences = [genome.fetch(coord[0], coord[1], coord[2]) for coord in coords]
    except ValueError as e:
        raise ValueError("Sequence extraction failed!") from e
    assert len(sequences) == len(coords)  # should match
    return "N".join(sequences)  # N is sequence separator


def design_guides(
    pamfile: str, genome: str, bedfile: str, mm_max: int, outname: str
) -> None:
    pam = read_pam(pamfile)  # read PAM
    pam_bits = encode_pam(pam.pam), encode_pam(reverse_complement(pam.pam))
    sequence = read_genome(genome)
    if bedfile:  # extract subsequences
        coordinates = read_coordinates(bedfile)
        sequence_bits = {
            f"{coords[0]}:{coords[1]}-{coords[2]}": encode_genome(
                sequence.fetch(coords[0], coords[1], coords[2]).upper()
            )
            for coords in coordinates
        }
    else:
        sequence_bits = {
            seqname: encode_genome(sequence[seqname].upper())
            for seqname in sequence.references
        }
    # search PAM occurrences within the query sequences
    pam_positions = {
        seqname: match_genome(
            sequence_bits[seqname],
            pam_bits,
            pam.length,
            pam.full_length,
            mm_max,
            pam.pam_at_beginning,
        )
        for seqname in sequence_bits
    }

    # print(sequence["chr_test"][88:112].upper())
    # print(sequence["chr_test"][95:119].upper())
    # print(sequence["chr_test"][105:129].upper())
    # print(sequence["chr_test"][127:151].upper())
    # print(sequence["chr_test"][150:174].upper())

    # write report and guides file
    recover_guides(pam_positions, sequence, pam, outname)


def extract_guides_from_genome(
    positions: Tuple, genome: str, guide_len: int, pam_len: int, pam_in_start: bool
) -> list:
    """_summary_

    Args:
        positions (Tuple): _description_
        genome (str): _description_
        guide_len (int): _description_
        pam_len (int): _description_
        pam_in_start (bool): _description_

    Returns:
        list: _description_
    """
    guides_list = list()
    ##pam at start == true
    if pam_in_start:
        for pos in positions[0]:
            guide = genome[pos + pam_len : pos + pam_len + guide_len]  # type: ignore
            if "N" in guide:
                continue
            guide = "N" * pam_len + guide
            pam = genome[pos : pos + pam_len]
            if "N" in pam:
                continue
            guides_list.append([guide, pam, "forward"])

        for pos in positions[1]:
            guide = genome[pos : pos + guide_len]
            if "N" in guide:
                continue
            guide = reverse_complement(guide)
            guide = "N" * pam_len + guide
            pam = genome[pos + guide_len : pos + pam_len + guide_len]
            if "N" in pam:
                continue
            guides_list.append([guide, pam, "reverse"])

        return guides_list
    else:
        ##pam at start == false
        for pos in positions[0]:
            guide = genome[pos : pos + guide_len]
            if "N" in guide:
                continue
            guide = guide + "N" * pam_len
            pam = genome[pos + guide_len : pos + pam_len + guide_len]
            if "N" in pam:
                continue
            guides_list.append([guide, pam, "forward"])

        for pos in positions[1]:
            guide = genome[pos + pam_len : pos + pam_len + guide_len]
            if "N" in guide:
                continue
            guide = reverse_complement(guide)
            guide = guide + "N" * pam_len
            pam = genome[pos : pos + pam_len]
            if "N" in pam:
                continue
            guides_list.append([guide, pam, "reverse"])

    return guides_list
