"""
The `reporter` module provides functions for processing and reporting guide sequences.

Functions:
- `is_valid_sequence(sequence: str) -> bool`: Checks if a sequence is valid by ensuring it does not contain 'N'.
- `process_seqname(coords: str) -> Tuple[str, int]`: Processes sequence coordinates and returns the sequence name and relative start position.
- `expand_iupac(nt: str, i: int) -> str`: Expands an input nucleotide into its IUPAC variants.
- `explode_iupac_sequence(sequence: str) -> List[str]`: Expands an IUPAC sequence into all possible combinations of nucleotides.
- `pad_guide(guide: str, pam_length: int, pam_at_beginning: bool) -> str`: Pads a guide sequence with 'N' nucleotides based on the PAM position.
- `extract_guide_pam(seqname: str, positions: Tuple[List[int], List[int]], genome: pysam.FastaFile, pam: PAM) -> List[List[str]]`: Extracts guide sequences and their corresponding PAM sequences from a genome based on given positions.
- `write_report(report: List[List[str]], outname: str) -> None`: Writes the guide report and summary to files.
- `recover_guides(pam_positions: Dict[str, Tuple[List[int], List[int]]], genome: pysam.FastaFile, pam: PAM, outname: str) -> None`: Recovers guide sequences and writes the report.
- `write_results(guide_list: list, guide_file: str, positions: Tuple, chr_name: str) -> None`: Writes the guide list and summary to files.
"""


from utils import IUPACTABLE, reverse_complement
from pam import PAM

from typing import Tuple, Dict, List

import pysam

REPORTCOLS = ["sgRNA", "PAM", "sgRNA_length", "PAM_length", "chr", "pos", "strand"]


def is_valid_sequence(sequence: str) -> bool:
    """
    Checks if a sequence is valid by ensuring it does not contain 'N'.

    Args:
        sequence (str): The sequence to be checked.

    Returns:
        bool: True if the sequence is valid, False otherwise.

    Example:
        ```python
        sequence = "ATGC"
        is_valid = is_valid_sequence(sequence)
        print(is_valid)  # Output: True
        ```
    """

    return "N" not in sequence


def process_seqname(coords: str) -> Tuple[str, int]:
    """
    Processes sequence coordinates and returns the sequence name and relative start position.

    Args:
        coords (str): The sequence coordinates.

    Returns:
        Tuple[str, int]: A tuple containing the sequence name and relative start position.

    Example:
        ```python
        coords = "chr1:100-200"
        seqname, start = process_seqname(coords)
        print(seqname, start)  # Output: ('chr1', 100)
        ```
    """

    fields = coords.strip().split(":")
    if len(fields) > 1:  # coordinates given as input
        return fields[0], int(fields[1].split("-")[0])
    return fields[0], 0


def expand_iupac(nt: str, i: int) -> str:
    """
    Expands an input nucleotide into its IUPAC variants.

    Args:
        nt (str): The input nucleotide.
        i (int): The position of the nucleotide.

    Returns:
        str: The expanded IUPAC variants of the nucleotide.

    Raises:
        ValueError: If the nucleotide does not belong to the IUPAC alphabet.

    Example:
        ```python
        nt = "R"
        i = 2
        expanded_nt = expand_iupac(nt, i)
        print(expanded_nt)  # Output: 'AG'
        ```
    """

    if nt in IUPACTABLE:
        return IUPACTABLE[nt]
    raise ValueError(f"Forbidden nucleotide ({nt}) at position {i}")


def explode_iupac_sequence(sequence: str) -> List[str]:
    """
    Expands an IUPAC sequence into all possible combinations of nucleotides.

    Args:
        sequence (str): The IUPAC sequence to be expanded.

    Returns:
        List[str]: A list of all possible combinations of nucleotides.

    Example:
        ```python
        sequence = "RY"
        expanded_sequences = explode_iupac_sequence(sequence)
        print(expanded_sequences)  # Output: ['AC', 'AT', 'GC', 'GT']
        ```
    """

    xpanded_sequences = [""]  # expanded IUPAC sequences
    for i, nt in enumerate(sequence):
        xpanded_nt = expand_iupac(nt, i)  # expand input nt in its IUPAC variants
        xpanded_sequences = [s + b for s in xpanded_sequences for b in xpanded_nt]
    assert len(xpanded_sequences) >= 1
    return xpanded_sequences


def pad_guide(guide: str, pam_length: int, pam_at_beginning: bool):
    """
    Pads a guide sequence with 'N' nucleotides based on the PAM position.

    Args:
        guide (str): The guide sequence.
        pam_length (int): The length of the PAM sequence.
        pam_at_beginning (bool): Flag indicating whether the PAM occurs at the beginning.

    Returns:
        str: The padded guide sequence.

    Example:
        ```python
        guide = "ATGC"
        pam_length = 3
        pam_at_beginning = True
        padded_guide = pad_guide(guide, pam_length, pam_at_beginning)
        print(padded_guide)  # Output: 'NNNATGC'
        ```
    """

    if pam_at_beginning:
        return "N" * pam_length + guide
    return guide + "N" * pam_length


# TODO: improve
def extract_guide_pam(
    seqname: str,
    positions: Tuple[List[int], List[int]],
    genome: pysam.FastaFile,
    pam: PAM,
) -> List[List[str]]:
    report = []
    seqname, pos = process_seqname(
        seqname
    )  # recover seqname and relative start position
    for p in positions[0]:  # forward
        p += pos
        start = p + len(pam) if pam.pam_at_beginning else p
        end = (
            p + len(pam) + pam.guide_length
            if pam.pam_at_beginning
            else p + pam.guide_length
        )
        guide_iupac = genome.fetch(seqname, start, end).upper()  # allow IUPAC chars
        if is_valid_sequence(guide_iupac):
            guides = explode_iupac_sequence(guide_iupac.upper())  # replace IUPAC chars
            start = p if pam.pam_at_beginning else p + pam.guide_length
            end = (
                p + pam.guide_length
                if pam.pam_at_beginning
                else p + len(pam) + pam.guide_length
            )
            pam_seq = genome.fetch(seqname, start, end)
            if is_valid_sequence(pam_seq):
                report.extend(
                    list(
                        map(
                            str,
                            [
                                pad_guide(guide, len(pam), pam.pam_at_beginning),
                                pam_seq,
                                pam.guide_length,
                                len(pam),
                                seqname,
                                p,
                                "+",
                            ],
                        )
                    )
                    for guide in guides
                )
    for p in positions[1]:  # reverse
        p += pos
        start = p if pam.pam_at_beginning else p + len(pam)
        end = (
            p + pam.guide_length
            if pam.pam_at_beginning
            else p + len(pam) + pam.guide_length
        )
        print(p, start, end)
        guide_iupac = reverse_complement(genome.fetch(seqname, start, end).upper())
        if is_valid_sequence(guide_iupac):
            guides = explode_iupac_sequence(guide_iupac)  # replace IUPAC chars
            start = p + pam.guide_length if pam.pam_at_beginning else p
            end = (
                p + len(pam) + pam.guide_length
                if pam.pam_at_beginning
                else p + len(pam)
            )
            pam_seq = genome.fetch(seqname, start, end)
            if is_valid_sequence(pam_seq):
                report.extend(
                    list(
                        map(
                            str,
                            [
                                pad_guide(guide, len(pam), pam.pam_at_beginning),
                                pam_seq,
                                pam.guide_length,
                                len(pam),
                                seqname,
                                p,
                                "-",
                            ],
                        )
                    )
                    for guide in guides
                )
    return report


def write_report(report: List[List[str]], outname: str) -> None:
    """
    Writes the guide report and summary to files.

    Args:
        report (List[List[str]]): The report data to be written.
        outname (str): The base name for the output files.

    Returns:
        None

    Raises:
        IOError: If there is an error while writing the report.

    Example:
        ```python
        report = [
            ["ATGC", "GGG", 4, 3, "chr1", 100, "+"],
            ["TACG", "CCC", 4, 3, "chr1", 200, "-"],
        ]
        outname = "output"
        write_report(report, outname)
        # Creates output.guides.txt and output.summary.tsv files
        ```
    """

    try:
        guides_file = f"{outname}.guides.txt"
        report_file = f"{outname}.summary.tsv"
        with open(guides_file, mode="w") as outguide, open(
            report_file, mode="w"
        ) as outreport:
            # write report header
            header = "\t".join(REPORTCOLS)
            outreport.write(f"{header}\n")
            for line in report:
                rline = "\t".join(line)  # report line
                outguide.write(f"{line[0]}\n")
                outreport.write(f"{rline}\n")
    except IOError as e:
        raise IOError("Report writing failed!") from e


def recover_guides(
    pam_positions: Dict[str, Tuple[List[int], List[int]]],
    genome: pysam.FastaFile,
    pam: PAM,
    outname: str,
) -> None:
    """
    Recovers guide sequences based on PAM positions and writes the report.

    Args:
        pam_positions (Dict[str, Tuple[List[int], List[int]]]): A dictionary mapping sequence names to PAM positions.
        genome (pysam.FastaFile): The genome file.
        pam (PAM): The PAM object.
        outname (str): The base name for the output files.

    Returns:
        None

    Example:
        ```python
        pam_positions = {
            "seq1": ([100, 200], [300, 400]),
            "seq2": ([500, 600], [700, 800]),
        }
        genome = pysam.FastaFile("genome.fasta")
        pam = PAM("path/to/pam.txt")
        outname = "output"
        recover_guides(pam_positions, genome, pam, outname)
        # Creates output.guides.txt and output.summary.tsv files
        ```
    """

    report = []
    for seqname in pam_positions:
        report += extract_guide_pam(seqname, pam_positions[seqname], genome, pam)
    write_report(report, outname)
