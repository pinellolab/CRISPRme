"""
This file contains functions for designing guides for CRISPR/Cas9 genome editing.

The file includes the following functions:

- `read_pam(pamfile: str) -> PAM`: Reads the PAM (Protospacer Adjacent Motif) sequence from a file and returns a `PAM` object.
- `read_genome(genome: str) -> pysam.FastaFile`: Reads the genome sequence from a file and returns a `pysam.FastaFile` object.
- `read_coordinates(bedfile: str) -> List[Tuple[str, int, int]]`: Reads genomic coordinates from a BED file and returns a list of tuples.
- `design_guides(pamfile: str, genome: str, bedfile: str, mm_max: int, outname: str) -> None`: Designs guides for CRISPR/Cas9 genome editing based on the PAM sequence, genome sequence, BED file, maximum number of mismatches, and output file name.
- `extract_guides_from_genome(positions: Tuple, genome: str, guide_len: int, pam_len: int, pam_in_start: bool) -> list`: Extracts guides from the genome based on the given positions, guide length, PAM length, and PAM position.

Please refer to the code for more details on each function.
"""

from utils import reverse_complement
from encoder import encode_pam, encode_genome
from search import match_genome
from pam import PAM
from reporter import recover_guides

from typing import Tuple, List

import pysam


def read_pam(pamfile: str) -> PAM:
    """
    Reads the PAM (Protospacer Adjacent Motif) sequence from a file and returns a PAM object.

    Args:
        pamfile (str): The path to the PAM file.

    Returns:
        PAM: A PAM object representing the PAM sequence.

    Example:
        ```python
        pamfile = "path/to/pam.txt"
        pam = read_pam(pamfile)
        ```
    """

    return PAM(pamfile)


def read_genome(genome: str) -> pysam.FastaFile:
    """
    Reads the genome sequence from a file and returns a `pysam.FastaFile` object.

    Args:
        genome (str): The path to the genome file.

    Returns:
        pysam.FastaFile: A `pysam.FastaFile` object representing the genome sequence.

    Example:
        ```python
        genome_file = "path/to/genome.fasta"
        genome = read_genome(genome_file)
        ```
    """

    return pysam.FastaFile(genome)  # load FastaFile object


def read_coordinates(bedfile: str) -> List[Tuple[str, int, int]]:
    """
    Reads genomic coordinates from a BED file and returns a list of tuples.

    Args:
        bedfile (str): The path to the BED file.

    Returns:
        List[Tuple[str, int, int]]: A list of tuples representing the genomic coordinates.

    Raises:
        IOError: If the BED file parsing fails.

    Example:
        ```python
        bedfile = "path/to/coordinates.bed"
        coordinates = read_coordinates(bedfile)
        ```
    """

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


def design_guides(
    pamfile: str, genome: str, bedfile: str, mm_max: int, outname: str
) -> None:
    """
    Designs guides for CRISPR/Cas9 genome editing based on the given PAM sequence, genome sequence, BED file, maximum number of mismatches, and output file name.

    Args:
        pamfile (str): The path to the PAM file.
        genome (str): The path to the genome file.
        bedfile (str): The path to the BED file (optional).
        mm_max (int): The maximum number of mismatches allowed.
        outname (str): The name of the output file.

    Returns:
        None

    Example:
        ```python
        pamfile = "path/to/pam.txt"
        genome = "path/to/genome.fasta"
        bedfile = "path/to/bedfile.bed"
        mm_max = 2
        outname = "output.txt"
        design_guides(pamfile, genome, bedfile, mm_max, outname)
        ```
    """

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
    # write report and guides file
    recover_guides(pam_positions, sequence, pam, outname)
