"""
This file contains utility functions for DNA sequence manipulation and file reading.

Functions:
- `complement(nt: str) -> str`: Returns the complement of a given nucleotide.
- `reverse_complement(sequence: str) -> str`: Returns the reverse complement of a given DNA sequence.
- `read_pam(pamfile: str) -> PAM`: Reads the PAM sequence from a file and returns a PAM object.
- `read_genome(genome: str) -> pysam.FastaFile`: Reads the genome sequence from a file and returns a `pysam.FastaFile` object.
- `read_coordinates(bedfile: str) -> List[Tuple[str, int, int]]`: Reads genomic coordinates from a BED file and returns a list of tuples.

Constants:
- `IUPAC`: A list of IUPAC characters representing nucleotides and their degenerate bases.
- `RC`: A dictionary mapping nucleotides to their complement.
- `IUPACTABLE`: A dictionary mapping IUPAC characters to their corresponding nucleotides.

Note: The `complement` function raises a `KeyError` if a forbidden IUPAC character is encountered.

Example:
    ```python
    nt = 'A'
    complement = complement(nt)
    print(complement)  # Output: 'T'

    sequence = 'ATCG'
    reverse_complement = reverse_complement(sequence)
    print(reverse_complement)  # Output: 'CGAT'

    pamfile = "path/to/pam.txt"
    pam = read_pam(pamfile)

    genome_file = "path/to/genome.fasta"
    genome = read_genome(genome_file)

    bedfile = "path/to/coordinates.bed"
    coordinates = read_coordinates(bedfile)
    ```
"""


from pam import PAM

from typing import List, Tuple

import pysam

# constant variables
IUPAC = ["A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]
RC = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "R": "Y",
    "Y": "R",
    "M": "K",
    "K": "M",
    "H": "D",
    "D": "H",
    "B": "V",
    "V": "B",
    "N": "N",
}
IUPACTABLE = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "M": "AC",
    "K": "GT",
    "S": "CG",
    "W": "AT",
    "H": "ACT",
    "B": "CGT",
    "V": "ACG",
    "D": "AGT",
    "N": "ACGT",
}


def complement(nt: str) -> str:
    """
    Returns the complement of a given nucleotide.

    Args:
        nt (str): The nucleotide to find the complement of.

    Returns:
        str: The complement of the given nucleotide.

    Raises:
        KeyError: Raised when a forbidden IUPAC character is encountered.

    Example:
        ```python
        nt = 'A'
        complement = complement(nt)
        print(complement)  # Output: 'T'
        ```
    """
    try:
        return RC[nt]
    except KeyError as e:
        raise KeyError(f"Forbidden IUPAC character encountered ({nt})") from e


def reverse_complement(sequence: str) -> str:
    """
    Returns the reverse complement of a given DNA sequence.

    Args:
        sequence (str): The DNA sequence to find the reverse complement of.

    Returns:
        str: The reverse complement of the given DNA sequence.

    Example:
        ```python
        sequence = 'ATCG'
        reverse_complement = reverse_complement(sequence)
        print(reverse_complement)  # Output: 'CGAT'
        ```
    """
    return "".join([complement(nt) for nt in sequence[::-1]])


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
