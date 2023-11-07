"""
This file contains utility functions for DNA sequence manipulation.

Functions:
- `complement(nt: str) -> str`: Returns the complement of a given nucleotide.
- `reverse_complement(sequence: str) -> str`: Returns the reverse complement of a given DNA sequence.

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
    ```
"""

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
