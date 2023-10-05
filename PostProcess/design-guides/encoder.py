"""
"""

from bitset import Bitset, SIZE
from utils import IUPAC

from typing import List

import numpy as np


def encoder(nt: str, position: int) -> Bitset:  # type: ignore
    """
    Encodes a nucleotide into a Bitset based on the IUPAC encoding scheme.

    Args:
        nt (str): The nucleotide to encode.
        position (int): The position of the nucleotide.

    Raises:
        TypeError: If nt is not a string.
        ValueError: If nt is not a valid IUPAC character.

    Returns:
        Bitset: The encoded Bitset.

    Example:
        ```python
        encoded = encoder("A", 0)
        print(encoded)  # Output: '0001'
        ```
    """
    if not isinstance(nt, str):
        raise TypeError(f"Expected {str.__name__}, got {type(nt).__name__}")
    if nt not in IUPAC:
        raise ValueError(
            f"The nucleotide {nt} (position {position}) is not a IUPAC character"
        )
    bitset = Bitset(SIZE)  # 4 - bits encoder
    if nt == IUPAC[0]:  # A - 0001
        bitset.set(0)
    elif nt == IUPAC[1]:  # C - 0010
        bitset.set(1)
    elif nt == IUPAC[2]:  # G - 0100
        bitset.set(2)
    elif nt == IUPAC[3]:  # T - 1000
        bitset.set(3)
    elif nt == IUPAC[4]:  # R - 0101
        bitset.set_bits("0101")
    elif nt == IUPAC[5]:  # Y - 1010
        bitset.set_bits("1010")
    elif nt == IUPAC[6]:  # S - 0110
        bitset.set_bits("0110")
    elif nt == IUPAC[7]:  # W  - 1001
        bitset.set_bits("1001")
    elif nt == IUPAC[8]:  # K - 1100
        bitset.set_bits("1100")
    elif nt == IUPAC[9]:  # M - 0011
        bitset.set_bits("0011")
    elif nt == IUPAC[10]:  # B - 1110
        bitset.set_bits("1110")
    elif nt == IUPAC[11]:  # D - 1101
        bitset.set_bits("1101")
    elif nt == IUPAC[12]:  # H - 1011
        bitset.set_bits("1011")
    elif nt == IUPAC[13]:  # V - 0111
        bitset.set_bits("0111")
    elif nt == IUPAC[14]:  # N - 1111
        bitset.set_bits("1111")
    # else:  # not a valid IUPAC character
    #     assert str(bitset) == "0000"
    return bitset


def encode_pam(pam: str) -> List[Bitset]:  # type: ignore
    """
    Encodes a PAM sequence into a list of Bitsets using the IUPAC encoding scheme.

    Args:
        pam (str): The PAM sequence to encode.

    Raises:
        TypeError: If pam is not a string.

    Returns:
        List[Bitset]: The list of encoded Bitsets.

    Example:
        ```python
        encoded_pam = encode_pam("NGG")
        print(encoded_pam)  # Output: [Bitset('0000'), Bitset('0100'), Bitset('0100')]
        ```
    """
    if not isinstance(pam, str):
        raise TypeError(f"Expected {str.__name__}, got {type(pam).__name__}")
    return [encoder(nt, i) for i, nt in enumerate(pam)]


def encode_genome(genome: str) -> List[Bitset]:  # type: ignore
    """
    Encodes a given genome into a list of bitsets.

    Args:
        genome (str): The genome to encode.

    Returns:
        List[Bitset]: A list of bitsets representing the encoded genome.

    Raises:
        TypeError: Raised when the input genome is not a string.

    Example:
        ```python
        genome = 'ATCG'
        encoded_genome = encode_genome(genome)
        print(encoded_genome)  # Output: [Bitset('0001'), Bitset('1000'), Bitset('0010'), Bitset('0100')]
        ```
    """
    if not isinstance(genome, str):
        raise TypeError(f"Expected {str.__name__}, got {type(genome).__name__}")
    return [encoder(nt, i) for i, nt in enumerate(genome)]
