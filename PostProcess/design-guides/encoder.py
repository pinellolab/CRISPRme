"""
"""

from bitset import Bitset, SIZE

from typing import List

import numpy as np

IUPAC = ["A", "C", "G", "T", "R", "Y", "S", "W", "M", "B", "D", "H", "V", "N"]


def encoder(nt: str, position: int) -> Bitset:  # type: ignore
    if nt not in IUPAC:
        raise ValueError(
            f"The nucleotide {nt} (position {position}) is not a IUPAC character"
        )
    bitset = Bitset(SIZE)  # 4 - bits encoder
    if nt == IUPAC[0]:  # A
        bitset.set(0)
    elif nt == IUPAC[1]:  # C
        bitset.set(1)
    elif nt == IUPAC[2]:  # G
        bitset.set(2)
    elif nt == IUPAC[3]:  # T
        bitset.set(3)
    elif nt == IUPAC[4]:  # R
        bitset.set(0), bitset.set(2)
    elif nt == IUPAC[5]:  # Y
        bitset.set(1), bitset.set(3)
    elif nt == IUPAC[6]:  # S
        bitset.set(1), bitset.set(2)
    elif nt == IUPAC[7]:  # W
        bitset.set(0), bitset.set(3)
    elif nt == IUPAC[8]:
        bitset.set(2), bitset.set(3)
    elif nt == IUPAC[9]:
        bitset.set(0), bitset.set(1)
    elif nt == IUPAC[10]:
        bitset.set(0), bitset.set(3)
    return bitset


def encode(sequence: str) -> List[Bitset]:
    pass
