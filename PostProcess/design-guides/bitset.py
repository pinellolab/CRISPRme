"""
"""

SIZE = 4


class Bitset(object):
    def __init__(self, size: int) -> None:
        self._size = size
        self._bits = 0

    def __str__(self) -> str:
        return bin(self._bits)[2:].zfill(self._size)

    def set(self, index: int) -> None:
        self._bits |= 1 << index

    def reset(self, index: int) -> None:
        if index >= self._size:
            raise IndexError("Bitset out of bounds")
        self._bits &= ~(1 << index)

    def test(self, index: int):
        if index >= self._size:
            raise IndexError("Bitset out of bounds")
        return bool(self._bits & (1 << index))
