"""
"""

SIZE = 4


class Bitset(object):
    """
    A class representing a Bitset.

    Args:
        size (int): The size of the Bitset.

    Raises:
        TypeError: If size is not an integer.
        ValueError: If size is less than 1.

    Example:
        ```python
        bitset = Bitset(8)
        print(bitset)  # Output: '00000000'
        ```
    """

    def __init__(self, size: int) -> None:
        """
        Initializes a Bitset instance with the given size.

        Args:
            size (int): The size of the Bitset.

        Raises:
            TypeError: If size is not an integer.
            ValueError: If size is less than 1.

        Returns:
            None
        """
        if not isinstance(size, int):
            raise TypeError(f"Expected {int.__name__}, got {type(size).__name__}")
        if size < 1:
            raise ValueError("Forbidden Bitset size")
        self._size = size
        self._bits = 0

    def __str__(self) -> str:
        """
        Returns a binary string representation of the Bitset.

        Returns:
            str: The binary string representation of the Bitset.
        """
        return bin(self._bits)[2:].zfill(self._size)

    def set(self, index: int) -> None:
        """
        Sets the bit at the given index in the Bitset.

        Args:
            index (int): The index of the bit to set.

        Raises:
            TypeError: If index is not an integer.
            IndexError: If index is out of bounds.

        Returns:
            None
        """
        if not isinstance(index, int):
            raise TypeError(f"Expected {int.__name__}, got {type(index).__name__}")
        if index >= self._size:
            raise IndexError("Index out of bounds")
        self._bits |= 1 << index

    def reset(self, index: int) -> None:
        """
        Resets the bit at the given index in the Bitset.

        Args:
            index (int): The index of the bit to reset.

        Raises:
            TypeError: If index is not an integer.
            IndexError: If index is out of bounds.

        Returns:
            None
        """
        if not isinstance(index, int):
            raise TypeError(f"Expected {int.__name__}, got {type(index).__name__}")
        if index >= self._size:
            raise IndexError("Index out of bounds")
        self._bits &= ~(1 << index)

    def set_bits(self, bits: str) -> None:
        """
        Sets the bits in the Bitset based on the given bit string.

        Args:
            bits (str): The bit string to set in the Bitset.

        Raises:
            TypeError: If bits is not a string.
            ValueError: If bits contains non-binary digits.

        Returns:
            None

        Example:
            ```python
            bitset = Bitset(8)
            bitset.set_bits("10101010")
            print(bitset)  # Output: '10101010'
            ```
        """
        if not isinstance(bits, str):
            raise TypeError(f"Expected {str.__name__}, got {type(bits).__name__}")
        if any(bit not in "01" for bit in bits):
            raise ValueError(f"{bits} is not a bit string")
        bitstring_size = len(bits)
        for i, bit in enumerate(bits):
            if bit == "0":  # reset bit
                self.reset(bitstring_size - 1 - i)
            else:  # set bit (bit == 1)
                self.set(bitstring_size - 1 - i)
        assert str(self) == bits

    def test(self, index: int) -> bool:
        """
        Tests the value of the bit at the given index in the Bitset.

        Args:
            index (int): The index of the bit to test.

        Raises:
            IndexError: If index is out of bounds.

        Returns:
            bool: True if the bit is set, False otherwise.
        """
        if index >= self._size:
            raise IndexError("Bitset out of bounds")
        return bool(self._bits & (1 << index))
