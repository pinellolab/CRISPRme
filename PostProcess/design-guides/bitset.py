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

    def __repr__(self) -> str:
        """
        Returns a string representation of the Bitset object.

        Returns:
            str: A string representation of the Bitset object.

        Example:
            ```python
            bitset = Bitset(2)
            print(repr(bitset))  # Output: "<Bitset object> <value: 00, size: 2>"
            ```
        """
        return f"<{self.__class__.__name__} object> <value: {str(self)}, size: {self._size}>"

    def __and__(self, bitset: "Bitset") -> "Bitset":
        """
        Performs a bitwise AND operation between two Bitset objects.

        Args:
            bitset (Bitset): The Bitset object to perform the AND operation with.

        Returns:
            Bitset: A new Bitset object representing the result of the AND operation.

        Raises:
            ValueError: Raised when the Bitsets have different sizes.

        Example:
            ```python
            bitset1 = Bitset(4)
            bitset2 = Bitset(4)
            bitset1.set_bits("1001")
            bitset2.set_bits("1110")
            result = bitset1 & bitset2
            print(result)  # Output: 1000
            ```
        """
        if self._size != bitset.size:
            raise ValueError("Bitsets must have the same size for AND operation")
        result = Bitset(self._size)
        result._bits = self._bits & bitset.bits
        return result

    def _getsize(self) -> int:
        """
        Returns the size of the Bitset object.

        Returns:
            int: The size of the Bitset object.

        Example:
            ```python
            bitset = Bitset()
            size = bitset._getsize()
            print(size)  # Output: 0
            ```
        """
        return self._size

    @property
    def size(self) -> int:
        """
        Returns the size of the Bitset object.

        Returns:
            int: The size of the Bitset object.

        Example:
            ```python
            bitset = Bitset()
            size = bitset.size
            print(size)  # Output: 0
            ```
        """
        return self._getsize()

    def _getbits(self) -> int:
        """
        Returns the bits of the Bitset object.

        Returns:
            int: The bits of the Bitset object.

        Example:
            ```python
            bitset = Bitset()
            bits = bitset._getbits()
            print(bits)  # Output: 0
            ```
        """
        return self._bits

    @property
    def bits(self) -> int:
        """
        Returns the bits of the Bitset object.

        Returns:
            int: The bits of the Bitset object.

        Example:
            ```python
            bitset = Bitset()
            bits = bitset._getbits()
            print(bits)  # Output: 0
            ```
        """
        return self._getbits()

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

    def to_bool(self) -> bool:
        """
        Converts the Bitset object to a boolean value.

        Returns:
            bool: The boolean representation of the Bitset object.

        Example:
            ```python
            bitset = Bitset(4)
            is_true = bitset.to_bool()
            print(is_true)  # Output: False
            ```
        """
        return bool(self._bits)

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
