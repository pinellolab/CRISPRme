"""
The `PAM` class represents a Protospacer Adjacent Motif sequence.

Args:
    pamfile (str): The path to the PAM file.

Example:
    ```python
    pamfile = "path/to/pam.txt"
    pam = PAM(pamfile)
    ```

---

The `PAM` class provides the following methods and properties:

- `__len__(self) -> int`: Returns the length of the PAM sequence.
- `__str__(self) -> str`: Returns the PAM sequence as a string.
- `file(self) -> str`: Returns the path to the PAM file.
- `pam(self) -> str`: Returns the PAM sequence.
- `length(self) -> int`: Returns the length of the PAM sequence.
- `pam_at_beginning(self) -> bool`: Returns a flag indicating if the PAM occurs at the beginning of the sequence.
- `guide_length(self) -> int`: Returns the length of the guide sequence.
- `full_length(self) -> int`: Returns the full length of the PAM sequence including the guide.

Please refer to the code for more details on each method and property.
"""

from typing import Tuple

import numpy as np


class PAM(object):
    """
    The `PAM` class represents a Protospacer Adjacent Motif sequence.

    Args:
        pamfile (str): The path to the PAM file.

    Example:
        ```python
        pamfile = "path/to/pam.txt"
        pam = PAM(pamfile)
        ```

    ---

    The `PAM` class provides methods and properties for working with PAM sequences. It allows for reading the PAM file, retrieving the PAM sequence, determining the length of the PAM, and more.

    Please refer to the code for more details on each method and property.
    """

    def __init__(self, pamfile: str) -> None:
        """
        Initializes a `PAM` object with the given PAM file.

        Args:
            pamfile (str): The path to the PAM file.

        Example:
            ```python
            pamfile = "path/to/pam.txt"
            pam = PAM(pamfile)
            ```

        ---

        The `__init__` method initializes a `PAM` object by reading the PAM file, establishing the PAM position, recovering the guide length, and recovering the PAM sequence.

        Please refer to the code for more details on the implementation.
        """

        self._pamfile = pamfile
        self._pam_full, self._pam_length = self._parse_pamfile()  # read PAM
        self._pam_at_beginning = self._pam_length < 0  # establish where PAM occurs
        self._guide_length = len(self._pam_full) - self.length  # recover guide length
        self._pam = self._recover_pam()

    def __len__(self) -> int:
        """
        Returns the length of the PAM sequence.

        Returns:
            int: The length of the PAM sequence.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            length = len(pam)
            print(length)  # Output: 3
            ```
        """

        return self.length

    def __str__(self) -> str:
        """
        Returns the PAM sequence as a string.

        Returns:
            str: The PAM sequence.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            sequence = str(pam)
            print(sequence)  # Output: 'NGG'
            ```
        """

        return self._pam

    def _parse_pamfile(self) -> Tuple[str, int]:
        """
        Parses the PAM file and returns the PAM sequence and its length as a tuple.

        Returns:
            Tuple[str, int]: A tuple containing the PAM sequence and its length.

        Raises:
            IOError: If the PAM file parsing fails or the file is empty.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            pam_sequence, pam_length = pam._parse_pamfile()
            print(pam_sequence, pam_length)  # Output: ('NGG', 3)
            ```
        """

        try:
            with open(self._pamfile, mode="r") as infile:
                while True:  # parse PAM file
                    line = infile.readline().strip()
                    if line:  # found PAM
                        pam, pam_length = line.split()[:2]
                        return pam, int(pam_length)
                raise ValueError(f"PAM file {pamfile} is empty!")
        except IOError as e:
            raise IOError("PAM file parsing failed!") from e

    def _recover_pam(self) -> str:
        """
        Recovers the PAM sequence based on the PAM position and the length of the PAM object.

        Returns:
            str: The recovered PAM sequence.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            recovered_pam = pam._recover_pam()
            print(recovered_pam)  # Output: 'NGG'
            ```
        """

        pam = (
            self._pam_full[: len(self)]
            if self.pam_at_beginning
            else self._pam_full[-len(self) :]
        )
        assert len(pam) == len(self)
        return pam

    def _get_pamfile(self) -> str:
        """
        Returns the path to the PAM file.

        Returns:
            str: The path to the PAM file.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            pamfile = pam._get_pamfile()
            print(pamfile)  # Output: 'path/to/pam.txt'
            ```
        """

        return self._pamfile

    @property
    def file(self) -> str:
        """
        Returns the path to the PAM file.

        Returns:
            str: The path to the PAM file.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            pamfile = pam.file
            print(pamfile)  # Output: 'path/to/pam.txt'
            ```
        """

        return self._get_pamfile()

    def _get_pam(self) -> str:
        """
        Returns the PAM sequence.

        Returns:
            str: The PAM sequence.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            pam_sequence = pam._get_pam()
            print(pam_sequence)  # Output: 'NGG'
            ```
        """

        return self._pam

    @property
    def pam(self) -> str:
        """
        Returns the PAM sequence.

        Returns:
            str: The PAM sequence.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            pam_sequence = pam.pam
            print(pam_sequence)  # Output: 'NGG'
            ```
        """

        return self._get_pam()

    def _get_length(self) -> int:
        """
        Returns the length of the PAM sequence.

        Returns:
            int: The length of the PAM sequence.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            length = pam._get_length()
            print(length)  # Output: 3
            ```
        """

        return np.abs(self._pam_length)

    @property
    def length(self) -> int:
        """
        Returns the length of the PAM sequence.

        Returns:
            int: The length of the PAM sequence.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            length = pam.length
            print(length)  # Output: 3
            ```
        """

        return self._get_length()

    def _get_pam_at_beginning(self) -> bool:
        """
        Returns a flag indicating whether the PAM occurs at the beginning of the sequence.

        Returns:
            bool: True if the PAM occurs at the beginning, False otherwise.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            at_beginning = pam._get_pam_at_beginning()
            print(at_beginning)  # Output: True
            ```
        """

        return self._pam_at_beginning

    @property
    def pam_at_beginning(self) -> bool:
        """
        Returns a flag indicating whether the PAM occurs at the beginning of the sequence.

        Returns:
            bool: True if the PAM occurs at the beginning, False otherwise.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            at_beginning = pam.pam_at_beginning
            print(at_beginning)  # Output: True
            ```
        """

        return self._get_pam_at_beginning()

    def _get_guide_length(self) -> int:
        """
        Returns the length of the guide sequence.

        Returns:
            int: The length of the guide sequence.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            guide_length = pam._get_guide_length()
            print(guide_length)  # Output: 0
            ```
        """

        return self._guide_length

    @property
    def guide_length(self) -> int:
        """
        Returns the length of the guide sequence.

        Returns:
            int: The length of the guide sequence.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            guide_length = pam.guide_length
            print(guide_length)  # Output: 0
            ```
        """

        return self._get_guide_length()

    def _get_full_length(self) -> int:
        """
        Returns the full length of the PAM sequence including the guide.

        Returns:
            int: The full length of the PAM sequence.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            full_length = pam._get_full_length()
            print(full_length)  # Output: 3
            ```
        """

        return self.length + self.guide_length

    @property
    def full_length(self) -> int:
        """
        Returns the full length of the PAM sequence including the guide.

        Returns:
            int: The full length of the PAM sequence.

        Example:
            ```python
            pam = PAM("path/to/pam.txt")
            full_length = pam.full_length
            print(full_length)  # Output: 3
            ```
        """

        return self._get_full_length()
