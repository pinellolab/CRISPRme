"""
"""

from typing import Tuple

import numpy as np


class PAM(object):
    def __init__(self, pamfile: str) -> None:
        self._pamfile = pamfile
        self._pam_full, self._pam_length = self._parse_pamfile()  # read PAM
        # establish where PAM occurs
        self._pam_at_beginning = True if self._pam_length < 0 else False
        self._guide_length = len(self._pam_full) - self.length  # recover guide length
        self._pam = self._recover_pam()

    def __len__(self) -> int:
        return self.length

    def __str__(self) -> str:
        return self._pam

    def _parse_pamfile(self) -> Tuple[str, int]:
        try:
            with open(self._pamfile, mode="r") as infile:
                while True:  # parse PAM file
                    line = infile.readline().strip()
                    if line:  # found PAM
                        pam, pam_length = line.split()[:2]
                        return pam, int(pam_length)
                raise ValueError(f"PAM file {pamfile} is empty!")
        except IOError as e:
            raise IOError(f"PAM file parsing failed!") from e

    def _recover_pam(self) -> str:
        pam = (
            self._pam_full[: len(self)]
            if self.pam_at_beginning
            else self._pam_full[-len(self) :]
        )
        assert len(pam) == len(self)
        return pam

    def _get_pamfile(self) -> str:
        return self._pamfile

    @property
    def file(self) -> str:
        return self._get_pamfile()

    def _get_pam(self) -> str:
        return self._pam

    @property
    def pam(self) -> str:
        return self._get_pam()

    def _get_length(self) -> int:
        return np.abs(self._pam_length)

    @property
    def length(self) -> int:
        return self._get_length()

    def _get_pam_at_beginning(self) -> bool:
        return self._pam_at_beginning

    @property
    def pam_at_beginning(self) -> bool:
        return self._get_pam_at_beginning()

    def _get_guide_length(self) -> int:
        return self._guide_length

    @property
    def guide_length(self) -> int:
        return self._get_guide_length()

    def _get_full_length(self) -> int:
        return self.length + self.guide_length

    @property
    def full_length(self) -> int:
        return self._get_full_length()
