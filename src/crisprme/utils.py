"""Define static variables and utilities functions used throughout CRISPRme.
"""

import sys
import os

__version__ = "2.1.1"

# --- static variables 
CRISPRME_PATH = os.path.dirname(os.path.abspath(__file__))
CRISPRME_DIRS = [
    "Genomes", "Results", "Dictionaries", "VCFs", "Annotations", "PAMs", "samplesIDs"
]
CRISPRME_COMMANDS = [
    "complete-search", "targets-integration", "gnomAD-converter", "generate-personal-card", "web-interface"
]
IUPAC_DNA = {
    "a",
    "A",
    "t",
    "T",
    "c",
    "C",
    "g",
    "G",
    "R",
    "Y",
    "S",
    "W",
    "K",
    "M",
    "B",
    "D",
    "H",
    "V",
    "r",
    "y",
    "s",
    "w",
    "k",
    "m",
    "b",
    "d",
    "h",
    "v",
}

# --- utils functions
def check_directories(basedir: str) -> None:
    """The function checks the consistency of CRISPRme's directory tree.
    If a directory is not found in the tree, it will be created.
    
    ...

    Parameters
    ----------
    basedir : str
        Base directory

    Returns
    -------
    None
    """

    if not isinstance(basedir, str):
        raise TypeError(f"Expected {str.__name__}, got {type(basedir).__name__}")
    if not os.path.exists(basedir):
        raise FileNotFoundError(f"Unable to locate {basedir}")
    for d in CRISPRME_DIRS:
        if not os.path.exists(os.path.join(basedir, d)):
            os.makedirs(os.path.join(basedir, d))


