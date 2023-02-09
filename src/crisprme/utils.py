"""Define static variables and utilities functions used throughout CRISPRme.
"""

from colorama import Fore, init
from typing import NoReturn

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
def exception_handler(exception_type: Exception, exception: str, debug: bool) -> NoReturn:
    """Trigger a runtime error that halts execution, and optionally display a 
    full error stack trace if debugging is enabled

    :param exception_type: exception type
    :type exception_type: Exception
    :param exception: exception message
    :type exception: str
    :param debug: debug mode
    :type debug: bool
    :return: trigger runtime error
    :rtype: NoReturn
    """
    init()
    if debug:  # display full error stack
        raise exception_type(f"\n\n{exception}")
    # gracefully trigger runtime error and exit
    sys.stderr.write(Fore.RED + f"\n\nERROR: {exception}" + Fore.RESET)
    sys.exit(1)

def raise_warning(message: str) -> None:
    """Emit a yellow-colored warning message without interrupting program 
    execution

    :param message: warning message
    :type message: str
    :raises TypeError: check message type
    """
    if not isinstance(message, str):  # trace this error
        raise TypeError(f"Expected {str.__name__}, got {type(message).__name__}")
    message = Fore.YELLOW + f"Warning: {message}" + Fore.RESET
    sys.stderr.write(f"\n{message}\n")

def check_directories(basedir: str) -> None:
    if not isinstance(basedir, str):
        raise TypeError(f"Expected {str.__name__}, got {type(basedir).__name__}")
    if not os.path.exists(basedir):
        raise FileNotFoundError(f"Unable to locate {basedir}")
    for d in CRISPRME_DIRS:
        if not os.path.exists(os.path.join(basedir, d)):
            os.makedirs(os.path.join(basedir, d))


