"""Define static variables and utilities functions used throughout CRISPRme.
"""

from colorama import Fore, init
from typing import NoReturn, Optional

import tempfile
import sys
import os

__version__ = "2.1.1"

# --- static variables 
BUFSIZE = 1024 * 1024
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

def process_personal_annotation_line(line: str, debug: bool) -> str:
    """Append the '_personal' suffix to the fourth column (feature field) of the 
    personal annotation file in its corresponding BED format.

    :param line: BED file line
    :type line: str
    :param debug: debug mode
    :type debug: bool
    :raises TypeError: check line type
    :return: processed line
    :rtype: str
    """
    if not isinstance(line, str):
        exception_handler(
            TypeError, f"Expected {str.__name__}, got {type(line).__name__}", debug
        )
    line_split = line.strip().split()
    line_split[3] = f"{line_split[3]}_personal"  # modify 4th column
    line_processed = "\t".join(line_split)
    line_processed = line_processed.replace(" ", "\t")  # replace blank with tabs
    return f"{line_processed}\n"

def process_personal_annotation(personal_annotation: str, annotation: str, debug: bool, onlypann: Optional[bool] = False) -> str:
    """Combine the contents of the annotation file and the personal annotation 
    file while preserving the origin of each line, with lines from the personal 
    annotation file clearly identified

    :param personal_annotation: personal annotation file
    :type personal_annotation: str
    :param annotation: annotation file
    :type annotation: str
    :param debug: debug mode
    :type debug: bool
    :param onlypann: only personal annotation available, defaults to False
    :type onlypann: Optional[bool], optional
    :return: _description_
    :rtype: str
    """
    try:
        with open(personal_annotation, mode="r") as infile:
            # create temporary personal annotation file in /tmp
            personal_annotation_tmp = tempfile.NamedTemporaryFile().name
            with open(personal_annotation_tmp, mode="w") as outfile:
                while True:
                    lines = infile.readlines(BUFSIZE)  # read file in chunks of BUFSIZE
                    if not lines:
                        break
                    outfile.writelines(
                        process_personal_annotation_line(line, debug) for line in lines
                    )
    except OSError:
        exception_handler(
            OSError, "A problem occurred while processing personal annotation file", debug
        )
    # merge the two files into a single annotation file
    try:
        with open(personal_annotation_tmp, mode="r") as pann_infile:
            with open(annotation, mode="r") as ann_infile:
                with open(f"{annotation}+personal.bed", mode="w") as outfile:
                    outfile.write(pann_infile.read())
                    if not onlypann:  # use annotation only if available
                        outfile.write(ann_infile.read())
    except OSError:
        exception_handler(
            OSError, "A problem occurred while processing personal annotation file", debug
        )
    os.remove(personal_annotation_tmp)  # delete the tmp personal annotation file
    assert not os.path.isfile(personal_annotation_tmp)
    assert os.path.isfile(f"{annotation}+personal.bed")
    return f"{annotation}+personal.bed"

def check_directories(basedir: str) -> None:
    if not isinstance(basedir, str):
        raise TypeError(f"Expected {str.__name__}, got {type(basedir).__name__}")
    if not os.path.exists(basedir):
        raise FileNotFoundError(f"Unable to locate {basedir}")
    for d in CRISPRME_DIRS:
        if not os.path.exists(os.path.join(basedir, d)):
            os.makedirs(os.path.join(basedir, d))


