"""Define static variables and utilities functions used throughout CRISPRme.
"""

from version import __version__

from argparse import Namespace
from colorama import Fore, init
from typing import NoReturn, Optional, Tuple
from time import ctime

import traceback
import resource
import tempfile
import logging
import shutil
import sys
import os

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
PAM_DICT = {
    "A":  "ARWMDHV",
    "C":  "CYSMBHV",
    "G":  "GRSKBDV",
    "T":  "TYWKBDH",
    "R":  "ARWMDHVSKBG",
    "Y":  "CYSMBHVWKDT",
    "S":  "CYSMBHVKDRG",
    "W":  "ARWMDHVYKBT",
    "K":  "GRSKBDVYWHT",
    "M":  "ARWMDHVYSBC",
    "B":  "CYSMBHVRKDGWT",
    "D":  "ARWMDHVSKBGYT",
    "H":  "ARWMDHVYSBCKT",
    "V":  "ARWMDHVYSBCKG",
    "N":  "ACGTRYSWKMBDHV",
}
LOG = "log.txt"
ERRORLOG = "crisprme_error.log"

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

def write_logerror() -> None:
    """Redirect the full error stack to the log error file for debugging purposes
    """
    errorlog = open(ERRORLOG, mode="w")  # open log error file
    # write run info
    errorlog.write(f"CRISPRme v{__version__}\n\n")
    errorlog.write(f"Command line:\n{' '.join(sys.argv[:])}\n\n")
    errorlog.write(f"Time: {ctime()}\n")
    errorlog.write(f"Memory usage: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024 / 1024} GB\n\n")
    traceback.print_exc(file=errorlog)  # trace the error stack

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

def add_n(guide: str, pam_len: str, pam_at_beginning: bool) -> str:
    """Add PAM_LEN Ns zs prefix or suffix to the guide sequence

    :param guide: guide sequence
    :type guide: str
    :param pam_len: PAM length
    :type pam_len: str
    :param pam_at_beginning: PAM occurs at the guide beginning
    :type pam_at_beginning: bool
    :return: complete guide sequence 
    :rtype: str
    """
    ns = "N" * pam_len
    if pam_at_beginning:
        return ns + guide
    return guide + ns

def write(message: str) -> None:
    """Write the message to stderr

    :param message: message
    :type message: str
    """
    sys.stderr.write(f"{message}\n")

def move(source: str, dest: str, debug: bool) -> None:
    """wrapper function to call the mv shell command, which allows you to move 
    or rename files in Unix-like operating systems

    :param source: source file/directory
    :type source: str
    :param dest: destination file/directory
    :type dest: str
    :param debug: debug mode
    :type debug: bool
    :raises OSError: an error occurs while moving the source to destination
    """
    try:
        os.system(f"mv {source} {dest}")
    except OSError:
        exception_handler(OSError, f"An error occurred while renaming/moving {source}", debug)

def remove_dir(folder: str, debug: bool) -> None:
    """Delete the input directory, regardless of wheter it's empty or not

    :param folder: folder to delete
    :type folder: str
    :param debug: debug mode
    :type debug: bool
    :raises ValueError: the input folder is not a directory
    :raises OSError: an error occurs while removing the directory
    """
    if not os.path.isdir(folder):
        exception_handler(ValueError, f"{folder} is not a directory", debug)
    try:
        shutil.rmtree(folder)
    except OSError:
        exception_handler(OSError, f"An error occurred while deleting {folder}", debug)
        
def check_directories(basedir: str) -> None:
    if not isinstance(basedir, str):
        raise TypeError(f"Expected {str.__name__}, got {type(basedir).__name__}")
    if not os.path.exists(basedir):
        raise FileNotFoundError(f"Unable to locate {basedir}")
    for d in CRISPRME_DIRS:
        if not os.path.exists(os.path.join(basedir, d)):
            os.makedirs(os.path.join(basedir, d))


