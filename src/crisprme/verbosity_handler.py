"""
"""

from crisprme_errors import VerbosityHandlerError
from utils import exception_handler, write

from time import time


def write_verbosity(message: str, verbosity: int, thresh: int, debug: bool) -> None:
    """Write the input message to stderr if the verbosity threshold is met

    :param message: message
    :type message: str
    :param verbosity: verbosity level
    :type verbosity: int
    :param thresh: verbosity threshold
    :type thresh: int
    :param debug: debug mode
    :type debug: bool
    :raises VerbosityHandlerError: raise on write() error
    """
    try:
        if verbosity > thresh:
            write(message)
    except Exception:
        exception_handler(
            VerbosityHandlerError,
            "An error occurred while handling verbosity write",
            debug,
        )


def time_verbosity(verbosity: int, thresh: int, debug: bool) -> float:
    """Retrieve current time if verbosity threshold is met

    :param verbosity: verbosity level
    :type verbosity: int
    :param thresh: verbosity threshold
    :type thresh: int
    :param debug: debug mode
    :type debug: bool
    :return: current time
    :rtype: float
    :raises VerbosityHandlerError: raise on time() error
    """
    try:
        if verbosity > thresh:
            return time()
    except Exception:
        exception_handler(
            VerbosityHandlerError,
            "An error occurred while handling verbosity time measuring",
            debug,
        )
    return .0
