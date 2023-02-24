"""
"""

from version import __version__
from utils import CRISPRME_COMMANDS

from argparse import SUPPRESS, ArgumentParser, HelpFormatter
from typing import Dict, NoReturn, Optional, Tuple
from colorama import Fore

import sys


class CRISPRmeArgumentParser(ArgumentParser):
    """Class extending the original :class:`ArgumentParser` printing custom
    help messages for CRISPRme.

    :param ArgumentParser: object of type :class:`ArgumentParser`
    :type ArgumentParser: ArgumentParser
    """

    class CRISPRmeHelpFormatter(HelpFormatter):
        """This class extends :class:`HelpFormatter` from `argparse` package.
        This extension allows to print custom help messages.

        :param HelpFormatter: object of type :class:`HelpFormatter`
        :type HelpFormatter: HelpFormatter
        """

        def add_usage(
            self, usage: str, actions: str, groups: str, prefix: Optional[str] = "None"
        ) -> None:
            """Add the usage message to the HelpFormatter class object.

            Args:
                usage (str): usage description
                actions (str): actions
                groups (str): arguments groups
                prefix (Optional[str], optional): arguments prefixes. Defaults to "None".
            """
            if usage is not SUPPRESS:
                args = usage, actions, groups, ""
                self._add_item(self._format_usage, args)

    def __init__(self, *args: Tuple, **kwargs: Dict) -> None:
        """Constructor method"""
        kwargs["formatter_class"] = self.CRISPRmeHelpFormatter
        kwargs["usage"] = kwargs["usage"].replace("{version}", __version__)
        super().__init__(*args, **kwargs)

    def error(self, message: str) -> None:
        """Prints the error message and exits CRISPRme execution

        Args:
            message (str): error message
        """
        # recover command raising error
        command = sys.argv[1] if sys.argv[1] in CRISPRME_COMMANDS else ""
        if "invalid choice:" in message:
            message = message.replace("choice", "command")  # for clarity
        message = Fore.RED + "\nERROR: " + f"{message}." + Fore.RESET  # error in red
        message = message + f"\n\nRun crisprme.py {command} -h for usage\n\n"
        sys.stderr.write(message)  # write the error message
        sys.exit(2)

    def error_noargs(self) -> None:
        """When no arguments is parsed from the command line, prints the help
        message and exits.
        """
        self.print_help()
        sys.exit(2)
