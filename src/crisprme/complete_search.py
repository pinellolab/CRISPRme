"""
"""

from crisprme_argparse import CRISPRmeArgumentParser

from typing import List


def parseargs_completesearch() -> CRISPRmeArgumentParser:
    """Parse the command line arguments (complete-search options)

    :return: :class:`CRISPRmeArgumentParser` object storing the command line arguments
    :rtype: CRISPRmeArgumentParser
    """
    parser = CRISPRmeArgumentParser(
        description="CRISPRme version {version}\n\nCopyright (C) 2023 Pinello "
                    "lab <lpinello@mgh.harvard.edu>\n\nIMPORTANT NOTE:\n\nALL "
                    "FASTA FILEs USED BY THE SOFTWARE MUST BE UNZIPPED AND "
                    "CHROMOSOME SEPARATED, ALL VCFs USED BY THE SOFTWARE MUST "
                    "BE ZIPPED AND CHROMOSOME SEPARATED",
        usage="crisprme.py complete-search"
    )
    group = parser.add_argument_group("Complete-search options")
    # group.add_argument("-h", "--help", action="help", help="Show this message and exit")

def complete_search(args: List) -> None:
    parser = parseargs_completesearch()
    args = parser.parse_args(args)




