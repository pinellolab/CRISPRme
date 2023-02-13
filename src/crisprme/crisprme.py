#!/usr/bin/env python

"""
CRISPRme version {version}

Copyright (C) 2023 Pinello lab <lpinello@mgh.harvard.edu>

CRISPRme is a tool for comprehensive off-target assessment, integrating human 
genetic variant datasets with orthogonal genomic annotations to predict and 
prioritize CRISPR-Cas off-target sites at scale.

usage: crisprme.py <command> [options]

IMPORTANT NOTE:

ALL FASTA FILEs USED BY THE SOFTWARE MUST BE UNZIPPED AND CHROMOSOME SEPARATED, 
ALL VCFs USED BY THE SOFTWARE MUST BE ZIPPED AND CHROMOSOME SEPARATED
"""

from crisprme_argparse import CRISPRmeArgumentParser
from subcommands import complete_search
from utils import CRISPRME_COMMANDS, CRISPRME_PATH, IUPAC_DNA, __version__

from Bio.Seq import Seq

import subprocess
import itertools
import sys
import os
import re



def parseargs_crisprme() -> CRISPRmeArgumentParser:
    """Parse the command line arguments (basic options)

    :return: :class:`CRISPRmeArgumentParser` object storing the command line arguments
    :rtype: CRISPRmeArgumentParser
    """
    parser = CRISPRmeArgumentParser(usage=__doc__, add_help=False)
    group = parser.add_argument_group("Basic options")
    group.add_argument("-h", "--help", action="help", help="Show this message and exit")
    group.add_argument("-v", "--version", action="version", help="Show CRISPRme version", version=__version__)
    # create different parsers for each command
    subparsers = parser.add_subparsers(
        title="Commands",
        description="Type 'crisprme.py <command> -h' to view the help of each command"
    )
    # complete-search arguments
    parser_completesearch = subparsers.add_parser(
        f"{CRISPRME_COMMANDS[0]}", 
        description="Automated off-targets search process, covering all analysis steps",
        usage=f"crisprme.py {CRISPRME_COMMANDS[0]}",
        help="Search the whole reference genome and variant genome (if requested) "
             "for off-targets occurrences, perform CFD, CRISTA and Mismatches+Bulges "
             "analyses, and targets selection",

    )
    group = parser_completesearch.add_argument_group("Options")
    group.add_argument(
        "--genome", 
        type=str, 
        required=True,
        metavar="GENOME-DIR", 
        help="Reference genome directory",
    )
    group.add_argument(
        "--vcf", 
        type=str, 
        metavar="VCF-FILE", 
        default="",
        nargs="?",  # optional arg
        help="File storing the list of VCF folders"
    )
    # --guide and --sequence options mutually exclude each other
    exclusive_group = group.add_mutually_exclusive_group()
    exclusive_group.add_argument(
        "--guide",
        type=str,
        metavar="GUIDE-FILE",
        help="File storing the input guides (excludes --sequence option)"
    )
    exclusive_group.add_argument(
        "--sequence",
        type=str,
        metavar="SEQUENCE-FILE",
        help="File storing the DNA sequences or genomic coordinates in BED "
             " (excludes --guide option)"
    )
    group.add_argument(
        "--pam",
        type=str,
        required=True,
        metavar="PAM-FILE",
        help="File storing the PAM sequences"
    )
    group.add_argument(
        "--be-window",
        type=str,
        metavar="START,STOP",
        dest="be_window",
        nargs="?",
        default="",
        help="Specify the window to search for base editing susceptibilty "
             "(example --be-window 4,8)"
    )
    group.add_argument(
        "--be-base",
        type=str,
        metavar="NUCLEOTIDE",
        dest="be_base",
        nargs="?",
        default="",
        help="Base(s) to check for the chosen base editor (example --be-base A,C)"
    )
    group.add_argument(
        "--annotation",
        type=str,
        metavar="ANNOTATION-FILE",
        help="File storing annotations for the reference genome"
    )
    group.add_argument(
        "--personal-annotation",
        type=str,
        metavar="PERSONAL-ANNOTATION",
        dest="personal_annotation",
        nargs="?",
        default="",
        help="specify the file that holds personal annotations for the "
             "reference genome"
    )
    group.add_argument(
        "--samples-id",
        type=str,
        metavar="SAMPLES-ID",
        dest="samples_id",
        nargs="?",
        default="",
        help="specify a file that lists, one per line, the information about "
             "samples present in the VCF files"
    )
    group.add_argument(
        "--gene-annotation",
        type=str, 
        metavar="GENE-ANNOTATION",
        dest="gene_annotation",
        nargs="?",
        default="",
        help="specify a gencode or similar annotation file, allowing the nearest "
             "gene to be determined for each target found"
    )
    group.add_argument(
        "--mm",
        type=int,
        required=True,
        metavar="MM",
        help="Number of allowed mismatches during the search"
    )
    group.add_argument(
        "--bDNA",
        type=int,
        metavar="BULGES-DNA",
        nargs="?",
        default=0,
        help="Number of DNA bulges allowed during the search"
    )
    group.add_argument(
        "--bRNA",
        type=int,
        metavar="BULGES-RNA",
        nargs="?",
        default=0,
        help="Number of RNA bulges allowed during the search"
    )
    group.add_argument(
        "--genome-index",
        type=str,
        metavar="GENOME-INDEX",
        dest="genome_index",
        nargs="?",
        default="",
        help="Folder containing a precomputed genome index (If provided, "
             "skip genome index construction)"
    )
    group.add_argument(
        "--merge",
        type=int,
        metavar="NUC-WINDOW",
        nargs="?",
        default=3,
        help="Specify the window (# of nucleotides) within which to merge "
             "candidate off-targets, using the off-target with the highest score "
             "as pivot"
    )
    group.add_argument(
        "--output",
        type=str,
        required=True,
        metavar="OUTPUT-NAME",
        help="Specify output name (Results written in CRISPRme/Results/<OUTPUT-NAME>)"
    )
    group.add_argument(
        "-t",
        "--threads",
        type=int,
        metavar="THREADS",
        nargs="?",
        default=8,
        help="Number of threads to use during the search (use 0 to autodetect "
             "and use all the available resources)"
    )
    group.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Run CRISPRme complete-search in debug mode"
    )
    # web-interface arguments
    parser_website = subparsers.add_parser(
        f"{CRISPRME_COMMANDS[4]}", 
        description="Starts a local web-server to use CRISPRme via a local web browser",
        usage=f"crisprme.py {CRISPRME_COMMANDS[4]}",
        help="Create and start a local server to run CRISPRme's web-interface. "
             "The website can be accessed locally through a web browser."
    )
    return parser

def main():
    parser = parseargs_crisprme()
    args = sys.argv[1:]  # read command or help options
    if len(args) == 0:
        parser.error_noargs()  # no arg seen, print help
        sys.exit(2)
    args = parser.parse_args(args)
    command = sys.argv[1]  # recover the called command
    assert command in CRISPRME_COMMANDS
    if command == CRISPRME_COMMANDS[0]:  # complete-search
        complete_search(parser, args)

if __name__ == "__main__":
    main()




