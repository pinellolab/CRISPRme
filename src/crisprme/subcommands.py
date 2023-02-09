"""
"""

from crisprme_argparse import CRISPRmeArgumentParser
from utils import IUPAC_DNA, raise_warning

from argparse import Namespace
from colorama import Fore

import multiprocessing
import os

def complete_search(parser: CRISPRmeArgumentParser, args: Namespace) -> None:
    """_summary_

    :param args: _description_
    :type args: Namespace
    """
    # check base editing arguments consistency
    if args.be_window and not args.be_base:
        parser.error(
            "Base(s) editor missing, when specifying the base editing window the " 
            "base(s) editor is required"
        )
    if args.be_base and not args.be_window:
        parser.error(
            "Base editing window missing, when specifying the base editor the "
            "base editing window is required"
        )
    if args.be_window:  # check base-editing window
        start, stop = (1, 0)
        try:
            window = list(map(int, args.be_window.split(",")))
        except ValueError:
            parser.error(
                f"Invalid base editing window ({args.be_window}), example --be-window 4,8"
            )
    if args.be_base:
        unknown_nts = set(args.be_base.split(",")).difference(IUPAC_DNA)
        if bool(unknown_nts):
            parser.error(f"Unknown bases ({','.join(list(unknown_nts))})")
    # check guide and sequence arguments consistency
    if not args.guide and not args.sequence:  # one of the two is required
        parser.error(
            "Either a guide or sequence file must be specified, but both are missing"
        ) 
    useguide = True if args.guide else False 
    usesequence = True if args.sequence else False
    assert useguide or usesequence  # at least one of the two should be used
    assert (useguide + usesequence) == 1  # only one of the two
    if useguide:
        if not os.path.isfile(args.guide):  # check guide file existence
            parser.error(f"Unable to locate {args.sequence}")
    if usesequence:
        if not os.path.isfile(args.sequence):  # check sequence file existence
            parser.error(f"Unable to locate {args.sequence}")
    # check genome argument consistency
    if not os.path.isdir(args.genome):
        parser.error(f"Unable to locate {args.genome}")
    # check threads number consistency
    assert isinstance(args.threads, int)
    if args.threads < 0:
        warnmessage = str(
            f"forbidden number of threads requested ({args.threads}), reverting "
            f"to default ({8})" 
        ) 
        raise_warning(warnmessage)
        args.threads = 8
    if args.threads == 0:  # reserve all resources 
        args.threads = multiprocessing.cpu_count()
    # check VCF consistency
    usevariants = True if args.vcf else False  # wheter to use or not variants in the search
    if usevariants:
        assert bool(args.vcf)
        if not os.path.isfile(args.vcf):
            parser.error(f"Unable to locate {args.vcf}")
    # check gene-annotation consistency
    usegeneannotation = True if args.gene_annotation else False
    if usegeneannotation: 
        if not os.path.isfile(args.gen_annotation):
            parser.error(f"Unable to locate {args.gene_annotation}")
    # check pam consistency
    if not os.path.isfile(args.pam):
        parser.error(f"Uanble to locate {args.pam}")