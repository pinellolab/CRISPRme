"""
"""

from crisprme_errors import CompleteSearchError
from crisprme_argparse import CRISPRmeArgumentParser
from parsers import parse_guide, parse_pam, parse_sequence
from utils import ERRORLOG, IUPAC_DNA, exception_handler, process_personal_annotation, raise_warning, write_logerror
from complete_search import run_complete_search

from argparse import Namespace
from colorama import Fore

import multiprocessing
import traceback
import logging
import os

logging.basicConfig(
    filename=ERRORLOG,
    level=logging.ERROR,
    format='%(asctime)s %(levelname)s %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p'
)

def complete_search(parser: CRISPRmeArgumentParser, args: Namespace) -> None:
    """_summary_

    :param parser: _description_
    :type parser: CRISPRmeArgumentParser
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
        if not os.path.exists(args.guide):  # check guide file existence
            parser.error(f"Unable to locate {args.sequence}")
        if not os.path.isfile(args.guide):
            parser.error(f"{args.guide} is not a file")
    if usesequence:
        if not os.path.exists(args.sequence):  # check sequence file existence
            parser.error(f"Unable to locate {args.sequence}")
        if not os.path.isfile(args.sequence):
            parser.error(f"{args.sequence} is not a file")
    # check genome argument consistency
    if not os.path.exists(args.genome):
        parser.error(f"Unable to locate {args.genome}")
    if not os.path.isdir(args.genome):
        parser.error(f"{args.genome} is not a directory")
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
        if not os.path.exists(args.vcf):
            parser.error(f"Unable to locate {args.vcf}")
        if not os.path.isfile(args.vcf):
            parser.error(f"{args.vcf} is not a file")
    # check gene-annotation consistency
    usegeneannotation = True if args.gene_annotation else False
    if usegeneannotation: 
        if not os.path.isfile(args.gene_annotation):
            parser.error(f"Unable to locate {args.gene_annotation}")
    # check pam consistency
    if not os.path.isfile(args.pam):
        parser.error(f"Unable to locate {args.pam}")
    # check functional annotation consistency
    if not args.annotation:  # TODO: should it be linked to empty.txt 
        if args.personal_annotation:  # use only personal annotation
            if not os.path.isfile(args.personal_annotation):
                parser.error(f"Unable to locate {args.personal_annotation}")
            raise_warning("only personal annotation provided")  # warn the user that only personal annotation data have been provided
            args.annotation = process_personal_annotation(
                args.personal_annotation, args.personal_annotation, args.debug, onlypann=True
            )
        else:
            raise_warning("annotation file not provided, skipping annotation")
    else:
        if not os.path.isfile(args.annotation):
            parser.error(f"Unable to locate {args.annotation}")
        # check if any personal-annotation file was provided
        if args.personal_annotation:
            if not os.path.isfile(args.personal_annotation):
                parser.error(f"Unable to locate {args.personal_annotation}")
            # combine annotation and personal annotation data
            args.annotation = process_personal_annotation(
                args.personal_annotation, args.annotation, args.debug
            )
    # check samples-id consistency
    if usevariants and not args.samples_id:
        parser.error("Sample IDs required when when searching with variants (--vcf option)")
    if not usevariants and bool(args.samples_id):
        parser.error("Sample IDs supplied, but the VCF files are missing")
    if args.samples_id:
        if not os.path.exists(args.samples_id):
            parser.error(f"Unable to locate {args.samples_id}")
        if not os.path.isfile(args.samples_id):
            parser.error(f"{args.samples_id} is not a file")
    # check mismatch argument consistency
    if args.mm < 0:
        parser.error(f"Forbidden number of mismatches ({args.mm})")
    # check DNA bulges consistency
    if args.bDNA < 0:
        parser.error(f"Forbidden number of DNA bulges ({args.bDNA})")
    # check RNA bulges consistency
    if args.bRNA < 0:
        parser.error(f"Forbidden number of RNA bulges ({args.bRNA})")
    bmax = max(args.bDNA, args.bRNA)  # number of bulges used during indexing
    # check genome-index consistency
    if args.genome_index:
        if not os.path.isfile(args.genome_index):
            parser.error(f"Unable to locate {args.genome_index}")
    # check merge consistency
    if args.merge < 0:
        parser.error(f"Forbidden merging window width ({args.merge})")
    # check output consistency
    if os.path.isdir(args.output):
        raise_warning(f"{args.output} already exists, some content could be lost")
    if not os.path.exists(args.output):  # create the directory
        os.makedirs(args.output)
    assert os.path.isdir(args.output)
    # check verbosity consistency
    if args.verbosity < 0 or args.verbosity > 3:
        parser.error(f"Forbidden verbosity value ({args.verbosity}). Allowed verbosity values: < 0|1|2|3 >")
    # recover genome directory name
    ref_genome = args.genome[:-1] if args.genome.endswith("/") else args.genome
    ref_genome = os.path.basename(ref_genome)
    # recover the PAM sequence
    pam, guide_expected_len, pam_at_beginning = parse_pam(args.pam, args.debug) 
    search_index = True if bmax !=0 else False  # construct the genome index
    compare_ref_genome = False  # compare results with those recovered on the ref
    if usevariants:  # enrich the genome
        assert args.vcf
        # list each VCF file contained in the input file (--vcf)
        genome_indexes = []
        with open(args.vcf, mode="r") as infile:
            for line in infile:
                line = line.strip()
                if line:
                    if line.endswith("/"):  # to recover basename
                        line = line[:-1]
                    vcfbasename = os.path.basename(line)
                    genome_indexes.append(f"{pam}_{bmax}_{ref_genome}+{vcfbasename}")
        genome_index = ",".join(genome_indexes)
        compare_ref_genome = True
    else:
        genome_index = f"{pam}_{bmax}_{ref_genome}"
    # TODO: when verbosity is high print this info 
    # TODO: for web-site write the parameters
    # TODO: recover nuclease for web-site
    if usesequence:  # sequence provided in input
        assert not useguide
        guides = parse_sequence(
            args.sequence, args.genome, pam, guide_expected_len, pam_at_beginning, args.debug
        )
    if useguide:  # guides provided in input
        assert not usesequence
        guides = parse_guide(args.guide, args.debug)
    # write the guides to a file stored in the output directory
    try:
        with open(os.path.join(args.output, "guides.txt"), mode="w") as outfile:
            for guide in guides:
                outfile.write(f"{guide}\n")
    except OSError:
        exception_handler(
            OSError, "An error occurred while writing the guides file", args.debug
        )
    # TODO: assign mail for web-site
    # TODO: verbosity level to launch the job
    # TODO: launch the job
    try:
        run_complete_search(args.genome, args.vcf, pam, args.pam, bmax, guides, args.output, args.threads, args.verbosity, args.debug)
    except Exception:
        write_logerror()
        exception_handler(
            CompleteSearchError, "An error occurred while running complete-search", True
        )
    
    
    




        