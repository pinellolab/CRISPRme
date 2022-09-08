"""TODO: write docstring

Functions required to run CRISPRme complete search
"""


from typing import Optional
from crisprme.crisprme_commands import CompleteSearch
from crisprme.utils import (
    exception_handler, 
    redir_stdout,
    redir_stderr,
    LOG_FILE, 
    LOG_VERBOSE,
    LOG_ERROR, 
    QUEUE_FILE,
    CRISPRme_DIRS, 
    CURRENT_WORKING_DIRECTORY,
)

from time import time
from glob import glob

import sys
import os


def run_complete_search(
    args: CompleteSearch, online: Optional[bool] = False
) -> None:
    """Run CRISPRme complete search using the input parameters selected by the 
    user.

    ...

    Parameters
    ----------
    args : CompleteSearch
        Parsed commandline arguments
    online : bool 
        Flag value to establish if the function was called by the online version
        of CRISPRme

    Returns
    -------
    None
    """

    if not isinstance(args, CompleteSearch):
        exception_handler(
            TypeError,
            f"Expected {CompleteSearch.__name__}, got {type(args).__name__}",
            args.debug
        )
    # create log.txt file storing run information
    log_file = os.path.join(
        CURRENT_WORKING_DIRECTORY, CRISPRme_DIRS[1], args.outname, LOG_FILE
    )
    if not os.path.exists(log_file):  # create log.txt if it doesn't exist
        try:
            open(log_file, mode="w")
        except:
            exception_handler(
                IOError, f"An error occurred while creating {LOG_FILE}", args.debug
            )
    if online:  
        # if the function was called form the online version then redirect
        # stdout and sterr
        # redirect stdout to log file called log_verbose.txt
        log_verbose = os.path.join(
                CURRENT_WORKING_DIRECTORY, CRISPRme_DIRS[1], args.outname, LOG_FILE
            )
        # create log_verbose.txt if it doesn't exist and write stdout to it
        try:
            sys.stdout = open(log_verbose, mode="w")
        except:
            exception_handler(
                IOError, f"An error occurred while creating {LOG_VERBOSE}", args.debug
            )
        # redirect stderr to log file called log_error.txt
        log_error = os.path.join(
            CURRENT_WORKING_DIRECTORY, CRISPRme_DIRS[1], args.outname, LOG_ERROR
        )
        # create log_error.txt if it doesn't exist and write stderr to it
        try:
            sys.stderr = open(log_error, mode="w")
        except:
            exception_handler(
                IOError, f"An error occurred while creating {LOG_ERROR}", args.debug
            )
    # remove the queue file and start the current job
    queue_file = os.path.join(
        CURRENT_WORKING_DIRECTORY, CRISPRme_DIRS[1], args.outname, QUEUE_FILE
    )
    # if not found could be ignored -> CRISPRme run from command line
    if os.path.isfile(queue_file):
        os.remove(queue_file)
    # keep track of used VCFs
    if args.vcf == "_":  # VCF not provided
        tmp_vcf_list = os.path.join(
                CURRENT_WORKING_DIRECTORY, 
                CRISPRme_DIRS[1], 
                args.outname, 
                "tmp_list_vcf.txt"
            )
        # add dummy name to list of employed VCFs
        try:
            handle = open(tmp_vcf_list, mode="w")
            handle.write("_\n")
        except IOError as ioerr:
            exception_handler(
                IOError, f"Unable to write to {tmp_vcf_list}", args.debug
            )
        finally:
            handle.close()
        handle = open(tmp_vcf_list, mode="r")
    else:  # recover VCF files
        handle = open(args.vcf, mode="r")  # assume one directory per line
    vcf_list = [line.strip() for line in handle if line.strip()]
    handle.close()
    # recover FASTA files within the reference genome directory
    chroms_fasta = glob(os.path.join(args.genome, "*.fa"))
    # recover the chromosome names list
    chroms_list = [
        os.path.basename(chrom).split(".")[0] for chrom in chroms_fasta
    ]
    # perform the analysis for each VCF in the list
    for vcf in vcf_list:
        vcf_dir = os.path.join(CURRENT_WORKING_DIRECTORY, CRISPRme_DIRS[3], vcf)
        start_analysis = time()
        sys.stderr.write(f"Starting analysis for {vcf}\n")
        
        # TODO: see fake chromosomes creation

        # recover PAM data
        pam_seq, pam_full, pam_start = args.pam
        

        
        



        stop_analysis = time()
