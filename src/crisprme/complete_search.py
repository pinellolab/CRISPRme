"""
"""

from crispritz import crispritz_add_variants
from utils import LOG, exception_handler, write

from typing import List
from glob import glob
from time import ctime, time

import os

def run_complete_search(genome: str, vcf: str, pam: str, guides: List[str], output: str, verbosity: int, debug: bool) -> None:
    """_summary_

    :param genome: _description_
    :type genome: str
    :param vcf: _description_
    :type vcf: str
    :param guides: _description_
    :type guides: List[str]
    :param output: _description_
    :type output: str
    :param verbosity: _description_
    :type verbosity: int
    :param debug: _description_
    :type debug: bool
    """
    # TODO: base editor check in resultintegration
    # open the log file
    try:
        logfile = open(os.path.join(output, LOG), mode="w")
    except OSError:
        exception_handler(
            OSError, f"An error occurred while writing to {os.path.join(output, LOG)}", debug
        )
    if verbosity > 0:  # write job start time
        write(f"Job start {ctime()}")
    # TODO: remove queue file
    # parse VCF file
    if verbosity > 2:
        write(f"Opening VCF list file {vcf}")
    if verbosity > 2:
        start_vcf_parsing = time()
        write(f"Start parsing {vcf} content")
    try:
        with open(vcf, mode="r") as vcfinfile:
            for line in vcfinfile:
                vcf_dataset = line.strip()
                if not vcf_dataset:  # skip empty lines
                    continue 
                if verbosity > 0:
                    write(f"Starting analyis on {vcf_dataset}")
                # TODO: write start time to log file
                # recover chromosomes FASTA files
                chromosomes = [
                    os.path.splitext(os.path.basename(chrom_fasta))[0]
                    for chrom_fasta in glob(os.path.join(genome, "*.fa"))
                ]
                if verbosity > 2:
                    write(f"Chromosomes read: {', '.join(chromosomes)}")
                if vcf:  # chromosomes with indels
                    if verbosity > 2:
                        write("Reading VCFs associated to each chromosome")
                    chromosomes_indels = []
                    for chrom_vcf in glob(os.path.join(vcf_dataset, "*.vcf.gz")):
                        chrom = [
                            p for p in os.path.basename(chrom_vcf).split(".") if "chr" in p
                        ][0]
                        chromosomes_indels.append(f"fake{chrom}")
                # TODO: create directories structure for web-site?
                if vcf:  # enrich the genome with crispritz
                    if verbosity > 0:
                        write(f"Started adding variants to the genome at {ctime}")
                        if verbosity > 1:
                            start_addvariants = time()
                    # add variants to the input genome with crispritz
                    crispritz_add_variants(vcf_dataset, genome, debug)
                    # TODO: move output directory for SNPs
                    if verbosity > 0:
                        write(f"Ended adding variants at {ctime}")
                        if verbosity > 1:
                            write("Adding variants took %.2fs" % (time() - start_addvariants))
                    # start indexing indels
                    if verbosity > 0:
                        write(f"Started indexing indels at {ctime}")
                        if verbosity > 1:
                            start_indel_index = time()
                    
                    
                    
                

    except OSError:
        exception_handler(OSError, f"An error occurred while reading {vcf}", debug)
    



