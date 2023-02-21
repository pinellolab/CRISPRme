"""
"""

from crispritz import ENRICHED_GENOME, GENOME_SNPS, crispritz_add_variants
from indels_index import index_indels
from utils import CRISPRME_DIRS, LOG, exception_handler, move, remove_dir, write

from typing import List
from glob import glob
from time import ctime, time

import sys
import os

def run_complete_search(genome: str, vcf: str, pam: str, pam_file: str, bmax: int, guides: List[str], output: str, threads: int, verbosity: int, debug: bool) -> None:
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
                    write(f"Chromosomes found in {genome}: {','.join(chromosomes)}")
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
                # TODO: check existence of enriched genome
                if vcf:  # enrich the genome with crispritz
                    if verbosity > 0:
                        write(f"Started adding variants to the genome at {ctime()}")
                        if verbosity > 1:
                            start_addvariants = time()
                    # add variants to the input genome with crispritz
                    crispritz_add_variants(vcf_dataset, genome, verbosity, debug)
                    if verbosity > 0:
                        write(f"Ended adding variants at {ctime()}")
                        if verbosity > 1:
                            write("Adding variants took %.2fs" % (time() - start_addvariants))
                    # change name to variants genome
                    gname = os.path.basename(genome[:-1]) if genome.endswith("/") else os.path.basename(genome)
                    vname = os.path.basename(vcf_dataset[:-1]) if vcf_dataset.endswith("/") else os.path.basename(vcf_dataset)
                    source = os.path.join(GENOME_SNPS, f"{gname}_enriched")
                    variants_genome = f"{gname}+{vname}"
                    move(source, variants_genome, debug)
                    # store snps dictionaries
                    variants_dictionaries = os.path.join(CRISPRME_DIRS[2], f"dictionaries_{vname}")
                    if not os.path.isdir(CRISPRME_DIRS[2]):  # Dictionaries folder
                        os.mkdir(CRISPRME_DIRS[2])
                    if not os.path.isdir(variants_dictionaries):
                        os.mkdir(variants_dictionaries)
                    move(os.path.join(GENOME_SNPS, "*.json"), variants_dictionaries, debug)
                    # store indels dictionaries
                    assert os.path.isdir(CRISPRME_DIRS[2])
                    indels_dictionaries = os.path.join(CRISPRME_DIRS[2], f"log_indels_{vname}")
                    if not os.path.isdir(indels_dictionaries):
                        os.mkdir(indels_dictionaries)
                    move(os.path.join(GENOME_SNPS, "log*.txt"), indels_dictionaries, debug)
                    # start indexing indels
                    if verbosity > 0:
                        write(f"Started indels indexing at {ctime()}")
                        if verbosity > 1:
                            start_indel_index = time()
                    # genome library stores the results of indexed chromosomes
                    genome_library = "genome_library"  # genome library directory
                    if not os.path.isdir(genome_library):
                        os.mkdir(genome_library)
                    indels_folder = os.path.join(genome_library, f"{pam}_{bmax}_{gname}+{vname}_INDELS")
                    if not os.path.isdir(indels_folder):
                        os.mkdir(indels_folder)
                    # index indels by chromosome using crispritz
                    index_indels(ENRICHED_GENOME, pam_file, pam, gname, vname, bmax, threads, verbosity, debug)
                    if verbosity > 0:
                        write(f"Ended indels indexing at {ctime()}")
                        if verbosity > 1:
                            write("Indels indexing took %.2fs" % (time() - start_indel_index))
                    # store indels genome
                    variants_genome_indels = f"{variants_genome}_INDELS"
                    if not os.path.isdir(variants_genome_indels):
                        os.mkdir(variants_genome_indels)
                    move(os.path.join(ENRICHED_GENOME, "fake*"), variants_genome_indels, debug)
                    # delete temporary directories storing enriched genomes
                    remove_dir(ENRICHED_GENOME)

                    
                    
                

    except OSError:
        exception_handler(OSError, f"An error occurred while reading {vcf}", debug)
    



