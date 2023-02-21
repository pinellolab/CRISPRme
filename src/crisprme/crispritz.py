"""
"""

from utils import exception_handler, write
from crisprme_errors import CRISPRitzError

import subprocess
import os

# crispritz commands
CRISPRITZ = "crispritz.py"
ADDVARIANTS = "add-variants"
INDEXGENOME = "index-genome"
# crispritz outputs
ENRICHED_GENOME = "variants_genome"
GENOME_SNPS = os.path.join(ENRICHED_GENOME, "SNPs_genome")

def run_crispritz(run: str, command: str, verbosity: int, debug: bool) -> None:
    """Run the selected crispritz command option

    :param run: crispritz run 
    :type run: str
    :param command: selected command
    :type command: str
    :param verbosity: verbosity level
    :type verbosity: int
    :param debug: debug mode
    :type debug: bool
    """
    try:
        if verbosity > 2:
            code = subprocess.call(run, shell=True)
        else:
            code = subprocess.call(run, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        if code != 0:
            exception_handler(
                CRISPRitzError, f"An error occurred while running {CRISPRITZ} {command}", debug
            )
    except Exception:
        exception_handler(
            CRISPRitzError, "An error occurred during CRISPRitz analysis", debug
        )

def crispritz_add_variants(vcf_folder: str, genome_folder: str, verbosity: int, debug: bool) -> None:
    """The wrapper function invokes crispritz to generate a genome sequence that 
    enriched with the input genetic variants

    :param vcf_folder: VCFs folder
    :type vcf_folder: str
    :param genome_folder: reference genome folder
    :type genome_folder: str
    :param verbosity: verbosity level
    :type verbosity: int
    :param debug: debug mode
    :type debug: bool
    """
    if not isinstance(vcf_folder, str):
        exception_handler(
            TypeError, f"Expected {str.__name__}, got {type(vcf_folder).__name__}", debug
        )
    if not isinstance(genome_folder, str):
        exception_handler(
            TypeError, f"Expected {str.__name__}, got {type(genome_folder)}", debug
        )
    # build crispritz command 
    crispritz_run = f"{CRISPRITZ} {ADDVARIANTS} {vcf_folder} {genome_folder} true"
    run_crispritz(crispritz_run, ADDVARIANTS, verbosity, debug)  # run crispritz add-variants

def crispritz_index_genome(genome: str, genome_folder: str, pam: str, bulges: int, verbosity: int, debug: bool) -> None:
    """The wrapper function invokes crispritz to generate a genome index using  
    the input PAM and number of bulges

    :param genome: genome name
    :type genome: str
    :param genome_folder: genome folder
    :type genome_folder: str
    :param pam: pam file
    :type pam: str
    :param bulges: number of bulges 
    :type bulges: int
    :param verbosity: verbosity level
    :type verbosity: int
    :param debug: debug mode
    :type debug: bool
    """
    if not isinstance(genome, str):
        exception_handler(
            TypeError, f"Expected {str.__name__}, got {type(genome).__name__}", debug
        )
    if not os.path.isdir(genome_folder):
        exception_handler(FileNotFoundError, f"Unable to locate {genome_folder}", debug)
    if not os.path.isfile(pam):
        exception_handler(FileNotFoundError, f"Unable to locate {pam}", debug)
    if not isinstance(bulges, int):
        exception_handler(
            TypeError, f"Expected {int.__name__}, got {type(bulges).__name__}", debug
        )
    # build crispritz command
    crispritz_run = f"{CRISPRITZ} {INDEXGENOME} {genome} {genome_folder} {pam} -bMax {bulges} -th 1"  # force one thread 
    run_crispritz(crispritz_run, INDEXGENOME, verbosity, debug)  # run crispritz index-genome