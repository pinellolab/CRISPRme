"""
"""

from utils import exception_handler, write
from crisprme_errors import CRISPRitzError

import subprocess
import os

CRISPRITZ = "crispritz.py"
ADDVARIANTS = "add-variants"

def crispritz_add_variants(vcf_folder: str, genome_folder: str, debug: bool) -> None:
    if not isinstance(vcf_folder, str):
        exception_handler(
            TypeError, f"Expected {str.__name__}, got {type(vcf_folder).__name__}", debug
        )
    if not isinstance(genome_folder, str):
        exception_handler(
            TypeError, f"Expected {str.__name__}, got {type(genome_folder)}", debug
        )
    crispritz_run = f"{CRISPRITZ} {ADDVARIANTS} {vcf_folder} {genome_folder} true"
    try:
        code = subprocess.run(crispritz_run, shell=True)
        if code != 0:
            exception_handler(
                CRISPRitzError, f"An error occurred while running {CRISPRITZ} {ADDVARIANTS}", debug
            )
    except Exception:
        exception_handler(
            CRISPRitzError, f"An error occurred during CRISPRitz analysis", debug
        )
    