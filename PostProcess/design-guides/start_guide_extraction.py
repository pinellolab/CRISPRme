"""
"""

import argparse

from utils import reverse_complement
from encoder import encode_pam, encode_genome
from search import match_genome
from design_guides import design_guides
from reporter import write_report

parser = argparse.ArgumentParser(
    description="Given input fasta file, chromosome coordinates file (BED format [optional]) and PAM sequence, find all possible sgRNA guides in the fasta file (SUPPORTS IUPAC NOTATIONS)"
)
parser.add_argument(
    "-fa",
    "--fasta",
    type=str,
    help="fasta file of single chromosome OR target sequence extracted from genome",
    dest="sequence_file",
    required=True,
)
parser.add_argument(
    "-p",
    "--pam",
    type=str,
    help="PAM file (CRISPRitz/CRISPRme format)",
    dest="pam_file",
    required=True,
)
parser.add_argument(
    "-coo",
    "--coordinates",
    type=str,
    help="coordinates file (BED format, only one chromosome allowed at the same time)",
    dest="coordinate_file",
    required=False,
    default="",
)
parser.add_argument(
    "-mm",
    "--mismatches",
    type=int,
    help="allowed mismatches in the PAM sequence",
    dest="mm_max",
    required=False,
    default=0,
)
parser.add_argument(
    "-o",
    "--output_file",
    type=str,
    help="output file name, this file contains all the guides generated",
    dest="out_file_name",
    required=False,
    default="guides_out.txt",
)

args = parser.parse_args()

# start guide search, given the input PAM sequence
design_guides(
    args.pam_file,
    args.sequence_file,
    args.coordinate_file,
    args.mm_max,
    args.out_file_name,
)
