"""
"""

import argparse

from read_input import read_sequence,read_pam,read_coordinate
# from .read_input import read_pam,read_sequence,read_coordinate


parser = argparse.ArgumentParser(description='Given input fasta file, chromosome coordinates file (BED format [optional]) and PAM sequence, find all possible sgRNA guides in the fasta file (SUPPORTS IUPAC NOTATIONS)')
parser.add_argument('--fasta', type=str, help='fasta file of single chromosome OR target sequence extracted from genome',dest='sequence_file',required=True)
parser.add_argument('--pam', type=str, help='PAM file (CRISPRitz/CRISPRme format)',dest='pam_file',required=True)
parser.add_argument('--coordinate', type=str, help='coordinate file (BED format, only one chromosome allowed at the same time)',dest='coordinate_file',required=False)

#parser.add_argument('--output', type=str, help='output file',dest='output_file',default=None)
args=parser.parse_args()
# print(args)

pam=read_pam(args.pam_file)
print(pam)

seq=read_sequence(args.sequence_file)
print(seq[0:30])
