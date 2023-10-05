"""
"""

import argparse

from read_input import read_sequence,read_pam,read_coordinate
from utils import reverse_complement
from encoder import encode_pam, encode_genome
from search import match_genome
from design_guides import extract_guides_from_genome
from output import write_results
# from .read_input import read_pam,read_sequence,read_coordinate


parser = argparse.ArgumentParser(description='Given input fasta file, chromosome coordinates file (BED format [optional]) and PAM sequence, find all possible sgRNA guides in the fasta file (SUPPORTS IUPAC NOTATIONS)')
parser.add_argument('-fa','--fasta', type=str, help='fasta file of single chromosome OR target sequence extracted from genome',dest='sequence_file',required=True)
parser.add_argument('-p','--pam', type=str, help='PAM file (CRISPRitz/CRISPRme format)',dest='pam_file',required=True)
parser.add_argument('-coo','--coordinates', type=str, help='coordinates file (BED format, only one chromosome allowed at the same time)',dest='coordinate_file',required=False,default=None)
parser.add_argument('-mm','--mismatches', type=int, help='allowed mismatches in the PAM sequence',dest='mm_max',required=False,default=0)
parser.add_argument('-o','--output_file', type=str, help='output file name, this file contains all the guides generated',dest='out_file_name',required=False,default="guides_out.txt")

args=parser.parse_args()

pam=read_pam(args.pam_file) ##pamseq,pam_len,guide_len,pam_in_start
seq=read_sequence(args.sequence_file) ##sequence,chr_name

if args.coordinate_file is not None:
    coordinates=read_coordinate(args.coordinate_file) ##list_of_coordinates (chr,start,end)
    seq_list=list()
    for cord in coordinates:
        seq_list.append(seq[0][cord[1]:cord[2]])
    seq[0]="N".join(seq_list)

# print("seq-coordinates-filtered: ",seq[0])

pam_bits = encode_pam(pam[0]), encode_pam(reverse_complement(pam[0]))
# print("pam_bits: ",pam_bits)
genome_bits = encode_genome(seq[0])

pam_positions = match_genome(genome_bits, pam_bits, pam[1], pam[1]+pam[2], args.mm_max, pam[3])
# print("pam_positions: ",pam_positions)

guides=extract_guides_from_genome(pam_positions,seq[0],pam[2],pam[1],pam[3])

write_results(guides,args.out_file_name,pam_positions,seq[1])
