#!/usr/bin/env python

from multiprocessing import Pool

import os
import sys


output_folder = sys.argv[1]
ref_folder = sys.argv[2]
vcf_name = sys.argv[3]
guide_file = sys.argv[4]
mm = sys.argv[5]
bDNA = sys.argv[6]
bRNA = sys.argv[7]
# annotation_file = sys.argv[8]
annotations = sys.argv[8]
annotation_name = sys.argv[9]
pam_file = sys.argv[10]
# sampleID=sys.argv[10]
dict_folder = sys.argv[11]
final_res = sys.argv[12]
final_res_alt = sys.argv[13]
ncpus = int(sys.argv[14])


def start_analysis(f):
    # splitted = f.split('.')
    chrom = str(f).replace(".fa", "")
    # for elem in splitted:
    # if "chr" in elem:
    #     chrom = elem
    # os.system(
    #     f'./post_analisi_snp.sh "{output_folder}" "{ref_folder}" "{vcf_name}" "{guide_file}" "{mm}" "{bDNA}" "{bRNA}" {annotation_file} {pam_file} {dict_folder} {final_res} {final_res_alt} {chrom}'
    # )
    os.system(
        f'./post_analisi_snp.sh "{output_folder}" "{ref_folder}" "{vcf_name}" "{guide_file}" "{mm}" "{bDNA}" "{bRNA}" {annotations} {annotation_name} {pam_file} {dict_folder} {final_res} {final_res_alt} {chrom}'
    )


chrs = []
for f in os.listdir(ref_folder):
    if ".fa" in f and ".fai" not in f:
        chrs.append(f)

# t = 6
# if ncpus < 6:
#     t = ncpus
with Pool(processes=ncpus) as pool:
    pool.map(start_analysis, chrs)
