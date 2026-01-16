#!/usr/bin/env python

from multiprocessing import Pool

import subprocess
import os
import sys

output_folder = sys.argv[1]
ref_folder = sys.argv[2]
vcf_name = sys.argv[3]
guide_file = sys.argv[4]
mm = sys.argv[5]
bDNA = sys.argv[6]
bRNA = sys.argv[7]
annotation_file = sys.argv[8]
pam_file = sys.argv[9]
# sampleID=sys.argv[10]
dict_folder = sys.argv[10]
final_res = sys.argv[11]
final_res_alt = sys.argv[12]
ncpus = int(sys.argv[13])


def start_analysis(chrom):
    cmd = f'./post_analisi_snp.sh "{output_folder}" "{ref_folder}" "{vcf_name}" "{guide_file}" "{mm}" "{bDNA}" "{bRNA}" {annotation_file} {pam_file} {dict_folder} {final_res} {final_res_alt} {chrom}'
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError(f"Post-analysis SNP failed on chromsomes {chrom}")

        

chroms = [
    os.path.splitext(os.path.basename(f))[0]
    for f in os.listdir(ref_folder)
    if f.endswith(".fa") and not f.endswith(".fai") 
]

with Pool(processes=ncpus) as pool:
    pool.map(start_analysis, chroms)
