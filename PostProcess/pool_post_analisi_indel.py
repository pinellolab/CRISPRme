#!/usr/bin/env python

from multiprocessing import Pool

import subprocess
import sys
import os

# post-analysis script name
POSTANALYSIS = "./post_analisi_indel.sh"

# read input arguments
output_folder = sys.argv[1]
ref_folder = sys.argv[2]
vcf_folder = sys.argv[3]
guide_file = sys.argv[4]
mm = sys.argv[5]
bDNA = sys.argv[6]
bRNA = sys.argv[7]
annotation_file = sys.argv[8]
pam_file = sys.argv[9]
dict_folder = sys.argv[10]
final_res = sys.argv[11]
final_res_alt = sys.argv[12]
ncpus = int(sys.argv[13])


def start_analysis(fname: str) -> None:
    chrom = next((e for e in fname.split(".") if e.startswith("chr")), None)
    assert chrom
    code = subprocess.call(
        f"{POSTANALYSIS} {output_folder} {ref_folder} {vcf_folder} {guide_file} "
        f"{mm} {bDNA} {bRNA} {annotation_file} {pam_file} {dict_folder} "
        f"{final_res} {final_res_alt} {chrom}",
        shell=True
    )
    if code != 0:
        raise subprocess.SubprocessError(
            f"Post-analysis on indels failed on chromosome {chrom}"
        )


# chromosome-wise vcfs list
chrs = [f for f in os.listdir(vcf_folder) if "vcf.gz" in f]
with Pool(processes=ncpus) as pool:  # run chrom-wise post-analysis in parallel
    pool.map(start_analysis, chrs)
