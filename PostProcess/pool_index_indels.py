#!/usr/bin/env python

import os
import sys
from multiprocessing import Pool
from datetime import datetime

indels_folder = sys.argv[1]
pam_file = sys.argv[2]
true_pam = sys.argv[3]
ref_name = sys.argv[4]
vcf_name = sys.argv[5]
ncpus = sys.argv[6]


def index_indels(chrom):
    print("Indexing INDELs in", chrom)
    os.system(
        f"crispritz.py index-genome {ref_name}+{vcf_name}_INDELS/{true_pam}_2_fake{chrom} {indels_folder}/fake_{vcf_name}_{chrom} {pam_file} -bMax 2 -th 1  >/dev/null"
    )  # {indels_folder}/fake_{vcf_name}_{chrom}
    print("Indexing ended for INDELs in", chrom)


chrs = []
for f in os.listdir(indels_folder):
    if "chr" in f:
        chrs.append(f.split("_")[-1])

with Pool(processes=int(ncpus)) as pool:
    pool.map(index_indels, chrs)


# os.system('echo "Indexing INDELs End: '+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'" >> '+output_folder+'/../log.txt')
