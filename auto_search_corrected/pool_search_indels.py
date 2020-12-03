#!/usr/bin/env python

import os
import sys
from multiprocessing import Pool
from datetime import datetime

ref_folder = sys.argv[1]
ref_name = os.path.basename(sys.argv[1])
vcf_dir = sys.argv[2]
vcf_name = sys.argv[3]
guide_file = sys.argv[4]
guide_name = os.path.basename(sys.argv[4])
pam_file = sys.argv[5]
pam_name = os.path.basename(sys.argv[5])
bMax = sys.argv[6]
mm = sys.argv[7]
bDNA = sys.argv[8]
bRNA = sys.argv[9]
output_folder = sys.argv[10]
true_pam = sys.argv[11]
cpus=int(sys.argv[12])

def search_indels(f):
    global use_thread
    splitted = f.split('.')
    chrom = splitted[1]
    print("Searching for INDELs in", chrom)
    os.system(f"crispritz.py search genome_library/{true_pam}_{bMax}_fake{chrom}/ {pam_file} {guide_file} fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA} -index -mm {mm} -bDNA {bDNA} -bRNA {bRNA} -t  -th 1 >/dev/null")
    print("Search ended for INDELs in", chrom)

chrs = []
for f in os.listdir(vcf_dir):
    if 'vcf.gz' in f:
        chrs.append(f)

#cpus = len(os.sched_getaffinity(0))
if cpus - 1 < 10:
    if cpus - 1 < 0:
        t = 1
    else:
        t = cpus - 1
else:
    t = 10

os.chdir(output_folder)
with Pool(processes = t) as pool:
    pool.map(search_indels, chrs)


for key in chrs:
    chrom = key.split('.')[1]
    os.system(f"tail -n +2 {output_folder}/fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt >> {output_folder}/indels_{ref_name}+{vcf_name}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt")
    header = os.popen(f"head -1 {output_folder}/fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt").read()
    os.system(f"rm {output_folder}/fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt")

os.system(f"sed -i 1i\"{header}\" {output_folder}/indels_{ref_name}+{vcf_name}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt")


#os.system('echo "Search INDELs End: '+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'" >> '+output_folder+'/../log.txt')