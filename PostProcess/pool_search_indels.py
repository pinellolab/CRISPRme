#!/usr/bin/env python

import gzip
import os
import sys
from multiprocessing import Pool
from datetime import datetime


def _chrom_from_vcf(vcf_path: str) -> str:
    """Return the chromosome name from the first data line of a VCF (plain or gzipped)."""
    open_fn = gzip.open if vcf_path.endswith(".gz") else open
    with open_fn(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                return line.split("\t")[0]
    raise ValueError(f"No data lines found in VCF: {vcf_path}")


def _normalize_chrom(chrom: str) -> str:
    """Ensure chromosome name has chr prefix (e.g. '22' -> 'chr22').
    Some VCF datasets (e.g. 1000G GRCh38) store chromosomes without the
    prefix in the CHROM field, while the genome indices use 'chr'-prefixed names.
    """
    return chrom if chrom.startswith("chr") else "chr" + chrom

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
current_working_directory = sys.argv[12]
threads = int(sys.argv[13])


def search_indels(f):
    chrom = _normalize_chrom(_chrom_from_vcf(os.path.join(vcf_dir, f)))
    print("Searching for INDELs in", chrom)
    if bDNA != "0" or bRNA != "0":
        os.system(
            f"crispritz.py search {current_working_directory}/genome_library/{true_pam}_{bMax}_{ref_name}+{vcf_name}_INDELS/{true_pam}_{bMax}_fake{chrom}/ {pam_file} {guide_file} fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA} -index -mm {mm} -bDNA {bDNA} -bRNA {bRNA} -t -th 1 >/dev/null"
        )
    else:
        os.system(
            f"crispritz.py search {current_working_directory}/Genomes/{ref_name}+{vcf_name}_INDELS/fake_{vcf_name}_{chrom}/ {pam_file} {guide_file} fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA} -mm {mm} -t -th 1 >/dev/null"
        )
    print("Search ended for INDELs in", chrom)


chrs = []
for f in os.listdir(vcf_dir):
    if "vcf.gz" == f[-6:]:
        chrs.append(f)

# cpus = len(os.sched_getaffinity(0))
# if cpus - 3 < 10:
#     if cpus - 3 < 0:
#         t = 1
#     else:
#         t = cpus - 3
# else:
#     t = 10

os.chdir(output_folder)
# with Pool(processes=t) as pool:
with Pool(processes=threads) as pool:
    pool.map(search_indels, chrs)


for key in chrs:
    chrom = _normalize_chrom(_chrom_from_vcf(os.path.join(vcf_dir, key)))
    os.system(
        f"tail -n +2 {output_folder}/fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt >> {output_folder}/indels_{ref_name}+{vcf_name}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt"
    )
    header = os.popen(
        f"head -1 {output_folder}/fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt"
    ).read()
    os.system(
        f"rm {output_folder}/fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt"
    )

os.system(
    f'sed -i 1i"{header}" {output_folder}/indels_{ref_name}+{vcf_name}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt'
)


# os.system('echo "Search INDELs End: '+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'" >> '+output_folder+'/../log.txt')
