#!/usr/bin/env python

from multiprocessing import Pool

import gzip
import subprocess
import sys
import os


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
    chrom = _normalize_chrom(_chrom_from_vcf(os.path.join(vcf_folder, fname)))
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
chrs = [f for f in os.listdir(vcf_folder) if f.endswith(".vcf.gz")]
with Pool(processes=ncpus) as pool:  # run chrom-wise post-analysis in parallel
    pool.map(start_analysis, chrs)
