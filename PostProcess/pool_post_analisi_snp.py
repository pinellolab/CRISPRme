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


def memory_capped_workers(requested, n_tasks):
    """Bound the number of concurrent post-analysis workers to a memory budget.

    Each per-chromosome worker loads that chromosome's SNP targets into memory,
    so running one worker per CPU makes peak RAM scale with the core count. On a
    genome-wide 1000G run this spiked to ~100 GB. We therefore cap concurrency so
    the estimated peak stays within a memory budget (default 64 GB). Both the
    budget and the per-worker estimate are overridable via environment variables
    (`CRISPRME_MAX_MEM_GB`, `CRISPRME_POSTPROC_WORKER_GB`) for large-memory or
    memory-constrained machines.
    """
    try:
        budget_gb = float(os.environ.get("CRISPRME_MAX_MEM_GB", "64"))
    except ValueError:
        budget_gb = 64.0
    try:
        per_worker_gb = float(os.environ.get("CRISPRME_POSTPROC_WORKER_GB", "4"))
    except ValueError:
        per_worker_gb = 4.0
    if per_worker_gb <= 0:
        per_worker_gb = 4.0
    cap = max(1, int(budget_gb // per_worker_gb))
    return max(1, min(requested, cap, n_tasks))


chroms = [
    os.path.splitext(os.path.basename(f))[0]
    for f in os.listdir(ref_folder)
    if f.endswith(".fa") and not f.endswith(".fai")
]

workers = memory_capped_workers(ncpus, len(chroms))
# NOTE: write this diagnostic to STDOUT (log_verbose.txt), never STDERR.
# The caller (submit_job_automated_new_multiple_vcfs.sh) treats a non-empty
# stderr log (`[ -s $logerror ]`) as a fatal post-analysis failure, so any
# informational text on stderr here would abort the run with a false error.
sys.stdout.write(
    f"Post-analysis SNPs: {workers} concurrent worker(s) "
    f"(cores={ncpus}, memory budget "
    f"{os.environ.get('CRISPRME_MAX_MEM_GB', '64')} GB)\n"
)
sys.stdout.flush()
with Pool(processes=workers) as pool:
    pool.map(start_analysis, chroms)
