""" """

from complete_test import COMPLETETESTRESDIR
from utils import download, compute_md5, BFTARGETSMD5, CRISPRME_DIRS

from typing import List, Tuple

import pandas as pd

import sys
import os

# define crispritz targets folder
TARGETSDIR = f"{CRISPRME_DIRS[1]}/{COMPLETETESTRESDIR}/crispritz_targets"

# define crispritz targets reports (reference/alternative)
TARGETSREPORT_REF = "hg38_20bp-NGG-SpCas9.txt_guides.txt_4_1_1.targets.txt"
TARGETSREPORT_1000G = "hg38+hg38_1000G_20bp-NGG-SpCas9.txt_guides.txt_4_1_1.targets.txt"

# define brute-force targets folder
BFDIR = f"{CRISPRME_DIRS[1]}/{COMPLETETESTRESDIR}/benchmark/"

# define brute-force targets url and report
BFTARGETSURL = (
    "https://raw.githubusercontent.com/pinellolab/CRISPRme/refs/heads/v2.1.9/"
    "test/benchmark/brute-force-1000G/"
)
BFTARGETS = "brute_force_1000G.tsv"

# define allowed 
CHROMS = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', "chrX"]


def check_crispritz_targets() -> Tuple[str, str]:
    crispritz_targets = {
        "reference": os.path.join(TARGETSDIR, TARGETSREPORT_REF),
        "alternative (1000G)": os.path.join(TARGETSDIR, TARGETSREPORT_1000G),
    }
    missing_files = [f"{dataset_type}: {targets_path}" for dataset_type, targets_path in crispritz_targets.items() if not os.path.isfile(targets_path)]
    if missing_files:
        raise FileNotFoundError(
            "Missing required CRISPRitz target report file(s) for validation:\n"
            + "\n".join(f"  - {entry}" for entry in missing_files)
            + "\n\nPlease ensure that `complete-test` was executed successfully "
            "using 1000 Genomes variant data before running `validate-test`."
        )
    return crispritz_targets["reference"], crispritz_targets["alternative (1000G)"]
    
def check_variant_dataset() -> None:
    # read which dataset has been used for test in vcf config
    vcf_config_file = "vcf.config.test.txt"
    if not os.path.isfile(vcf_config_file):
        raise FileNotFoundError(
            f"Variant dataset configuration file not found: {vcf_config_file}"
        )
    with open(vcf_config_file, mode="r") as fin:
        datasets = [d for line in fin if (d := line.strip())]
    # check variant datasetn consistency (expected 1000G)
    if len(datasets) != 1:
        raise ValueError(
            "Invalid variant dataset configuration: exactly one dataset must be specified "
            f"(found {len(datasets)})."
        )
    if datasets[0] != "hg38_1000G":
        raise ValueError(
            "Unsupported variant dataset for validation.\n"
            "The `validate-test` functionality requires the variant dataset to be "
            "`hg38_1000G`.\n"
            f"Found: '{datasets[0]}'."
        )
    
def download_brute_force_targets() -> str:
    os.makedirs(BFDIR, exist_ok=True)  # create brute-force targets directory
    # construct download URL and download file
    bf_targets_url = f"{BFTARGETSURL.rstrip('/')}/{BFTARGETS}"
    bf_targets_path = download(BFDIR, http_url=bf_targets_url)
    # verify file existence
    if not bf_targets_path or not os.path.isfile(bf_targets_path):
        raise FileNotFoundError(
            "Failed to download brute-force targets file.\n"
            f"Expected file not found at: {bf_targets_path}"
        )
    # verify MD5 checksum
    md5 = compute_md5(bf_targets_path)
    if md5 != BFTARGETSMD5:
        raise ValueError(
            "MD5 checksum mismatch for brute-force targets file.\n"
            f"Expected: {BFTARGETSMD5}\n"
            f"Found:    {md5}\n"
            "The downloaded file may be corrupted or incomplete."
        )
    return bf_targets_path

def _load_targets(fname: str, chrom: str, crisprme: bool = False) -> pd.DataFrame:
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"Targets file not found: {fname}")
    df = pd.read_csv(fname, sep="\t")
    chrom_col = "Chromosome" if crisprme else "CHR"
    if chrom_col not in df.columns:
        raise ValueError(
            f"Invalid targets file format: missing required column '{chrom_col}'."
        )
    df = df[df[chrom_col].isin(CHROMS)]   # Keep only supported chromosomes
    if chrom != "all":  # filter by specific chromosome if requested
        if chrom not in CHROMS:
            raise ValueError(
                f"Invalid chromosome '{chrom}'. "
                f"Allowed values are 'all' or one of: {', '.join(CHROMS)}"
            )
        df = df[df[chrom_col] == chrom]
    return df  # type: ignore


def load_targets(crisprme_targets_ref_fname: str, crisprme_targets_alt_fname: str, bf_targets_fname: str, chrom: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # read crisprme and brute force targets
    crisprme_targets = pd.concat(
        [
            _load_targets(crisprme_targets_ref_fname, chrom, crisprme=True),
            _load_targets(crisprme_targets_alt_fname, chrom, crisprme=True),
        ]
    )
    bf_targets = _load_targets(bf_targets_fname, chrom)
    return crisprme_targets, bf_targets

def _compute_site(chrom: str, pos: int, strand: str) -> str:
    return "_".join(list(map(str, [chrom, pos, strand]))) 

def compute_sites(targets: pd.DataFrame, crisprme: bool = False) -> List[str]:
    if crisprme:
        targets["site"] = targets.apply(lambda x: _compute_site(x[3], x[4], x[6]), axis=1)
    else:
        targets["site"] = targets.apply(lambda x: _compute_site(x[0], x[4], x[3]), axis=1)
    return targets["site"].tolist()

def validate(crisprme_targets: pd.DataFrame, bf_targets: pd.DataFrame) -> None:
    # compute off-target site identifiers
    crisprme_sites = set(compute_sites(crisprme_targets, crisprme=True))
    bf_sites = set(compute_sites(bf_targets))
    # fast path: perfect match
    if crisprme_sites == bf_sites:
        sys.stderr.write(
            "Validation test passed: CRISPRme and brute-force off-target sets match.\n"
        )
        sys.stderr.write("Enjoy CRISPRme!\n")
        sys.exit(0)
    # detailed diagnostics on mismatch
    missing_in_crisprme = bf_sites - crisprme_sites
    extra_in_crisprme = crisprme_sites - bf_sites
    sys.stderr.write(
        "Validation test failed: mismatch detected between CRISPRme and "
        "brute-force off-target sites.\n"
    )
    sys.stderr.write(
        f"  CRISPRme sites:    {len(crisprme_sites)}\n"
        f"  Brute-force sites: {len(bf_sites)}\n"
        f"  Missing in CRISPRme: {len(missing_in_crisprme)}\n"
        f"  Extra in CRISPRme:   {len(extra_in_crisprme)}\n"
    )
    # print a small sample to help debugging
    max_examples = 5
    if missing_in_crisprme:
        sys.stderr.write(
            "  Example sites missing in CRISPRme:\n"
            + "\n".join(f"    - {s}" for s in list(missing_in_crisprme)[:max_examples])
            + "\n"
        )
    if extra_in_crisprme:
        sys.stderr.write(
            "  Example extra sites in CRISPRme:\n"
            + "\n".join(f"    - {s}" for s in list(extra_in_crisprme)[:max_examples])
            + "\n"
        )
    sys.exit(0)

def run_test_validation(chrom: str) -> None:
    # check crispritz targets existence and dataset used
    crisprme_targets_ref, crisprme_targets_alt = check_crispritz_targets()
    check_variant_dataset()
    # download brute-force targets on 1000G
    bf_targets_fname = download_brute_force_targets()
    # read and load crisprme and brute-force targets
    crisprme_targets, bf_targets = load_targets(crisprme_targets_ref, crisprme_targets_alt, bf_targets_fname, chrom)
    # validate offtarget sites
    validate(crisprme_targets, bf_targets)

def main():
    chrom = sys.argv[1]  # read command line args
    run_test_validation(chrom)  # validate results from complete-test 

if __name__ == "__main__":
    main()