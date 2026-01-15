""" """

from .utils import compute_md5, BFTARGETSMD5

from typing import List, NoReturn

import pandas as pd

import sys
import os

# define crispritz targets folder
TARGETSDIR = "Results/crisprme-test-out/crispritz_targets"

# define crispritz targets reports (reference/alternative)
TARGETSREPORT_REF = "hg38_20bp-NGG-SpCas9.txt_guides.txt_4_1_1.targets.txt"
TARGETSREPORT_1000G = "hg38+hg38_1000G_20bp-NGG-SpCas9.txt_guides.txt_4_1_1.targets.txt"

# define brute-force targets folder
BFDIR = "crisprme-test-out/benchmark/brute-force-1000G"

# define brute-force targets report
BFTARGETS = "brute_force_1000G.tsv"

# define allowed 
CHROMS = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', "chrX"]


def check_crispritz_targets() -> None:
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


def load_targets(fname: str, crisprme: bool = False) -> pd.DataFrame:
    df = pd.read_csv(fname, sep="\t")
    if crisprme:
        df = df[df["Chromosome"].isin(CHROMS)]
    else:
        df = df[df["CHR"].isin(CHROMS)]
    return df

def _compute_site(chrom: str, pos: int, strand: str) -> str:
    return "_".join(list(map(str, [chrom, pos, strand]))) 

def compute_sites(targets: pd.DataFrame, crisprme: bool = False) -> List[str]:
    if crisprme:
        targets["site"] = targets.apply(lambda x: _compute_site(x[3], x[4], x[6]), axis=1)
    else:
        targets["site"] = targets.apply(lambda x: _compute_site(x[0], x[4], x[3]), axis=1)
    return targets["site"].tolist()

def run_test_validation(chrom: str) -> None:
    # check xrispritz targets existence and dataset used
    check_crispritz_targets()
    check_variant_dataset()

def main():
    chrom = sys.argv[1]  # read command line args
    run_test_validation(chrom)  # validate results from complete-test 




    # read crisprme and brute force targets
    # crisprme_targets = pd.concat(
    #     [
    #         load_targets(os.path.join(TARGETSDIR, TARGETSREPORT_REF), crisprme=True),
    #         load_targets(os.path.join(TARGETSDIR, TARGETSREPORT_1000G), crisprme=True),
    #     ]
    # )
    # bf_targets = load_targets(os.path.join(BFDIR, BFTARGETS))
    # # compute offtarget sites ids
    # crisprme_sites = set(compute_sites(crisprme_targets, crisprme=True))
    # bf_sites = set(compute_sites(bf_targets))
    # if crisprme_sites != bf_sites:
    #     x = list(bf_sites.difference(crisprme_sites))[:20]
    #     print(bf_targets[bf_targets["site"].isin(x)].head())
    #     raise ValueError(f"Mismatch found between CRISPRme and brute force sites: CRISPRme: {len(crisprme_sites)} - Brute force: {len(bf_sites)}")
    # print("Test passed!")

if __name__ == "__main__":
    main()