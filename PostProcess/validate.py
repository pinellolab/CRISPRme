""" """

from .utils import compute_md5, BFTARGETSMD5

from typing import List, NoReturn

import pandas as pd

import os

TARGETSDIR = "Results/crisprme-test-out/crispritz_targets"
TARGETSREPORT_REF = "hg38_20bp-NGG-SpCas9.txt_guides.txt_4_1_1.targets.txt"
TARGETSREPORT_1000G = "hg38+hg38_1000G_20bp-NGG-SpCas9.txt_guides.txt_4_1_1.targets.txt"
BFDIR = "crisprme-test-out/benchmark/brute-force-1000G"
BFTARGETS = "brute_force_1000G.tsv"
CHROMS = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9']


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

def main():






    # read crisprme and brute force targets
    crisprme_targets = pd.concat(
        [
            load_targets(os.path.join(TARGETSDIR, TARGETSREPORT_REF), crisprme=True),
            load_targets(os.path.join(TARGETSDIR, TARGETSREPORT_1000G), crisprme=True),
        ]
    )
    bf_targets = load_targets(os.path.join(BFDIR, BFTARGETS))
    # compute offtarget sites ids
    crisprme_sites = set(compute_sites(crisprme_targets, crisprme=True))
    bf_sites = set(compute_sites(bf_targets))
    if crisprme_sites != bf_sites:
        x = list(bf_sites.difference(crisprme_sites))[:20]
        print(bf_targets[bf_targets["site"].isin(x)].head())
        raise ValueError(f"Mismatch found between CRISPRme and brute force sites: CRISPRme: {len(crisprme_sites)} - Brute force: {len(bf_sites)}")
    print("Test passed!")

if __name__ == "__main__":
    main()