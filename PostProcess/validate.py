"""
This module validates CRISPRme off-target predictions against brute-force
benchmarks to ensure the correctness of complete-test results.

For every benchmark registered in ``test/benchmark/benchmarks.json`` it loads the
CRISPRme and brute-force target reports, harmonizes and compares off-target site
identifiers, and reports per-benchmark pass/fail. Adding a new benchmark (e.g. a
new nuclease/guide) requires only a registry entry plus its precomputed reference
TSV -- no code change here. `validate-test` exits non-zero if any present
benchmark fails.
"""

from complete_test import COMPLETETESTRESDIR
from utils import download, compute_md5, CRISPRME_DIRS

from typing import Dict, List, Optional, Tuple

import pandas as pd

import json
import sys
import os

# crispritz targets folder produced by complete-test
TARGETSDIR = f"{CRISPRME_DIRS[1]}/{COMPLETETESTRESDIR}/crispritz_targets"
# local (cached) folder for downloaded brute-force references
BFDIR = f"{CRISPRME_DIRS[1]}/{COMPLETETESTRESDIR}/benchmark/"
# benchmark registry (test/benchmark/benchmarks.json), resolved relative to repo
BENCHMARKS_JSON = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, "test", "benchmark",
    "benchmarks.json",
)
# fallback (canonical Cas9 benchmark) if the registry is unavailable
_FALLBACK_BENCHMARKS = {
    "reference_base_url": (
        "https://raw.githubusercontent.com/pinellolab/CRISPRme/refs/heads/main/"
        "test/benchmark/"
    ),
    "thresholds": {"mm": 4, "bDNA": 1, "bRNA": 1},
    "benchmarks": [
        {
            "name": "cas9_sg1617", "pam_name": "20bp-NGG-SpCas9.txt",
            "reference": "brute-force-1000G/brute_force_1000G.tsv",
            "md5": "2c1721f6c4586698bfd70b6ddbe34b8b",
        }
    ],
}

CHROMS = [
    "chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
    "chr17", "chr18", "chr19", "chr2", "chr20", "chr21", "chr22", "chr3",
    "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX",
]


def load_benchmarks() -> Dict:
    """Load the benchmark registry, falling back to the canonical Cas9 case."""
    try:
        with open(BENCHMARKS_JSON) as fin:
            return json.load(fin)
    except (OSError, ValueError):
        sys.stderr.write(
            f"WARNING: benchmark registry not found/parseable at {BENCHMARKS_JSON}; "
            "validating the canonical Cas9 benchmark only.\n"
        )
        return _FALLBACK_BENCHMARKS


def target_report_names(pam_name: str, mm: int, bdna: int, brna: int) -> Tuple[str, str]:
    """Reference and 1000G crispritz target report filenames for a benchmark."""
    suffix = f"{pam_name}_guides.txt_{mm}_{bdna}_{brna}.targets.txt"
    return f"hg38_{suffix}", f"hg38+hg38_1000G_{suffix}"


def find_crispritz_targets(pam_name: str, mm: int, bdna: int, brna: int) -> Optional[Tuple[str, str]]:
    """Return (ref, 1000G) target report paths if both exist, else None."""
    ref_name, alt_name = target_report_names(pam_name, mm, bdna, brna)
    ref_path = os.path.join(TARGETSDIR, ref_name)
    alt_path = os.path.join(TARGETSDIR, alt_name)
    if os.path.isfile(ref_path) and os.path.isfile(alt_path):
        return ref_path, alt_path
    return None


def check_variant_dataset() -> None:
    """Ensure complete-test used the 1000G dataset (required for validation)."""
    vcf_config_file = "vcf.config.test.txt"
    if not os.path.isfile(vcf_config_file):
        raise FileNotFoundError(
            f"Variant dataset configuration file not found: {vcf_config_file}"
        )
    with open(vcf_config_file, mode="r") as fin:
        datasets = [d for line in fin if (d := line.strip())]
    if len(datasets) != 1 or datasets[0] != "hg38_1000G":
        raise ValueError(
            "The `validate-test` functionality requires the variant dataset to be "
            f"`hg38_1000G`. Found: {datasets}."
        )


def download_brute_force_targets(base_url: str, reference: str, expected_md5: str) -> str:
    """Download (and cache) a brute-force reference, verifying its MD5."""
    os.makedirs(BFDIR, exist_ok=True)
    local = os.path.join(BFDIR, os.path.basename(reference))
    # cache: reuse a previously downloaded file if its md5 already matches
    if os.path.isfile(local) and compute_md5(local) == expected_md5:
        sys.stderr.write(f"Using cached brute-force reference: {local}\n")
        return local
    url = f"{base_url.rstrip('/')}/{reference}"
    sys.stderr.write(f"Downloading brute-force reference: {url}\n")
    path = download(BFDIR, http_url=url)
    if not path or not os.path.isfile(path):
        raise FileNotFoundError(f"Failed to download brute-force reference: {url}")
    md5 = compute_md5(path)
    if expected_md5 and expected_md5 != "TODO-after-generation" and md5 != expected_md5:
        raise ValueError(
            f"MD5 mismatch for {reference}. Expected {expected_md5}, found {md5}."
        )
    return path


def _load_targets(fname: str, chrom: str, crisprme: bool = False) -> pd.DataFrame:
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"Targets file not found: {fname}")
    df = pd.read_csv(fname, sep="\t")
    chrom_col = "Chromosome" if crisprme else "CHR"
    if chrom_col not in df.columns:
        raise ValueError(f"Invalid targets file: missing column '{chrom_col}'.")
    df = df[df[chrom_col].isin(CHROMS)]
    if chrom != "all":
        if chrom not in CHROMS:
            raise ValueError(f"Invalid chromosome '{chrom}'.")
        df = df[df[chrom_col] == chrom]
    return df  # type: ignore


def load_targets(ref_fname: str, alt_fname: str, bf_fname: str, chrom: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    crisprme_targets = pd.concat([
        _load_targets(ref_fname, chrom, crisprme=True),
        _load_targets(alt_fname, chrom, crisprme=True),
    ])
    return crisprme_targets, _load_targets(bf_fname, chrom)


def _compute_site(chrom: str, pos: int, strand: str) -> str:
    return "_".join(list(map(str, [chrom, pos, strand])))


def compute_sites(targets: pd.DataFrame, crisprme: bool = False) -> List[str]:
    if crisprme:
        targets["site"] = targets.apply(lambda x: _compute_site(x[3], x[4], x[6]), axis=1)
    else:
        targets["site"] = targets.apply(lambda x: _compute_site(x[0], x[4], x[3]), axis=1)
    return targets["site"].tolist()


def validate(name: str, crisprme_targets: pd.DataFrame, bf_targets: pd.DataFrame) -> bool:
    """Compare CRISPRme vs brute-force sites for one benchmark. Returns True if they match."""
    crisprme_sites = set(compute_sites(crisprme_targets, crisprme=True))
    bf_sites = set(compute_sites(bf_targets))
    if crisprme_sites == bf_sites:
        sys.stderr.write(f"[{name}] PASSED: CRISPRme and brute-force off-target sets match.\n")
        return True
    missing = bf_sites - crisprme_sites
    extra = crisprme_sites - bf_sites
    sys.stderr.write(
        f"[{name}] FAILED: mismatch.\n"
        f"  CRISPRme sites: {len(crisprme_sites)}  Brute-force sites: {len(bf_sites)}\n"
        f"  Missing in CRISPRme: {len(missing)}  Extra in CRISPRme: {len(extra)}\n"
    )
    for tag, s in (("missing in CRISPRme", missing), ("extra in CRISPRme", extra)):
        if s:
            sys.stderr.write(f"  Example {tag}:\n" +
                             "\n".join(f"    - {x}" for x in list(s)[:5]) + "\n")
    return False


def run_test_validation(chrom: str) -> None:
    """Validate every registered benchmark that has complete-test output present."""
    check_variant_dataset()
    registry = load_benchmarks()
    th = registry.get("thresholds", {"mm": 4, "bDNA": 1, "bRNA": 1})
    base_url = registry["reference_base_url"]
    validated = failed = skipped = 0
    for bench in registry["benchmarks"]:
        name = bench["name"]
        targets = find_crispritz_targets(bench["pam_name"], th["mm"], th["bDNA"], th["bRNA"])
        if targets is None:
            sys.stderr.write(f"[{name}] SKIPPED: complete-test targets not found "
                             "(this benchmark was not run).\n")
            skipped += 1
            continue
        ref_path, alt_path = targets
        bf_path = download_brute_force_targets(base_url, bench["reference"], bench.get("md5", ""))
        crisprme_targets, bf_targets = load_targets(ref_path, alt_path, bf_path, chrom)
        if validate(name, crisprme_targets, bf_targets):
            validated += 1
        else:
            failed += 1
    sys.stderr.write(
        f"\nValidation summary: {validated} passed, {failed} failed, {skipped} skipped.\n"
    )
    if validated == 0 and failed == 0:
        sys.stderr.write("ERROR: no benchmarks could be validated (no targets found).\n")
        sys.exit(1)
    if failed:
        sys.exit(1)
    sys.stderr.write("Enjoy CRISPRme!\n")
    sys.exit(0)


def main():
    run_test_validation(sys.argv[1])


if __name__ == "__main__":
    main()
