"""
This module validates CRISPRme off-target predictions against a brute-force
benchmark to ensure the correctness of complete-test results.

It loads CRISPRme and brute-force target reports, harmonizes and compares
off-target site identifiers across datasets, and reports detailed diagnostics
for any mismatches found during validation.
"""

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
    "https://raw.githubusercontent.com/pinellolab/CRISPRme/refs/heads/main/"
    "test/benchmark/brute-force-1000G/"
)
BFTARGETS = "brute_force_1000G.tsv"

# define allowed
CHROMS = [
    "chr1",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr2",
    "chr20",
    "chr21",
    "chr22",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chrX",
]


def check_crispritz_targets() -> Tuple[str, str]:
    """
    Locate and verify the CRISPRitz target report files required for validation.

    This function checks for the presence of both reference and alternative
    (1000G) CRISPRitz targets produced by complete-test and returns their paths
    if they are available.

    Returns:
        Tuple[str, str]: The file paths for the reference and alternative
            (1000G) CRISPRitz target reports.

    Raises:
        FileNotFoundError: If one or both required CRISPRitz target report
            files are missing.
    """
    crispritz_targets = {
        "reference": os.path.join(TARGETSDIR, TARGETSREPORT_REF),
        "alternative (1000G)": os.path.join(TARGETSDIR, TARGETSREPORT_1000G),
    }
    missing_files = [
        f"{dataset_type}: {targets_path}"
        for dataset_type, targets_path in crispritz_targets.items()
        if not os.path.isfile(targets_path)
    ]
    if missing_files:
        raise FileNotFoundError(
            "Missing required CRISPRitz target report file(s) for validation:\n"
            + "\n".join(f"  - {entry}" for entry in missing_files)
            + "\n\nPlease ensure that `complete-test` was executed successfully "
            "using 1000 Genomes variant data before running `validate-test`."
        )
    return crispritz_targets["reference"], crispritz_targets["alternative (1000G)"]


def check_variant_dataset() -> None:
    """
    Verify that the variant dataset used for complete-test is compatible with
    validation requirements.

    This function reads the test VCF configuration, ensures that exactly one
    dataset is specified, and checks that it corresponds to the expected
    `hg38_1000G` dataset.

    Raises:
        FileNotFoundError: If the variant dataset configuration file cannot be
            found.
        ValueError: If the configuration specifies an invalid number of
            datasets or a dataset other than `hg38_1000G`.
    """
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
    """
    Download the brute-force off-target benchmark file required for validation.

    This function fetches the benchmark targets from the remote repository into
    the local benchmark directory, verifies that the file exists, and checks its
    integrity using an MD5 checksum before returning the local path.

    Returns:
        str: The local filesystem path to the downloaded brute-force targets
            file.

    Raises:
        FileNotFoundError: If the download fails or the expected file cannot be
            found on disk.
        ValueError: If the downloaded file does not match the expected MD5
            checksum, indicating possible corruption or incompleteness.
    """
    sys.stderr.write("Downloading brute-force targets\n")
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
    sys.stderr.write(
        f"Brute-force targets successfully downloaded: {bf_targets_path}\n"
    )
    return bf_targets_path


def _load_targets(fname: str, chrom: str, crisprme: bool = False) -> pd.DataFrame:
    """
    Load and filter target records from a CRISPRme or brute-force results file.

    This function reads a tab-separated targets file, restricts entries to
    supported chromosomes, and optionally filters to a specific chromosome while
    adapting to the appropriate chromosome column for CRISPRme or brute-force
    formats.

    Args:
        fname (str): Path to the targets file to load.
        chrom (str): Chromosome to retain (e.g., "chr1") or "all" to keep all
            supported chromosomes.
        crisprme (bool): Whether the input file is a CRISPRme output
            (uses the "Chromosome" column instead of "CHR").

    Returns:
        pd.DataFrame: A DataFrame containing only targets on allowed
            chromosomes, optionally restricted to a single chromosome.

    Raises:
        FileNotFoundError: If the targets file does not exist.
        ValueError: If the expected chromosome column is missing, or if an
            invalid chromosome value is provided.
    """
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"Targets file not found: {fname}")
    df = pd.read_csv(fname, sep="\t")
    chrom_col = "Chromosome" if crisprme else "CHR"
    if chrom_col not in df.columns:
        raise ValueError(
            f"Invalid targets file format: missing required column '{chrom_col}'."
        )
    df = df[df[chrom_col].isin(CHROMS)]  # Keep only supported chromosomes
    if chrom != "all":  # filter by specific chromosome if requested
        if chrom not in CHROMS:
            raise ValueError(
                f"Invalid chromosome '{chrom}'. "
                f"Allowed values are 'all' or one of: {', '.join(CHROMS)}"
            )
        df = df[df[chrom_col] == chrom]
    return df  # type: ignore


def load_targets(
    crisprme_targets_ref_fname: str,
    crisprme_targets_alt_fname: str,
    bf_targets_fname: str,
    chrom: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load CRISPRme and brute-force target reports for a given chromosome to
    prepare them for validation.

    This function reads reference and alternative CRISPRme targets together with
    brute-force targets, filters them to the requested chromosome, and returns
    harmonized DataFrames for downstream comparison.

    Args:
        crisprme_targets_ref_fname (str): Path to the CRISPRme reference
            targets report file.
        crisprme_targets_alt_fname (str): Path to the CRISPRme alternative
            (1000G) targets report file.
        bf_targets_fname (str): Path to the brute-force targets file.
        chrom (str): Chromosome to filter on (e.g., "chr1") or "all" to include
            all supported chromosomes.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing the concatenated
            CRISPRme targets DataFrame and the brute-force targets DataFrame.

    Raises:
        FileNotFoundError: If any of the specified target files cannot be found.
        ValueError: If an invalid chromosome or malformed targets file is
            encountered while loading data.
    """
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
    """
    Build a compact site identifier string from genomic coordinates and strand.

    This function concatenates chromosome, position, and strand information into
    a single underscore-separated token suitable for comparing off-target sites
    across datasets.

    Args:
        chrom (str): Chromosome name in UCSC format (e.g., "chr1").
        pos (int): Genomic 1-based position of the site.
        strand (str): Strand symbol, typically "+" or "-".

    Returns:
        str: An underscore-separated site identifier of the form
            "<chrom>_<pos>_<strand>".
    """
    return "_".join(list(map(str, [chrom, pos, strand])))


def compute_sites(targets: pd.DataFrame, crisprme: bool = False) -> List[str]:
    """
    Compute compact site identifiers for all target records in a results table.

    This function derives an underscore-separated site string for each row,
    using the appropriate chromosome, position, and strand columns for either
    CRISPRme or brute-force formats, and returns the collection of identifiers.

    Args:
        targets (pd.DataFrame): DataFrame containing target records from either
            CRISPRme or brute-force results.
        crisprme (bool): Whether the input table is from CRISPRme output, which
            uses CRISPRme-specific column ordering.

    Returns:
        List[str]: A list of site identifier strings, one per target row.
    """
    if crisprme:
        targets["site"] = targets.apply(
            lambda x: _compute_site(x[3], x[4], x[6]), axis=1
        )
    else:
        targets["site"] = targets.apply(
            lambda x: _compute_site(x[0], x[4], x[3]), axis=1
        )
    return targets["site"].tolist()


def validate(crisprme_targets: pd.DataFrame, bf_targets: pd.DataFrame) -> None:
    """
    Validate CRISPRme off-target predictions against a brute-force benchmark set.

    This function compares site identifiers derived from CRISPRme and brute-force
    targets, reports whether they match, and prints detailed diagnostics for any
    discrepancies to help debugging.

    Args:
        crisprme_targets (pd.DataFrame): DataFrame containing CRISPRme-derived
            off-target records.
        bf_targets (pd.DataFrame): DataFrame containing brute-force-derived
            off-target records.

    Raises:
        SystemExit: Always raised at the end of validation, with exit status 0
            regardless of whether the comparison passes or fails.
    """
    sys.stderr.write("Running off-target sites validation\n")
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
    """
    Run the complete validation workflow for CRISPRme off-target predictions on
    a specified chromosome.

    This function checks that required CRISPRitz outputs and variant datasets
    are available, downloads the brute-force benchmark, loads all targets for
    the requested chromosome, and invokes the validator to compare prediction
    sets.

    Args:
        chrom (str): Chromosome to validate (e.g., "chr1") or "all" to validate
            across all supported chromosomes.
    """
    # check crispritz targets existence and dataset used
    crisprme_targets_ref, crisprme_targets_alt = check_crispritz_targets()
    check_variant_dataset()
    # download brute-force targets on 1000G
    bf_targets_fname = download_brute_force_targets()
    # read and load crisprme and brute-force targets
    crisprme_targets, bf_targets = load_targets(
        crisprme_targets_ref, crisprme_targets_alt, bf_targets_fname, chrom
    )
    # validate offtarget sites
    validate(crisprme_targets, bf_targets)


def main():
    chrom = sys.argv[1]  # read command line args
    run_test_validation(chrom)  # validate results from complete-test


if __name__ == "__main__":
    main()
