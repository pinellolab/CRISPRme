"""
"""

from .utils import (
    check_crisprme_directory_tree,
    download,
    gunzip,
    untar,
    rename,
    CHROMS,
    CRISPRME_DIRS,
)

import sys
import os

HG38URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38"


def ensure_hg38_directory(dest: str) -> str:
    """
    Ensure the existence of the 'hg38' directory within the specified destination
    directory.

    Args:
        dest (str): The destination directory where the 'hg38' directory should
            be created.

    Returns:
        str: The path to the 'hg38' directory.
    """

    hg38_dir = os.path.join(dest, "hg38")
    if not os.path.exists(hg38_dir):  # create hg38 directory within Genomes
        os.mkdir(hg38_dir)
    return hg38_dir


def download_genome_data(chrom: str, dest: str) -> None:
    """
    Download genome data for a specified chromosome to the destination directory.

    Args:
        chrom (str): The chromosome identifier in UCSC format.
        dest (str): The destination directory to save the downloaded genome data.

    Returns:
        None

    Raises:
        ValueError: If the input chromosome is not valid.
        FileExistsError: If the destination directory does not exist.
    """

    # assume chromosomes given in UCSC format (chr1, chr2, etc.)
    if chrom not in CHROMS + ["all"]:
        raise ValueError(f"Forbidden input chromosome ({chrom})")
    if not os.path.isdir(dest):  # check dest directory existence
        raise FileExistsError(f"Unable to locate {dest}")
    if chrom == "all":
        chromstar = download(f"{HG38URL}/bigZips/hg38.chromFa.tar.gz", dest)
        chromsdir = untar(chromstar, dest, "chroms")  # decompress archive
        # rename chroms dir to hg38
        chromsdir = rename(chromsdir, os.path.join(os.path.dirname(chromsdir), "hg38"))
        assert os.path.isdir(chromsdir)
    else:
        chromgz = download(f"{HG38URL}/chromosomes/{chrom}.fa.gz", dest)
        dest = ensure_hg38_directory(dest)  # create hg38 directory
        chromfa = gunzip(
            chromgz,
            os.path.join(dest, f"{os.path.splitext(os.path.basename(chromgz))[0]}"),
        )  # decompress chrom FASTA
        assert os.path.isfile(chromfa)


def run_crisprme_test(chrom: str) -> None:
    # check crisprme directory tree
    check_crisprme_directory_tree(os.getcwd())
    # download genome data
    download_genome_data(chrom, CRISPRME_DIRS[0])


if __name__ == "__main__":
    run_crisprme_test("all")
