"""
"""

from utils import (
    check_crisprme_directory_tree,
    download,
    gunzip,
    untar,
    rename,
    CHROMS,
    CRISPRME_DIRS,
)

from typing import Tuple

import subprocess
import sys
import os

HG38URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38"
VCF1000GPURL = (
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/"
    "release/20190312_biallelic_SNV_and_INDEL/"
    "ALL.{}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
)
VCFHGDPURL = (
    "ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/"
    "hgdp_wgs.20190516.full.{}.vcf.gz"
)
SAMPLESIDURL = (
    "https://raw.githubusercontent.com/pinellolab/CRISPRme/test-function/download_data"
)
ANNOTATIONURL = "https://github.com/pinellolab/CRISPRme/raw/test-function/download_data"


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


def ensure_vcf_dataset_directory(dest: str, dataset: str) -> str:
    """
    Ensure the existence of a directory for a specific VCF dataset within the
    specified destination directory.

    Args:
        dest (str): The destination directory where the VCF dataset directory
            should be created.
        dataset (str): The name or identifier of the VCF dataset.

    Returns:
        str: The path to the VCF dataset directory.
    """

    vcf_dataset_dir = os.path.join(dest, f"hg38_{dataset}")
    if not os.path.exists(vcf_dataset_dir):
        os.mkdir(vcf_dataset_dir)
    return vcf_dataset_dir


def download_vcf_data(chrom: str, dest: str, dataset: str) -> None:
    """
    Download VCF data for a specific chromosome and variant dataset to the destination
    directory.

    Args:
        chrom (str): The chromosome identifier in UCSC format.
        dest (str): The destination directory to save the downloaded VCF data.
        dataset (str): The name or identifier of the variant dataset (e.g., "1000G",
            "HGDP").

    Returns:
        None

    Raises:
        ValueError: If the input chromosome or variant dataset is invalid.
        FileExistsError: If the destination directory does not exist.
    """

    # assume chromosomes given in UCSC format (chr1, chr2, etc.)
    if chrom not in CHROMS + ["all"]:
        raise ValueError(f"Forbidden input chromosome ({chrom})")
    if not os.path.isdir(dest):  # check dest directory existence
        raise FileExistsError(f"Unable to locate {dest}")
    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # create VCF dataset directory within VCFs folder
    vcf_dataset_dir = ensure_vcf_dataset_directory(dest, dataset)
    vcf_url = VCF1000GPURL if dataset == "1000G" else VCFHGDPURL
    chroms = CHROMS if chrom == "all" else [chrom]
    for c in chroms:
        download(vcf_url.format(c), vcf_dataset_dir)


def ensure_samplesids_directory(dest: str) -> str:
    """
    Ensure the existence of the 'samplesIDs' directory within the specified
    destination directory.

    Args:
        dest (str): The destination directory where the 'samplesIDs' directory
            should be created.

    Returns:
        str: The path to the 'samplesIDs' directory.
    """

    samplesids_dir = os.path.join(dest, CRISPRME_DIRS[6])
    if not os.path.exists(samplesids_dir):
        os.mkdir(samplesids_dir)
    return samplesids_dir


def download_samples_ids_data(dataset: str) -> None:
    """
    Download samples IDs data for a specific variant dataset.

    Args:
        dataset (str): The name or identifier of the variant dataset (e.g.,
            "1000G", "HGDP").

    Returns:
        None

    Raises:
        ValueError: If the variant dataset is unknown.
    """

    if dataset not in ["1000G", "HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # samples ids folder must be located within current directory
    # -- see check_crisprme_directory_tree() for details
    samplesids_dir = ensure_samplesids_directory(os.getcwd())
    samplesid_fname = (
        "hg38_1000G.samplesID.txt" if dataset == "1000G" else "hg38_HGDP.samplesID.txt"
    )
    download(f"{SAMPLESIDURL}/{samplesid_fname}", samplesids_dir)


def ensure_annotation_directory(dest: str) -> str:
    """
    Ensure the existence of the 'annotation' directory within the specified
    destination directory.

    Args:
        dest (str): The destination directory where the 'annotation' directory
            should be created.

    Returns:
        str: The path to the 'annotation' directory.
    """

    annotation_dir = os.path.join(dest, CRISPRME_DIRS[4])
    if not os.path.exists(annotation_dir):
        os.mkdir(annotation_dir)
    return annotation_dir


def download_annotation_data() -> Tuple[str, str]:
    """
    Download gencode and encode annotation data to the 'annotation' directory
    within the current working directory.

    Returns:
        Tuple[str, str]: Paths to the downloaded gencode and encode annotation
            files.
    """

    annotation_dir = ensure_annotation_directory(os.getcwd())
    # download gencode annotation
    gencodetar = download(
        f"{ANNOTATIONURL}/gencode.protein_coding.tar.gz", annotation_dir
    )
    gencode = os.path.join(
        untar(gencodetar, annotation_dir), "gencode.protein_coding.bed"
    )
    # download encode annotation
    encodetar = download(f"{ANNOTATIONURL}/encode+gencode.hg38.tar.gz", annotation_dir)
    encode = os.path.join(untar(encodetar, annotation_dir), "encode+gencode.hg38.bed")
    return gencode, encode


def ensure_pams_directory(dest: str) -> str:
    """
    Ensure the existence of the 'PAMs' directory within the specified destination
    directory.

    Args:
        dest (str): The destination directory where the 'PAMs' directory should
            be created.

    Returns:
        str: The path to the 'PAMs' directory.
    """

    pams_dir = os.path.join(dest, CRISPRME_DIRS[5])
    if not os.path.exists(pams_dir):
        os.mkdir(pams_dir)
    return pams_dir


def write_ngg_pamfile() -> str:
    """
    Write a test PAM file containing the NGG sequence to the 'PAMs' directory
    within the current working directory.

    Returns:
        str: The path to the created test PAM file.
    """

    pams_dir = ensure_pams_directory(
        os.getcwd()
    )  # PAMs directory must be in current working dir
    pamfile = os.path.join(pams_dir, "20bp-NGG-SpCas9.txt")
    try:
        with open(pamfile, mode="w") as outfile:
            outfile.write("NNNNNNNNNNNNNNNNNNNNNGG 3\n")  # 20 + 3 bp (NGG)
    except IOError as e:
        raise IOError("An error occurred while writing the test PAM file") from e
    return pamfile


def write_sg1617_guidefile() -> str:
    """
    Write a test guide file containing the sg1617 guide sequence.

    Returns:
        str: The path to the created test guide file.
    """

    guidefile = "sg1617_test_guide.txt"
    try:
        with open(guidefile, mode="w") as outfile:
            outfile.write("CTAACAGTTGCTTTTATCACNNN\n")  # sg1617 guide
    except IOError as e:
        raise IOError("An error occerred while writing the test guide file") from e
    return guidefile


def write_vcf_list(dataset: str) -> str:
    """
    Write a test VCF list file for a specific variant dataset.

    Args:
        dataset (str): The name or identifier of the variant dataset (e.g., "1000G",
            "HGDP").

    Returns:
        str: The path to the created test VCF list file.
    """

    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # select the test vcf list
    vcflist = f"hg38_{dataset}"
    vcflistfile = "vcf_list_test.txt"
    try:
        with open(vcflistfile, mode="w") as outfile:
            outfile.write(f"{vcflist}\n")
    except IOError as e:
        raise IOError("An error occurred while writing the test VCF list") from e
    return vcflistfile


def write_samplesids_list(dataset: str) -> str:
    """
    Write a test samples ID list file for a specific variant dataset.

    Args:
        dataset (str): The name or identifier of the variant dataset (e.g., "1000G",
            "HGDP").

    Returns:
        str: The path to the created test samples ID list file.
    """

    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # select the test vcf list
    samplesidslist = (
        "hg38_1000G.samplesID.txt" if dataset == "1000G" else "hg38_HGDP.samplesID.txt"
    )
    samplesidslistfile = "samplesID_list_test.txt"
    try:
        with open(samplesidslistfile, mode="w") as outfile:
            outfile.write(f"{samplesidslist}\n")
    except IOError as e:
        raise IOError("An error occurred while writing the test VCF list") from e
    return samplesidslistfile


def run_crisprme_test(chrom: str, dataset: str, debug: bool) -> None:
    """
    Run CRISPRme test on specified chromosome and dataset.

    Args:
        chrom (str): The chromosome to run the test on.
        dataset (str): The dataset to use for the test.
        debug (bool): Flag to enable debug mode.

    Returns:
        None

    Raises:
        None
    """

    check_crisprme_directory_tree(os.getcwd())  # check crisprme directory tree
    download_genome_data(chrom, CRISPRME_DIRS[0])  # download genome data
    download_vcf_data(chrom, CRISPRME_DIRS[3], dataset)  # download vcf data
    download_samples_ids_data(dataset)  # download vcf dataset samples ids
    gencode, encode = (
        download_annotation_data()
    )  # download gencode and encode annotation data
    pam = write_ngg_pamfile()  # write test NGG PAM file
    guide = write_sg1617_guidefile()  # write test sg1617 guide
    vcf = write_vcf_list(dataset)  # write test vcf list
    samplesids = write_samplesids_list(dataset)  # write test samples ids list
    debug_arg = "--debug" if debug else ""
    crisprme_cmd = (
        f"crisprme.py complete-search --genome {CRISPRME_DIRS[0]}/hg38 "
        f"--thread 4 --bmax 1 --mm 4 --bDNA 1 --bRNA 1 --merge 3 --pam {pam} "
        f"--guide {guide} --vcf {vcf} --samplesID {samplesids} --annotation {encode} "
        f"--gene_annotation {gencode} --output crisprme-test-out {debug_arg}"
    )
    subprocess.call(crisprme_cmd, shell=True)  # run crisprme test


def main():
    chrom, dataset, debug = sys.argv[1:]  # read commandline args
    debug = debug == "True"
    run_crisprme_test(chrom, dataset, debug)  # run crisprme test


if __name__ == "__main__":
    main()
