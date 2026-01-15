"""
This module provides functionality to execute the CRISPRme test workflow, including 
downloading genomic and variant data, preparing input files, and running the CRISPRme 
command-line tool. It includes functions for managing directories, downloading data 
from various sources, and configuring test parameters.

Key functions include:
- `ensure_hg38_directory`: Ensures the existence of the 'hg38' directory within the 
    specified destination.
- `download_genome_data`: Downloads genome data for a specified chromosome to the 
    destination directory.
- `ensure_vcf_dataset_directory`: Ensures the existence of a directory for a specific 
    VCF dataset.
- `download_vcf_data`: Downloads VCF data for a specific chromosome and variant dataset.
- `ensure_samplesids_directory`: Ensures the existence of the 'samplesIDs' directory.
- `download_samples_ids_data`: Downloads samples IDs data for a specific variant dataset.
- `ensure_annotation_directory`: Ensures the existence of the 'annotation' directory.
- `download_annotation_data`: Downloads gencode and encode annotation data to the 
    'annotation' directory.
- `write_ngg_pamfile`: Writes a test PAM file containing the NGG sequence.
- `write_sg1617_guidefile`: Writes a test guide file containing the sg1617 guide 
    sequence.
- `write_vcf_config`: Writes a test VCF list file for a specific variant dataset.
- `write_samplesids_config`: Writes a test samples ID list file for a specific 
    variant dataset.
- `run_crisprme_test`: Executes the CRISPRme test workflow for a specified chromosome 
    and dataset.
- `main`: The entry point of the module that orchestrates the test execution.

This module is designed to facilitate the testing and validation of the CRISPRme 
tool, ensuring that all necessary data and configurations are correctly handled 
before running the analysis.
"""

from utils import (
    check_crisprme_directory_tree,
    download,
    gunzip,
    untar,
    rename,
    compute_md5,
    CHROMS,
    CRISPRME_DIRS,
    MD5GENOME,
    MD51000G,
    MD5HGDP,
    MD5SAMPLES,
    MD5ANNOTATION,
)

from typing import Tuple

import subprocess
import sys
import os

HG38URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38"
VCF1000GSERVER = "ftp.1000genomes.ebi.ac.uk"
VCF1000GURL = (
    "/vol1/ftp/data_collections/1000_genomes_project/release/"
    "20190312_biallelic_SNV_and_INDEL/"
    "ALL.{}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
)
VCFHGDPSERVER = "ngs.sanger.ac.uk"
VCFHGDPURL = "/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.{}.vcf.gz"
TESTDATAURL = "https://raw.githubusercontent.com/pinellolab/CRISPRme/refs/heads/main/test/data/"


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
    sys.stderr.write(f"Downloading fasta file for chromosome(s) {chrom}\n")
    if chrom == "all":
        chromstar = download(dest, http_url=f"{HG38URL}/bigZips/hg38.chromFa.tar.gz")
        # check genome md5
        if MD5GENOME[os.path.basename(chromstar)] != compute_md5(chromstar):
            raise ValueError(f"Download for {os.path.basename(chromstar)} failed")
        chromsdir = untar(chromstar, dest, "chroms")  # decompress archive
        # rename chroms dir to hg38
        chromsdir = rename(chromsdir, os.path.join(os.path.dirname(chromsdir), "hg38"))
        assert os.path.isdir(chromsdir)
    else:
        chromgz = download(dest, http_url=f"{HG38URL}/chromosomes/{chrom}.fa.gz")
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
            "HGDP", or 1000G+HGDP).

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
    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # create VCF dataset directory within VCFs folder
    sys.stderr.write(f"Downloading VCF data for chromsome(s) {chrom}\n")
    for ds in dataset.split("+"):
        vcf_dataset_dir = ensure_vcf_dataset_directory(dest, ds)
        ftp_server = VCF1000GSERVER if ds == "1000G" else VCFHGDPSERVER
        vcf_url = VCF1000GURL if ds == "1000G" else VCFHGDPURL
        chroms = CHROMS if chrom == "all" else [chrom]
        for c in chroms:  # request FTP connection
            vcf = download(
                vcf_dataset_dir,
                ftp_conn=True,
                ftp_server=ftp_server,
                ftp_path=vcf_url.format(c),
            )
            md5data = MD51000G if ds == "1000G" else MD5HGDP
            if md5data[os.path.basename(vcf)] != compute_md5(vcf):
                raise ValueError(f"Download for {os.path.basename(vcf)} failed")


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

    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # samples ids folder must be located within current directory
    # -- see check_crisprme_directory_tree() for details
    sys.stderr.write(f"Downloading sample ids for dataset(s) {dataset}\n")
    samplesids_dir = ensure_samplesids_directory(os.getcwd())
    for ds in dataset.split("+"):
        samplesid_fname = (
            "samplesIDs.1000G.txt" if ds == "1000G" else "samplesIDs.HGDP.txt"
        )
        samplesids = download(
            samplesids_dir, http_url=f"{TESTDATAURL}/samplesIDs/{samplesid_fname}"
        )
        if MD5SAMPLES[os.path.basename(samplesids)] != compute_md5(samplesids):
            raise ValueError(f"Download for {os.path.basename(samplesids)} failed")


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

    sys.stderr.write("Downloading ENCODE and GENCODE annotation data\n")
    annotation_dir = ensure_annotation_directory(os.getcwd())
    gencode = _retrieve_ann_data(
        annotation_dir,
        "Annotations/gencode.protein_coding.bed.tar.gz",
        "gencode.protein_coding.bed",
    )
    encode = _retrieve_ann_data(
        annotation_dir,
        "Annotations/dhs+encode+gencode.hg38.bed.tar.gz",
        "dhs+encode+gencode.hg38.bed",
    )
    return gencode, encode

def _bgzip_ann_data(ann_fname: str) -> str:
    """
    Compress an annotation file using bgzip and verify the output.

    This function calls the bgzip utility to compress the specified annotation file,
    checks that the compressed file exists, and returns its path. It raises an error
    if the compression fails.

    Args:
        ann_fname (str): The path to the annotation file to be compressed.

    Returns:
        str: The path to the compressed annotation file.

    Raises:
        subprocess.SubprocessError: If bgzip compression fails.
    """
    
    try: 
        subprocess.call(f"bgzip -f {ann_fname}", shell=True)
    except (subprocess.SubprocessError, Exception) as e:
        raise subprocess.SubprocessError(f"Bgzip compression failed on {ann_fname}") from e
    ann_fname_gz = f"{ann_fname}.gz"
    assert os.path.isfile(ann_fname_gz)  # check that the bgzipped bed exists
    return ann_fname_gz
    


def _retrieve_ann_data(annotation_dir: str, url: str, fname: str) -> str:
    """
    Download and extract an annotation file, then compress it with bgzip.

    This function downloads an annotation archive, verifies its integrity, extracts 
    the specified file, and compresses it using bgzip. It returns the path to the 
    compressed annotation file.

    Args:
        annotation_dir (str): The directory to store the annotation data.
        url (str): The URL of the annotation archive to download.
        fname (str): The name of the file to extract and compress.

    Returns:
        str: The path to the compressed annotation file.

    Raises:
        ValueError: If the downloaded file fails the MD5 check.
    """

    # download gencode annotation
    annfile_tar = download(annotation_dir, http_url=os.path.join(TESTDATAURL, url))
    if MD5ANNOTATION[os.path.basename(annfile_tar)] != compute_md5(annfile_tar):
        raise ValueError(f"Download for {os.path.basename(annfile_tar)} failed")
    return _bgzip_ann_data(os.path.join(untar(annfile_tar, annotation_dir), fname))


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

    sys.stderr.write("Creating PAM file\n")
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

    sys.stderr.write("Creating guide file\n")
    guidefile = "sg1617_test_guide.txt"
    try:
        with open(guidefile, mode="w") as outfile:
            outfile.write("CTAACAGTTGCTTTTATCACNNN\n")  # sg1617 guide
    except IOError as e:
        raise IOError("An error occerred while writing the test guide file") from e
    return guidefile


def write_vcf_config(dataset: str) -> str:
    """
    Write a test VCF list file for a specific variant dataset.

    Args:
        dataset (str): The name or identifier of the variant dataset (e.g., "1000G",
            "HGDP").

    Returns:
        str: The path to the created test VCF list file.
    """

    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # config vcf list file
    sys.stderr.write(f"Creating VCF config file for dataset(s) {dataset}\n")
    vcf_config = "vcf.config.test.txt"
    try:
        with open(vcf_config, mode="w") as outfile:
            for ds in dataset.split("+"):
                outfile.write(f"hg38_{ds}\n")
    except IOError as e:
        raise IOError("An error occurred while writing the test VCF list") from e
    return vcf_config


def write_samplesids_config(dataset: str) -> str:
    """
    Write a test samples ID list file for a specific variant dataset.

    Args:
        dataset (str): The name or identifier of the variant dataset (e.g., "1000G",
            "HGDP").

    Returns:
        str: The path to the created test samples ID list file.
    """

    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # configure sample ids list
    sys.stderr.write(f"Creating samples config file for dataset(s) {dataset}\n")
    samples_config = "samplesIDs.config.test.txt"
    try:
        with open(samples_config, mode="w") as outfile:
            for ds in dataset.split("+"):
                samplesidslist = (
                    "samplesIDs.1000G.txt" if ds == "1000G" else "samplesIDs.HGDP.txt"
                )
                outfile.write(f"{samplesidslist}\n")
    except IOError as e:
        raise IOError("An error occurred while writing the test VCF list") from e
    return samples_config


def run_crisprme_test(chrom: str, dataset: str, threads: int, debug: bool) -> None:
    """Execute the CRISPRme test workflow for a specified chromosome and dataset.

    This function orchestrates the downloading of necessary genomic and VCF data,
    prepares input files, and runs the CRISPRme command-line tool to perform a
    complete search.

    Args:
        chrom (str): The chromosome to be analyzed.
        dataset (str): The dataset identifier for VCF data.
        threads (int): The number of threads to use for processing.
        debug (bool): A flag indicating whether to run in debug mode.

    Raises:
        Any exceptions raised by the called functions or subprocess.
    """

    check_crisprme_directory_tree(os.getcwd())  # check crisprme directory tree
    download_genome_data(chrom, CRISPRME_DIRS[0])  # download genome data
    download_vcf_data(chrom, CRISPRME_DIRS[3], dataset)  # download vcf data
    vcf = write_vcf_config(dataset)  # write test vcf list
    download_samples_ids_data(dataset)  # download vcf dataset samples ids
    samplesids = write_samplesids_config(dataset)  # write test samples ids list
    # download gencode and encode annotation data
    gencode, encode = download_annotation_data()
    pam = write_ngg_pamfile()  # write test NGG PAM file
    guide = write_sg1617_guidefile()  # write test sg1617 guide
    debug_arg = "--debug" if debug else ""
    # TODO: replace call to local crisprme
    crisprme_cmd = (
        f"crisprme.py complete-search --genome {CRISPRME_DIRS[0]}/hg38 "
        f"--thread 4 --bmax 1 --mm 4 --bDNA 1 --bRNA 1 --merge 3 --pam {pam} "
        f"--guide {guide} --vcf {vcf} --samplesID {samplesids} --annotation {encode} "
        f"--gene_annotation {gencode} --output crisprme-test-out --thread {threads} "
        f"{debug_arg} --ci-cd-test"
    )
    subprocess.call(crisprme_cmd, shell=True)  # run crisprme test


def main():
    chrom, dataset, threads, debug = sys.argv[1:]  # read commandline args
    debug = debug == "True"
    run_crisprme_test(chrom, dataset, int(threads), debug)  # run crisprme test


if __name__ == "__main__":
    main()
