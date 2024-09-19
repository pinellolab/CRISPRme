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
VCF1000GSERVER = "ftp.1000genomes.ebi.ac.uk"
VCF1000GURL = (
    "/vol1/ftp/data_collections/1000_genomes_project/release/"
    "20190312_biallelic_SNV_and_INDEL/"
    "ALL.{}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
)
VCFHGDPSERVER = "ngs.sanger.ac.uk"
VCFHGDPURL = "/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.{}.vcf.gz"
SAMPLESIDURL = "https://raw.githubusercontent.com/pinellolab/CRISPRme/gnomad-4.1-converter/download_data"
ANNOTATIONURL = "https://raw.githubusercontent.com/pinellolab/CRISPRme/gnomad-4.1-converter/download_data"


# TODO: before PR fix urls and script call


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
    Download genome data for a specified chromosome or all chromosomes.

    This function retrieves FASTA files for the specified chromosome in UCSC format
    and saves them to the designated destination directory. It handles both individual
    chromosome downloads and bulk downloads for all chromosomes, ensuring that the
    destination directory exists and is valid.

    Args:
        chrom (str): The chromosome to download, specified in UCSC format (e.g., "chr1").
            Use "all" to download data for all chromosomes.
        dest (str): The destination directory where the downloaded files will be saved.

    Raises:
        ValueError: If the specified chromosome is not valid.
        FileExistsError: If the destination directory does not exist.
        RuntimeError: If the genome data download fails.
    """

    # assume chromosomes given in UCSC format (chr1, chr2, etc.)
    if chrom not in CHROMS + ["all"]:
        raise ValueError(f"Forbidden input chromosome ({chrom})")
    if not os.path.isdir(dest):  # check dest directory existence
        raise FileExistsError(f"Unable to locate {dest}")
    sys.stderr.write(f"Downloading fasta file for chromosome(s) {chrom}\n")
    try:
        if chrom == "all":
            chromstar = download(
                dest, http_url=f"{HG38URL}/bigZips/hg38.chromFa.tar.gz", resume=True
            )
            chromsdir = untar(chromstar, dest, "chroms")  # decompress archive
            # rename chroms dir to hg38
            chromsdir = rename(
                chromsdir, os.path.join(os.path.dirname(chromsdir), "hg38")
            )
            assert os.path.isdir(chromsdir)
        else:
            chromgz = download(
                dest, http_url=f"{HG38URL}/chromosomes/{chrom}.fa.gz", resume=True
            )
            dest = ensure_hg38_directory(dest)  # create hg38 directory
            chromfa = gunzip(
                chromgz,
                os.path.join(dest, f"{os.path.splitext(os.path.basename(chromgz))[0]}"),
            )  # decompress chrom FASTA
            assert os.path.isfile(chromfa)
    except RuntimeError as e:
        raise RuntimeError("Genome data download failed!") from e


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
    Download VCF data for a specified chromosome or all chromosomes from a given dataset.

    This function retrieves Variant Call Format (VCF) files for the specified chromosome
    or all chromosomes in UCSC format from the selected dataset. It ensures that the
    destination directory exists and that the dataset is valid before proceeding with the
    download.

    Args:
        chrom (str): The chromosome to download, specified in UCSC format (e.g., "chr1").
            Use "all" to download data for all chromosomes.
        dest (str): The destination directory where the downloaded VCF files will be saved.
        dataset (str): The dataset to download from, which can be "1000G", "HGDP", or
            "1000G+HGDP".

    Raises:
        ValueError: If the specified chromosome or dataset is not valid.
        FileExistsError: If the destination directory does not exist.
        RuntimeError: If the VCF dataset download fails.
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
    try:
        for ds in dataset.split("+"):
            vcf_dataset_dir = ensure_vcf_dataset_directory(dest, ds)
            ftp_server = VCF1000GSERVER if ds == "1000G" else VCFHGDPSERVER
            vcf_url = VCF1000GURL if ds == "1000G" else VCFHGDPURL
            chroms = CHROMS if chrom == "all" else [chrom]
            for c in chroms:  # request FTP connection
                download(
                    vcf_dataset_dir,
                    ftp_conn=True,
                    ftp_server=ftp_server,
                    ftp_path=vcf_url.format(c),
                )
    except RuntimeError as e:
        raise RuntimeError(f"VCF dataset {ds} download failed!") from e


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
    Download sample IDs for the specified variant dataset.

    This function retrieves sample ID files for the given dataset, which can include
    "1000G", "HGDP", or both. It ensures that the dataset is valid before proceeding
    with the download and saves the files in the appropriate directory.

    Args:
        dataset (str): The dataset for which to download sample IDs. Valid options are
            "1000G", "HGDP", or "1000G+HGDP".

    Raises:
        ValueError: If the specified dataset is not valid.
        RuntimeError: If the download of sample IDs fails.
    """

    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # samples ids folder must be located within current directory
    # -- see check_crisprme_directory_tree() for details
    sys.stderr.write(f"Downloading sample ids for dataset(s) {dataset}\n")
    try:
        samplesids_dir = ensure_samplesids_directory(os.getcwd())
        for ds in dataset.split("+"):
            samplesid_fname = (
                "hg38_1000G.samplesID.txt"
                if ds == "1000G"
                else "hg38_HGDP.samplesID.txt"
            )
            download(samplesids_dir, http_url=f"{SAMPLESIDURL}/{samplesid_fname}")
    except RuntimeError as e:
        raise RuntimeError(f"Download samples ids for dataset {ds} failed!") from e


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
    Download ENCODE and GENCODE annotation data files.

    This function retrieves the annotation data for protein-coding genes from the
    GENCODE and ENCODE projects and saves them in the appropriate directory. It
    ensures that the necessary directories exist before downloading the data.

    Returns:
        Tuple[str, str]: A tuple containing the file paths of the downloaded
            GENCODE and ENCODE annotation data files.

    Raises:
        RuntimeError: If the download of annotation data fails.
    """

    sys.stderr.write("Downloading ENCODE and GENCODE annotation data\n")
    try:
        annotation_dir = ensure_annotation_directory(os.getcwd())
        # download gencode annotation
        gencodetar = download(
            annotation_dir,
            http_url=f"{ANNOTATIONURL}/gencode.protein_coding.bed.tar.gz",
        )
        gencode = os.path.join(
            untar(gencodetar, annotation_dir), "gencode.protein_coding.bed"
        )
        # download encode annotation
        encodetar = download(
            annotation_dir,
            http_url=f"{ANNOTATIONURL}/dhs+encode+gencode.hg38.bed.tar.gz",
        )
        encode = os.path.join(
            untar(encodetar, annotation_dir), "dhs+encode+gencode.hg38.bed"
        )
    except RuntimeError as e:
        raise RuntimeError("Annotation data download failed!") from e
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
    Create a PAM (Protospacer Adjacent Motif) file for SpCas9.

    This function generates a PAM file containing the sequence for the NGG motif,
    which is essential for the SpCas9 gene-editing system. It ensures that the
    necessary directory exists before writing the PAM sequence to a file.

    Returns:
        str: The path to the created PAM file.

    Raises:
        IOError: If an error occurs while writing the PAM file.
        RuntimeError: If the PAM file generation fails.
    """

    sys.stderr.write("Creating PAM file\n")
    try:
        pams_dir = ensure_pams_directory(
            os.getcwd()
        )  # PAMs directory must be in current working dir
        pamfile = os.path.join(pams_dir, "20bp-NGG-SpCas9.txt")
        try:
            with open(pamfile, mode="w") as outfile:
                outfile.write("NNNNNNNNNNNNNNNNNNNNNGG 3\n")  # 20 + 3 bp (NGG)
        except IOError as e:
            raise IOError("An error occurred while writing the test PAM file") from e
    except RuntimeError as e:
        raise RuntimeError("PAM file generation failed!") from e
    return pamfile


def write_sg1617_guidefile() -> str:
    """
    Create a guide file for the sg1617 sequence.

    This function generates a guide file containing the sg1617 sequence, which is
    used in gene editing applications. It ensures that the file is created successfully
    and handles any potential errors during the writing process.

    Returns:
        str: The path to the created guide file.

    Raises:
        IOError: If an error occurs while writing the guide file.
        RuntimeError: If the guide file generation fails.
    """

    sys.stderr.write("Creating guide file\n")
    try:
        guidefile = "sg1617_test_guide.txt"
        try:
            with open(guidefile, mode="w") as outfile:
                outfile.write("CTAACAGTTGCTTTTATCACNNN\n")  # sg1617 guide
        except IOError as e:
            raise IOError("An error occerred while writing the test guide file") from e
    except RuntimeError as e:
        raise RuntimeError("Guide file generation failed!") from e
    return guidefile


def write_vcf_list(dataset: str) -> str:
    """
    Create a VCF configuration file for the specified variant dataset.

    This function generates a VCF list file based on the provided dataset, which can
    include "1000G", "HGDP", or both. It ensures that the dataset is valid before
    writing the corresponding entries to the file.

    Args:
        dataset (str): The variant dataset for which to create the VCF list.
            Valid options are "1000G", "HGDP", or "1000G+HGDP".

    Returns:
        str: The path to the created VCF list file.

    Raises:
        ValueError: If the specified dataset is not valid.
        IOError: If an error occurs while writing the VCF list file.
        RuntimeError: If the VCF configuration file generation fails.
    """

    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # select the test vcf list
    sys.stderr.write(f"Creating VCF config file for dataset(s) {dataset}\n")
    try:
        vcflistfile = "vcf_list_test.txt"
        try:
            with open(vcflistfile, mode="w") as outfile:
                for ds in dataset.split("+"):
                    outfile.write(f"hg38_{ds}\n")
        except IOError as e:
            raise IOError("An error occurred while writing the test VCF list") from e
    except RuntimeError as e:
        raise RuntimeError("VCF configuration file generation failed!") from e
    return vcflistfile


def write_samplesids_list(dataset: str) -> str:
    """
    Create a samples ID configuration file for the specified variant dataset.

    This function generates a list of sample ID files based on the provided dataset,
    which can include "1000G", "HGDP", or both. It ensures that the dataset is valid
    before writing the corresponding sample ID entries to the file.

    Args:
        dataset (str): The variant dataset for which to create the samples ID list.
            Valid options are "1000G", "HGDP", or "1000G+HGDP".

    Returns:
        str: The path to the created samples ID list file.

    Raises:
        ValueError: If the specified dataset is not valid.
        IOError: If an error occurs while writing the samples ID list file.
        RuntimeError: If the samples ID configuration file generation fails.
    """

    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # select the test vcf list
    sys.stderr.write(f"Creating samples config file for dataset(s) {dataset}\n")
    try:
        samplesidslistfile = "samplesID_list_test.txt"
        try:
            with open(samplesidslistfile, mode="w") as outfile:
                for ds in dataset.split("+"):
                    samplesidslist = (
                        "hg38_1000G.samplesID.txt"
                        if ds == "1000G"
                        else "hg38_HGDP.samplesID.txt"
                    )
                    outfile.write(f"{samplesidslist}\n")
        except OSError as e:
            raise IOError("An error occurred while writing the test VCF list") from e
    except RuntimeError as e:
        raise RuntimeError("Samples ids configuration file generation failed!") from e
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
    # download gencode and encode annotation data
    gencode, encode = download_annotation_data()
    pam = write_ngg_pamfile()  # write test NGG PAM file
    guide = write_sg1617_guidefile()  # write test sg1617 guide
    vcf = write_vcf_list(dataset)  # write test vcf list
    samplesids = write_samplesids_list(dataset)  # write test samples ids list
    debug_arg = "--debug" if debug else ""
    crisprme_cmd = (
        f"python ./crisprme.py complete-search --genome {CRISPRME_DIRS[0]}/hg38 "
        f"--thread 4 --bmax 1 --mm 4 --bDNA 1 --bRNA 1 --merge 3 --pam {pam} "
        f"--guide {guide} --vcf {vcf} --samplesID {samplesids} --annotation {encode} "
        f"--gene_annotation {gencode} --output crisprme-test-out {debug_arg}"
    )
    code = subprocess.call(crisprme_cmd, shell=True)  # run crisprme test
    if code != 0:
        raise subprocess.SubprocessError(
            "Test search failed! See Results/crisprme-test-out/log_error.txt for details\n"
        )


def main():
    chrom, dataset, debug = sys.argv[1:]  # read commandline args
    debug = debug == "True"  # debug mode
    run_crisprme_test(chrom, dataset, debug)  # run crisprme test


if __name__ == "__main__":
    main()
