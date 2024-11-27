"""
This module provides utility functions for file handling, downloading data from 
FTP and HTTP sources, and managing directory structures for the CRISPRme project. 
It includes functions for checking directory structures, downloading files, 
decompressing archives, and computing file hashes.

Key functions include:
- `check_crisprme_directory_tree`: Checks and creates the necessary CRISPRme 
    directory structure.
- `ftp_download`: Downloads a file from an FTP server to a specified local destination.
- `http_download`: Downloads a file from an HTTP or HTTPS URL to a specified 
    local destination.
- `download`: Downloads a file from either an FTP server or an HTTP/HTTPS URL 
    based on provided parameters.
- `remove`: Removes a file or directory specified by the given path.
- `rename`: Renames a file or directory from the original path to the new 
    specified name.
- `untar`: Decompresses and extracts the contents of a tar.gz file to a specified 
    destination.
- `gunzip`: Decompresses a gzip file to the specified output file.
- `compute_md5`: Computes the MD5 hash of a file for integrity verification.

This module is designed to facilitate data management and processing for genomic
analysis workflows, ensuring that necessary files and directories are correctly 
handled and maintained.
"""

from io import TextIOWrapper
from typing import Optional, Union
from ftplib import FTP

import subprocess
import requests
import hashlib
import tarfile
import gzip
import sys
import os

CRISPRME_DIRS = [
    "Genomes",
    "Results",
    "Dictionaries",
    "VCFs",
    "Annotations",
    "PAMs",
    "samplesIDs",
]
CHROMS = [f"chr{i}" for i in list(range(1, 23)) + ["X"]]  # canonical chroms
# md5 hashes stored in dictionary, md5s are used to check test files consistency
MD5GENOME = {"hg38.chromFa.tar.gz": "a5aa5da14ccf3d259c4308f7b2c18cb0"}
MD51000G = {
    "ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "77f154e53c2b7c36b04d03bab3af8b74",
    "ALL.chr2.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "f9d29c4935e591b2b269eed7cd7e35d8",
    "ALL.chr3.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "6e59d00235de71562b4199e09b7e5934",
    "ALL.chr4.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "70a2c1ede97eceb7baeea06c8e46cf3c",
    "ALL.chr5.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "74d5486c0fd29b0e6add24d3740fc3b4",
    "ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "8c5d83c1a9253058120368af39baf0c8",
    "ALL.chr7.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "dfaa282712fc1292146173dd2ffeb1d9",
    "ALL.chr8.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "ddf7b370fcee63462037c237f12b4444",
    "ALL.chr9.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "5ade69521dc50d88ad7c91bf4ec6fcd8",
    "ALL.chr10.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "1c409a674426eda2fd29b49078137c5d",
    "ALL.chr11.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "65339bffc61bc97f2130832fe9f84d7c",
    "ALL.chr12.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "9a1bda389121140d30c768ef6a1b1370",
    "ALL.chr13.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "47b0463541be137a8bbfe40f6aade864",
    "ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "241aedf0792c45d5345d421105c782af",
    "ALL.chr15.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "b48e7c64e35b727d34786faa76467f94",
    "ALL.chr16.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "1ce7d66799cab6718852d78dd2aab765",
    "ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "ecc22783fd1ee7a1c66b053491873192",
    "ALL.chr18.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "fdf3e460e91cd955a9e8cebf01b5d815",
    "ALL.chr19.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "a2f17e4ec552fc07cbd05c1eac0cf7ec",
    "ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "155c3b440d7990630132e4756f7fcc85",
    "ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "52882490028507e5d4e606b0905072b1",
    "ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "57a1722e6ed7d9df08cb3c0e42b62d53",
    "ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz": "e6a3d41811faee60de177061edcd6fe6",
}
MD5HGDP = {
    "hgdp_wgs.20190516.full.chr1.vcf.gz": "70d82ae3ae65cb73858738f547f64e93",
    "hgdp_wgs.20190516.full.chr2.vcf.gz": "539d4eb31355b90f0262453fa1349ae6",
    "hgdp_wgs.20190516.full.chr3.vcf.gz": "0d37ba60afd5ff092cf1bc75bde3588e",
    "hgdp_wgs.20190516.full.chr4.vcf.gz": "68d57e5c2129bbafa1a9dd75f630cf89",
    "hgdp_wgs.20190516.full.chr5.vcf.gz": "929eb66a26e9679320bcc26df0bd4116",
    "hgdp_wgs.20190516.full.chr6.vcf.gz": "28c9ad734025e7292bde533da908cf68",
    "hgdp_wgs.20190516.full.chr7.vcf.gz": "ed7eaf339cd7964b9f1e7581de5bdeb1",
    "hgdp_wgs.20190516.full.chr8.vcf.gz": "2d47d60ff6b63e1163d219d964999ee3",
    "hgdp_wgs.20190516.full.chr9.vcf.gz": "52f917fc3068eff76f0ba8bde0c59292",
    "hgdp_wgs.20190516.full.chr10.vcf.gz": "7131981641e886173da90d215346e857",
    "hgdp_wgs.20190516.full.chr11.vcf.gz": "2f127b0006cbc36fb32c66860d4b31d9",
    "hgdp_wgs.20190516.full.chr12.vcf.gz": "023d0e2d852c167490d4578f814d043d",
    "hgdp_wgs.20190516.full.chr13.vcf.gz": "afcfba8b01258e418f5fb230b14daa02",
    "hgdp_wgs.20190516.full.chr14.vcf.gz": "90b0c15b61fd9c47a9751495f2b784ce",
    "hgdp_wgs.20190516.full.chr15.vcf.gz": "665e844d7e2e85e226d25827ea8014be",
    "hgdp_wgs.20190516.full.chr16.vcf.gz": "0d6f1b6141c78489a2b2e27eeec848dd",
    "hgdp_wgs.20190516.full.chr17.vcf.gz": "d53421438b3bc3c5ce5ab51b90578182",
    "hgdp_wgs.20190516.full.chr18.vcf.gz": "6351d9b20995cf500ac4b11490ff31c7",
    "hgdp_wgs.20190516.full.chr19.vcf.gz": "167ce7a43876b32e586978a75f3b0d39",
    "hgdp_wgs.20190516.full.chr20.vcf.gz": "d90130b11620378bed7c2cc43be94b7e",
    "hgdp_wgs.20190516.full.chr21.vcf.gz": "8f44e4daa3952cd73751141f66b6e5ae",
    "hgdp_wgs.20190516.full.chr22.vcf.gz": "84f4a1d86f54bdc0cd9b19502ff8d2c2",
    "hgdp_wgs.20190516.full.chrX.vcf.gz": "8d0e4e178fdfa07db76d0218a9b2ceab",
    "hgdp_wgs.20190516.full.chrY.vcf.gz": "54b3aba28600c8d0d8a695c8dcfdc4cd",
}
MD5ANNOTATION = {
    "dhs+encode+gencode.hg38.bed.tar.gz": "4f5eb631af903d4091bb2f57558c7b46",
    "gencode.protein_coding.bed.tar.gz": "04297ade436db70784733a5b13d42723",
}
MD5SAMPLES = {
    "samplesIDs.1000G.txt": "720af666c9a938de74a2808033aa4509",
    "samplesIDs.HGDP.txt": "f92e14e5317221486f20597560ca3a31",
}


def check_crisprme_directory_tree(basedir: str) -> None:
    """
    Check if the CRISPRme directory tree exists in the specified base directory
    and create it if not.

    Args:
        basedir (str): The base directory path to check and create the CRISPRme
            directory tree.

    Returns:
        None

    Raises:
        TypeError: If basedir is not a string.
        FileExistsError: If basedir does not exist.

    Examples:
        check_crisprme_directory_tree('/path/to/base/directory')
    """

    if not isinstance(basedir, str):
        raise TypeError(f"Expected {str.__name__}, got {type(basedir).__name__}")
    if not os.path.exists(basedir):
        raise FileExistsError(f"Unable to locate {basedir}")
    for d in CRISPRME_DIRS:  # create CRISPRme directory tree
        # if directory not found, create it
        if not os.path.isdir(os.path.join(basedir, d)):
            os.makedirs(os.path.join(basedir, d))


def ftp_download(
    ftp_server: Union[str, None],
    ftp_path: Union[str, None],
    dest: str,
    fname: Optional[str] = None,
) -> str:
    """Download a file from an FTP server to a specified destination.

    This function connects to the given FTP server and retrieves a file from the
    specified path, saving it to the local destination. It handles errors related
    to FTP connection and file writing.

    Args:
        ftp_server (Union[str, None]): The address of the FTP server.
        ftp_path (Union[str, None]): The path of the file on the FTP server.
        dest (str): The local destination directory where the file will be saved.
        fname (Optional[str]): The name to save the file as locally. If not provided,
            the base name of the FTP path will be used.

    Returns:
        str: The path to the downloaded file.

    Raises:
        ValueError: If the FTP server or path is not provided.
        TypeError: If the FTP server or path is not a string.
        OSError: If an error occurs while saving the file.
        FileNotFoundError: If the downloaded file is not created.
    """

    if ftp_server is None or ftp_path is None:
        raise ValueError(
            "FTP server and path must be provided if FTP connection is requested"
        )
    if not isinstance(ftp_server, str):
        raise TypeError(f"Expected {str.__name__}, got {type(ftp_server).__name__}")
    if not isinstance(ftp_path, str):
        raise TypeError(f"Expected {str.__name__}, got {type(ftp_path).__name__}")
    fname = os.path.join(dest, fname or os.path.basename(ftp_path))
    try:
        with FTP(ftp_server) as ftp:  # initialize ftp server
            ftp.login()  # open connection to server
            with open(fname, mode="wb") as outfile:  # write binary data
                ftp.retrbinary(f"RETR {ftp_path}", outfile.write)
    except IOError as e:
        raise OSError(f"An error occurred while saving {fname}") from e
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"{fname} not created")
    return fname


def http_download(
    http_url: Union[str, None], dest: str, fname: Optional[str] = None
) -> str:
    """Download a file from an HTTP or HTTPS URL to a specified destination.

    This function retrieves a file from the provided HTTP URL and saves it to the
    specified local destination. It ensures that the URL is valid and handles
    errors related to the HTTP request and file writing.

    Args:
        http_url (Union[str, None]): The URL of the file to download.
        dest (str): The local destination directory where the file will be saved.
        fname (Optional[str]): The name to save the file as locally. If not provided,
            the base name of the URL will be used.

    Returns:
        str: The path to the downloaded file.

    Raises:
        ValueError: If the HTTP URL is not provided or is invalid.
        TypeError: If the HTTP URL is not a string.
        FileNotFoundError: If the downloaded file is not created.
    """

    if http_url is None:
        raise ValueError("HTTP URL must be provided if HTTP connection is requested")
    if not isinstance(http_url, str):
        raise TypeError(f"Expected {str.__name__}, got {type(http_url).__name__}")
    if not (http_url.startswith("http://") or http_url.startswith("https://")):
        raise ValueError(
            "Invalid HTTP URL. It must start with 'http://' or 'https://'."
        )
    fname = os.path.join(dest, fname or os.path.basename(http_url))
    response = requests.get(http_url, stream=True)  # download data from http
    response.raise_for_status()  # ensure the request was successful
    with open(fname, mode="wb") as outfile:
        # write downloaded data in fixed size chunks
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                outfile.write(chunk)
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"{fname} not created")
    return fname


def download(
    dest: str,
    fname: Optional[str] = None,
    ftp_conn: Optional[bool] = False,
    ftp_server: Optional[str] = None,
    ftp_path: Optional[str] = None,
    http_url: Optional[str] = None,
) -> str:
    """Download a file from either an FTP server or an HTTP/HTTPS URL.

    This function determines the appropriate method to download a file based on the
    provided parameters, either using FTP or HTTP. It validates the input parameters
    and ensures that the destination directory exists before proceeding with the
    download.

    Args:
        dest (str): The local destination directory where the file will be saved.
        fname (Optional[str]): The name to save the file as locally. If not provided,
            the base name of the source will be used.
        ftp_conn (Optional[bool]): Flag indicating whether to use FTP for the download.
        ftp_server (Optional[str]): The FTP server address if using FTP.
        ftp_path (Optional[str]): The path of the file on the FTP server.
        http_url (Optional[str]): The URL of the file to download via HTTP/HTTPS.

    Returns:
        str: The path to the downloaded file.

    Raises:
        TypeError: If dest or fname is not a string.
        ValueError: If the destination directory does not exist, or if both FTP and
            HTTP connections are requested simultaneously.
    """

    if not isinstance(dest, str):
        raise TypeError(f"Expected {str.__name__}, got {type(dest).__name__}")
    if not os.path.isdir(dest):
        raise ValueError(
            f"Destination directory {dest} does not exist or is not a directory"
        )
    if fname is not None and not isinstance(fname, str):
        raise TypeError(f"Expected {str.__name__}, got {type(fname).__name__}")
    if ftp_conn and http_url:
        raise ValueError(
            "Both ftp and http connection cannot be requested at the same time"
        )
    if ftp_conn:  # ftp connection requested
        return ftp_download(ftp_server, ftp_path, dest, fname)
    else:  # http/https connection requested
        return http_download(http_url, dest, fname)


def remove(fname: str) -> None:
    """
    Remove a file or directory specified by the given path.

    Args:
        fname (str): The path of the file or directory to be removed.

    Returns:
        None

    Raises:
        subprocess.SubprocessError: If an error occurs during the removal process.
    """

    code = subprocess.call(f"rm -rf {fname}", shell=True)
    if code != 0:
        raise subprocess.SubprocessError(f"An error occurred while removing {fname}")


def rename(orig: str, newname: str) -> str:
    """
    Rename a file or directory from the original path to the new specified name.

    Args:
        orig (str): The original path of the file or directory to be renamed.
        newname (str): The new name or path to rename the file or directory to.

    Returns:
        str: The new path or name after renaming.

    Raises:
        subprocess.SubprocessError: If the renaming process fails.
    """

    code = subprocess.call(f"mv {orig} {newname}", shell=True)
    if code != 0:
        subprocess.SubprocessError(f"Renaming {orig} failed")
    assert os.path.exists(newname)
    return newname


def untar(fname_tar_gz: str, dest: str, outdir: Optional[str] = "") -> str:
    """
    Decompress and extract the contents of a tar.gz file to the specified
    destination directory.

    Args:
        fname_tar_gz (str): The path to the tar.gz file to decompress and extract.
        dest (str): The destination directory to extract the contents to.
        outdir (Optional[str]): Optional subdirectory within the destination to
            extract the contents to.

    Returns:
        str: The path to the directory where the contents were extracted.

    Raises:
        IOError: If an error occurs during the decompression process.
    """

    try:
        with gzip.open(fname_tar_gz, mode="rb") as fin:
            with tarfile.open(fileobj=fin, mode="r") as tar:  # type: ignore
                tar.extractall(dest)
    except IOError as e:
        raise IOError(f"An error occurred while decompressing {fname_tar_gz}") from e
    outdir = os.path.join(dest, outdir) if outdir else dest
    assert os.path.isdir(outdir)
    remove(fname_tar_gz)  # delete compressed archive
    return outdir


def gunzip(fname_gz: str, fname_out: str) -> str:
    """
    Decompress a gzip file to the specified output file.

    Args:
        fname_gz (str): The path to the gzip file to decompress.
        fname_out (str): The path to save the decompressed output.

    Returns:
        str: The path to the decompressed output file.

    Raises:
        IOError: If an error occurs during the decompression process.
    """

    try:
        with gzip.open(fname_gz, mode="rb") as fin:
            with open(fname_out, mode="wb") as fout:
                fout.write(fin.read())
    except IOError as e:
        raise IOError(f"An error occurred while decompressing {fname_gz}") from e
    assert os.stat(fname_out).st_size > 0
    remove(fname_gz)  # delete compressed archive
    return fname_out


def compute_md5(fname: str) -> str:
    """Compute the MD5 hash of a file.

    This function reads a file in binary mode and calculates its MD5 hash, returning
    the hash as a hexadecimal string. It handles file access errors gracefully.

    Args:
        fname (str): The path to the file for which to compute the MD5 hash.

    Returns:
        str: The computed MD5 hash of the file as a hexadecimal string.

    Raises:
        FileNotFoundError: If the specified file cannot be found.
        PermissionError: If there is a permission issue accessing the file.
        Exception: For any unexpected errors encountered during the process.
    """

    hashmd5 = hashlib.md5()  # initialize md5
    try:  # open file and compute its md5
        with open(fname, mode="rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):  # read 4096 bytes per chunk
                hashmd5.update(chunk)
        return hashmd5.hexdigest()  # return md5 as string
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Unable to locate {fname}") from e
    except PermissionError as e:
        raise PermissionError(f"Permission denied when accessing {fname}") from e
    except Exception as e:
        # sourcery skip: raise-specific-error
        raise Exception(
            f"Unexpected error encountered while computing md5 on {fname}"
        ) from e
