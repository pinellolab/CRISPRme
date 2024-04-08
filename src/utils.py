"""
"""

from typing import Optional

import subprocess
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
CHROMS = [f"chr{i}" for i in list(range(1, 23)) + ["X"]]


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


def download(url: str, dest: str) -> str:
    """
    Download a file from the specified URL to the destination directory.

    Args:
        url (str): The URL of the file to download.
        dest (str): The destination directory to save the downloaded file.

    Returns:
        str: The file path of the downloaded file.

    Raises:
        TypeError: If url is not a string.
        subprocess.SubprocessError: If the download process fails.

    Examples:
        download('http://example.com/file.txt', '/path/to/destination/directory')
    """

    if not isinstance(url, str):
        raise TypeError(f"Expected {str.__name__}, got {type(url).__name__}")
    fname = os.path.join(dest, os.path.basename(url))
    code = subprocess.call(f"wget {url} -O {fname}", shell=True)  # download using wget
    if code != 0:
        raise subprocess.SubprocessError(f"Download from {url} failed")
    assert os.path.isfile(fname)
    return fname


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
            with tarfile.open(fileobj=fin, mode="r") as tar:
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
