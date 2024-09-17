"""
"""

from io import TextIOWrapper
from typing import Optional, Union
from ftplib import FTP

import subprocess
import requests
import tarfile
import gzip
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


def ftp_download(
    ftp_server: Union[str, None],
    ftp_path: Union[str, None],
    dest: str,
    fname: Optional[str] = None,
) -> str:
    """
    Download a file from an FTP server and save it to a specified destination.

    This function connects to an FTP server to retrieve a file located at a specified
    path and saves it to a designated directory. It ensures that the necessary
    parameters are provided and that the file is successfully created after the
    download.

    Args:
        ftp_server (Union[str, None]): The address of the FTP server. Must be provided.
        ftp_path (Union[str, None]): The path of the file on the FTP server.
            Must be provided.
        dest (str): The destination directory where the file will be saved.
        fname (Optional[str]): The name of the file to save as. If not provided,
            the base name of the FTP path will be used.

    Returns:
        str: The path to the downloaded file.

    Raises:
        ValueError: If the FTP server or path is not provided.
        TypeError: If the FTP server or path is not a string.
        FileNotFoundError: If the file is not created after the download attempt.
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
    """
    Download a file from a specified HTTP or HTTPS URL and save it to a destination.

    This function retrieves a file from the provided HTTP URL and saves it to a specified
    directory. It ensures that the URL is valid and that the file is successfully created
    after the download.

    Args:
        http_url (Union[str, None]): The URL of the file to download. Must be provided.
        dest (str): The destination directory where the file will be saved.
        fname (Optional[str]): The name of the file to save as. If not provided,
            the base name of the URL will be used.

    Returns:
        str: The path to the downloaded file.

    Raises:
        ValueError: If the HTTP URL is not provided or is invalid.
        TypeError: If the HTTP URL is not a string.
        FileNotFoundError: If the file is not created after the download attempt.
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
    """
    Download a file from either an FTP server or an HTTP/HTTPS URL and save it to
    a destination.

    This function determines the appropriate method to download a file based on the
    provided parameters, either using FTP or HTTP. It validates the input parameters
    and ensures that the file is saved correctly in the specified directory.

    Args:
        dest (str): The destination directory where the file will be saved.
        fname (Optional[str]): The name of the file to save as. If not provided,
            the base name of the URL or FTP path will be used.
        ftp_conn (Optional[bool]): A flag indicating whether to use FTP for the download.
        ftp_server (Optional[str]): The FTP server address, required if ftp_conn is True.
        ftp_path (Optional[str]): The path on the FTP server to retrieve the file from,
            required if ftp_conn is True.
        http_url (Optional[str]): The URL of the file to download via HTTP/HTTPS.

    Returns:
        str: The path to the downloaded file.

    Raises:
        TypeError: If the destination or file name is not a string.
        ValueError: If the destination directory does not exist, or if both FTP and HTTP
            connections are requested simultaneously.
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
