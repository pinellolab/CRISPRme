#!/usr/bin/env python
"""Processes a report file to replace 'n' and '.' values with 'NA'.

This module reads a tab-separated report file in chunks, cleans specific values,
and overwrites the original file with the cleaned data. It is intended for use
in post-processing steps where standardized missing value representation is required.
"""

from pandas.io.parsers import TextFileReader

import pandas as pd

import subprocess
import warnings
import sys
import os

warnings.simplefilter(action="ignore", category=FutureWarning)

# dataframe read chunksize
CHUNKSIZE = 50000


def read_report_chunks(report_fname: str) -> TextFileReader:
    """Reads a tab-separated report file in chunks for efficient processing.

    This function returns an iterator that yields DataFrame chunks from the
    specified file.

    Args:
        report_fname: Path to the tab-separated report file.

    Returns:
        TextFileReader: An iterator over DataFrame chunks.
    """
    return pd.read_csv(report_fname, sep="\t", chunksize=CHUNKSIZE)


def polish_report(report_chunks: TextFileReader, report_fname: str) -> None:
    """Cleans and writes report data by replacing specific values with 'NA'.

    This function processes each chunk of the report, replacing 'n' and '.' values
    with 'NA', and writes the cleaned data to a temporary file.

    Args:
        report_chunks: An iterator over DataFrame chunks from the report.
        report_fname: Path to the original report file.

    Returns:
        None
    """
    header = True  # only for first iteration
    for chunk in report_chunks:
        assert "rsID" in chunk.columns.tolist()
        chunk: pd.DataFrame = chunk.replace("n", "NA")  # replace ns with NAs
        chunk["rsID"] = chunk["rsID"].str.replace(
            ".", "NA"
        )  # replace . in rsids with NAs
        chunk.to_csv(
            f"{report_fname}.tmp", header=header, mode="a", sep="\t", index=False
        )
        header = False  # not required anymore


def remove_n_dots() -> None:
    """Replaces 'n' and '.' values in a report file with 'NA' and overwrites the
    original file.

    This function reads a report file, cleans specific values, and saves the cleaned
    data back to the original file.

    Returns:
        None
    """
    report_fname = sys.argv[1]  # read input report filename
    if not os.path.isfile(report_fname) or os.stat(report_fname).st_size <= 0:
        raise FileNotFoundError(f"Report {report_fname} not found or empty")
    # remove ns and dots and replace them with NAs
    polish_report(read_report_chunks(report_fname), report_fname)
    if not os.path.isfile(f"{report_fname}.tmp"):
        raise FileNotFoundError(f"Polished report {report_fname}.tmp not created?")
    code = subprocess.call(f"mv {report_fname}.tmp {report_fname}")
    if code != 0:
        raise subprocess.SubprocessError(f"Failed renaming {report_fname}.tmp")


remove_n_dots()
