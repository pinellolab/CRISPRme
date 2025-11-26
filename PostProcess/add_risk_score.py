#!/usr/bin/env python
"""Adds risk score columns to a report file and computes risk scores for each entry.

This module processes a tab-separated report file, calculates risk score differences,
and writes the results to a new file with updated columns. It supports both primary
and alternative target scenarios based on user input.
"""

from typing import List, Tuple

import sys

# risk score column names
RISKCOLNAMES = ["Highest_CFD_Risk_Score", "Highest_CFD_Absolute_Risk_Score"]


def add_risk_score_columns(header: List[str], alternative: bool) -> List[str]:
    """Appends risk score columns and optionally a cluster ID to the header list.

    This function updates the provided header list by adding predefined risk score
    columns, and adds a cluster ID column if the alternative flag is set.

    Args:
        header: The list of column names to be updated.
        alternative: Whether to add the cluster ID column.

    Returns:
        List[str]: The updated list of column names including risk score columns.
    """
    header.extend(iter(RISKCOLNAMES))
    if alternative:
        header.append("CLUSTER_ID")
    return header


def _compute_risk_score(line: str) -> Tuple[float, float]:
    """Calculates the risk score difference and its absolute value from a line of
    report data.

    This function extracts two score values from the input line, computes their
    difference, and returns both the difference and its absolute value.

    Args:
        line: A tab-separated string containing report data.

    Returns:
        Tuple[float, float]: The risk score difference and its absolute value.
    """
    fields = line.strip().split("\t")
    score, score_alt = float(fields[20]), float(fields[21])
    score_diff = score - score_alt
    return score_diff, abs(score_diff)


def compute_risk_score(report_fname: str, report_outfname: str, alternative: bool):
    """Reads a report file, computes risk scores for each line, and writes the
    results to a new file.

    This function processes the input report, appends risk score columns, and
    outputs the updated data. It handles both primary and alternative target
    scenarios based on the provided flag.

    Args:
        report_fname: Path to the input report file.
        report_outfname: Path to the output report file.
        alternative: Whether to process as alternative targets.

    Returns:
        None

    Raises:
        OSError: If an error occurs while reading or writing files.
    """
    try:
        with open(report_fname, mode="r") as fin, open(
            report_outfname, mode="w"
        ) as fout:
            header = fin.readline().strip().split("\t")  # read file header
            header = add_risk_score_columns(
                header, alternative
            )  # append risk score columns
            fout.write("\t".join(header) + "\n")  # write header to out put report
            for line in fin:
                score_diff, score_diff_abs = _compute_risk_score(line)
                fout.write(f"{line.strip()}\t{score_diff}\t{score_diff_abs}\n")
    except (IOError, Exception) as e:
        raise OSError(
            f"An error occurred while computing risk scores for {report_fname}"
        ) from e


def risk_score() -> None:
    """Parses command-line arguments and computes risk scores for a report file.

    This function reads input and output file paths and a flag from the command
    line, then processes the report to add risk score columns accordingly.

    Returns:
        None
    """
    report_fname, report_outfname, alternative = sys.argv[1:4]  # read input args
    alternative = alternative == "True"  # are primary or alternative targets?
    # compute risk score on CFD/CRISTA
    compute_risk_score(report_fname, report_outfname, alternative)


risk_score()
