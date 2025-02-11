"""Provides functions for displaying and managing CRISPRme job history.

This module implements the history page of the CRISPRme web application,
allowing users to view and interact with past job results. It includes
functions for retrieving job parameters, displaying a history table, and
handling user interactions such as row selection and filtering.
"""

from .pages_utils import RESULTS_DIR, PARAMS_FILE, LOG_FILE, GUIDES_FILE

from typing import Dict, List, Optional, Tuple
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from app import app, current_working_directory, URL

import dash_html_components as html
import pandas as pd
import numpy as np

import math
import os

SUMMARYTABCOLS = [
    "Job", "Genome", "Variants", "Mismatches", "DNA bulge", "RNA bulge", "PAM", "Number of Guides", "Start"
]


def support_filter_history(
    result_df: pd.DataFrame, genome_filter: str, pam_filter: str
) -> Tuple[pd.DataFrame, int]:
    """Filter the results history dataframe based on genome and PAM criteria.

    This function filters the input dataframe to include only jobs matching the
    specified genome build and PAM sequence(s). It also calculates the maximum
    number of pages for displaying the filtered results.

    Args:
        result_df: The dataframe containing the job history.
        genome_filter: The genome build to filter by.
        pam_filter: The PAM sequence(s) to filter by.

    Returns:
        A tuple containing the filtered dataframe and the maximum number of pages.
    """
    if not isinstance(result_df, pd.DataFrame):
        raise TypeError(
            f"Expected {type(pd.DataFrame).__name__}, got {type(result_df).__name__}"
        )
    if not isinstance(genome_filter, str):
        raise TypeError(f"Expected {str.__name__}, got {type(genome_filter).__name__}")
    if not isinstance(pam_filter, str):
        raise TypeError(f"Exepcted {str.__name__}, got {type(pam_filter).__name__}")
    if genome_filter is not None:
        # keep only rows related to the requested genome type
        drop_rows = result_df[result_df["Genome"] != genome_filter].index
        result_df.drop(drop_rows, axis=0, inplace=True)
    if pam_filter is not None:
        drop_rows = result_df[~result_df["PAM"].isin(pam_filter)].index
        result_df.drop(drop_rows, axis=0, inplace=True)
    # allow also empty results history?
    max_page = result_df.shape[0]
    max_page = math.floor(max_page / 1000000) + 1  # max size
    return result_df, max_page


# trigger row highlighting
@app.callback(
    Output("results-table", "style_data_conditional"),
    [Input("results-table", "selected_cells")],
    [State("results-table", "data")],
)
def highlight_row(sel_cel: List, all_guides: List) -> List[Dict[str, str]]:
    """Highlight the selected row in the results table.

    This callback function dynamically styles the results table to highlight the
    row corresponding to the selected cell.

    Args:
        sel_cel: A list representing the selected cell(s).
        all_guides: The data of the results table.

    Returns:
        A list of dictionaries specifying the styling for the highlighted row.

    Raises:
        PreventUpdate: If no cell is selected or the table data is empty.
    """
    if sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate
    job_name = all_guides[int(sel_cel[0]["row"])]["Job"]
    assert isinstance(job_name, str)
    return [
        {
            "if": {"filter_query": "".join(["{Job} eq ", f'"{job_name}"'])},
            "background-color": "rgba(0, 0, 255,0.15)",  # highlighting color
        }
    ]


def retrieve_resultsdirs(resultsdir: str) -> List[str]:
    """Retrieve a list of valid result directories.

    This function scans the specified results directory and returns a list of
    subdirectories that represent valid CRISPRme job results. A directory is
    considered valid if it contains a PARAMS_FILE.

    Args:
        resultsdir: The path to the results directory.

    Returns:
        A list of strings, where each string is the path to a valid result
        directory.
    """
    return [
        d
        for d in os.listdir(resultsdir)
        if os.path.isdir(os.path.join(resultsdir, d)) and os.path.isfile(os.path.join(resultsdir, d, PARAMS_FILE))
    ]

def read_params(paramsfile: str, jobid: str) -> Dict[str, str]:
    """Read job parameters from a file.

    This function reads the parameters used for a specific CRISPRme job from
    the given parameters file.

    Args:
        paramsfile: The path to the parameters file.
        jobid: The ID of the job.

    Returns:
        A dictionary containing the job parameters, where keys are parameter
        names and values are parameter values.

    Raises:
        IOError: If an error occurs while reading the parameters file.
    """
    try:
        with open(paramsfile, mode="r") as infile:
            params = {
                fields[0]: fields[1] 
                for line in infile 
                for fields in [line.strip().split()]
            }  # read search parameters
    except OSError as e:
        raise IOError(f"An error occurred while collecting results in {jobid}") from e
    return params

def read_job_info(logfile: str, jobid: str) -> str:
    """Read job start time from the log file.

    This function extracts the job start time from the specified log file.

    Args:
        logfile: The path to the log file.
        jobid: The ID of the job.

    Returns:
        A string representing the job start time.

    Raises:
        IOError: If an error occurs while reading the log file.
    """
    try:
        with open(logfile, mode="r") as infile:
            log = infile.read()  # read job log data
    except OSError as e:
        raise IOError(f"An error occurred while collecting results in {jobid}") from e
    return (next(s for s in log.split("\n") if "Job\tStart" in s)).split("\t")[-1]

def count_guides(guidesfile: str, jobid: str) -> int:
    """Count the number of guides in a guides file.

    This function reads the specified guides file and returns the number of
    guides found within it.

    Args:
        guidesfile: The path to the guides file.
        jobid: The ID of the job.

    Returns:
        The number of guides found in the file.

    Raises:
        IOError: If an error occurs while reading the guides file.
    """
    try:
        with open(guidesfile, mode="r") as infile:
            return len(infile.read().strip().split("\n"))
    except OSError as e:
        raise IOError(f"An error occurred while collecting results in {jobid}") from e

def process_genome(genome_selected: str, genome_idx: str) -> Tuple[str, str]:
    """Process genome information from job parameters.

    This function extracts the genome build and variant information from the
    given genome selection and index parameters.

    Args:
        genome_selected: The selected genome string.
        genome_idx: The genome index string.

    Returns:
        A tuple containing the genome build and a comma-separated string of
        variants.
    """
    genome = genome_selected.split("+")[0] if "+" in genome_selected else genome_selected
    variants = "Reference"
    if "+" in genome_idx:
        variants = ",".join([e.split("+")[-1] for e in genome_idx.split(",")])
    return genome, variants


def construct_history_summary(results: List[str]) -> pd.DataFrame:
    """Construct a summary dataframe of CRISPRme job history.

    This function reads job information from parameter and log files for each
    completed job ID and compiles a summary dataframe. The dataframe includes
    job parameters, start time, and number of guides.

    Args:
        results: A list of job IDs.

    Returns:
        A pandas DataFrame summarizing the job history.
    """
    summary = {c: [] for c in SUMMARYTABCOLS}  # initialize table
    for jobid in results:
        params = read_params(os.path.join(current_working_directory, RESULTS_DIR, jobid, PARAMS_FILE), jobid)
        jobinfo = read_job_info(os.path.join(current_working_directory, RESULTS_DIR, jobid, LOG_FILE), jobid)
        guidesnum = count_guides(os.path.join(current_working_directory, RESULTS_DIR, jobid, GUIDES_FILE), jobid)
        genome, variants = process_genome(params["Genome_selected"], params["Genome_idx"])
        summary[SUMMARYTABCOLS[0]].append(jobid)  # job id
        summary[SUMMARYTABCOLS[1]].append(genome)  # genome 
        summary[SUMMARYTABCOLS[2]].append(variants)  # variants
        summary[SUMMARYTABCOLS[3]].append(int(params["Mismatches"]))  # mms
        summary[SUMMARYTABCOLS[4]].append(int(params["DNA"]))  # dna bulges
        summary[SUMMARYTABCOLS[5]].append(int(params["RNA"]))  # rna bulges
        summary[SUMMARYTABCOLS[6]].append(params["Pam"])  # pam
        summary[SUMMARYTABCOLS[7]].append(guidesnum)  # number of guides
        summary[SUMMARYTABCOLS[8]].append(jobinfo)  # job start time
        print(jobinfo)
    summary = pd.DataFrame(summary)
    summary[SUMMARYTABCOLS[8]] = pd.to_datetime(summary[SUMMARYTABCOLS[8]])
    summary = summary.sort_values([SUMMARYTABCOLS[8]], ascending=False)
    return summary


def get_available_results() -> pd.DataFrame:
    """Retrieve available CRISPRme job results.

    This function scans the results directory for valid job result directories
    and constructs a summary dataframe of the available results.

    Returns:
        A pandas DataFrame summarizing the available job results.
    """
    # retrieve results stored in Results folder
    results_directory = os.path.join(current_working_directory, RESULTS_DIR)
    results_dirs = retrieve_resultsdirs(results_directory) 
    assert len(results_dirs) <= len(os.listdir(results_directory))  # empty results are skipped
    return construct_history_summary(results_dirs)


def table_header_(results: pd.DataFrame) -> html.Thead:
    return html.Thead(
        html.Tr(
            [
                html.Th(c, style={"vertical-align": "middle", "text-align": "center"})
                if c not in ["Load", "Delete"]
                else html.Th("", style={"vertical-align": "middle", "text-align": "center"})
                for c in results.columns.tolist()
            ]
        )
    )

def history_table_body_(results: pd.DataFrame, page: int, max_rows: int, remaining_rows: int) -> List[html.Tr]:
    """Generate the body of the history table.

    This function creates the rows for the history table, displaying job
    information for each result. It handles pagination by displaying only a
    subset of rows based on the current page and maximum rows per page.

    Args:
        results: The DataFrame containing job history data.
        page: The current page number.
        max_rows: The maximum number of rows to display per page.
        remaining_rows: The number of rows remaining to be displayed.

    Returns:
        A list of html.Tr elements representing the table rows.
    """
    return [
        html.Tr([
            html.Td(
                html.A(
                    str(results.at[i + (page - 1) * max_rows, "Job"]),
                    target="_blank",
                    href=os.path.join(URL, f"load?job={results.at[i + (page - 1) * max_rows, 'Job']}")
                ),
                style={"vertical-align": "middle", "text-align": "center"},
            ) if col == "Job" else 
            html.Td(
                results.at[i + (page - 1) * max_rows, col],
                style={"vertical-align": "middle", "text-align": "center"},
            )
            for col in results.columns
        ])
        for i in range(min(remaining_rows, max_rows))
    ]


def display_history_table(
    results: pd.DataFrame, page: int, max_rows: Optional[int] = 1000000
) -> html.Div:
    """Display the history table with pagination.

    This function creates the HTML representation of the history table, including
    pagination controls. It generates the table header and body, and handles
    displaying only a subset of rows based on the current page and maximum rows
    per page.

    Args:
        results: The DataFrame containing job history data.
        page: The current page number.
        max_rows: The maximum number of rows to display per page.

    Returns:
        An html.Div element containing the history table.
    """
    remaining_rows = results.shape[0] - (page - 1) * max_rows
    hist_tab_body = history_table_body_(results, page, max_rows, remaining_rows)
    return [
        html.Table(
            [table_header_(results), html.Tbody(hist_tab_body)],
            style={"display": "inline-block"},
        ),
    ] + [
        html.Button(
            f"{i}",
            id=f"button-delete-history-{i}",
            **{"data-jobid": "None"},
            style={"display": "none"},
        )
        for i in range(min(remaining_rows, max_rows) - 1, 10)
    ]


def retrieve_mode() -> str:
    """Retrieve the CRISPRme running mode.

    This function reads the running mode (server or local) from the mode type
    file.

    Returns:
        A string representing the running mode ("server" or "local").
    """
    with open(os.path.join(current_working_directory, ".mode_type.txt"), mode="r") as infile:
        modetype = infile.readlines()
    return modetype[0]

def history_header() -> html.Div:
    """Create the header for the history page.

    This function generates the HTML header section for the history page,
    including the title and a brief description.

    Returns:
        An html.Div element containing the header content.
    """
    return html.Div(
        [
            html.H3("Results History"),
            html.P(
                "List of available results. Click on a link to open the " 
                "corresponding results page in a new tab."
            )
        ]
    )

def history_table(results: pd.DataFrame, mode: str) -> html.Div:
    """Create the history table or display a message if not available.

    This function generates the HTML for the history table if the CRISPRme
    running mode allows it (local mode). If running in server mode, it displays
    a message indicating that history is not available.

    Args:
        results: The DataFrame containing job history data.
        mode: The CRISPRme running mode ("server" or "local").

    Returns:
        An html.Div element containing the history table or a message.
    """
    if mode != "server":
        return html.Div(
            display_history_table(results, 1), id="div-history-table", style={"text-align": "center"}
        )
    return html.Div("History is not available while using website mode")


def history_page() -> html.Div:
    """Create the layout for the history page.

    This function generates the HTML structure for the history page, displaying
    a table of previous CRISPRme job results if available. If running in server
    mode, the history table is not displayed.

    Returns:
        An html.Div element containing the history page layout.
    """
    results = get_available_results()  # retrieve available results 
    if retrieve_mode() == "server":  # running online -> do not disply history
        results = pd.DataFrame()
    html_divs = [
        history_header(), 
        html.Div("None,None", id="div-history-filter-query", style={"display": "none"}),
        history_table(results, retrieve_mode()),
        html.Div(id="div-remove-jobid", style={"display": "none"}),
        html.Div(html.Br(), style={"text-align": "center"})
    ]
    max_page = results.shape[0]
    max_page = np.floor(max_page / 1000000) + 1
    return html.Div(html_divs, style={"margin": "1%"})
