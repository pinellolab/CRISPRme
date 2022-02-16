"""Script to manage and print the results page. The results page allows the user 
to navigate through the CRISPRme analysis results, providing several different
filtering options to visualize the results.

The result page consists of 2 (?) layers the ... on the top of the page and 
tab selection layer to select the way to visualize the results. 
Options (tabs) to visualize the results:
    - Custom Ranking
        Display the best targets for each input guide, according to the scoring
        criterion selected by the user.
    - Summary by Mismatch/Bulges
    - Summary by Sample
        Display the best targets for each individual sample considered during 
        CRISPRme analysis.
    - Query Genomic Region
        Display the best targets for each input guide in a specific genomic
        region. The results are sorted according to the scoring criterion 
        selected by the user.
    - Graphical Reports
        Plot information regarding the target guides (plots computed at execution
        time -> could require few secs to complete).
    - Personal Risk Card (displayed only if used individual data)

The results could be sorted and filtered according to 3 criteria:
    - CFD score
    - CRISTA score
    - Number of mismacthes and bulges

TODO: complete doc string with missing info --> read paper carefully
"""


from .results_page_utils import (
    PAGE_SIZE,
    BARPLOT_LEN,
    COL_REF,
    COL_REF_TYPE,
    COL_REF_RENAME,
    COL_BOTH,
    COL_BOTH_TYPE,
    COL_BOTH_RENAME,
    GENOME_DATABASE,
    GUIDE_COLUMN,
    CHR_COLUMN,
    POS_COLUMN,
    MM_COLUMN,
    BLG_COLUMN,
    TOTAL_COLUMN,
    TOTAL_FEWEST_COLUMN,
    BLG_T_COLUMN,
    CFD_COLUMN,
    CRISTA_COLUMN,
    RISK_COLUMN,
    SAMPLES_COLUMN,
    SAMPLES_CRISTA_COLUMN,
    SAMPLES_FEWEST_COLUMN,
    RESULTS_DIR,
    DATA_DIR,
    FILTERING_CRITERIA,
    PARAMS_FILE,
    drop_columns,
    write_json,
    read_json,
    get_query_column,
    split_filter_part
)

from typing import Dict, List, Tuple, Type, final
from glob import glob

import os

from sqlite3.dbapi2 import Row
import sys
from dash.exceptions import PreventUpdate
from dash.dependencies import Input, Output, State
from numpy.lib.function_base import _diff_dispatcher
from app import URL, app

# from app import app
import pandas as pd

# from datatable import dt, f, sort
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_table
from app import current_working_directory, cache, app_main_directory, operators
from PostProcess import CFDGraph
from PostProcess.supportFunctions.loadSample import associateSample
from os.path import isfile, isdir, join  # for getting directories

from os import listdir
import subprocess
import math
import base64  # for decoding upload content
import time
import re
import webbrowser as wb
import sqlite3
from PostProcess import query_manager
import flask


#-------------------------------------------------------------------------------
# Result page layout 
#

def result_page(job_id: str) -> html.Div:
    """Print the results page layout (guides table + images).
    The guides table contains the research profile found during 
    target search. Creates 10 buttons (mismatch number + 2), the 
    remaining ones are set to style = {"display":None}, in order 
    to have the right number of buttons, based on mismatches required 
    in input during the target search. This choice solves some 
    callback issues that have in input elements not created. In this
    case, all the possible buttons are created, but are shown only 
    those correct based on the selected number of mismatches.

    ...

    Parameters
    ----------
    job_id : str
        Unique job identifier

    Returns
    -------
    html.Div
        Result page layout
    """

    # check input function arguments
    if not isinstance(job_id, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(job_id).__name__}")
    # start result page creation code
    value = job_id
    job_directory = os.path.join(
        current_working_directory, "Results", f"{job_id}")
    integrated_file_name = glob(
        os.path.join(current_working_directory, "Results",
                     f"{job_id}", "*integrated*")
    )[
        0
    ]  # take the first list element
    assert isinstance(integrated_file_name, str)
    # integrated_file_name = str(integrated_file_name)
    integrated_file_name_zip = integrated_file_name.replace("tsv", "zip")
    if not os.path.isdir(job_directory):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))
    count_guides = 0
    guides_file = os.path.join(
        current_working_directory, "Results", f"{value}", ".guides.txt"
    )
    assert os.path.isfile(guides_file)
    try:
        with open(guides_file) as handle:
            for line in handle:
                count_guides += 1
    except:
        raise IOError(f"Unable to read {guides_file}.")
    finally:
        handle.close()
    # Load mismatches
    try:
        with open(
            os.path.join(
                current_working_directory, RESULTS_DIR, value, PARAMS_FILE
            )
        ) as p:
            all_params = p.read()
            real_genome_name = (
                next(
                    s for s in all_params.split("\n") if "Genome_idx" in s
                )
            ).split("\t")[-1]
            mms = (
                next(
                    s for s in all_params.split("\n") if "Mismatches" in s
                )
            ).split("\t")[-1]
            bulge_dna = (
                next(
                    s for s in all_params.split("\n") if "DNA" in s
                )
            ).split("\t")[-1]
            bulge_rna = (
                next(
                    s for s in all_params.split("\n") if "RNA" in s
                )
            ).split("\t")[-1]
            genome_type_f = (
                next(
                    s for s in all_params.split("\n") if "Genome_selected" in s
                )
            ).split("\t")[-1]
            ref_comp = (
                next(
                    s for s in all_params.split("\n") if "Ref_comp" in s
                )
            ).split("\t")[-1]
            max_bulges = (
                next(
                    s for s in all_params.split("\n") if "Max_bulges" in s
                )
            ).split("\t")[-1]
            pam_name = (
                next(
                    s for s in all_params.split("\n") if "Pam" in s
                )
            ).split("\t")[-1]
    except OSError as e:
        raise e
    finally:
        p.close()
    # recover genome name
    genome_name = genome_type_f
    if "+" in real_genome_name:
        genome_name = [genome_name] + [
            name.split("+")[1] for name in real_genome_name.strip().split(",")
        ]
        genome_name = "+".join(genome_name)
    if "True" in ref_comp:
        genome_type = "both"
    else:
        genome_type = "ref"
    mms = int(mms[0])
    # load acfd for each guide
    acfd_file = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        "".join([".", job_id, ".acfd_CFD.txt"]),
    )
    if not os.path.isfile(acfd_file):
        # something went wrong
        raise FileNotFoundError(f"Unable to locate {acfd_file}")
    try:
        with open(acfd_file) as handle:
            all_scores = handle.read().strip().split("\n")
    except OSError as e:
        raise e
    finally:
        handle.close()

    guides_error_file = os.path.join(
        current_working_directory, RESULTS_DIR, job_id, "guides_error.txt"
    )
    list_error_guides = []
    if os.path.exists(guides_error_file):
        try:
            with open(guides_error_file) as handle_error_g:
                for e_g in handle_error_g:
                    list_error_guides.append(e_g.strip())
        except OSError as e:
            raise e
        finally:
            handle_error_g.close()
    col_targetfor = "("
    for i in range(1, (mms + int(max_bulges))):
        col_targetfor = "".join([col_targetfor, str(i), "-"])
    col_targetfor = "".join([col_targetfor, str(mms + int(max_bulges))])
    col_targetfor = " ".join([col_targetfor, "Mismatches + Bulges)"])
    # Column of headers. Remove the entries accordingly when checking genome type
    columns_profile_table = [
        {"name":["", "gRNA (spacer+PAM)"], "id":"Guide","type": "text"},
        {"name":["", "Nuclease", ""], "id":"Nuclease", "type":"text"},
        {
            "name":["", "Aggregated Specificity Score (0-100)"],
            "id":"CFD",
            "type":"text",
        },
        {
            "name":["Off-targets for Mismatch (MM) and Bulge (B) Value", "Total"],
            "id":"Total",
            "type":"text",
        },
    ]
    columns_profile_table.append(
        {
            "name":["Off-targets for Mismatch (MM) and Bulge (B) Value", "# Bulges"],
            "id":"# Bulges",
            "type":"text",
        }
    )
    for i in range(mms + 1):
        columns_profile_table.append(
            {
                "name":[
                    "Off-targets for Mismatch (MM) and Bulge (B) Value",
                    "".join([str(i), "MM"]),
                ],
                "id":"".join([str(i), "MM"]),
                "type":"text",
            }
        )
    remove_indices = set()
    if "NO SCORES" in all_scores:
        # remove CFD and Doench header from table
        remove_indices.add("CFD", "Doench 2016", "Reference", "Enriched")
    if genome_type == "ref":
        # remove reference header 
        remove_indices.update(["Reference", "Enriched"])
    else:
        # remove reference and reference target headers
        remove_indices.update(
            [
                "Reference",
                "Enriched",
                "On-Targets Reference",
                "Samples in Class 0 - 0+ - 1 - 1+",
            ]
        )
    # Remove headers not used in the selected search results
    columns_profile_table = [
        i
        for j, i in enumerate(columns_profile_table)
        if columns_profile_table[j]["id"] not in remove_indices
    ]
    # build final result page with corresponding fields
    final_list = []
    if list_error_guides:
        final_list.append(
            dbc.Alert(
                [
                    "Warning: Some guides have too many targets! ",
                    html.A(
                        "Click here",
                        href=os.path.join(URL, DATA_DIR, job_id, "guides_error.txt"),
                        className="alert-link",
                    ),
                    " to view them",
                ],
                color="warning",
            )
        )
    final_list.append(
        html.H3(
            " ".join(
                [
                    "Result Summary",
                    "-",
                    genome_name,
                    "-",
                    pam_name,
                    "-",
                    "Mismatches",
                    str(mms),
                    "-",
                    "DNA bulges",
                    bulge_dna,
                    "-",
                    "RNA bulges",
                    bulge_rna,
                ]
            )
        )
    )
    # short description
    if genome_type == "both":
        add_to_description = html.P(
            [
                str(
                    "General summary for input guides. For each guide, is "
                    "reported the count of targets in reference and variant "
                    "genome grouped by mismatches count and bulge size."
                ),
            ]
        )
    else:
        add_to_description = html.P(
            str(
                "General summary for input guides. For each guide, is reported the "
                "count of targets in reference and variant genome grouped by "
                "mismatches count and bulge size."
            )
        )
    final_list.append(add_to_description)  # add description line to page layout
    # define upper page box
    final_list.append(
        html.Div(
            dbc.Row(
                dbc.Col(
                    [
                        html.Div(
                            [
                                html.P(
                                    "Generating download link, Please wait...",
                                    id="download-link-general-table",
                                ),
                                dcc.Interval(
                                    interval=1 * 1000, id="interval-general-table"
                                ),
                                html.Div(
                                    os.path.join(
                                        current_working_directory,
                                        RESULTS_DIR,
                                        job_id,
                                        ".".join([job_id, "general_table.txt"])
                                    ),
                                    style={"display":"none"},
                                    id="div-info-general-table"
                                ),
                            ]
                        ),
                        html.Div(
                            [
                                html.P(
                                    "Generating download link, Please wait...",
                                    id="download-link-integrated-results",
                                ),
                                dcc.Interval(
                                    interval=1 * 1000, id="interval-integrated-results"
                                ),
                                html.Div(
                                    integrated_file_name_zip,
                                    style={"display": "none"},
                                    id="div-info-integrated-results",
                                ),
                            ]
                        ),
                    ]
                )
            )
        )
    )
    # results table (middle of page layout)
    final_list.append(
        html.Div(
            html.Div(
                dash_table.DataTable(
                    id="general-profile-table",
                    # page_size=PAGE_SIZE,
                    columns=columns_profile_table,
                    merge_duplicate_headers=True,
                    # fixed_rows={ 'headers': True, 'data': 0 },
                    # data = profile.to_dict('records'),
                    selected_cells=[{"row":0, "column":0}],
                    # layout CSS style
                    css=[
                        {
                            "selector":".row",
                            "rule":"margin: 0",
                            "selector":"td.cell--selected, td.focused",
                            "rule":"background-color: rgba(0, 0, 255,0.15) !important;",
                        },
                        {
                            "selector":"td.cell--selected *, td.focused *",
                            "rule":"background-color: rgba(0, 0, 255,0.15) !important;",
                        },
                    ],
                    page_current=0,
                    page_size=10,
                    page_action="custom",
                    # virtualization = True,
                    filter_action="custom",
                    filter_query="",
                    sort_action="custom",
                    sort_mode="multi",
                    sort_by=[],
                    style_table={
                        # 'margin-left': "10%",
                        "max-height":"260px",
                        "overflowY":"scroll",
                        # 'overflowX': 'hidden',
                    },
                    style_data={
                        "whiteSpace":"pre",
                        "height":"auto",
                        "font-size":"1.30rem",
                    },
                    # style_cell={
                    #    'width':f'{1/len(columns_profile_table)*100}%'
                    # },
                    style_data_conditional=[
                        {
                            "if":{"column_id":"Genome"},
                            "font-weight":"bold",
                            "textAlign":"center",
                        },
                        # {'if': {'column_id': 'Guide'},
                        #                    'width': '10%',
                        #                    }
                    ],
                    style_cell_conditional=[
                        {
                            "if":{"column_id": "Guide"},
                            "width":"20%",
                        },
                        {
                            "if":{"column_id": "Total"},
                            "width":"15%",
                        },
                        {
                            "if":{"column_id": "Doench 2016"},
                            "width":"5%",
                        },
                        {
                            "if":{"column_id": "# Bulges"},
                            "width":"5%",
                        },
                        {
                            "if":{"column_id": "Nuclease"},
                            "width":"5%",
                        },
                    ],
                    #                        {'if': {'column_id': 'Reference'},
                    #                        'width': '10%',
                    #                        }],
                ),
                id="div-general-profile-table",
                style={"margin-left":"5%", "margin-right":"5%"},
            )
        )
    )
    final_list.append(html.Br()) # add space between HTML lines
    # drop-down bar (filetring criterion selection)
    final_list.append(
        html.Div(
            dbc.Row(
                dbc.Col(
                    html.Div(
                        [
                            html.H4("Select filter criteria for targets"),
                            dcc.Dropdown(
                                options=[
                                    {"label":"CFD score", "value":"CFD"},
                                    {"label":"CRISTA Score", "value":"CRISTA"},
                                    {
                                        "label":"Fewest Mismatches and Bulges",
                                        "value":"fewest",
                                    },
                                ],
                                value="CFD",
                                id="target_filter_dropdown",
                            ),
                            dcc.Store(id="store"),
                        ]
                    )
                )
            ),
        )
    )
    final_list.append(html.Br())  # add space between HTML lines
    if genome_type == "ref":
        final_list.append(
            dcc.Tabs(
                id="tabs-reports",
                value="tab-query-table",
                children=[
                    dcc.Tab(label="Custom Ranking", value="tab-query-table"),
                    dcc.Tab(
                        label="Summary by Mismatches/Bulges",
                        value="tab-summary-by-guide",
                    ),
                    dcc.Tab(
                        label="Query Genomic Region", value="tab-summary-by-position"
                    ),
                    dcc.Tab(
                        label="Graphical Reports",
                        value="tab-summary-graphical"
                    ),
                ],
            )
        )
    else:
        # Barplot for population distributions
        final_list.append(
            html.Div(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Button(
                                    "Show/Hide Target Distribution in SuperPopulations",
                                    id="btn-collapse-populations",
                                )
                            ),
                        ]
                    ),
                    dbc.Collapse(
                        dbc.Card(
                            dbc.CardBody(
                                html.Div(id="content-collapse-population")
                            )
                        ),
                        id="collapse-populations",
                    ),
                ],
                hidden=True,
            )
        )
        final_list.append(html.Br())  # add space between HTML lines
        # define results page tabs
        final_list.append(
            dcc.Tabs(
                id="tabs-reports",
                value="tab-query-table",
                children=[
                    dcc.Tab(label="Custom Ranking", value="tab-query-table"),
                    dcc.Tab(
                        label="Summary by Mismatches/Bulges",
                        value="tab-summary-by-guide",
                    ),
                    dcc.Tab(
                        label="Summary by Sample",
                        value="tab-summary-by-sample"
                    ),
                    dcc.Tab(
                        label="Query Genomic Region", value="tab-summary-by-position"
                    ),
                    dcc.Tab(
                        label="Graphical Reports",
                        value="tab-summary-graphical"
                    ),
                    dcc.Tab(
                        label="Personal Risk Cards", value="tab-graphical-sample-card"
                    ),
                ],
            )
        )
    final_list.append(html.Div(id="div-tab-content"))

    final_list.append(
        html.Div(genome_type, style={"display":"none"}, id="div-genome-type")
    )
    result_page = html.Div(final_list, style={"margin":"1%"})
    return result_page


# store drop-down value in auxiliary file
@app.callback(
    Output("store", "data"),
    [Input("target_filter_dropdown", "value")],
    [State("url", "search")]
)
def sendto_write_json(filter_criterion: str, search: str) -> None:
    """Write auxiliary file to store the table filtering criterion
    (received from the drop-down) and filter the tables displayed in 
    Summary by Mismatches/Bulges accordingly.

    The function is triggered by the user, when choosing the filtering
    criterion from the drop-down bar.
    
    ...

    Parameters
    ----------
    filter_criterion : str
        Table filtering criterion
    search : str
        Target search name

    Returns 
    -------
    None
    """
    if not isinstance(filter_criterion, str):
        raise TypeError(f"Expected {str.__name__}, got {type(filter_criterion).__name__}")
    if not filter_criterion in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filtering criterion ({filter_criterion})")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    job_id = search.split("=")[-1]
    write_json(filter_criterion, job_id)


#-------------------------------------------------------------------------------
# Download links generation and actions definition
#

# Generate download link summary_by_sample
@app.callback(
    [
        Output("download-link-summary_by_sample", "children"),
        Output("interval-summary_by_sample", "disabled"),
    ],
    [Input("interval-summary_by_sample", "n_intervals")],
    [State("div-info-summary_by_sample", "children"), State("url", "search")],
)
def download_link_sample(
    n: int, file_to_load: str, search: str
) -> Tuple[str, bool]:  # file to load =
    """Create the link to download CRISPRme result files.
    
    ...

    Parameters
    ----------
    n : int
    file_to_load : str
        File to download
    search : str
        Target search name

    Returns
    -------
    str 
    bool
    """

    if not isinstance(file_to_load, str):
        raise TypeError(f"Expected {str.__name__}, got {type(file_to_load).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if n is None:
        raise PreventUpdate  # nothing to do
    job_id = search.split("=")[-1]
    file_to_load = ".".join([file_to_load, "txt"])
    file_to_load = file_to_load.strip().split("/")[-1]
    # print(file_to_load)
    if os.path.exists(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, file_to_load)
    ):
        return (
            html.A(
                "Download file",
                href=os.path.join(URL, RESULTS_DIR, job_id, file_to_load),
                target="_blank",
            ),
            True,
        )
    return "Generating download link, Please wait...", False


# download summary result table
@app.callback(
    [
        Output("download-link-general-table", "children"),
        Output("interval-general-table", "disabled"),
    ],
    [Input("interval-general-table", "n_intervals")],
    [State("div-info-general-table", "children"), State("url", "search")],
)
def download_general_table(
    n: int, file_to_load: str, search: str
) -> Tuple[str, bool]:  # file to load =
    """Create the link to download CRISPRme result summary table.
    
    ...

    Parameters
    ----------
    n : int
    file_to_load : str
        File to download
    search : str
        Target search name

    Returns
    -------
    str 
    bool
    """

    if not isinstance(file_to_load, str):
        raise TypeError(f"Expected {str.__name__}, got {type(file_to_load).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if n is None:
        raise PreventUpdate
    job_id = search.split("=")[-1]
    file_to_load = file_to_load.split("/")[-1]
    # print(file_to_load)
    if os.path.exists(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, file_to_load)
    ):
        return (
            html.A(
                "Download General Table",
                href=os.path.join(URL, RESULTS_DIR, job_id, file_to_load),
                target="_blank",
            ),
            True,
        )
    return "Generating download link, Please wait...", False


# download integrated results
@app.callback(
    [
        Output("download-link-integrated-results", "children"),
        Output("interval-integrated-results", "disabled"),
    ],
    [Input("interval-integrated-results", "n_intervals")],
    [State("div-info-integrated-results", "children"), State("url", "search")],
)
def download_general_table(
    n: int, file_to_load: str, search: str
) -> Tuple[str, bool]:  # file to load =
    """Create the link to download CRISPRme integrated result table.
    
    ...

    Parameters
    ----------
    n : int
    file_to_load : str
        File to download
    search : str
        Target search name

    Returns
    -------
    str 
    bool
    """

    if not isinstance(file_to_load, str):
        raise TypeError(f"Expected {str.__name__}, got {type(file_to_load).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if n is None:
        raise PreventUpdate
    job_id = search.split("=")[-1]
    file_to_load = file_to_load.split("/")[-1]
    # print(file_to_load)
    if os.path.exists(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, file_to_load)
    ):
        return (
            html.A(
                "Download Integrated Results",
                href=os.path.join(URL, RESULTS_DIR, job_id, file_to_load),
                target="_blank",
            ),
            True,
        )
    return "Generating download link, Please wait...", False


# Generate download link sumbysample
@app.callback(
    [
        Output("download-link-sumbysample", "children"),
        Output("interval-sumbysample", "disabled"),
    ],
    [Input("interval-sumbysample", "n_intervals")],
    [State("div-info-sumbysample-targets", "children"), State("url", "search")],
)
def download_link_sample(
    n: int, file_to_load: str, search: str
) -> Tuple[str, bool]:  # file to load = job_id.HG001.guide
    """Create the link to download CRISPRme results by sample table.
    
    ...

    Parameters
    ----------
    n : int
    file_to_load : str
        File to download
    search : str
        Target search name

    Returns
    -------
    str 
    bool
    """

    if not isinstance(file_to_load, str):
        raise TypeError(f"Expected {str.__name__}, got {type(file_to_load).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if n is None:
        raise PreventUpdate
    job_id = search.split("=")[-1]
    file_to_load = ".".join([file_to_load, "zip"])
    if os.path.exists(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, file_to_load)
    ):
        return (
            html.A(
                "Download zip",
                href=os.path.join(URL, RESULTS_DIR, job_id, file_to_load),
                target="_blank",
            ),
            True,
        )
    return "Generating download link, Please wait...", False


# Generate download link sumbyguide
@app.callback(
    [
        Output("download-link-sumbyguide", "children"),
        Output("interval-sumbyguide", "disabled"),
    ],
    [Input("interval-sumbyguide", "n_intervals")],
    [State("div-info-sumbyguide-targets", "children"), State("url", "search")],
)
def downloadLinkGuide(
    n: int, file_to_load: str, search: str
) -> Tuple[str, bool]:  # file to load = job_id.RNA.1.0.guide
    """Create the link to download CRISPRme results by sample table.
    
    ...

    Parameters
    ----------
    n : int
    file_to_load : str
        File to download
    search : str
        Target search name

    Returns
    -------
    str 
    bool
    """

    if not isinstance(file_to_load, str):
        raise TypeError(f"Expected {str.__name__}, got {type(file_to_load).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if n is None:
        raise PreventUpdate
    job_id = search.split("=")[-1]
    file_to_load = ".".join([file_to_load, "zip"])
    if os.path.exists(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, file_to_load)
    ):
        return (
            html.A(
                "Download zip",
                href=os.path.join(URL, RESULTS_DIR, job_id, file_to_load),
                target="_blank",
            ),
            True,
        )
    return "Generating download link, Please wait...", False


# trigger file download
@app.server.route("/Results/<path:path>")
def download_file(path: str) -> flask.Response:
    """Download the chosen file.
    
    ...

    Parameters
    ----------
    path : str
        Path to file location

    Returns
    -------
    flask.Response
    """
    
    if not isinstance(path, str):
        raise TypeError(f"Expected {str.__name__}, got {type(path).__name__}")
    # print(current_working_directory)
    # print('test', path)
    return flask.send_from_directory(
        os.path.join(current_working_directory, "Results/"), path, as_attachment=True
    )


# Filter/sort IUPAC decomposition table for cluster page
@app.callback(
    Output("table-scomposition-cluster", "data"),
    [
        Input("table-scomposition-cluster", "page_current"),
        Input("table-scomposition-cluster", "page_size"),
        Input("table-scomposition-cluster", "sort_by"),
        Input("table-scomposition-cluster", "filter_query"),
    ],
    [State("url", "search"), State("url", "hash")],
)
def update_iupac_decomposition_table_cluster(
    page_current: int, 
    page_size: int, 
    filter_criterion: str, 
    search: str, 
    hash_term: str
) -> Dict[str, str]:
    """

    ...

    Parameters
    ----------
    page_current : int
        Current page
    page_size : int
        Page size
    filter_criterion : str
        Data table filter
    search : str
        Unique search ID
    hash_term : str
        Hashing

    Returns 
    -------
    Dict[str, str]
    """

    if not isinstance(page_current, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_current).__name__}")
    if not isinstance(page_size, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_size).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(f"Expected {str.__name__}, got {type(filter_criterion).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if not isinstance(hash_term, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_term).__name__}")
    job_id = search.split("=")[-1]
    hash_term = hash_term.split("#")[1]
    guide = hash_term[:hash_term.find("-Pos-")]
    chr_pos = hash_term[(hash_term.find("-Pos-") + 5):]
    chromosome = chr_pos.split("-")[0]
    position = chr_pos.split("-")[1]
    try:
        with open(
            os.path.join(current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE)
        ) as handle:
            all_params = handle.read()
            genome_type_f = (
                next(s for s in all_params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (
                next(s for s in all_params.split("\n") if "Ref_comp" in s)
            ).split("\t")[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
    if genome_type == "ref":
        raise PreventUpdate
    filtering_expressions = filter_criterion.split(" && ")
    decomp_fname = (
        job_id + "." + chromosome + "_" + position + "." + guide + ".scomposition.txt"
    )
    # load data and cache the data table (in pd.DataFrame)
    df_cached = global_store_general(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, decomp_fname)
    )
    if df_cached is None:  #  nothing to display and do not update the page
        raise PreventUpdate
    df_cached.rename(columns=COL_BOTH_RENAME, inplace=True)
    # filter data table
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)
        if operator in ("eq", "ne", "lt", "le", "gt", "ge"):
            # these operators match pandas series operator method names
            df_cached = df_cached.loc[
                getattr(df_cached[col_name], operator)(filter_value)
            ]
        elif operator == "contains":
            df_cached = df_cached.loc[
                df_cached[col_name].str.contains(filter_value)
            ]
        elif operator == "datestartswith":
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            df_cached = df_cached.loc[
                df_cached[col_name].str.startswith(filter_value)
            ]
    # Calculate sample count
    data_to_send = df_cached.iloc[
        page_current * page_size: (page_current + 1) * page_size
    ].to_dict("records")
    return data_to_send


@app.callback(
    Output("table-position-target", "data"),
    [
        Input("table-position-target", "page_current"),
        Input("table-position-target", "page_size"),
        Input("table-position-target", "sort_by"),
        Input("table-position-target", "filter_query"),
        Input("hide-reference-targets", "value"),
    ],
    [State("url", "search"), State("url", "hash")],
)
def update_table_cluster(
    page_current: int, 
    page_size: int, 
    sort_by: List[str], 
    filter_criterion: str, 
    hide_reference: str, 
    search: str, 
    hash_term: str
) -> Dict[str, str]:
    """

    ...

    Parameters
    ----------
    page_current : int
        Current page
    page_size : int
        Page size
    sort_by : List[str]
        Columns used while sorting the data table
    filter_criterion : str
        Data table filter
    hide_reference : str
        Hide reference data
    search : str
        Unique search ID
    has_term : str
        Hashing

    Returns
    -------
    Dict[str, str]
    """
    
    if not isinstance(page_current, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_current).__name__}")
    if not isinstance(page_size, int):
        raise TypeError(f"Exepcted {int.__name__}, got {type(page_size).__name__}")
    if not isinstance(sort_by, list):
        raise TypeError(f"Expected {list.__name__}, got {type(sort_by).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(f"Exepcted {str.__name__}, got {type(filter_criterion).__name__}")
    if not isinstance(hide_reference, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hide_reference).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if not isinstance(hash_term, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_term).__name__}")
    job_id = search.split("=")[-1]
    job_directory = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    hash_term = hash_term.split("#")[1]
    guide = hash_term[: hash_term.find("-Pos-")]
    chr_pos = hash_term[hash_term.find("-Pos-") + 5:]
    chromosome = chr_pos.split("-")[0]
    position = chr_pos.split("-")[1]
    try:
        with open(
            os.path.join(
                current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE
            )
        ) as handle:
            all_params = handle.read()
            genome_type_f = (
                next(s for s in all_params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (
                next(s for s in all_params.split("\n") if "Ref_comp" in s)
            ).split("\t")[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
    filtering_expressions = filter_criterion.split(" && ")
    guide_fname = job_id + "." + chromosome + "_" + position + "." + guide + ".txt"
    # cache guide data table
    df_cached = global_store_general(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, guide_fname)
    )
    if df_cached is None:  # empty file -> nothing cached and nothing to do
        raise PreventUpdate
    if genome_type == "ref":
        df_cached.rename(columns=COL_BOTH_RENAME, inplace=True)
    else:
        df_cached.rename(columns=COL_BOTH_RENAME, inplace=True)
    # drop unused columns
    if "hide-ref" in hide_reference or genome_type == "var":
        df_cached.drop(
            df_cached[(df_cached["Samples"] == "n")].index, inplace=True
        )
    # hide reference data
    if "hide-cluster" in hide_reference:
        df_cached = df_cached.head(1)
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)
        if operator in ("eq", "ne", "lt", "le", "gt", "ge"):
            # these operators match pandas series operator method names
            df_cached = df_cached.loc[
                getattr(df_cached[col_name], operator)(filter_value)
            ]
        elif operator == "contains":
            df_cached = df_cached.loc[
                df_cached[col_name].str.contains(filter_value)
            ]
        elif operator == "datestartswith":
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            df_cached = df_cached.loc[
                df_cached[col_name].str.startswith(filter_value)
            ]
    # sort data table by the defined columns
    if bool(sort_by):
        df_cached = df_cached.sort_values(
            [
                "Samples" if col["column_id"] == "Samples Summary" else col["column_id"]
                for col in sort_by
            ],
            ascending=[col["direction"] == "asc" for col in sort_by],
            inplace=False,
        )
    # Calculate sample count
    data_to_send = df_cached.iloc[
        (page_current * page_size):((page_current + 1) * page_size)
    ].to_dict("records")
    if genome_type != "ref":
        (
            dict_sample_to_pop,
            dict_pop_to_superpop,
        ) = associateSample.loadSampleAssociation(
            job_directory + ".sampleID.txt"
        )[:2]
        for row in data_to_send:
            summarized_sample_cell = dict()
            for s in row["Samples"].split(","):
                if s == "n":
                    break
                try:
                    summarized_sample_cell[
                        dict_pop_to_superpop[dict_sample_to_pop[s]]
                    ] += 1
                except:
                    summarized_sample_cell[
                        dict_pop_to_superpop[dict_sample_to_pop[s]]
                    ] = 1
            if summarized_sample_cell:
                row["Samples Summary"] = ", ".join(
                    [
                        str(summarized_sample_cell[sp]) + " " + sp
                        for sp in summarized_sample_cell
                    ]
                )
            else:
                row["Samples Summary"] = "n"
    return data_to_send


def cluster_page(job_id: str, hash_term: str) -> html.Div:
    """Recover CRISPR targets for the selected cluster.

    ...

    Parameters
    ----------
    job_id : str
        Unique job identifier
    hash_term : str
        Hashing

    Returns 
    -------
    html.Div
        Sample page layout
    """

    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if not isinstance(hash_term, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_term).__name__}")
    guide = hash_term[:hash_term.find("-Pos-")]
    chr_pos = hash_term[(hash_term.find("-Pos-") + 5):]
    chromosome = chr_pos.split("-")[0]
    position = chr_pos.split("-")[1]
    if not os.path.isdir(
        os.path.join(current_working_directory, RESULTS_DIR, job_id)
    ):
        return html.Div(
            dbc.Alert("The selected result does not exist", color="danger")
        )
    try:
        with open(
            os.path.join(
                current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE
            )
        ) as handle_params:
            params = handle_params.read()
            genome_type_f = (
                next(s for s in params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (
                next(s for s in params.split("\n") if "Ref_comp" in s)
            ).split("\t")[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    style_hide_reference = {"display":"none"}  # display reference data
    value_hide_reference = []
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
        style_hide_reference = {}
        value_hide_reference = ["hide-ref", "hide-cluster"]  # hide reference data
    # begin page body construction
    final_list = []  # HTML page handler
    assert isinstance(chromosome, str)
    assert isinstance(position, str)
    final_list.append(
        html.H3(f"Selected Position: {chromosome} - {position}")
    )
    if genome_type == "ref":
        cols = [
            {"name":i, "id":i, "type":t, "hideable":True}
            for i, t in zip(COL_BOTH, COL_BOTH_TYPE)
        ]
        file_to_grep = ".bestMerge.txt"
    else:
        cols = [
            {"name":i, "id":i, "type":t, "hideable":True}
            for i, t in zip(COL_BOTH, COL_BOTH_TYPE)
        ]
        file_to_grep = ".bestMerge.txt"
    cluster_grep_result = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        ".".join([job_id, f"{chromosome}_{position}", guide, "txt"])
    )
    put_header_cmd = " ".join(
        [
            "head -1",
            os.path.join(
                current_working_directory, 
                RESULTS_DIR, 
                job_id,
                f".{job_id}{file_to_grep}"
            ),
            f"> {cluster_grep_result} ; "
        ]
    )
    # Example    job_id.chr3_100.guide.txt
    if not os.path.exists(cluster_grep_result):
        # os.system(f'touch {cluster_grep_result}')
        # Grep annotation for ref
        cmd = f"head -1 {file_to_grep} > {cluster_grep_result}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise ValueError(f"An error occurred while running {cmd}")
        if genome_type == "ref":  # NOTE HEADER NOT SAVED
            cmd = " ".join(
                [
                    "grep -F",
                    guide,
                    os.path.join(
                        current_working_directory,
                        RESULTS_DIR,
                        job_id,
                        f"{job_id}.Annotation.targets.txt"
                    ),
                    "|",
                    f"awk '$6=={position} && $4==\"{chromosome}\"'"
                ]
            )
            get_annotation = subprocess.Popen(
                [cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
            out, err = get_annotation.communicate()
            annotation_type = out.decode("UTF-8").strip().split("\t")[-1]
            os.popen(
                put_header_cmd
                + " grep -F "
                + guide
                + " "
                + current_working_directory
                + "Results/"
                + job_id
                + "/"
                + job_id
                + file_to_grep
                + " | awk '$6=="
                + position
                + ' && $4=="'
                + chromosome
                + '" {###print $0"\\t'
                + annotation_type
                + "\"}' >> "
                + cluster_grep_result
            ).read()
        else:  # NOTE HEADER NOT SAVED
            os.popen(
                " ".join(
                    [
                        put_header_cmd,
                        "grep -F",
                        guide,
                        os.path.join(
                            current_working_directory,
                            RESULTS_DIR,
                            job_id,
                            f"{job_id}{file_to_grep}"
                        ),
                        "|",
                        f"awk '$6=={position} && $4==\"{chromosome}\"'",
                        ">>",
                        cluster_grep_result
                    ]
                )
            ).read()  
            # NOTE top1 will have sample and annotation, other targets will 
            # have '.'-> 18/03 all samples and annotation are already writter 
            # for all targets

        # TODO: review this part    
        os.system(
            f"python {app_main_directory}/PostProcess/change_headers_bestMerge.py {cluster_grep_result} {cluster_grep_result}.tmp"
        )
        os.system(
            f"mv -f {cluster_grep_result}.tmp {cluster_grep_result} > /dev/null 2>&1"
        )
        # zip cluster results
        cmd = f"zip -j {cluster_grep_result.replace('txt', 'zip')} {cluster_grep_result} &"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise ValueError(f"An error occurred while running {cmd}")
    final_list.append(
        html.Div(
            f"{job_id}.{chromosome}_{position}.{guide}",
            style={"display":"none"},
            id="div-info-sumbyposition-targets",
        )
    )
    decomp_fname = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        f"{job_id}.{chromosome}_{position}.{guide}.scomposition.txt"
    )
    iupac_decomp_visibility = {"display":"none"}
    if genome_type != "ref":
        iupac_decomp_visibility = {}
        # Example    job_id.chr_pos.guide.scomposition.txt
        # if not os.path.exists(scomposition_file):
        # os.system(f'touch {scomposition_file}')
        cmd = " ".join(
            [
                "grep -F",
                guide,
                os.path.join(
                    current_working_directory,
                    RESULTS_DIR,
                    job_id,
                    f".{job_id}{file_to_grep}"
                ),
                "|",
                f"awk '$6=={position} && $4==\"{chromosome}\" && $13!=\"n\"'",
                ">",
                decomp_fname
            ]
        )
        os.popen(cmd).read()
    final_list.append(
        html.P(
            [
                html.P(
                    "List of all the configurations for the target in the selected position.",
                    style=iupac_decomp_visibility,
                ),
                dcc.Checklist(
                    options=[
                        {"label":"Hide Reference Targets", "value":"hide-ref"},
                        {"label":"Show only TOP1 Target", "value":"hide-cluster"},
                    ],
                    id="hide-reference-targets",
                    value=value_hide_reference,
                    style=style_hide_reference,
                ),
                html.Div(
                    [
                        html.P(
                            "Generating download link, Please wait...",
                            id="download-link-sumbyposition",
                        ),
                        dcc.Interval(
                            interval=5 * 1000, id="interval-sumbyposition"
                        ),
                    ]
                ),
            ]
        )
    )
    cols_for_decomp = cols.copy()
    cols_for_decomp.append(
        {"name":"Samples", "id":"Samples", "type":"text", "hideable":True}
    )
    final_list.append(
        html.Div(
            dash_table.DataTable(
                # Table storing IUPAC decomposition of the selected target
                # rows are recovered from top1.samples.txt
                id="table-scomposition-cluster",
                columns=cols_for_decomp,
                virtualization=True,
                fixed_rows={"headers":True, "data":0},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action="custom",
                sort_action="custom",
                sort_mode="multi",
                sort_by=[],
                filter_action="custom",
                filter_query="",
                style_table={"max-height": "600px"},
                style_cell_conditional=[
                    {
                        "if":{"column_id":"Variant_samples_(highest_CFD)"},
                        "textAlign":"left",
                        "minWidth":"180px",
                        "width":"180px",
                        "maxWidth":"180px",
                        "overflow":"hidden",
                    },
                    {
                        "if":{"column_id":"Variant_samples_(fewest_mm+b)"},
                        "textAlign":"left",
                        "minWidth":"180px",
                        "width":"180px",
                        "maxWidth":"180px",
                        "overflow":"hidden",
                    },
                    {
                        "if": {"column_id":"Variant_samples_(highest_CRISTA)"},
                        "textAlign":"left",
                        "minWidth":"180px",
                        "width":"180px",
                        "maxWidth":"180px",
                        "overflow":"hidden",
                    },
                ],
                css=[
                    {"selector":".row", "rule":"margin: 0"},
                    {
                        "selector":"td.cell--selected, td.focused",
                        "rule":"background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                    {
                        "selector":"td.cell--selected *, td.focused *",
                        "rule":"background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                ],
            ),
            style=iupac_decomp_visibility,
        )
    )
    final_list.append(html.Hr())
    # Build cluster Table
    final_list.append(
        # if rows are highlighted in red, the target was found only in 
        # non-reference genome (enriched with variants)
        str(
            "List of Targets found for the selected position. Other possible "
            "configurations of the target are listed in the table above, along " 
            "with the corresponding samples list."
        ),
    )
    final_list.append(
        html.Div(
            dash_table.DataTable(
                id="table-position-target",
                columns=cols,
                virtualization=True,
                fixed_rows={"headers":True, "data":0},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action="custom",
                sort_action="custom",
                sort_mode="multi",
                sort_by=[],
                filter_action="custom",
                filter_query="",
                style_table={
                    "max-height":"600px",
                    "overflowY":"scroll",
                },
                style_cell_conditional=[
                    {
                        "if":{"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign":"left",
                        "minWidth":"180px",
                        "width":"180px",
                        "maxWidth":"180px",
                        "overflow":"hidden",
                    },
                    {
                        "if": {"column_id":"Variant_samples_(fewest_mm+b)"},
                        "textAlign":"left",
                        "minWidth":"180px",
                        "width":"180px",
                        "maxWidth":"180px",
                        "overflow":"hidden",
                    },
                    {
                        "if":{"column_id":"Variant_samples_(highest_CRISTA)"},
                        "textAlign":"left",
                        "minWidth":"180px",
                        "width":"180px",
                        "maxWidth":"180px",
                        "overflow":"hidden",
                    },
                ],
                css=[
                    {"selector":".row", "rule":"margin: 0"},
                    {
                        "selector":"td.cell--selected, td.focused",
                        "rule":"background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                    {
                        "selector":"td.cell--selected *, td.focused *",
                        "rule":"background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                ],
            ),
            id="div-result-table",
        )
    )
    return html.Div(final_list, style={"margin":"1%"})


# Filter and sorting sample targets

#-------------------------------------------------------------------------------
# Summary by Sample tab
#
def global_get_sample_targets(
    job_id: str, sample: str, guide: str, page: int
) -> pd.DataFrame:
    """Recover CRISPRme analysis report regarding the selected sample.
    The sample related report can be filtered using the criteria available
    in the drop-down bar, above the report tabs:
    - CFD score
    - CRISTA score
    - Fewest Mismatches and Bulges

    ...

    Parameters
    ----------
    job_id : str
        Unique job identifier
    sample : str
        Sample identifier
    guide : str
        CRISPR guide
    page : int
        Current page 
    
    Returns
    -------
    pd.DataFrame
        Data table reporting CRISPRme analysis results related to the
        selected sample
    """

    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if not isinstance(sample, str):
        raise TypeError(f"Expected {str.__name__}, got {type(sample).__name__}")
    if not isinstance(guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide).__name__}")
    if not isinstance(page, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page).__name__}")
    if job_id is None:
        return ""
    db_path = glob(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, ".*.db")
    )[0]
    assert isinstance(db_path, str)
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    filter_criterion = read_json(job_id)  # recover filter criterion selected
    query_cols = get_query_column(filter_criterion)
    # query the db
    result = pd.read_sql_query(
        "SELECT * FROM final_table WHERE \"{}\"='{}' AND \"{}\" LIKE '%{}%' LIMIT {} OFFSET {}".format(
            GUIDE_COLUMN, 
            guide, 
            query_cols["samples"], 
            sample, 
            PAGE_SIZE, 
            page * PAGE_SIZE
        ),
        conn,
    )
    return result


# callback to update the samples table 
@app.callback(
    [Output("table-sample-target", "data"),
     Output("table-sample-target", "columns")],
    [
        Input("table-sample-target", "page_current"),
        Input("table-sample-target", "page_size"),
        Input("table-sample-target", "sort_by"),
        Input("table-sample-target", "filter_query"),
    ],
    [State("url", "search"), State("url", "hash")],
)
def update_table_sample(
    page_current: int, 
    page_size: int, 
    sort_by: str, 
    filter_criterion: str, 
    search: str, 
    hash_term : str
) -> Tuple[Dict[str, str], pd.DataFrame]:
    """Update the sample table accordingly to the filtering criterion 
    selected in the drop-down bar.

    ...

    Parameters
    ----------
    page_current : int
        Current webpage
    page_size : int
        Webpage size
    sort_by : str
        Data table sorting criterion
    filter_criterion : str
        Data table filtering criterion
    search : str
        Search identifier
    hash_term : str
        Hashing term

    Returns
    -------
    Tuple[Dict[str, str], pd.DataFrame]
    """

    if not isinstance(page_current, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_current).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if not isinstance(hash_term, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_term).__name__}")
    job_id = search.split("=")[-1]
    filter_criterion = read_json(job_id)  # recover filter criterion
    assert isinstance(filter_criterion, str)
    assert filter_criterion in FILTERING_CRITERIA
    hash_term = hash_term.split("#")[1]
    guide = hash_term[:hash_term.find("-Sample-")]
    sample = str(hash_term[hash_term.rfind("-") + 1:])
    try:
        with open(
            os.path.join(
                current_working_directory,
                RESULTS_DIR,
                job_id,
                ".Params.txt"
            )
        ) as handle:
            all_params = handle.read()
            genome_type_f = (
                next(s for s in all_params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (
                next(s for s in all_params.split("\n") if "Ref_comp" in s)
            ).split("\t")[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
    # populate the sample table
    sample_df = global_get_sample_targets(job_id, sample, guide, page_current)
    # filter the sample table
    drop_cols = drop_columns(sample_df, filter_criterion)
    sample_df.drop(drop_cols, inplace=True, axis=1)
    # personal targets report filename
    integrated_sample_personal_fname = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        ".".join(
            [
                job_id,
                sample,
                guide,
                "personal_targets.tsv"
            ]
        )
    )
    # store sample table to personal targets file
    sample_df.to_csv(
        integrated_sample_personal_fname, sep='\t', na_rep='NA', index=False
    )
    # personal targets report ZIP
    integrated_sample_personal_zip_fname = integrated_sample_personal_fname.replace(
        "tsv", "zip"
    )
    # zip operation, non blocking
    cmd =  f"zip -j {integrated_sample_personal_zip_fname} {integrated_sample_personal_fname} &"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise ValueError(f"An error occurred while running \"{cmd}\"")
    columns_df = [
        {"name":i, "id":i, "hideable":True} 
        for col, i in enumerate(sample_df.columns)
    ]
    return sample_df.to_dict("records"), columns_df


# Return the targets found for the selected sample


def samplePage(job_id, hash):
    # ###print("SAMPLE PAGE LOADED FOR", job_id, hash)
    guide = hash[: hash.find("-Sample-")]
    sample = str(hash[hash.rfind("-") + 1:])
    if not isdir(current_working_directory + "Results/" + job_id):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))

    with open(current_working_directory + "Results/" + job_id + "/.Params.txt") as p:
        all_params = p.read()
        genome_type_f = (
            next(s for s in all_params.split("\n") if "Genome_selected" in s)
        ).split("\t")[-1]
        ref_comp = (next(s for s in all_params.split("\n") if "Ref_comp" in s)).split(
            "\t"
        )[-1]

    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"

    final_list = []
    final_list.append(
        # html.P('List of Targets found for the selected Sample - ' + sample + ' - and guide - ' + guide + ' -')
        html.H3("Selected Sample: " + sample)
    )
    final_list.append(
        html.P(
            [
                # 'The rows highlighted in red indicates that the target was found only in the genome with variants.',
                "List of Targets found for the selected sample.",
                html.Div(
                    [
                        html.P(
                            "Generating download link, Please wait...",
                            id="download-link-sumbysample",
                        ),
                        dcc.Interval(interval=5 * 1000,
                                     id="interval-sumbysample"),
                    ]
                ),
            ]
        )
    )

    header = current_working_directory + "Results/" + job_id + "/header.txt"

    # file_to_grep = current_working_directory + 'Results/' + \
    #     job_id + '/.' + job_id + '.bestMerge.txt'
    integrated_file_name = glob(
        current_working_directory + "Results/" + job_id + "/" + "*integrated*"
    )[0]
    integrated_file_name = str(integrated_file_name)
    file_to_grep = (
        current_working_directory
        + "Results/"
        + job_id
        + "/"
        + job_id
        + ".bestMerge.txt.integrated_results.tsv"
    )
    sample_grep_result = (
        current_working_directory
        + "Results/"
        + job_id
        + "/"
        + job_id
        + "."
        + sample
        + "."
        + guide
        + ".txt"
    )
    integrated_sample_personal = (
        current_working_directory
        + "Results/"
        + job_id
        + "/"
        + job_id
        + "."
        + sample
        + "."
        + guide
        + ".personal_targets.tsv"
    )
    integrated_sample_personal_zip = integrated_sample_personal.replace(
        "tsv", "zip")
    final_list.append(
        html.Div(
            job_id + "." + sample + "." + guide + ".personal_targets",
            style={"display": "none"},
            id="div-info-sumbysample-targets",
        )
    )

    path_db = glob(current_working_directory +
                   "Results/" + job_id + "/.*.db")[0]
    path_db = str(path_db)
    conn = sqlite3.connect(path_db)
    c = conn.cursor()
    total_private_sample = f"SELECT * FROM final_table LIMIT 1"
    rows = c.execute(total_private_sample)
    header = [description[0] for description in rows.description]
    conn.commit()
    conn.close()

    cols = [{"name": i, "id": i, "hideable": True} for i in header]

    final_list.append(
        html.Div(
            dash_table.DataTable(
                id="table-sample-target",
                columns=cols,
                style_cell={"textAlign": "left"},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action="custom",
                style_table={
                    "overflowX": "scroll",
                    "overflowY": "scroll",
                    "max-height": "300px",
                },
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
                css=[
                    {"selector": ".row", "rule": "margin: 0"},
                    {
                        "selector": "td.cell--selected, td.focused",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                    {
                        "selector": "td.cell--selected *, td.focused *",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                ],
            ),
            id="div-result-table",
        )
    )
    return html.Div(final_list, style={"margin": "1%"})


@cache.memoize()
def global_store_general(path_file_to_load: str) -> pd.DataFrame:
    """Cache target files to improve results visualization and get better
    performances.

    ...

    Parameters
    ----------
    path_file_to_load : str
        Path to file to cache

    Returns
    -------
    pandas.DataFrame
        Results table
    """

    if not isinstance(path_file_to_load, str):
        raise TypeError(f"Expected {str.__name__}, got {type(path_file_to_load).__name__}")
    if path_file_to_load is not None and not os.path.isfile(path_file_to_load):
        raise FileNotFoundError(f"Unable to locate {path_file_to_load}")
    if path_file_to_load is None:
        return ""  # do not cache anything
    if "scomposition" in path_file_to_load:
        rows_to_skip = 1
    else:
        rows_to_skip = 1  # Skip header
    # make sure file to cache is not empty
    if os.path.getsize(path_file_to_load) > 0:  
        # TSV format -> sep="\t"
        df = pd.read_csv(
            path_file_to_load, sep="\t", index_col=False, na_filter=False
        )
    else:
        df = None  # empty file, no need for caching
    return df


# Filter etc for second table


# Update primary table of 'Show targets' of Summary by Guide
@app.callback(
    [
        Output("table-subset-target", "data"),
        Output("table-subset-target", "columns"),
    ],
    [
        Input("table-subset-target", "page_current"),
        Input("table-subset-target", "page_size"),
        Input("table-subset-target", "sort_by"),
        Input("table-subset-target", "filter_query"),
        Input("hide-reference-targets", "value"),
    ],
    [
        State("url", "search"),
        State("url", "hash"),
    ],
)
def update_table_subset(
    page_current,
    page_size,
    sort_by,
    filter,
    hide_reference,
    # filter_criterion,
    search,
    hash_guide,
):
    """
    La funzione ritorna uno split dei risultati in base ad un filtering o a un sort da parte dell'utente. Inoltre aggiorna i risultati
    visualizzati quando il bottone next page / prev page è cliccato. (Codice preso dalla pagina dash datatable sul sorting con python)
    Inoltre carica i file targets, o scores se presente, e lo trasforma in un dataframe, cambiando il nome delle colonne per farle corrispondere
    all'id delle colonne della tabella nella pagina.
    Se non ci sono targets ritorna un avviso di errore
    """

    job_id = search.split("=")[-1]
    filter_criterion = read_json(job_id)
    job_directory = current_working_directory + "Results/" + job_id + "/"
    with open(current_working_directory + "Results/" + job_id + "/.Params.txt") as p:
        all_params = p.read()
        genome_type_f = (
            next(s for s in all_params.split("\n") if "Genome_selected" in s)
        ).split("\t")[-1]
        ref_comp = (next(s for s in all_params.split("\n") if "Ref_comp" in s)).split(
            "\t"
        )[-1]

    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
    value = job_id
    if search is None:
        raise PreventUpdate
    if not (filter is None):
        filtering_expressions = filter.split(" && ")
    # filtering_expressions.append(['{crRNA} = ' + guide])
    guide = hash_guide[1: hash_guide.find("new")]
    mms = hash_guide[-1:]
    bulge_s = hash_guide[-2:-1]
    if "DNA" in hash_guide:
        bulge_t = "DNA"
    elif "RNA" in hash_guide:
        bulge_t = "RNA"
    else:
        bulge_t = "X"

    if "hide-ref" in hide_reference or genome_type == "var":
        result = global_store_subset_no_ref(
            value, bulge_t, bulge_s, mms, guide, page_current, job_id
        )
    else:
        result = global_store_subset(
            value, bulge_t, bulge_s, mms, guide, page_current, job_id
        )
    drop_cols = drop_columns(result, filter_criterion)
    result = result.drop(drop_cols, axis=1)

    # name of target file filtered with bul-type, mm and bul
    targets_with_mm_bul = (
        current_working_directory
        + "Results/"
        + job_id
        + "/"
        + job_id
        + "."
        + str(bulge_t)
        + "."
        + str(mms)
        + "."
        + str(bulge_s)
        + "."
        + guide
        + ".targets.tsv"
    )
    # save df to tsv with filtered data
    result.to_csv(targets_with_mm_bul, sep='\t', na_rep='NA', index=False)
    # change name to zip file
    targets_with_mm_bul_zip = targets_with_mm_bul.replace("tsv", "zip")
    # zip operation, non blocking
    os.system(f"zip -j {targets_with_mm_bul_zip} {targets_with_mm_bul} &")

    columns_result = [
        {"name": i, "id": i, "hideable": True} for col, i in enumerate(result.columns)
    ]
    data_to_send = result.to_dict("records")
    return [data_to_send, columns_result]


def guidePagev3(job_id, hash):
    guide = hash[: hash.find("new")]
    mms = hash[-1:]
    bulge_s = hash[-2:-1]
    if "DNA" in hash:
        bulge_t = "DNA"
    elif "RNA" in hash:
        bulge_t = "RNA"
    else:
        bulge_t = "X"
    add_header = " - Mismatches " + str(mms)
    if bulge_t != "X":
        add_header += " - " + str(bulge_t) + " " + str(bulge_s)
    value = job_id
    if not isdir(current_working_directory + "Results/" + job_id):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))
    with open(current_working_directory + "Results/" + value + "/.Params.txt") as p:
        all_params = p.read()
        genome_type_f = (
            next(s for s in all_params.split("\n") if "Genome_selected" in s)
        ).split("\t")[-1]
        ref_comp = (next(s for s in all_params.split("\n") if "Ref_comp" in s)).split(
            "\t"
        )[-1]
        pam = (next(s for s in all_params.split(
            "\n") if "Pam" in s)).split("\t")[-1]

    job_directory = current_working_directory + "Results/" + job_id + "/"
    genome_type = "ref"
    style_hide_reference = {"display": "none"}
    value_hide_reference = []
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
        style_hide_reference = {}
        value_hide_reference = ["hide-ref"]

    pam_at_start = False
    if str(guide)[0] == "N":
        pam_at_start = True

    final_list = []
    if pam_at_start:
        final_list.append(
            html.H3(
                "Selected Guide: " +
                str(pam) + str(guide).replace("N", "") + add_header
            )
        )
    else:
        final_list.append(
            html.H3(
                "Selected Guide: " +
                str(guide).replace("N", "") + str(pam) + add_header
            )
        )
    final_list.append(
        html.P(
            [
                # 'Select a row to view the target IUPAC character scomposition. The rows highlighted in red indicates that the target was found only in the genome with variants.',
                "List of Targets found for the selected guide.",
                dcc.Checklist(
                    options=[
                        {"label": "Hide Reference Targets", "value": "hide-ref"}],
                    id="hide-reference-targets",
                    value=value_hide_reference,
                    style=style_hide_reference,
                ),
                html.Div(
                    [
                        html.P(
                            "Generating download link, Please wait...",
                            id="download-link-sumbyguide",
                        ),
                        dcc.Interval(interval=5 * 1000,
                                     id="interval-sumbyguide"),
                    ]
                ),
            ]
        )
    )
    integrated_file_name = glob(
        current_working_directory + "Results/" + job_id + "/" + "*integrated*"
    )[0]
    integrated_file_name = str(integrated_file_name)
    file_to_grep = job_directory + job_id + ".bestMerge.txt.integrated_results.tsv"
    # file_to_grep_alt = job_directory + job_id + '.altMerge.txt'

    guide_grep_result = (
        job_directory
        + job_id
        + "."
        + bulge_t
        + "."
        + bulge_s
        + "."
        + mms
        + "."
        + guide
        + ".txt"
    )
    # put_header = 'head -1 ' + job_directory + job_id + file_to_grep + ' > ' + guide_grep_result + ' ; '

    final_list.append(
        html.Div(
            job_id
            + "."
            + str(bulge_t)
            + "."
            + str(mms)
            + "."
            + str(bulge_s)
            + "."
            + guide
            + ".targets",
            style={"display": "none"},
            id="div-info-sumbyguide-targets",
        )
    )

    path_db = glob(current_working_directory +
                   "Results/" + job_id + "/.*.db")[0]
    path_db = str(path_db)
    conn = sqlite3.connect(path_db)
    c = conn.cursor()
    total_private_sample = f"SELECT * FROM final_table LIMIT 1"
    rows = c.execute(total_private_sample)
    header = [description[0] for description in rows.description]
    conn.commit()
    conn.close()

    cols = [{"name": i, "id": i, "hideable": True} for i in header]
    final_list.append(
        html.Div(
            dash_table.DataTable(
                id="table-subset-target",
                # columns=cols,
                style_cell={"textAlign": "left"},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action="custom",
                style_table={"max-height": "600px", "overflowX": "scroll"},
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
                css=[
                    {"selector": ".row", "rule": "margin: 0"},
                    {
                        "selector": "td.cell--selected, td.focused",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                    {
                        "selector": "td.cell--selected *, td.focused *",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                ],
                # style_data_conditional=[ ]
            ),
            id="div-result-table",
        )
    )
    final_list.append(html.Br())

    return html.Div(final_list, style={"margin": "1%"})


# @cache.memoize()
def global_store_subset_no_ref(value, bulge_t, bulge_s, mms, guide, page, job_id):
    """
    Caching dei file targets per una miglior performance di visualizzazione
    """
    if value is None:
        return ""
    path_db = glob(current_working_directory +
                   "Results/" + value + "/.*.db")[0]
    path_db = str(path_db)
    conn = sqlite3.connect(path_db)
    c = conn.cursor()

    filter_criterion = read_json(job_id)
    query_cols = get_query_column(filter_criterion)

    result = pd.read_sql_query(
        'SELECT * FROM final_table WHERE "{}"=\'{}\' AND "{}"=\'{}\' AND "{}"={} AND "{}"={} AND "{}"<>\'NA\' LIMIT {} OFFSET {}'.format(
            GUIDE_COLUMN,
            guide,
            query_cols['bul_type'],
            bulge_t,
            query_cols['bul'],
            bulge_s,
            query_cols['mm'],
            mms,
            query_cols['samples'],
            PAGE_SIZE,
            page * PAGE_SIZE,
        ),
        conn,
    )

    return result


# @cache.memoize()
def global_store_subset(value, bulge_t, bulge_s, mms, guide, page, job_id):
    """
    Caching dei file targets per una miglior performance di visualizzazione
    """
    if value is None:
        return ""
    path_db = glob(current_working_directory +
                   "Results/" + value + "/.*.db")[0]
    path_db = str(path_db)
    conn = sqlite3.connect(path_db)
    c = conn.cursor()

    filter_criterion = read_json(job_id)
    query_cols = get_query_column(filter_criterion)

    result = pd.read_sql_query(
        'SELECT * FROM final_table WHERE "{}"=\'{}\' AND "{}"=\'{}\' AND "{}"={} AND "{}"={} LIMIT {} OFFSET {}'.format(
            GUIDE_COLUMN,
            guide,
            query_cols['bul_type'],
            bulge_t,
            query_cols['bul'],
            bulge_s,
            query_cols['mm'],
            mms,
            PAGE_SIZE,
            page * PAGE_SIZE,
        ),
        conn,
    )

    targets_with_mm_bul = (
        current_working_directory
        + "Results/"
        + job_id
        + "/"
        + job_id
        + "."
        + str(bulge_t)
        + "."
        + str(mms)
        + "."
        + str(bulge_s)
        + "."
        + guide
        + ".targets.tsv"
    )
    result.to_csv(targets_with_mm_bul, sep='\t', na_rep='NA', index=False)

    return result


# Load barplot of population distribution for selected guide


@app.callback(
    Output("content-collapse-population", "children"),
    [Input("general-profile-table", "selected_cells")],
    [State("general-profile-table", "data"), State("url", "search")],
)
def loadDistributionPopulations(sel_cel, all_guides, job_id):
    if sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    job_id = job_id.split("=")[-1]

    with open(current_working_directory + "Results/" + job_id + "/.Params.txt") as p:
        all_params = p.read()
        mms = int(
            (next(s for s in all_params.split("\n") if "Mismatches" in s)).split("\t")[
                -1
            ]
        )
        max_bulges = int(
            (next(s for s in all_params.split("\n") if "Max_bulges" in s)).split("\t")[
                -1
            ]
        )

    distributions = [
        dbc.Row(
            html.P(
                "On- and Off-Targets distributions in the Reference and Variant Genome. For the Variant Genome, the targets are divided into SuperPopulations.",
                style={"margin-left": "0.75rem"},
            )
        )
    ]

    for i in range(math.ceil((mms + max_bulges + 1) / BARPLOT_LEN)):
        all_images = []
        for mm in range(i * BARPLOT_LEN, (i + 1) * BARPLOT_LEN):
            if mm < (mms + max_bulges + 1):
                try:
                    all_images.append(
                        dbc.Col(
                            [
                                html.A(
                                    html.Img(
                                        src="data:image/png;base64,{}".format(
                                            base64.b64encode(
                                                open(
                                                    current_working_directory
                                                    + "Results/"
                                                    + job_id
                                                    + "/imgs/populations_distribution_"
                                                    + guide
                                                    + "_"
                                                    + str(mm)
                                                    + "total.png",
                                                    "rb",
                                                ).read()
                                            ).decode()
                                        ),
                                        id="distribution-population" + str(mm),
                                        width="100%",
                                        height="auto",
                                    ),
                                    target="_blank",
                                    href="/Results/"
                                    + job_id
                                    + "/imgs/"
                                    + "populations_distribution_"
                                    + guide
                                    + "_"
                                    + str(mm)
                                    + "total.png",
                                ),
                                html.Div(
                                    html.P(
                                        "Distribution "
                                        + str(mm)
                                        + " Mismatches + Bulges ",
                                        style={"display": "inline-block"},
                                    ),
                                    style={"text-align": "center"},
                                ),
                            ]
                        )
                    )
                except:
                    all_images.append(
                        dbc.Col(
                            [
                                html.Div(
                                    html.P(
                                        "No Targets found with "
                                        + str(mm)
                                        + " Mismatches + Bulges",
                                        style={"display": "inline-block"},
                                    ),
                                    style={"text-align": "center"},
                                ),
                                # html.Div(html.P('Distribution ' + str(mm) + ' Mismatches + Bulges ', style = {'display':'inline-block'} ),style = {'text-align':'center'})
                            ],
                            align="center",
                        )
                    )
            else:
                all_images.append(dbc.Col(html.P("")))

        distributions.append(html.Div([dbc.Row(all_images)]))
    return distributions


# Open/close barplot for population distribution
@app.callback(
    Output("collapse-populations", "is_open"),
    [Input("btn-collapse-populations", "n_clicks")],
    [State("collapse-populations", "is_open")],
)
def toggleCollapseDistributionPopulations(n, is_open):
    if n:
        return not is_open
    return is_open


# Filtering e sorting per la pagina principale delle guide


@app.callback(
    [
        Output("general-profile-table", "data"),
        Output("general-profile-table", "selected_cells"),
    ],
    [
        Input("general-profile-table", "page_current"),
        Input("general-profile-table", "page_size"),
        Input("general-profile-table", "sort_by"),
        Input("general-profile-table", "filter_query"),
        Input("target_filter_dropdown", "value"),
    ],
    [State("url", "search")],
)
def update_table_general_profile(
    page_current, page_size, sort_by, filter, filter_criterion, search
):
    job_id = search.split("=")[-1]

    with open(current_working_directory + "Results/" + job_id + "/.Params.txt") as p:
        all_params = p.read()
        genome_type_f = (
            next(s for s in all_params.split("\n") if "Genome_selected" in s)
        ).split("\t")[-1]
        ref_comp = (next(s for s in all_params.split("\n") if "Ref_comp" in s)).split(
            "\t"
        )[-1]
        mms = int(
            (next(s for s in all_params.split("\n") if "Mismatches" in s)).split("\t")[
                -1
            ]
        )
        max_bulges = int(
            (next(s for s in all_params.split("\n") if "Max_bulges" in s)).split("\t")[
                -1
            ]
        )
        nuclease = (next(s for s in all_params.split("\n") if "Nuclease" in s)).split(
            "\t"
        )[-1]

    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"

    filtering_expressions = filter.split(" && ")

    # Get error guides
    list_error_guides = []
    if os.path.exists(
        current_working_directory + "Results/" + job_id + "/guides_error.txt"
    ):
        with open(
            current_working_directory + "Results/" + job_id + "/guides_error.txt"
        ) as error_g:
            for e_g in error_g:
                list_error_guides.append(e_g.strip())

    # Get guide from guide.txt
    with open(current_working_directory + "Results/" + job_id + "/.guides.txt") as g:
        guides = g.read().strip().split("\n")
        guides.sort()

    # MT add
    acfd_file = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        "".join([".", job_id, ".acfd_", filter_criterion, ".txt"]),
    )
    if not os.path.isfile(acfd_file):
        raise FileNotFoundError(f"Unable to locate {acfd_file}")

    # load acfd for each guide
    with open(acfd_file) as a:
        all_scores = a.read().strip().split("\n")

    # Load scores
    if "NO SCORES" not in all_scores:
        all_scores.sort()
        acfd = [
            float(a.split("\t")[1])
            for a in all_scores
            if a.split("\t")[0] not in list_error_guides
        ]
        doench = [
            a.split("\t")[2]
            for a in all_scores
            if a.split("\t")[0] not in list_error_guides
        ]
        if genome_type == "both":
            doench_enr = [
                a.split("\t")[3]
                for a in all_scores
                if a.split("\t")[0] not in list_error_guides
            ]
        # acfd = [int(round((100/(100 + x))*100)) for x in acfd]
        acfd = [
            float("{:.3f}".format(x * 100))
            if x < 1 and x >= 0
            else "CFD score not available"
            for x in acfd
        ]

    # Get target counting from summary by guide
    column_on_target = []
    column_off_target_ref = []
    column_sample_class = []
    column_total = []

    df = []
    table_to_file = list()
    for x, g in enumerate(guides):
        table_to_file.append(g)  # append guide to table
        # append nuclease to table
        table_to_file.append("Nuclease: " + str(nuclease))
        data_general_count = pd.read_csv(
            current_working_directory
            + "Results/"
            + job_id
            + "/."
            + job_id
            + ".general_target_count."
            + g
            + "_"
            + filter_criterion
            + ".txt",
            sep="\t",
            na_filter=False,
        )

        data_guides = dict()
        data_guides["Guide"] = g
        data_guides["Nuclease"] = nuclease
        data_general_count_copy = data_general_count.copy()
        count_bulges = list()
        origin_ref = list()
        origin_var = list()
        for the_bulge in range(max_bulges + 1):
            origin_ref.append("REF")
            origin_var.append("VAR")
            count_bulges.append(the_bulge)

        count_bulges_concat = count_bulges + count_bulges
        origin_concat = origin_ref + origin_var

        data_general_count_copy.insert(0, "Genome", origin_concat, True)
        data_general_count_copy.insert(1, "Bulges", count_bulges_concat, True)

        if "NO SCORES" not in all_scores:
            data_guides["CFD"] = acfd[x]
            table_to_file.append("CFD: " + str(acfd[x]))  # append CFD to table
            table_to_file.append("\t\t\t\tMismatches")

            # table_to_file.append('IN THE FOLLOWING MATRIX, THE FIRST GROUP OF '+str(max_bulges)+' LINES, ARE REFERED TO REFERENCE TARGET, THE SECOND GROUP OF '+str(max_bulges)+' LINES ARE REFERED TO VARIANT GENOME')

            table_to_file.append(
                data_general_count_copy.to_string(index=False))

            if genome_type == "both":
                data_guides["Doench 2016"] = doench[x]
                # data_guides['Enriched'] = doench_enr[x]
            else:
                data_guides["Doench 2016"] = doench[x]

        if genome_type == "both":
            tmp = [str(i) for i in range(max_bulges + 1)] * 2
            tmp.insert(len(tmp) // 2, "")
            data_guides["# Bulges"] = "\n".join(tmp)
        else:
            tmp = [str(i) for i in range(max_bulges + 1)]
            data_guides["# Bulges"] = "\n".join(tmp)

        data_guides["Total"] = []
        if genome_type == "both":
            if max_bulges == 2:
                for i in range(len(data_guides["# Bulges"].split("\n")) - 1):
                    if i == 1:
                        data_guides["Total"].append(
                            "REFERENCE\t" +
                            str(sum(data_general_count.iloc[i, :]))
                        )
                    elif i == 2:
                        data_guides["Total"].append(
                            "\t" + str(sum(data_general_count.iloc[i, :]))
                        )
                        data_guides["Total"].append("\t")
                    elif i == 4:
                        data_guides["Total"].append(
                            "VARIANT\t\t" +
                            str(sum(data_general_count.iloc[i, :]))
                        )
                    else:
                        data_guides["Total"].append(
                            "\t" + str(sum(data_general_count.iloc[i, :]))
                        )
            elif max_bulges == 1:
                for i in range(len(data_guides["# Bulges"].split("\n")) - 1):
                    if i == 1:
                        data_guides["Total"].append(
                            "REFERENCE\t" +
                            str(sum(data_general_count.iloc[i, :]))
                        )
                        data_guides["Total"].append("\t")
                    elif i == 3:
                        data_guides["Total"].append(
                            "VARIANT\t\t" +
                            str(sum(data_general_count.iloc[i, :]))
                        )
                    else:
                        data_guides["Total"].append(
                            "\t" + str(sum(data_general_count.iloc[i, :]))
                        )
            else:
                for i in range(len(data_guides["# Bulges"].split("\n")) - 1):
                    if i == 0:
                        data_guides["Total"].append(
                            "REFERENCE\t" +
                            str(sum(data_general_count.iloc[i, :]))
                        )
                        data_guides["Total"].append("\t")
                    elif i == 1:
                        data_guides["Total"].append(
                            "VARIANT\t\t" +
                            str(sum(data_general_count.iloc[i, :]))
                        )
        else:
            for i in range(len(data_guides["# Bulges"].split("\n"))):
                if i == len(data_guides["# Bulges"].split("\n")) // 2:
                    data_guides["Total"].append(
                        "REFERENCE\t" + str(sum(data_general_count.iloc[i, :]))
                    )
                else:
                    data_guides["Total"].append(
                        "\t" + str(sum(data_general_count.iloc[i, :]))
                    )

        if genome_type == "both":
            for i in range(mms + 1):
                tmp = list(data_general_count.iloc[:, i].values.astype(str))
                tmp.insert(len(tmp) // 2, "")
                data_guides[str(i) + "MM"] = "\n".join(tmp)
        else:
            for i in range(mms + 1):
                tmp = list(
                    data_general_count.iloc[: max_bulges +
                                            1, i].values.astype(str)
                )
                # tmp.insert(len(tmp)//2, "")
                data_guides[str(i) + "MM"] = "\n".join(tmp)

        data_guides["Total"] = "\n".join(data_guides["Total"])

        df.append(data_guides)
    dff = pd.DataFrame(df)

    table_to_file_save_dest = (
        current_working_directory
        + "Results/"
        + job_id
        + "/"
        + job_id
        + ".general_table.txt"
    )

    outfile = open(table_to_file_save_dest, "w")
    for elem in table_to_file:
        outfile.write(elem + "\n")
    outfile.close()

    # zip integrated results
    integrated_file_name = glob(
        current_working_directory + "Results/" + job_id + "/" + "*integrated*"
    )[0]
    integrated_file_name = str(integrated_file_name)
    integrated_file = integrated_file_name

    integrated_to_zip = integrated_file_name.replace("tsv", "zip")
    if not os.path.exists(integrated_to_zip):
        os.system(f"zip -j {integrated_to_zip} {integrated_file} &")

    if "NO SCORES" not in all_scores:
        try:
            dff = dff.sort_values(["CFD", "Doench 2016"],
                                  ascending=[False, False])
        except:  # for BOTH
            dff = dff.sort_values(["CFD", "Enriched"],
                                  ascending=[False, False])
    else:
        try:
            dff = dff.sort_values("On-Targets Reference", ascending=True)
        except:
            dff = dff.sort_values("On-Targets Enriched", ascending=True)

    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

        if operator in ("eq", "ne", "lt", "le", "gt", "ge"):
            # these operators match pandas series operator method names
            dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
        elif operator == "contains":
            dff = dff.loc[dff[col_name].str.contains(filter_value)]
        elif operator == "datestartswith":
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            dff = dff.loc[dff[col_name].str.startswith(filter_value)]

    if len(sort_by):
        dff = dff.sort_values(
            [
                "Samples" if col["column_id"] == "Samples Summary" else col["column_id"]
                for col in sort_by
            ],
            ascending=[col["direction"] == "asc" for col in sort_by],
            inplace=False,
        )

    # Calculate sample count

    data_to_send = dff.iloc[
        page_current * page_size: (page_current + 1) * page_size
    ].to_dict("records")
    return data_to_send, [{"row": 0, "column": 0}]


# Update color on selected row


@app.callback(
    Output("general-profile-table", "style_data_conditional"),
    [Input("general-profile-table", "selected_cells")],
    [State("general-profile-table", "data")],
)
def colorSelectedRow(sel_cel, all_guides):
    if sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    return [
        {
            "if": {
                "filter_query": '{Guide} eq "' + guide + '"',
            },
            "background-color": "rgba(0, 0, 255,0.15)",  # 'rgb(255, 102, 102)'
        },
        {"if": {"column_id": "Genome"}, "font-weight": "bold", "textAlign": "center"},
    ]


# ------------------------------------------------------------------------------
# Query genomic region tab
#
@app.callback(
    [
        Output("div-table-position", "children"),
        Output("div-current-page-table-position", "children"),
    ],
    [Input("div-position-filter-query", "children")],
    [
        State("button-filter-position", "n_clicks_timestamp"),
        State("target_filter_dropdown", "value"),
        State("url", "search"),
        State("general-profile-table", "selected_cells"),
        State("general-profile-table", "data"),
        State("div-current-page-table-position", "children"),
        # State("div-mms-bulges-position", "children"),
    ],
)
def filterPositionTable(
    filter_q: List[str],
    n: int,
    filter_criterion: str,
    search: str,
    sel_cel: List[int],
    all_guides: List[int],
    current_page: str,
    # mms_bulge: str,
) -> Tuple[List[html.P], str]:
    """Filter result table by genomic region. The table is filtered in order to
    display only those targets falling within the genomic interval defined
    by the user.

    The results can be furtherly filtered by scoring criterion. The available
    criteria are CFD score, CRISTA score and the number of mismatches and bulges.

    ...

    Parameters
    ----------
    filter_q : List[str]
        Filtering query (coordinates, filtering criterion)
    n : int
        Input click listened
    filter_criterion : str
        result table filtering criterion
    search : str
        Search ID
    sel_cel : List[int]
    all_guides : List[int]
        List of the guides
    current_page : str
        Current table page number

    Returns
    -------
    List[html.P]
        HTML result page
    str
        Page numeration

    """
    if n is not None:
        if not isinstance(n, int):
            raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if not filter_criterion in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filtering criterion {filter_criterion}")
    if not isinstance(search, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(search).__name__}")
    if not isinstance(sel_cel, list):
        raise TypeError(
            f"Expected {list.__name__}, got {type(sel_cel).__name__}")
    if not bool(sel_cel):
        raise ValueError("Empty Dictionary, stopping execution.")
    if not isinstance(all_guides, list):
        raise TypeError(
            f"Expected {list.__name__}, got {type(all_guides).__name__}")
    if not bool(all_guides):
        raise ValueError("Empty Dictionary, stopping execution.")

    if sel_cel is None:
        raise PreventUpdate
    if n is None:
        raise PreventUpdate
    # recover filter query fields; if there are NULL fields -> prevent table update
    # query structure: chrom,start,stop
    if isinstance(filter_q, str):  # simple regular query
        filter_q = filter_q.split(",")
    elif isinstance(filter_q, list):  # updated by callback
        assert len(filter_q) == 2  # we should have just two elements
        filter_criterion = filter_q[1]  # recover table filtering criterion
        filter_q = filter_q[0].split(",")  # query genomic coordinates
    print(filter_criterion)
    assert filter_criterion in FILTERING_CRITERIA
    chrom = filter_q[0]
    if chrom == "None":
        raise PreventUpdate
    start = filter_q[1]
    if start == "None":
        raise PreventUpdate
    end = filter_q[2]
    if end == "None":
        raise PreventUpdate
    current_page = int(current_page.split("/")[0])
    job_id = search.split("=")[-1]
    job_directory = os.path.join(
        current_working_directory, RESULTS_DIR, job_id)
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    # recover db file
    db_path = glob(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, ".*.db")
    )[0]
    assert isinstance(db_path, str)
    assert os.path.isfile(db_path)
    # connect db with sqlite3
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    query = 'SELECT * FROM final_table WHERE "{}"=\'{}\' AND "{}">={} AND "{}"<={} AND "{}"=\'{}\''

    filter_criterion = read_json(job_id)
    query_cols = get_query_column(filter_criterion)

    result = pd.read_sql_query(
        query.format(
            GUIDE_COLUMN, guide, query_cols['start'], start, query_cols['start'], end, CHR_COLUMN, chrom
        ),
        conn,
    )
    conn.commit()
    conn.close()
    if result.shape[0] == 0:  # no guides found
        df_check = False
    else:  # check table fit to page
        df_check = True
    # filter diplayed column using filtering criterion
    drop_cols = drop_columns(result, filter_criterion)
    result.drop(drop_cols, inplace=True, axis=1)  # remove columns from table
    # check table characteristics to fit it into html page
    if df_check:
        out_1 = [
            dash_table.DataTable(
                css=[{"selector": ".row", "rule": "margin: 0"}],
                id="table-position",
                export_format="xlsx",
                columns=[
                    {"name": i, "id": i, "hideable": True}
                    for count, i in enumerate(result.columns)
                ],
                data=result.to_dict("records"),
                style_cell={"textAlign": "left"},
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
                style_table={
                    "overflowX": "scroll",
                },
                page_size=PAGE_SIZE,
            )
        ]
    else:
        out_1 = [html.P("No results found with this genomic coordinates")]

    return out_1, "/".join([str(1) + str(1)])


@app.callback(
    Output("div-position-filter-query", "children"),
    [Input("button-filter-position", "n_clicks")],
    [
        State("target_filter_dropdown", "value"),
        State("dropdown-chr-table-position", "value"),
        State("input-position-start", "value"),
        State("input-position-end", "value"),
    ],
)
def updatePositionFilter(
    n: int, filter_criterion: str, chrom: str, pos_start: str, pos_end: str
) -> Tuple[str, str, int]:
    """Callback to update the result table filtered by genomic location.

    ...

    Parameters
    ----------
    n : int
        Number of clicks listened
    filter_criterion : str
        Filtering citerion to apply to the table
    chrom : str
        Chromosome
    pos_start : str
        Start position
    pos_end : str
        Stop position

    Returns
    -------
    Tuple[str, str, int]
        New genomic locations and potential filtering criterion
    """

    if n is not None:
        if not isinstance(n, int):
            raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if not isinstance(chrom, str):
        raise TypeError(f"Expected {str.__name__}, got {type(chrom).__name__}")
    if not isinstance(pos_start, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(pos_start).__name__}")
    if not isinstance(pos_end, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(pos_end).__name__}")

    if n is None:  # no click -> no page update
        raise PreventUpdate
    if pos_start == "":
        pos_start = "None"
    if pos_end == "":
        pos_end = "None"
    coords = ",".join([chrom, pos_start, pos_end])
    return coords, filter_criterion


# Callback to view next/prev page on sample table


@app.callback(
    [
        Output("div-table-samples", "children"),
        Output("div-current-page-table-samples", "children"),
    ],
    [
        Input("prev-page-sample", "n_clicks_timestamp"),
        Input("next-page-sample", "n_clicks_timestamp"),
        Input("div-sample-filter-query", "children"),
        Input("target_filter_dropdown", "value"),
    ],
    [
        State("button-filter-population-sample", "n_clicks_timestamp"),
        State("url", "search"),
        State("general-profile-table", "selected_cells"),
        State("general-profile-table", "data"),
        State("div-current-page-table-samples", "children"),
    ],
)
def filterSampleTable(
    nPrev,
    nNext,
    filter_q,
    filter_criterion,
    n,
    search,
    sel_cel,
    all_guides,
    current_page,
):
    if sel_cel is None:
        raise PreventUpdate
    if nPrev is None and nNext is None and n is None:
        raise PreventUpdate

    if nPrev is None:
        nPrev = 0
    if nNext is None:
        nNext = 0
    if n is None:
        n = 0

    sup_pop = filter_q.split(",")[0]
    pop = filter_q.split(",")[1]
    samp = str(filter_q.split(",")[2])
    if sup_pop == "None":
        sup_pop = None
    if pop == "None":
        pop = None
    if samp == "None" or samp == "NONE":
        samp = None
    current_page = current_page.split("/")[0]
    current_page = int(current_page)
    btn_sample_section = []
    btn_sample_section.append(n)
    btn_sample_section.append(nPrev)
    btn_sample_section.append(nNext)
    job_id = search.split("=")[-1]
    job_directory = current_working_directory + "Results/" + job_id + "/"
    population_1000gp = associateSample.loadSampleAssociation(
        job_directory + ".sampleID.txt"
    )[2]
    with open(current_working_directory + "Results/" + job_id + "/.Params.txt") as p:
        all_params = p.read()
        genome_type_f = (
            next(s for s in all_params.split("\n") if "Genome_selected" in s)
        ).split("\t")[-1]
        ref_comp = (next(s for s in all_params.split("\n") if "Ref_comp" in s)).split(
            "\t"
        )[-1]

    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"

    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    if genome_type == "both":
        col_names_sample = [
            "Sample",
            "Sex",
            "Population",
            "Super Population",  # 'Targets in Reference',
            "Targets in Variant",
            "Targets in Population",
            "Targets in Super Population",
            "PAM Creation",
        ]  # , 'Class']
    else:
        col_names_sample = [
            "Sample",
            "Sex",
            "Population",
            "Super Population",  # 'Targets in Reference',
            "Targets in Variant",
            "Targets in Population",
            "Targets in Super Population",
            "PAM Creation",
        ]  # , 'Class']
    # Last button pressed is filtering, return the first page of the filtered table
    if max(btn_sample_section) == n:
        if genome_type == "both":
            df = pd.read_csv(
                job_directory
                + job_id
                + ".summary_by_samples."
                + guide
                + "_"
                + filter_criterion
                + ".txt",
                sep="\t",
                names=col_names_sample,
                skiprows=2,
                na_filter=False,
            )
            df = df.sort_values("Targets in Variant", ascending=False)
        else:
            df = pd.read_csv(
                job_directory
                + job_id
                + ".summary_by_samples."
                + guide
                + "_"
                + filter_criterion
                + ".txt",
                sep="\t",
                names=col_names_sample,
                skiprows=2,
                na_filter=False,
            )
            df = df.sort_values("Targets in Variant", ascending=False)

        more_info_col = []
        for i in range(df.shape[0]):
            more_info_col.append("Show Targets")
        df[""] = more_info_col
        if (
            (sup_pop is None or sup_pop == "")
            and (pop is None or pop == "")
            and (samp is None or samp == "")
        ):  # No filter value selected
            max_page = len(df.index)
            max_page = math.floor(max_page / 10) + 1
            return (
                generate_table_samples(df, "table-samples", 1, guide, job_id),
                "1/" + str(max_page),
            )
        if samp is None or samp == "":
            if pop is None or pop == "":
                df.drop(
                    df[(~(df["Population"].isin(population_1000gp[sup_pop])))].index,
                    inplace=True,
                )
            else:
                df.drop(df[(df["Population"] != pop)].index, inplace=True)
        else:
            df.drop(df[(df["Sample"] != samp)].index, inplace=True)
        max_page = len(df.index)
        max_page = math.floor(max_page / 10) + 1
        return (
            generate_table_samples(df, "table-samples", 1, guide, job_id),
            "1/" + str(max_page),
        )
    else:
        if max(btn_sample_section) == nNext:
            current_page = current_page + 1
            if genome_type == "both":
                df = pd.read_csv(
                    job_directory + job_id + ".summary_by_samples." + guide + ".txt",
                    sep="\t",
                    names=col_names_sample,
                    skiprows=2,
                    na_filter=False,
                )
                df = df.sort_values("Targets in Variant", ascending=False)
            else:
                df = pd.read_csv(
                    job_directory + job_id + ".summary_by_samples." + guide + ".txt",
                    sep="\t",
                    names=col_names_sample,
                    skiprows=2,
                    na_filter=False,
                )
                df = df.sort_values("Targets in Reference", ascending=False)

            more_info_col = []
            for i in range(df.shape[0]):
                more_info_col.append("Show Targets")
            df[""] = more_info_col
            # Active filter
            if pop or sup_pop or samp:
                if samp is None or samp == "":
                    if pop is None or pop == "":
                        df.drop(
                            df[
                                (~(df["Population"].isin(
                                    population_1000gp[sup_pop])))
                            ].index,
                            inplace=True,
                        )
                    else:
                        df.drop(df[(df["Population"] != pop)].index,
                                inplace=True)
                else:
                    df.drop(df[(df["Sample"] != samp)].index, inplace=True)

            if ((current_page - 1) * 10) > len(df):
                current_page = current_page - 1
                if current_page < 1:
                    current_page = 1
            max_page = len(df.index)
            max_page = math.floor(max_page / 10) + 1
            return (
                generate_table_samples(
                    df, "table-samples", current_page, guide, job_id
                ),
                str(current_page) + "/" + str(max_page),
            )

        else:  # Go to previous page
            current_page = current_page - 1
            if current_page < 1:
                current_page = 1
            if genome_type == "both":
                df = pd.read_csv(
                    job_directory + job_id + ".summary_by_samples." + guide + ".txt",
                    sep="\t",
                    names=col_names_sample,
                    skiprows=2,
                    na_filter=False,
                )
                df = df.sort_values("Targets in Variant", ascending=False)
            else:
                df = pd.read_csv(
                    job_directory + job_id + ".summary_by_samples." + guide + ".txt",
                    sep="\t",
                    names=col_names_sample,
                    skiprows=2,
                    na_filter=False,
                )
                df = df.sort_values("Targets in Variant", ascending=False)

            more_info_col = []
            for i in range(df.shape[0]):
                more_info_col.append("Show Targets")
            df[""] = more_info_col
            if pop or sup_pop or samp:
                if samp is None or samp == "":
                    if pop is None or pop == "":
                        df.drop(
                            df[
                                (~(df["Population"].isin(
                                    population_1000gp[sup_pop])))
                            ].index,
                            inplace=True,
                        )
                    else:
                        df.drop(df[(df["Population"] != pop)].index,
                                inplace=True)
                else:
                    df.drop(df[(df["Sample"] != samp)].index, inplace=True)
            max_page = len(df.index)
            max_page = math.floor(max_page / 10) + 1
            return (
                generate_table_samples(
                    df, "table-samples", current_page, guide, job_id
                ),
                str(current_page) + "/" + str(max_page),
            )
    raise PreventUpdate


# Callback to update the hidden div filter


@app.callback(
    Output("div-sample-filter-query", "children"),
    [Input("button-filter-population-sample", "n_clicks")],
    [
        State("dropdown-superpopulation-sample", "value"),
        State("dropdown-population-sample", "value"),
        State("input-sample", "value"),
    ],
)
def updateSampleFilter(n, superpopulation, population, sample):
    if n is None:
        raise PreventUpdate
    return (
        str(superpopulation)
        + ","
        + str(population)
        + ","
        + str(sample).replace(" ", "").upper()
    )


# Callback to update the sample based on population selected


@app.callback(
    [Output("dropdown-sample", "options"), Output("dropdown-sample", "value")],
    [Input("dropdown-population-sample", "value")],
    [State("url", "search")],
)
def updateSampleDrop(pop, search):
    if pop is None or pop == "":
        return [], None
    job_id = search.split("=")[-1]
    job_directory = current_working_directory + "Results/" + job_id + "/"
    dict_pop = associateSample.loadSampleAssociation(
        job_directory + ".sampleID.txt")[3]
    return [{"label": sam, "value": sam} for sam in dict_pop[pop]], None


# Callback to update the population tab based on superpopulation selected


@app.callback(
    [
        Output("dropdown-population-sample", "options"),
        Output("dropdown-population-sample", "value"),
    ],
    [Input("dropdown-superpopulation-sample", "value")],
    [State("url", "search")],
)
def updatePopulationDrop(superpop, search):
    if superpop is None or superpop == "":
        raise PreventUpdate
    job_id = search.split("=")[-1]
    job_directory = current_working_directory + "Results/" + job_id + "/"
    population_1000gp = associateSample.loadSampleAssociation(
        job_directory + ".sampleID.txt"
    )[2]
    return [{"label": i, "value": i} for i in population_1000gp[superpop]], None


def generate_table_position(
    dataframe, id_table, page, mms, bulges, guide="", job_id="", max_rows=10
):
    rows_remaining = len(dataframe) - (page - 1) * max_rows
    header = [
        html.Tr(
            [
                html.Th(
                    "Chromosome",
                    rowSpan="2",
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th(
                    "Position",
                    rowSpan="2",
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th(
                    "Best Target",
                    rowSpan="2",
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th(
                    "Min Mismatch",
                    rowSpan="2",
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th(
                    "Min Bulge",
                    rowSpan="2",
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th(
                    "Bulge",
                    rowSpan="2",
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th(
                    "Targets in Cluster by Mismatch Value",
                    colSpan=str(mms + 1),
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th("", rowSpan="2"),
            ]
        )
    ]
    mms_header = []
    for mm in range(mms + 1):
        mms_header.append(
            html.Th(
                str(mm) + " MM",
                style={"vertical-align": "middle", "text-align": "center"},
            )
        )
    header.append(html.Tr(mms_header))

    data = []
    for i in range(min(rows_remaining, max_rows)):
        first_cells = [
            html.Td(
                dataframe.iloc[i + (page - 1) * max_rows]["Chromosome"],
                rowSpan=str(bulges + 1),
                style={"vertical-align": "middle", "text-align": "center"},
            ),
            html.Td(
                dataframe.iloc[i + (page - 1) * max_rows]["Position"],
                rowSpan=str(bulges + 1),
                style={"vertical-align": "middle", "text-align": "center"},
            ),
            html.Td(
                dataframe.iloc[i + (page - 1) * max_rows]["Best Target"],
                rowSpan=str(bulges + 1),
                style={"vertical-align": "middle", "text-align": "center"},
            ),
            html.Td(
                dataframe.iloc[i + (page - 1) * max_rows]["Min Mismatch"],
                rowSpan=str(bulges + 1),
                style={"vertical-align": "middle", "text-align": "center"},
            ),
            html.Td(
                dataframe.iloc[i + (page - 1) * max_rows]["Min Bulge"],
                rowSpan=str(bulges + 1),
                style={"vertical-align": "middle", "text-align": "center"},
            ),
            html.Th(
                "0 Bulge",
                style={
                    "vertical-align": "middle",
                    "text-align": "center",
                    "padding-left": "0",
                },
            ),
        ]

        mm_cells = [
            html.Td(
                dataframe.iloc[i + (page - 1) * max_rows][col],
                style={"vertical-align": "middle", "text-align": "center"},
            )
            for col in dataframe.columns[5: 5 + mms + 1]
        ]
        data.append(
            html.Tr(
                first_cells
                + mm_cells
                + [
                    html.Td(
                        html.A(
                            "Show Targets",
                            href="result?job="
                            + job_id
                            + "#"
                            + guide
                            + "-Pos-"
                            + str(
                                dataframe.iloc[i + (page - 1)
                                               * max_rows]["Chromosome"]
                            )
                            + "-"
                            + str(
                                dataframe.iloc[i + (page - 1)
                                               * max_rows]["Position"]
                            ),
                            target="_blank",
                        ),
                        rowSpan=str(bulges + 1),
                        style={"vertical-align": "middle",
                               "text-align": "center"},
                    )
                ]
            )
        )
        for b in range(bulges):
            data.append(
                html.Tr(
                    [
                        html.Th(
                            str(b + 1) + " Bulge",
                            style={"vertical-align": "middle",
                                   "text-align": "center"},
                        )
                    ]
                    + [
                        html.Td(dataframe.iloc[i + (page - 1) * max_rows][col])
                        for col in dataframe.columns[
                            5 + (b + 1) * (mms + 1): 5 + (b + 1) * (mms + 1) + mms + 1
                        ]
                    ]
                )
            )

    return html.Table(header + data, style={"display": "inline-block"}, id=id_table)


def generate_table_samples(dataframe, id_table, page, guide="", job_id="", max_rows=10):
    """
    Per generare una html table. NOTE è diversa da una dash dataTable
    """
    dataframe = dataframe.astype(str)
    rows_remaining = len(dataframe) - (page - 1) * max_rows
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] +
        # Body
        [
            html.Tr(
                [
                    html.Td(
                        html.A(
                            dataframe.iloc[i + (page - 1) * max_rows][col],
                            href="result?job="
                            + job_id
                            + "#"
                            + guide
                            + "-Sample-"
                            + dataframe.iloc[i + \
                                             (page - 1) * max_rows]["Sample"],
                            target="_blank",
                        )
                    )
                    if col == ""
                    else html.Td(dataframe.iloc[i + (page - 1) * max_rows][col])
                    for col in dataframe.columns
                ]
            )
            for i in range(min(rows_remaining, max_rows))
        ],
        style={"display": "inline-block"},
        id=id_table,
    )


def generate_table(
    dataframe, id_table, genome_type, guide="", job_id="", max_rows=2600
):
    """
    Per generare una html table. NOTE è diversa da una dash dataTable
    """
    # if genome_type == 'both':
    header = [
        html.Tr(
            [
                html.Th(
                    "Bulge type",
                    rowSpan="2",
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th(
                    "Mismatches",
                    rowSpan="2",
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th(
                    "Bulge Size",
                    rowSpan="2",
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th(
                    "Targets found in Genome",
                    colSpan=str(3),
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th(
                    "PAM Creation",
                    rowSpan="2",
                    style={"vertical-align": "middle", "text-align": "center"},
                ),
                html.Th("", rowSpan="2"),
            ]
        )
    ]

    header.append(
        html.Tr(
            [
                html.Th(x, style={"vertical-align": "middle",
                        "text-align": "center"})
                for x in ["Reference", "Variant", "Combined"]
            ]
        )
    )
    return html.Table(
        header +
        # Body
        [
            html.Tr(
                [
                    html.Td(
                        html.A(
                            dataframe.iloc[i][col],
                            href="result?job="
                            + job_id
                            + "#"
                            + guide
                            + "new"
                            + dataframe.iloc[i]["Bulge Type"]
                            + str(dataframe.iloc[i]["Bulge Size"])
                            + str(dataframe.iloc[i]["Mismatches"]),
                            target="_blank",
                        ),
                        style={"vertical-align": "middle",
                               "text-align": "center"},
                    )
                    if col == ""
                    else html.Td(
                        dataframe.iloc[i][col],
                        style={"vertical-align": "middle",
                               "text-align": "center"},
                    )
                    for col in dataframe.columns
                ]
            )
            for i in range(min(len(dataframe), max_rows))
        ],
        style={"display": "inline-block"},
        id=id_table,
    )


def check_existance_sample(job_directory, job_id, sample):
    df = pd.read_csv(
        job_directory + job_id + ".sampleID.txt", sep="\t", na_filter=False
    )
    samples = df.iloc[:, 0]
    if sample in samples.values:
        return True
    else:
        return False


# Select figures on mms value, sample value


@app.callback(
    [
        Output("div-radar-chart-encode_gencode", "children"),
        Output("div-population-barplot", "children"),
        Output("div-sample-image", "children"),
        Output("row-radar-chart-sample", "children"),
    ],
    [
        Input("mm-dropdown", "value"),
        Input("general-profile-table", "selected_cells"),
        Input("target_filter_dropdown", "value"),
    ],
    [State("url", "search"), State("general-profile-table", "data")],
)
def updateImagesTabs(mm, sel_cel, filter_criterion, search, all_guides):
    bulge = 0
    job_id = search.split("=")[-1]
    job_directory = current_working_directory + "Results/" + job_id + "/"
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]

    # search for getting job id
    # get guide with sel_cel and all_data
    # radar_chart_images = list()
    radar_chart_encode_gencode = list()
    # radar_chart_gencode = list()
    population_barplots = list()
    guide_images = list()
    sample_images = list()

    try:
        population_barplots.extend(
            [
                html.A(
                    html.Img(
                        src="data:image/png;base64,{}".format(
                            base64.b64encode(
                                open(
                                    current_working_directory
                                    + "Results/"
                                    + job_id
                                    + "/imgs/populations_distribution_"
                                    + guide
                                    + "_"
                                    + str(int(mm) + int(bulge))
                                    + "total"
                                    + "_"
                                    + filter_criterion
                                    + ".png",
                                    "rb",
                                ).read()
                            ).decode()
                        ),
                        id="distribution-population" +
                        str(int(mm) + int(bulge)),
                        width="100%",
                        height="auto",
                    ),
                    target="_blank",
                    href="/Results/"
                    + job_id
                    + "/imgs/"
                    + "populations_distribution_"
                    + guide
                    + "_"
                    + str(int(mm) + int(bulge))
                    + "total"
                    + "_"
                    + filter_criterion
                    + ".png",
                ),
            ]
        )
    except:
        population_barplots = [
            html.Div(
                html.H2(
                    "No result found for this combination of mismatches and bulges")
            )
        ]

    radar_img_encode_gencode = (
        "/imgs/summary_single_guide_"
        + guide
        + "_"
        + str(mm)
        + "."
        + str(bulge)
        + "_TOTAL_"
        + filter_criterion
        + ".ENCODE+GENCODE.png"
    )
    os.system(
        f"python {app_main_directory}/PostProcess/generate_img_radar_chart.py {guide} {job_directory}/.guide_dict_{guide}_{filter_criterion}.json {job_directory}/.motif_dict_{guide}_{filter_criterion}.json {mm} {bulge} TOTAL_{filter_criterion} {job_directory}/imgs/"
    )

    img_found = False
    try:
        radar_src_encode_gencode = "data:image/png;base64,{}".format(
            base64.b64encode(
                open(
                    current_working_directory
                    + "Results/"
                    + job_id
                    + "/"
                    + radar_img_encode_gencode,
                    "rb",
                ).read()
            ).decode()
        )
        img_found = True
    except:
        pass

    try:
        radar_href_encode_gencode = (
            "/Results/" + job_id + "/" + radar_img_encode_gencode
        )
    except:
        radar_href = ""

    if img_found:
        radar_chart_encode_gencode.append(
            html.A(
                html.Img(
                    src=radar_src_encode_gencode,
                    id="radar-img-guide",
                    width="100%",
                    height="auto",
                ),
                target="_blank",
                href=radar_href_encode_gencode,
            )
        )

    if len(radar_chart_encode_gencode) == 0:
        radar_chart_encode_gencode.append(
            html.H2("No result found for this combination of mismatches and bulges")
        )

    # reverse list to print plots in correct order since they are append in reverse order into main sample_images list
    reversed_sample_images = sample_images[::-1]
    return (
        radar_chart_encode_gencode,
        population_barplots,
        guide_images,
        reversed_sample_images,
    )


@app.callback(
    [
        Output("download-link-personal-card", "children"),
        Output("download-link-personal-card", "hidden"),
        Output("div-personal-plot", "children"),
        Output("div-private-plot", "children"),
        Output("div-table-sample-card", "children"),
        Output("div-top-target-sample-card", "children"),
    ],
    [Input("button-sample-card", "n_clicks")],
    [
        State("target_filter_dropdown", "value"),
        State("dropdown-sample-card", "value"),
        State("general-profile-table", "selected_cells"),
        State("general-profile-table", "data"),
        State("url", "search"),
    ],
)
# FUNCTION TO GENERATE SAMPLE CARD, UPDATE WITH FILTER DROPDOWN
def generate_sample_card(n, filter_criterion, sample, sel_cel, all_guides, search):
    if n is None:
        raise PreventUpdate

    # convert sample to str to avoid concatenation errrors
    sample = str(sample)
    # print('leggo sample')
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    # print('leggo gen table')
    job_id = search.split("=")[-1]
    job_directory = current_working_directory + "Results/" + job_id + "/"
    file_to_grep = job_directory + "." + job_id + ".bestMerge.txt"
    sample_grep_result = (
        current_working_directory
        + "Results/"
        + job_id
        + "/"
        + job_id
        + "."
        + sample
        + "."
        + guide
        + ".private.txt"
    )

    # if not os.path.exists(current_working_directory + 'Results/' + job_id + '/' + job_id + '.' + sample + '.' + guide + '.sample_card.txt'):
    df = pd.read_csv(
        job_directory
        + job_id
        + ".summary_by_samples."
        + guide
        + "_"
        + filter_criterion
        + ".txt",
        sep="\t",
        skiprows=2,
        index_col=0,
        header=None,
        na_filter=False,
    )
    try:
        int_sample = int(sample)
    except:
        int_sample = sample

    personal = df.loc[int_sample, 4]
    pam_creation = df.loc[int_sample, 7]

    # file_to_grep = job_directory + '.' + job_id + '.bestMerge.txt'
    integrated_file_name = glob(job_directory + "*integrated*")[0]
    integrated_file_name = str(integrated_file_name)
    # integrated_to_grep = job_directory+job_id + \
    #     '.bestMerge.txt.integrated_results.tsv'
    integrated_to_grep = integrated_file_name
    integrated_personal = (
        job_directory + job_id + "." + sample + "." + guide + ".personal_targets.txt"
    )
    integrated_private = (
        job_directory + job_id + "." + sample + "." + guide + ".private_targets.tsv"
    )

    path_db = glob(current_working_directory +
                   "Results/" + job_id + "/.*.db")[0]
    path_db = str(path_db)
    conn = sqlite3.connect(path_db)
    c = conn.cursor()

    query_cols = get_query_column(filter_criterion)
    print(query_cols)

    result_personal = pd.read_sql_query(
        "SELECT * FROM final_table WHERE \"{}\"='{}' AND \"{}\" LIKE '%{}%'".format(
            GUIDE_COLUMN, guide, query_cols['samples'], sample
        ),
        conn,
    )
    # sort personal data targets
    order = False
    if filter_criterion == 'fewest':
        order = True
    result_personal = result_personal.sort_values(
        [query_cols['sort']], ascending=[order])
    # extract sample private targets
    result_private = result_personal[result_personal[query_cols['samples']] == sample]

    # print(result_personal)
    # print(result_private)

    conn.commit()
    conn.close()

    result_personal.to_csv(integrated_personal, sep="\t", index=False)
    result_private.to_csv(integrated_private, sep="\t", index=False)

    integrated_private_zip = integrated_private.replace("tsv", "zip")

    os.system(f"zip -j {integrated_private_zip} {integrated_private}")

    # plot for images in personal card
    # print('faccio personal')
    os.system(
        f"python {app_main_directory}/PostProcess/CRISPRme_plots_personal.py {integrated_personal} {current_working_directory}/Results/{job_id}/imgs/ {guide}.{sample}.personal > /dev/null 2>&1"
    )
    # print('faccio private')
    os.system(
        f"python {app_main_directory}/PostProcess/CRISPRme_plots_personal.py {integrated_private} {current_working_directory}/Results/{job_id}/imgs/ {guide}.{sample}.private > /dev/null 2>&1"
    )
    os.system(f"rm -f {integrated_personal}")

    private = result_private.shape[0]

    results_table = pd.DataFrame(
        [[len(result_personal.index), pam_creation, len(result_private.index)]],
        columns=["Personal", "PAM Creation", "Private"],
    ).astype(str)
    # else:
    #     pass

    try:  # to read the private targets file, if not created, pass
        ans = result_private
    except:
        pass

    # image for personal and private
    try:
        image_personal_top = "data:image/png;base64,{}".format(
            base64.b64encode(
                open(
                    current_working_directory
                    + "Results/"
                    + job_id
                    + f"/imgs/CRISPRme_{filter_criterion}_top_1000_log_for_main_text_{guide}.{sample}.personal.png",
                    "rb",
                ).read()
            ).decode()
        )
        image_private_top = "data:image/png;base64,{}".format(
            base64.b64encode(
                open(
                    current_working_directory
                    + "Results/"
                    + job_id
                    + f"/imgs/CRISPRme_{filter_criterion}_top_1000_log_for_main_text_{guide}.{sample}.private.png",
                    "rb",
                ).read()
            ).decode()
        )
    except:
        sys.stderr.write("PERSONAL AND PRIVATE LOLLIPOP PLOTS NOT GENERATED")

    filter_criterion = read_json(job_id)

    try:
        # file_to_load = job_id + '.' + sample + '.' + guide + '.private.zip'
        file_to_load = job_id + "." + sample + "." + guide + ".private_targets.zip"
        # #print(file_to_load)
        # ans = ans[COL_BOTH]
        out_1 = [
            html.A(
                "Download private targets",
                href=URL + "/Results/" + job_id + "/" + file_to_load,
                target="_blank",
            ),
            False,
            [
                html.P("Top 100 Personal Targets ordered by "+filter_criterion),
                html.A(
                    html.Img(
                        src=image_personal_top,
                        id="sample-personal-top",
                        width="100%",
                        height="auto",
                    ),
                    target="_blank",
                ),
            ],
            [
                html.P("Top 100 Private Targets ordered by "+filter_criterion),
                html.A(
                    html.Img(
                        src=image_private_top,
                        id="sample-private-top",
                        width="100%",
                        height="auto",
                    ),
                    target="_blank",
                ),
            ],
            dash_table.DataTable(
                css=[{"selector": ".row", "rule": "margin: 0"}],
                id="results-table",
                columns=[{"name": i, "id": i} for i in results_table.columns],
                data=results_table.to_dict("records"),
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
            ),
            dash_table.DataTable(
                css=[{"selector": ".row", "rule": "margin: 0"}],
                id="results-table-risk",
                # columns=[{"name": COL_BOTH[count], "id": i, 'hideable':True}
                #          for count, i in enumerate(ans.columns)],
                columns=[
                    {"name": i, "id": i, "hideable": True}
                    for count, i in enumerate(ans.columns)
                ],
                data=ans.to_dict("records"),
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
                style_table={
                    "overflowX": "scroll",
                    "overflowY": "scroll",
                    "max-height": "300px",
                },
            ),
        ]
    except:
        out_1 = [
            html.A(
                "Download private targets",
                href=URL + "/Results/" + job_id + "/" + file_to_load,
                target="_blank",
            ),
            True,
            [
                html.P("Top 100 Personal Targets ordered by "+filter_criterion),
                html.A(
                    html.Img(
                        src=image_personal_top,
                        id="sample-personal-top",
                        width="100%",
                        height="auto",
                    ),
                    target="_blank",
                ),
            ],
            [
                html.P("Top 100 Private Targets ordered by "+filter_criterion),
                html.A(
                    html.Img(
                        src=image_private_top,
                        id="sample-private-top",
                        width="100%",
                        height="auto",
                    ),
                    target="_blank",
                ),
            ],
            dash_table.DataTable(
                css=[{"selector": ".row", "rule": "margin: 0"}],
                id="results-table",
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
                columns=[{"name": i, "id": i} for i in results_table.columns],
                data=results_table.to_dict("records"),
            ),
            [],
        ]

    return list(out_1)


# Load the table/children under the tab value


@app.callback(
    Output("div-tab-content", "children"),
    [
        Input("tabs-reports", "value"),
        Input("general-profile-table", "selected_cells"),
        Input("target_filter_dropdown", "value"),
    ],
    [
        State("general-profile-table", "data"),
        State("url", "search"),
        State("div-genome-type", "children"),
    ],
)
def updateContentTab(
    value, sel_cel, filter_criterion, all_guides, search, genome_type
):
    if value is None or sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate

    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    job_id = search.split("=")[-1]
    job_directory = current_working_directory + "Results/" + job_id + "/"

    with open(current_working_directory + "Results/" + job_id + "/.Params.txt") as p:
        all_params = p.read()
        mms = (next(s for s in all_params.split("\n") if "Mismatches" in s)).split(
            "\t"
        )[-1]
        genome_selected = (
            next(s for s in all_params.split("\n") if "Genome_selected" in s)
        ).split("\t")[-1]
        max_bulges = (
            next(s for s in all_params.split("\n") if "Max_bulges" in s)
        ).split("\t")[-1]
        pam = (next(s for s in all_params.split(
            "\n") if "Pam" in s)).split("\t")[-1]
        nuclease = (next(s for s in all_params.split("\n") if "Nuclease" in s)).split(
            "\t"
        )[-1]

    fl = []
    fl.append(html.Br())

    if nuclease != "SpCas9":
        CFD_notification = html.Div(
            "CFD score is not calculated if the used nuclease is not SpCas9"
        )
        filter_criterion = "fewest"
    else:
        CFD_notification = html.Div("", hidden=True)

    pam_at_start = False
    if str(guide)[0] == "N":
        pam_at_start = True
    if pam_at_start:
        fl.append(html.H5("Focus on: " + str(pam) +
                  str(guide).replace("N", "")))
    else:
        fl.append(html.H5("Focus on: " + str(guide).replace("N", "") + str(pam)))

    if (
        value == "tab-summary-by-guide"
    ):  # BUG se cambio guida selezionata due volte mi cambia il mms mettendo a 0, provare con un div nascosto
        # Show Summary by Guide table
        fl.append(
            html.P(
                [
                    "Summary table counting the number of targets found in the Reference and Variant Genome for each combination of Bulge Type, Bulge Size and Mismatch. Select 'Show Targets' to view the corresponding list of targets. ",
                ]
            )
        )
        fl.append(html.Br())
        df = pd.read_csv(
            job_directory
            + job_id
            + ".summary_by_guide."
            + guide
            + "_"
            + filter_criterion
            + ".txt",
            sep="\t",
            na_filter=False,
        )
        more_info_col = []
        total_col = []
        for i in range(df.shape[0]):
            more_info_col.append("Show Targets")
            total_col.append(df["Bulge Size"])
        df[""] = more_info_col

        fl.append(
            html.Div(
                generate_table(
                    df, "table-summary-by-guide", genome_type, guide, job_id
                ),
                style={"text-align": "center"},
            )
        )
        return fl
    elif value == "tab-summary-by-sample":
        # Show Summary by Sample table
        fl.append(
            html.P(
                "Summary table counting the number of targets found in the Variant Genome for each sample. Filter the table by selecting the Population or Superpopulation desired from the dropdowns."
            )
        )
        if genome_type == "both":
            # col_names_sample = ['Sample', 'Sex', 'Population', 'Super Population',  'Targets in Reference', 'Targets in Enriched', 'Targets in Population', 'Targets in Super Population', 'PAM Creation', 'Class']
            col_names_sample = [
                "Sample",
                "Sex",
                "Population",
                "Super Population",  # 'Targets in Reference',
                "Targets in Variant",
                "Targets in Population",
                "Targets in Super Population",
                "PAM Creation",
            ]
            df = pd.read_csv(
                job_directory
                + job_id
                + ".summary_by_samples."
                + guide
                + "_"
                + filter_criterion
                + ".txt",
                sep="\t",
                names=col_names_sample,
                skiprows=2,
                na_filter=False,
            )
            df = df.sort_values("Targets in Variant", ascending=False)
        else:
            col_names_sample = [
                "Sample",
                "Sex",
                "Population",
                "Super Population",  # 'Targets in Reference',
                "Targets in Variant",
                "Targets in Population",
                "Targets in Super Population",
                "PAM Creation",
            ]
            df = pd.read_csv(
                job_directory
                + job_id
                + ".summary_by_samples."
                + guide
                + "_"
                + filter_criterion
                + ".txt",
                sep="\t",
                names=col_names_sample,
                skiprows=2,
                na_filter=False,
            )
            df = df.sort_values("Targets in Variant", ascending=False)

        more_info_col = []
        for i in range(df.shape[0]):
            more_info_col.append("Show Targets")
        df[""] = more_info_col

        population_1000gp = associateSample.loadSampleAssociation(
            job_directory + ".sampleID.txt"
        )[2]
        super_populations = [{"label": i, "value": i}
                             for i in population_1000gp.keys()]
        populations = []
        for k in population_1000gp.keys():
            for i in population_1000gp[k]:
                populations.append({"label": i, "value": i})
        fl.append(
            html.Div(
                [
                    html.Div(
                        job_directory + job_id + ".summary_by_samples." + guide,
                        style={"display": "none"},
                        id="div-info-summary_by_sample",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Div(
                                    dcc.Dropdown(
                                        options=super_populations,
                                        id="dropdown-superpopulation-sample",
                                        placeholder="Select a Super Population",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    dcc.Dropdown(
                                        options=populations,
                                        id="dropdown-population-sample",
                                        placeholder="Select a Population",
                                    )
                                )
                            ),
                            # dbc.Col(html.Div(dcc.Dropdown( id = 'dropdown-sample', placeholder = 'Select a Sample'))),
                            dbc.Col(
                                html.Div(
                                    dcc.Input(
                                        id="input-sample", placeholder="Select a Sample"
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    html.Button(
                                        "Filter", id="button-filter-population-sample"
                                    )
                                )
                            ),
                        ]
                    ),
                    dbc.Row(
                        dbc.Col(
                            html.Div(
                                [
                                    html.P(
                                        "Generating download link, Please wait...",
                                        id="download-link-summary_by_sample",
                                    ),
                                    dcc.Interval(
                                        interval=1 * 1000,
                                        id="interval-summary_by_sample",
                                    ),
                                ]
                            )
                        )
                    ),
                ],
                style={"width": "50%"},
            )
        )
        fl.append(
            html.Div(
                "None,None,None",
                id="div-sample-filter-query",
                style={"display": "none"},
            )
        )  # Folr keep current filter:  Superpop,Pop
        fl.append(
            html.Div(
                generate_table_samples(df, "table-samples", 1, guide, job_id),
                style={"text-align": "center"},
                id="div-table-samples",
            )
        )
        fl.append(
            html.Div(
                [
                    html.Button("Prev", id="prev-page-sample"),
                    html.Button("Next", id="next-page-sample"),
                ],
                style={"text-align": "center"},
            )
        )
        max_page = len(df.index)
        max_page = math.floor(max_page / 10) + 1
        fl.append(html.Div("1/" + str(max_page),
                  id="div-current-page-table-samples"))
        return fl
    elif value == "tab-summary-by-position":
        # Show Summary by position table
        fl.append(
            html.P(
                "Summary table containing all the targets found in a specific range of positions (chr, start, end) of the genome."
            )
        )

        fl.append(
            html.P(
                "Filter the table by selecting the chromosome of interest and writing the start and end position of the region to view."
            )
        )
        # Dropdown chromosomes
        try:
            onlyfile = [
                f
                for f in listdir(
                    current_working_directory + "Genomes/" + genome_selected
                )
                if (
                    isfile(
                        join(
                            current_working_directory + "Genomes/" + genome_selected, f
                        )
                    )
                    and (f.endswith(".fa") or f.endswith(".fasta"))
                )
            ]
        except:
            onlyfile = ["chr" + str(i) + ".fa" for i in range(1, 23)]
            onlyfile.append("chrX.fa")
            # NOTE in case no chr in GENOMES/ i put 22 chr + X Y M
            onlyfile.append("chrY.fa")
            onlyfile.append("chrM.fa")
        # removed .fa for better visualization
        onlyfile = [x[: x.rfind(".")] for x in onlyfile]
        chr_file = []
        chr_file_unset = []
        for chr_name in onlyfile:
            chr_name = chr_name.replace(".enriched", "")
            if "_" in chr_name:
                chr_file_unset.append(chr_name)
            else:
                chr_file.append(chr_name)
        chr_file.sort(
            key=lambda s: [
                int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)
            ]
        )
        chr_file_unset.sort(
            key=lambda s: [
                int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)
            ]
        )
        chr_file += chr_file_unset
        chr_file = [{"label": chr_name, "value": chr_name}
                    for chr_name in chr_file]

        # Colonne tabella: chr, pos, target migliore, min mm, min bulges, num target per ogni categoria di mm e bulge, show targets; ordine per total, poi mm e poi bulge
        # TODO inserire failsafe se non ci sono chr, esempio elenco chr da 1 a 22
        fl.append(
            html.Div(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Div(
                                    dcc.Dropdown(
                                        options=chr_file,
                                        id="dropdown-chr-table-position",
                                        placeholder="Select a chromosome",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    dcc.Input(
                                        placeholder="Start Position",
                                        id="input-position-start",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    dcc.Input(
                                        placeholder="End Position",
                                        id="input-position-end",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    html.Button(
                                        "Filter", id="button-filter-position")
                                )
                            ),
                            html.Br()
                            # )
                        ]
                    ),
                ],
                style={"width": "50%"},
            )
        )
        # ###print('Position dataframe ready', time.time() - start_time)
        # Folr keep current filter:  chr,pos_start,pos_end
        fl.append(
            html.Div(
                "None,None,None",
                id="div-position-filter-query",
                style={"display": "none"},
            )
        )
        # start_time = time.time()
        fl.append(html.Br())
        fl.append(
            html.Div(style={"text-align": "center"}, id="div-table-position"))
        max_page = 1
        fl.append(html.Div("1/" + str(max_page),
                  id="div-current-page-table-position"))
        fl.append(
            html.Div(
                mms + "-" + max_bulges,
                id="div-mms-bulges-position",
                style={"display": "none"},
            )
        )
        return fl
    elif value == "tab-graphical-sample-card":
        df = pd.read_csv(
            job_directory
            + job_id
            + ".summary_by_samples."
            + guide
            + "_"
            + filter_criterion
            + ".txt",
            skiprows=2,
            sep="\t",
            header=None,
            na_filter=False,
        )
        samples = df.iloc[:, 0]
        fl.append(
            html.P(
                "Summary page containing the single Personal Risk card to be inspected and downloaded"
            )
        )
        fl.append(
            html.Div(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Div(
                                    dcc.Dropdown(
                                        id="dropdown-sample-card",
                                        options=[
                                            {"label": sam, "value": sam}
                                            for sam in samples
                                        ],
                                        placeholder="Select a Sample",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    html.Button(
                                        "Generate", id="button-sample-card")
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    id="download-link-personal-card", hidden=True)
                            ),
                        ]
                    ),
                ],
                style={"width": "50%"},
            )
        )
        fl.append(
            html.Div(
                [
                    html.Br(),
                    dbc.Row(
                        [
                            dbc.Col(html.Div("", id="div-personal-plot")),
                            dbc.Col(html.Div("", id="div-private-plot")),
                        ]
                    ),
                ]
            )
        )
        fl.append(
            html.Div(
                "",
                id="div-table-sample-card",
                style={
                    "text-align": "center",
                    "margin-left": "1%",
                    "margin-right": "1%",
                },
            )
        )
        fl.append(
            html.Div(
                "",
                id="div-top-target-sample-card",
                style={
                    "text-align": "center",
                    "margin-left": "1%",
                    "margin-right": "1%",
                },
            )
        )
        fl.append(html.Div("", id="div-sample-card"))
        return fl
    elif value == "tab-query-table":
        fl.append(
            html.P(
                "Summary page to query the final result file selecting one/two column to group by the table and extract requested targets"
            )
        )
        all_value = {
            "Target1 :with highest CFD": [
                "Mismatches",
                "Bulges",
                "Mismatches+bulges",
                "CFD_score",
                "CFD_risk_Score",
            ],  # , 'Highest_CFD_Absolute_Risk_Score'
            "Target2 :with lowest Mismatches + Bulge Count": [
                "Mismatches",
                "Bulges",
                "Mismatches+bulges",
                "CFD_score",
                "CFD_risk_Score",
            ],
        }  # , 'CFD_Absolute_Risk_Score'
        all_options = {
            "Target1 :with highest CFD": [
                " Mismatches",
                " Bulges",
                " Mismatch+Bulges",
                " CFD",
                " Risk Score",
            ],  # , ' Absolute Risk Score'
            "Target2 :with lowest Mismatches + Bulges Count": [
                " Mismatches",
                " Bulges",
                " Mismatch+Bulges",
                " CFD",
                " Risk Score",
            ],
        }  # , ' Absolute Risk Score'

        label = [{"label": lab} for lab in all_options.keys()]
        value = [{"value": val} for val in all_value.keys()]
        target_opt = [label, value]

        query_tab_content = html.Div(
            [
                # dbc.Row(
                #     dbc.Col(
                #         html.Div(
                #             [
                #                 html.H4('Select filter criteria for targets'),
                #                 dcc.Dropdown(options=[
                #                     {'label': 'Fewest Mismatches and Bulges',
                #                         'value': 'fewest'},
                #                     {'label': 'CFD score', 'value': 'CFD'},
                #                     {'label': 'CRISTA Score', 'value': 'CRISTA'}
                #                 ], value='CFD',
                #                     id='target_filter_dropdown'
                #                 )
                #             ]))),
                dbc.Row(  # row with main group by, secondo group by and thresholds
                    [
                        dbc.Col(  # col0 phantom target select
                            [
                                html.Div(
                                    [
                                        html.H4("Order by"),
                                        dcc.RadioItems(
                                            id="target",
                                            options=target_opt,
                                            # options=[{'label': k, 'value': k} for k in all_options.keys()],
                                            value="Target1 :with highest CFD",
                                        ),
                                    ]
                                )
                            ],
                            style={"display": "none"},
                        ),
                        dbc.Col(  # col1 main group by
                            html.Div(
                                [
                                    html.H4("Group by"),
                                    dcc.RadioItems(
                                        id="order", value="CFD_score"),
                                ]
                            ),
                            width=3,
                        ),
                        dbc.Col(  # col2 second group by
                            html.Div(
                                [
                                    html.H4("And group by"),
                                    html.P(
                                        "First select the left group by value",
                                        id="secondtext",
                                    ),
                                    dcc.RadioItems(id="multiorder"),
                                ]
                            ),
                            width=3,
                        ),
                        dbc.Col(  # select threshold
                            html.Div(
                                [
                                    html.H4("Select thresholds"),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                [
                                                    html.Div(
                                                        [
                                                            html.H6("Min"),
                                                            dcc.Dropdown(
                                                                id="sholddrop"
                                                            ),
                                                        ]
                                                    ),
                                                ]
                                            ),
                                            dbc.Col(
                                                [
                                                    html.Div(
                                                        [
                                                            html.H6("Max"),
                                                            dcc.Dropdown(
                                                                id="maxdrop"),
                                                        ]
                                                    )
                                                ]
                                            ),
                                        ]
                                    ),
                                ]
                            ),
                            width=3,
                        ),
                        dbc.Col(
                            html.Div(
                                [
                                    html.H4("Select ordering"),
                                    dcc.RadioItems(
                                        id="Radio-asc-1",
                                        options=[
                                            {"label": " Ascending",
                                                "value": "ASC"},
                                            {"label": " Descending",
                                                "value": "DESC"},
                                        ],
                                        value="DESC",
                                        labelStyle={
                                            "display": "inline-block",
                                            "margin": "10px",
                                        },
                                    ),
                                ]
                            ),
                            width=3,
                        ),
                    ],
                    justify="center",
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                html.Div(
                                    html.Button(
                                        "Submit",
                                        id="submit-val",
                                        n_clicks=0,
                                        # style={
                                        #     'position': 'absolute', 'center': '50%'}
                                    ),
                                )
                            ],
                            width={"size": 1},
                        ),
                        dbc.Col(
                            [
                                html.Div(
                                    html.Button(
                                        "Reset",
                                        id="reset-val",
                                        n_clicks=0,
                                        # style={'position': 'absolute',
                                        #        'center': '50%'}
                                    )
                                )
                            ],
                            width={"size": 1, "offset": 1},
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.Br(),
                                        html.Hr(),
                                    ]
                                )
                            ]
                        ),
                    ],
                    justify="center",
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                CFD_notification,
                                # html.Br(),
                                html.P(
                                    "Export will download 1000 lines contained in the current view of the table"
                                ),
                                html.Div(
                                    dash_table.DataTable(
                                        css=[
                                            {
                                                "word-break": "break-all",
                                                "line-break": "anywhere",
                                                "overflow-wrap": "break-word",
                                                "selector": ".row",
                                                # 'rule': 'margin: 0',
                                                "rule": "margin: 0; overflow: inherit; word-break: break-all; overflow-wrap: break-word; line-break: anywhere;",
                                            }
                                        ],
                                        style_cell={
                                            "height": "auto",
                                            "textAlign": "left",
                                            # 'maxWidth': '500px'
                                        },
                                        export_format="xlsx",
                                        id="live_table",
                                        style_cell_conditional=[
                                            {
                                                "if": {
                                                    "column_id": "Variant_samples_(highest_CFD)"
                                                },
                                                "textAlign": "left",
                                                "minWidth": "180px",
                                                "width": "180px",
                                                "maxWidth": "180px",
                                                "overflow": "hidden",
                                            },
                                            {
                                                "if": {
                                                    "column_id": "Variant_samples_(fewest_mm+b)"
                                                },
                                                "textAlign": "left",
                                                "minWidth": "180px",
                                                "width": "180px",
                                                "maxWidth": "180px",
                                                "overflow": "hidden",
                                            },
                                            {
                                                "if": {
                                                    "column_id": "Variant_samples_(highest_CRISTA)"
                                                },
                                                "textAlign": "left",
                                                "minWidth": "180px",
                                                "width": "180px",
                                                "maxWidth": "180px",
                                                "overflow": "hidden",
                                            },
                                        ],
                                        style_table={
                                            "overflowX": "scroll",
                                            "overflowY": "scroll",
                                            "max-height": "300px",
                                        },
                                        page_current=0,
                                        page_size=1000,
                                        page_action="custom",
                                        tooltip_delay=0,
                                        tooltip_duration=None,
                                    ),
                                    id="div-query-table",
                                ),
                            ],
                        ),
                    ],
                ),
                html.Div(
                    [
                        dbc.Row(
                            dbc.Col(
                                [
                                    dbc.Alert(
                                        "Select a main order before submitting the query",
                                        id="message-alert",
                                        color="danger",
                                        dismissable=True,
                                        fade=True,
                                        is_open=False,
                                        duration=4000,
                                    ),
                                ]
                            ),
                        )
                    ],
                    style={"display": "inline-block"},
                ),
            ]
        )
        fl.append(query_tab_content)
        # ##print('table query', dff)
        # fl.append(

        return fl
    else:  # tab-graphical
        # Show Report images
        samp_style = {}
        if genome_type == "ref":
            samp_style = {"display": "none"}

        # fl.append(html.Br())
        fl.append(
            html.P(
                "Summary Graphical report collecting all the plots and images produced during the search"
            )
        )

        opt_mm = []
        total = int(mms) + int(max_bulges)
        for i in range(total + 1):
            opt_mm.append({"label": str(i), "value": str(i)})
        opt_blg = []
        for i in range(int(max_bulges) + 1):
            opt_blg.append({"label": str(i), "value": str(i)})

        if genome_type != "ref":
            population_1000gp = associateSample.loadSampleAssociation(
                job_directory + ".sampleID.txt"
            )[2]

            super_populations = [
                {"label": i, "value": i} for i in population_1000gp.keys()
            ]
            populations = []
            for k in population_1000gp.keys():
                for i in population_1000gp[k]:
                    populations.append({"label": i, "value": i})
        else:
            super_populations = []
            populations = []
        # CRISPRme_top_1000_log_for_main_text_{guide}_MMvBUL.png
        try:
            # if nuclease == 'SpCas9':  # IMPLEMENTARE DROPDOWN PER VISUALIZZARE CFD,CRISTA,MMBUL
            top1000_image = html.Div(
                html.A(
                    html.Img(
                        src="data:image/png;base64,{}".format(
                            base64.b64encode(
                                open(
                                    current_working_directory
                                    + "Results/"
                                    + job_id
                                    + f"/imgs/CRISPRme_{filter_criterion}_top_1000_log_for_main_text_{guide}.png",
                                    "rb",
                                ).read()
                            ).decode()
                        ),
                        id="top-1000-score",
                        width="80%",
                        height="auto",
                    ),
                    target="_blank",
                )
            )
            # else:
            #     top1000_image = html.Div(
            #         html.A(html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(
            #                         current_working_directory + 'Results/' + job_id + f'/imgs/CRISPRme_top_1000_log_for_main_text_{guide}_MMvBUL.png', 'rb').read()).decode()),
            #                         id='top-1000-score', width="80%", height="auto"),
            #                target="_blank")
            #     )
        except:
            top1000_image = html.Div("")

        total_buttons = [
            dbc.Col(
                html.Div(
                    [
                        html.P(
                            "Select total number of mismatches and/or bulges to consider, up to"
                        ),
                        dcc.Dropdown(
                            id="mm-dropdown",
                            options=opt_mm,
                            value="0",
                            clearable=False,
                        ),
                    ]
                )
            )
        ]
        sample_buttons = [
            dbc.Col(
                html.Div(
                    [
                        html.P("Select a Superpopulation", style=samp_style),
                        html.Div(
                            dcc.Dropdown(
                                options=super_populations,
                                id="dropdown-superpopulation-sample",
                                placeholder="SuperPopulation",
                                style=samp_style,
                            ),
                        ),
                    ]
                ),
                md=4,
            ),
            dbc.Col(
                html.Div(
                    [
                        html.P("Select a Population", style=samp_style),
                        html.Div(
                            dcc.Dropdown(
                                options=populations,
                                id="dropdown-population-sample",
                                placeholder="Population",
                                style=samp_style,
                            ),
                        ),
                    ]
                ),
                md=4,
            ),
            dbc.Col(
                html.Div(
                    [
                        html.P("Select a Sample", style=samp_style),
                        html.Div(
                            dcc.Dropdown(
                                id="dropdown-sample",
                                placeholder="Sample",
                                style=samp_style,
                            ),
                        ),
                    ]
                ),
                md=4,
            ),
        ]
        fl.append(
            html.Div(
                [
                    CFD_notification,
                    dbc.Row(dbc.Col(top1000_image, width={
                            "size": 10, "offset": 2})),
                    dbc.Row(total_buttons, justify="center"),
                    html.Br(),
                ]
            )
        )

        radar_chart_encode_gencode = dbc.Col(
            html.Div(id="div-radar-chart-encode_gencode")
        )
        populations_barplots = dbc.Col(html.Div(id="div-population-barplot"))
        radar_chart_sample_content = dbc.Row(id="row-radar-chart-sample")
        sample_image_content = html.Div(id="div-sample-image")

        if genome_type != "ref":
            graph_summary_both = [
                populations_barplots, radar_chart_encode_gencode]
        else:
            graph_summary_both = [radar_chart_encode_gencode]

        fl.append(html.Div([dbc.Row(graph_summary_both)]))

        fl.append(
            dbc.Row(
                dbc.Col(
                    html.Div(
                        [
                            "The GENCODE and ENCODE annotations are defined in detail ",
                            html.A(
                                "here",
                                target="_blank",
                                href="https://www.gencodegenes.org/human/",
                            ),
                            " and ",
                            html.A(
                                "here",
                                target="_blank",
                                href="https://screen.encodeproject.org/",
                            ),
                        ]
                    ),
                    width={"size": 6, "offset": 6},
                )
            )
        )

        cfd_path = job_directory + job_id + ".CFDGraph.txt"
        if not isfile(cfd_path):  # No file found
            return fl

        fl.extend(CFDGraph.CFDGraph(cfd_path))

        return fl
        # return guide + value
    raise PreventUpdate


# Read the uploaded file and converts into bit


def parse_contents(contents):
    content_type, content_string = contents.split(",")

    decoded = base64.b64decode(content_string)
    return decoded


# Perform expensive loading of a dataframe and save result into 'global store'
# Cache are in the Cache directory


@cache.memoize()
def global_store(value):
    """
    Caching dei file targets per una miglior performance di visualizzazione
    """
    if value is None:
        return ""
    target = [
        f
        for f in listdir(current_working_directory + "Results/" + value)
        if isfile(join(current_working_directory + "Results/" + value, f))
        and f.endswith("scores.txt")
    ]
    if not target:
        target = [
            f
            for f in listdir(current_working_directory + "Results/" + value)
            if isfile(join(current_working_directory + "Results/" + value, f))
            and f.endswith("targets.txt")
        ]

    df = pd.read_csv(
        current_working_directory + "Results/" + value + "/" + target[0],
        sep="\t",
        usecols=range(0, 38),
        na_filter=False,
    )
    df.rename(
        columns={
            "#Bulge type": "BulgeType",
            "#Bulge_type": "BulgeType",
            "Bulge Size": "BulgeSize",
            "Bulge_Size": "BulgeSize",
            "Doench 2016": "Doench2016",
            "Doench_2016": "Doench2016",
        },
        inplace=True,
    )
    return df


@app.callback(
    Output("result-table", "data"),
    [
        Input("result-table", "page_current"),
        Input("result-table", "page_size"),
        Input("result-table", "sort_by"),
        Input("result-table", "filter_query"),
    ],
    [State("url", "search"), State("url", "hash")],
)
def update_table(page_current, page_size, sort_by, filter, search, hash_guide):
    """
    La funzione ritorna uno split dei risultati in base ad un filtering o a un sort da parte dell'utente. Inoltre aggiorna i risultati
    visualizzati quando il bottone next page / prev page è cliccato. (Codice preso dalla pagina dash datatable sul sorting con python)
    Inoltre carica i file targets, o scores se presente, e lo trasforma in un dataframe, cambiando il nome delle colonne per farle corrispondere
    all'id delle colonne della tabella nella pagina.
    Se non ci sono targets ritorna un avviso di errore
    """
    job_id = search.split("=")[-1]
    job_directory = current_working_directory + "Results/" + job_id + "/"
    guide = hash_guide.split("#")[1]
    value = job_id
    if search is None:
        raise PreventUpdate

    filtering_expressions = filter.split(" && ")
    # filtering_expressions.append(['{crRNA} = ' + guide])
    df = global_store(value)
    dff = df[df["crRNA"] == guide]

    sort_by.insert(0, {"column_id": "Mismatches", "direction": "asc"})
    sort_by.insert(1, {"column_id": "BulgeSize", "direction": "asc"})
    # sort_by.insert(2, {'column_id': 'CFD', 'direction':'desc'})
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

        if operator in ("eq", "ne", "lt", "le", "gt", "ge"):
            # these operators match pandas series operator method names
            dff = dff.loc[getattr(dff[col_name], operator)(filter_value)].sort_values(
                [col["column_id"] for col in sort_by],
                ascending=[col["direction"] == "asc" for col in sort_by],
                inplace=False,
            )
        elif operator == "contains":
            dff = dff.loc[dff[col_name].str.contains(filter_value)]
        elif operator == "datestartswith":
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            dff = dff.loc[dff[col_name].str.startswith(filter_value)]

    if len(sort_by):
        dff = dff.sort_values(
            [
                "Samples" if col["column_id"] == "Samples Summary" else col["column_id"]
                for col in sort_by
            ],
            ascending=[col["direction"] == "asc" for col in sort_by],
            inplace=False,
        )

    # Check if results are not 0
    warning_no_res = ""
    with open(job_directory + job_id + ".targets.txt") as t:
        no_result = False
        t.readline()
        last_line = t.readline()
        if last_line == "" or last_line == "\n":
            no_result = True

    if no_result:
        warning_no_res = dbc.Alert(
            "No results were found with the given parameters", color="warning"
        )

    return dff.iloc[page_current * page_size: (page_current + 1) * page_size].to_dict(
        "records"
    )


# Callbacks for querying part--------------------------------------------------------------


# Return the table with the result of the query
@app.callback(
    # [Output('live_table', 'data'),
    [
        Output("live_table", "columns"),
        Output("live_table", "data"),
        Output("live_table", "tooltip_data"),
        Output("message-alert", "is_open"),
    ],
    [
        Input("submit-val", "n_clicks"),
        Input("live_table", "page_current"),
        Input("target_filter_dropdown", "value"),
    ],  # take this value (as state)
    [
        State("live_table", "page_size"),
        State("general-profile-table", "selected_cells"),
        State("target", "value"),
        State("order", "value"),
        State("general-profile-table", "data"),
        State("multiorder", "value"),
        State("sholddrop", "value"),
        State("Radio-asc-1", "value"),
        State("maxdrop", "value"),
        State("url", "search"),
        State("message-alert", "is_open"),
    ],
)
# see here
def update_output(
    n_clicks,
    page_current,
    filter_target_value,
    page_size,
    sel_cel,
    target,
    radio_order,
    all_guides,
    orderdrop,
    sholddrop,
    asc1,
    maxdrop,
    url,
    alert,
):
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]

    # target is the filter value to query on the db
    target = filter_target_value

    if n_clicks > 0:
        if radio_order == None:
            data = []
            tooltip_data = []
            return data, tooltip_data, not alert
        else:
            if sholddrop != None:
                alert = False
                data = query_manager.shold(
                    target,
                    n_clicks,
                    page_current,
                    page_size,
                    radio_order,
                    orderdrop,
                    sholddrop,
                    maxdrop,
                    asc1,
                    url,
                    guide,
                    current_working_directory,
                )
            else:
                data = query_manager.noshold(
                    target,
                    n_clicks,
                    page_current,
                    page_size,
                    radio_order,
                    orderdrop,
                    asc1,
                    url,
                    guide,
                    current_working_directory,
                )

            # find col to drop using the user filter
            # COPIARE PER FARE DROP COLONNE NON VOLUTE IN TARGET FILTER
            drop_col = list()
            for elem in list(data.columns):
                if filter_target_value == "fewest" and (
                    "highest_CFD" in elem or "highest_CRISTA" in elem
                ):
                    drop_col.append(elem)
                if filter_target_value == "CFD" and (
                    "fewest" in elem or "highest_CRISTA" in elem
                ):
                    drop_col.append(elem)
                if filter_target_value == "CRISTA" and (
                    "fewest" in elem or "highest_CFD" in elem
                ):
                    drop_col.append(elem)
            # drop column from datatable to show
            data.drop(drop_col, inplace=True, axis=1)
            # extract cols for datatable
            columns = [
                {"name": i, "id": i, "hideable": True}
                for count, i in enumerate(data.columns)
            ]

            # selct SNPs col to filter
            if filter_target_value == "fewest":
                snps = pd.DataFrame(data["Variant_info_genome_(fewest_mm+b)"]).to_dict(
                    "records"
                )
            if filter_target_value == "CFD":
                snps = pd.DataFrame(data["Variant_info_genome_(highest_CFD)"]).to_dict(
                    "records"
                )
            if filter_target_value == "CRISTA":
                snps = pd.DataFrame(
                    data["Variant_info_genome_(highest_CRISTA)"]
                ).to_dict("records")

            # extract data and list datas
            data = data.to_dict("records")
            tooltip_data = [
                {
                    column: {"value": str(value), "type": "markdown"}
                    for column, value in row.items()
                }
                for row in snps
            ]
    else:
        raise PreventUpdate
    # ##print('query table', data)
    return columns, data, tooltip_data, alert


# to get correct number of page
@app.callback(Output("live_table", "page_current"), [Input("submit-val", "n_clicks")])
def reset_pagenumber(n):
    if n > 0:
        a = 0
        return a
    else:
        raise PreventUpdate


@app.callback(Output("order", "options"), [Input("target", "value")])
def set_columns_options(selected_target):
    all_value = {
        "Target1 :with highest CFD": [
            "Mismatches",
            "Bulges",
            "Mismatches+bulges",
            "CFD_score",
            "CFD_risk_score",
        ],  # , 'Highest_CFD_Absolute_Risk_Score'
        "Target2 :with lowest Mismatches + Bulge Count": [
            "Mismatches",
            "Bulges",
            "Mismatches+bulges",
            "CFD_score",
            "CFD_risk_score",
        ],
    }  # , 'CFD_Absolute_Risk_Score'
    all_options = {
        "Target1 :with highest CFD": [
            " Mismatches",
            " Bulges",
            " Mismatch+Bulges",
            " Score",
            " Risk Score",
        ],  # , ' Absolute Risk Score'
        "Target2 :with lowest Mismatches + Bulges Count": [
            " Mismatches",
            " Bulges",
            " Mismatch+Bulges",
            " CFD",
            " Risk Score",
        ],
    }  # , ' Absolute Risk Score'
    gi = []
    for count in range(0, len(all_value[selected_target])):
        gi.append(
            {
                "label": all_options[selected_target][count],
                "value": all_value[selected_target][count],
            }
        )
    # return gi
    # ###print(main_order_dict)
    return gi


# callback to return the parameters in the various cases
@app.callback(
    [
        Output("multiorder", "options"),
        Output("sholddrop", "options"),
        Output(component_id="secondtext", component_property="style"),
    ],
    [Input("order", "value")],
)
def set_display_children(selected_order):
    target_value = {
        "Mismatches": ["Bulges", "Mismatches+bulges", "CFD"],
        "Bulges": ["Mismatches", "Mismatches+bulges", "CFD_score"],
        "Mismatches+bulges": ["Mismatches", "Bulges", "CFD_score"],
        "CFD_score": ["Mismatches", "Bulges", "Mismatches+bulges"],
        "CFD_risk_score": [],
    }
    target_label = {
        "Mismatches": [" Bulges", " Mismatch+Bulges", " Score"],
        "Bulges": [" Mismatches", " Mismatch+Bulges", " Score"],
        "Mismatches+bulges": [" Mismatches", " Bulges", " Score"],
        "CFD_score": [" Mismatches", " Bulges", " Mismatch+Bulges"],
        "Highest CFD Risk Score": [],
        "Highest CFD Absolute Risk Score": [],
        "CFD_risk_score": [],
        "CFD Absolute Risk Score": [],
    }

    gi = []
    if selected_order is not None:
        for count in range(0, len(target_value[selected_order])):
            gi.append(
                {
                    "label": target_label[selected_order][count],
                    "value": target_value[selected_order][count],
                }
            )

    if selected_order == None:
        return [], [], {"display": "block"}
    elif selected_order == "Mismatches":
        data = [
            {"label": "0", "value": "0"},
            {"label": "1", "value": "1"},
            {"label": "2", "value": "2"},
            {"label": "3", "value": "3"},
            {"label": "4", "value": "4"},
            {"label": "5", "value": "5"},
            {"label": "6", "value": "6"},
        ]
        # return [{'label': i, 'value': i} for i in target_options[selected_order]], data, {'display': 'none'}
        return gi, data, {"display": "none"}
    elif selected_order == "CFD_score":
        data = [
            {"label": "0.01", "value": "0.01"},
            {"label": "0.1", "value": "0.1"},
            {"label": "0.2", "value": "0.2"},
            {"label": "0.3", "value": "0.3"},
            {"label": "0.4", "value": "0.4"},
            {"label": "0.5", "value": "0.5"},
            {"label": "0.6", "value": "0.6"},
            {"label": "0.7", "value": "0.7"},
            {"label": "0.8", "value": "0.8"},
            {"label": "0.9", "value": "0.9"},
        ]
        # return [{'label': i, 'value': i} for i in target_options[selected_order]], data, {'display': 'none'}
        return gi, data, {"display": "none"}
    elif selected_order == "Mismatches+bulges":
        data = [
            {"label": "0", "value": "0"},
            {"label": "1", "value": "1"},
            {"label": "2", "value": "2"},
            {"label": "3", "value": "3"},
            {"label": "4", "value": "4"},
            {"label": "5", "value": "5"},
            {"label": "6", "value": "6"},
            {"label": "7", "value": "7"},
            {"label": "8", "value": "8"},
        ]
        # return [{'label': i, 'value': i} for i in target_options[selected_order]], data, {'display': 'none'}
        return gi, data, {"display": "none"}
    elif selected_order == "Bulges":
        data = [
            {"label": "0", "value": "0"},
            {"label": "1", "value": "1"},
            {"label": "2", "value": "2"},
        ]
        # return [{'label': i, 'value': i} for i in target_options[selected_order]], data, {'display': 'none'}
        return gi, data, {"display": "none"}
    else:
        return [], [], {"display": "none"}


@app.callback(
    Output("maxdrop", "options"), [
        Input("sholddrop", "value"), Input("order", "value")]
)
def maxdrop(sholddrop, order):
    if order == "Mismatches":
        if sholddrop:
            start_value = int(sholddrop)
            data = [{"label": str(i), "value": str(i)}
                    for i in range(start_value, 7)]
        else:
            data = []

    elif order == "CFD_score":
        if sholddrop:
            start_value = int(float(sholddrop) * 10)
            if start_value < 1:
                start_value = 1
                small = True
            else:
                small = False
            if start_value < 10:
                data = [
                    {"label": f"0.{i}", "value": f"0.{i}"}
                    for i in range(start_value, 10)
                ]
                data.append({"label": "1", "value": "1"})
                if small:
                    data.insert(0, {"label": "0.01", "value": "0.01"})
            else:
                data = []
        else:
            data = []

    elif order == "Bulges":
        if sholddrop:
            start_value = int(sholddrop)
            data = [{"label": str(i), "value": str(i)}
                    for i in range(start_value, 3)]
        else:
            data = []

    elif order == "Mismatches+bulges":
        if sholddrop:
            start_value = int(sholddrop)
            data = [{"label": str(i), "value": str(i)}
                    for i in range(start_value, 9)]
        else:
            data = []

    else:
        data = []
    return data


@app.callback(
    [
        Output("order", "value"),
        Output("multiorder", "value"),
        Output("maxdrop", "value"),
        Output("sholddrop", "value "),
        Output("Radio-asc-1", "value"),
    ],
    [Input("reset-val", "n_clicks")],
)
def resetbutton(n_clicks):
    if n_clicks > 0:
        return None, None, None, None, None
    else:
        return None, None, None, None, None
