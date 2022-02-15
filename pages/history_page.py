"""Construct the job history page. 

The user can navigate the CRISPRme job history page to retrieve the 
results obtained during past analyses.

The page dispalys a table with all the past CRISPRme jobs. The user can 
recover the results obtained from each analysis by clicking the link
on the job identifier.
"""


from .results_page_utils import RESULTS_DIR, PARAMS_FILE, LOG_FILE, GUIDES_FILE

from typing import Any, Dict, List, Optional, Tuple
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from app import app, current_working_directory, URL

import dash_html_components as html
import pandas as pd

import math
import os


def support_filter_history(
    result_df: pd.DataFrame, genome_filter: str, pam_filter: str
) -> Tuple[pd.DataFrame, int]:
    """Filter the results history and create the table to display on 
    the history page in CRISPRme webpage.

    ...

    Parameters
    ----------
    result_df : pd.DataFrame
        Results history
    genome_filter : str
        Genomes to exclude
    pam_filter : str
        PAMs to exclude

    Returns
    -------
    pd.DataFrame
        Filtered results history
    int
        Maximum page size
    """

    if not isinstance(result_df, pd.DataFrame):
        raise TypeError(f"Expected {type(pd.DataFrame).__name__}, got {type(result_df).__name__}")
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
    Output('results-table', 'style_data_conditional'),
    [Input('results-table', 'selected_cells')],
    [State('results-table', 'data')]
)
def highlight_row(sel_cel, all_guides) -> List[Dict[str, str]]:
    """
    """

    if sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate
    job_name = all_guides[int(sel_cel[0]["row"])]["Job"]
    assert isinstance(job_name, str)
    return [
        {
            "if":{
                "filter_query":"".join(["{Job} eq ", f"\"{job_name}\""])
            },
            "background-color":"rgba(0, 0, 255,0.15)"  # highlighting color
        }
    ]


def get_results() -> pd.DataFrame:
    """Retrieve the job history and create a table to display 
    past jobs along with the main info related to each job.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    pd.DataFrame
        Table of jobs history
    """
    
    # retrieve results stored in directories (in /Results)
    results_dirs = []
    for resdir in os.listdir(os.path.join(current_working_directory, RESULTS_DIR)):
        if (
            os.path.isdir(
                os.path.join(current_working_directory, RESULTS_DIR, resdir)
            ) and os.path.isfile(
                os.path.join(current_working_directory, RESULTS_DIR, resdir, PARAMS_FILE))
        ):
            results_dirs.append(resdir)
    # allow empty results?
    assert len(results_dirs) <= len(
        os.listdir(os.path.join(current_working_directory, RESULTS_DIR))
    )
    cols = [
            "Job", 
            "Genome", 
            "Variants", 
            "Mismatches", 
            "DNA bulge",
            "RNA bulge",
            "PAM",
            "Number of Guides"
    ]
    result_param_df = pd.DataFrame(columns=cols)
    for job_id in results_dirs:
        try:
            if os.path.exists(
                os.path.join(
                    current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE
                )
            ):
                try:
                    with open(
                        os.path.join(
                            current_working_directory, 
                            RESULTS_DIR, 
                            job_id,
                            PARAMS_FILE
                        )
                    ) as handle_params:
                        params = handle_params.read()
                        mms = (
                            next(
                                s for s in params.split("\n") if "Mismatches" in s
                            )
                        ).split("\t")[-1]
                        genome_selected = (
                            next(
                                s for s in params.split("\n") if "Genome_selected" in s
                            )
                        ).split("\t")[-1]
                        try:
                            with open(
                                os.path.join(
                                    current_working_directory, 
                                    RESULTS_DIR,
                                    job_id,
                                    LOG_FILE
                                )
                            ) as handle_log:
                                log = handle_log.read()
                        except OSError as e:
                            raise e
                        job_start = (
                            next(
                                s for s in log.split("\n") if "Job\tStart" in s
                            )
                        ).split("\t")[-1]
                        if "+" in genome_selected:
                            genome_selected = "".join(
                                [genome_selected.split("+")[0], "+"]
                            )
                        dna = (
                            next(
                                s for s in params.split("\n") if "DNA" in s
                            )
                        ).split("\t")[-1]
                        rna = (
                            next(
                                s for s in params.split("\n") if "RNA" in s
                            )
                        ).split("\t")[-1]
                        genome_idx = (
                            next(
                                s for s in params.split("\n") if "Genome_idx" in s
                            )
                        ).split("\t")[-1]
                        if "+" in genome_idx:
                            genome_idx_split = [
                                e.split("+")[-1] for e in genome_idx.split(",")
                            ]
                            genome_idx = ",".join(genome_idx_split)
                        else:
                            genome_idx = "Reference"
                        pam = (
                            next(
                                s for s in params.split("\n") if "Pam" in s
                            )
                        ).split("\t")[-1]
                        if os.path.exists(
                            os.path.join(
                                current_working_directory, 
                                RESULTS_DIR, 
                                job_id,
                                GUIDES_FILE
                            )
                        ):
                            try:
                                with open(
                                    os.path.join(
                                        current_working_directory, 
                                        RESULTS_DIR, 
                                        job_id,
                                        GUIDES_FILE
                                    )
                                ) as handle_guides:
                                    n_guides = str(
                                        len(
                                            handle_guides.read().strip().split("\n")
                                        )
                                    )
                            except OSError as e:
                                raise e
                        else:
                            n_guides = "NA"
                        result_param_df = result_param_df.append(
                            {
                                "Job":job_id, 
                                "Genome":genome_selected, 
                                "Variants":genome_idx, 
                                "Mismatches":mms, 
                                "DNA bulge":dna,
                                "RNA bulge":rna, 
                                "PAM":pam, 
                                "Number of Guides":n_guides, 
                                "Start":job_start
                            }, 
                            ignore_index=True
                        )
                except OSError as e:
                    raise e
        except:
            continue
    try:
        result_param_df["Start"] = pd.to_datetime(
            result_param_df["Start"]
        )
        result_param_df.sort_values(
            by=["Start"], ascending=False, inplace=True
        )   
    except:
        pass
    # resultParamDataframe = resultParamDataframe.sort_values(
    #     ['Mismatches', 'DNA_bulge', 'RNA_bulge'], ascending=[True, True, True])
    return result_param_df


def generate_table_results(
    results_df: pd.DataFrame, page: int, max_rows: Optional[int] = 1000000
) -> List[html.Table]:
    """Generate the table displaying the job history. The table is displayed
    in History page page of CRISPRme website.

    The user can access the results of past jobs and CRISPRme analysis
    selecting the requested job ID and retrieve the results going through
    the job's link.

    ...

    Parameters
    ----------
    results_df : pd.DataFrame
    page : int
    max_rows : int

    Returns
    -------
    List[html.Table]
        History page layout
    """

    if not isinstance(results_df, pd.DataFrame):
        raise TypeError(f"Expected {type(pd.DataFrame).__name__}, got {type(results_df).__name__}")
    if not isinstance(page, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page).__name__}")
    if not isinstance(max_rows, int):
        raise TypeError(f"Expected {int.__name__}, got {type(max_rows).__name__}")
    final_list = []  # page layout
    remaining_rows = results_df.shape[0] - (page - 1) * max_rows
    header = html.Thead(
        html.Tr(
            [
                html.Th(
                    col, 
                    style={"vertical-align":"middle", "text-align":"center"}
                ) if col != "Load" and col != "Delete" else html.Th(
                    "", 
                    style={"vertical-align":"middle", "text-align":"center"}
                ) for col in results_df.columns
            ]
        )
    )
    # build history page body
    body_history = []
    add_button = 0
    for i in range(min(remaining_rows, max_rows)):
        add_button += 1
        row_hist = []
        for col in results_df.columns:
            if col == "Job":
                job_id = str(results_df.iloc[i + (page - 1) * max_rows][col])
                row_hist.append(
                    html.Td(
                        html.A(
                            job_id, 
                            target="_blank", 
                            href=os.path.join(URL, f"load?job={job_id}"),
                        ), style={
                            "vertical-align":"middle", "text-align":"center"
                        }
                    )
                )
            else:
                row_hist.append(
                    html.Td(
                        results_df.iloc[i + (page - 1) * max_rows][col], 
                        style={
                            "vertical-align":"middle", "text-align":"center"
                        }
                    )
                )
        body_history.append(html.Tr(row_hist))
    final_list.append(
        html.Table([header,html.Tbody(body_history)], style={"display":"inline-block"},)
    )
    # Add hidden buttons for callback removeJobId compatibility
    for i in range(add_button, 10):  
        final_list.append(
            html.Button(
                str(i), 
                id="button-delete-history-" + str(i), **{"data-jobid":"None"}, 
                style={"display":"none"}
            )
        )
    return final_list


def historyPage():
    """Construct CRISPRme history webpage.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    int
        Current webpage
    """

    results = get_results()
    final_list = []
    final_list.append(
        html.Div(
            [
                html.H3("Results History"),
                html.P(
                    str(
                        "List of available results. Click on the link to open "
                        "the corresponding load page in a new tab."
                    )
                )
            ]
        )
    )
    final_list.append(
        html.Div(
            "None,None", 
            id="div-history-filter-query", 
            style={"display":"none"}
        )
    )
    final_list.append(
        html.Div(
            generate_table_results(results, 1),
            id="div-history-table",
            style={"text-align":"center"}
        ),
    )
    final_list.append(
        html.Div(id="div-remove-jobid", style={"display":"none"})
    )
    final_list.append(
        html.Div(
            [
                html.Br(),
            ],
            style={"text-align":"center"}
        )
    )
    max_page = len(results.index)
    max_page = math.floor(max_page / 1000000) + 1
    page = html.Div(final_list, style={"margin":"1%"})
    return page
