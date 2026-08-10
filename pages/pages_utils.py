"""Define static variables and utilities functions used throughout CRISPRme's
webpages.
"""

from app import operators, current_working_directory

from typing import Dict, List, Optional, Tuple
from glob import glob

from dash import html
import pandas as pd

import subprocess
import base64
import sys
import os


# Define DNA alphabet
DNA_ALPHABET = ["A", "C", "G", "T"]
# define IUPAC alphabet as valid characters for CRISPRme queries
VALID_CHARS = {
    "A",
    "T",
    "C",
    "G",
    "R",
    "Y",
    "S",
    "W",
    "K",
    "M",
    "B",
    "D",
    "H",
    "V",
    "a",
    "t",
    "c",
    "g",
    "r",
    "y",
    "s",
    "w",
    "k",
    "m",
    "b",
    "d",
    "h",
    "v",
}
# number of entries in report table (for each table page)
PAGE_SIZE = 10
# number of barplots in each row of Populations Distributions
BARPLOT_LEN = 4
# column names for reference report
COL_REF = [
    "Bulge Type",
    "crRNA",
    "Off target_motif",
    "Reference sequence",
    "Chromosome",
    "Position",
    "Direction",
    "Mismatches",
    "Bulge Size",
    "PAM gen",
    "Samples",
    "Variant",
    "CFD",
    "CFD ref",
    "Highest CFD Risk Score",
    "AF",
    "Annotation Type",
]
# reference column types
COL_REF_TYPE = [
    "text",
    "text",
    "text",
    "text",
    "text",
    "numeric",
    "numeric",
    "text",
    "numeric",
    "numeric",
    "text",
    "text",
    "text",
    "numeric",
    "numeric",
    "numeric",
    "numeric",
    "text",
]
# reference columns renaming
COL_REF_RENAME = {
    0: "Bulge Type",
    1: "crRNA",
    2: "Off target motif",
    3: "Reference sequence",
    4: "Chromosome",
    5: "Position",
    6: "Cluster Position",
    7: "Direction",
    8: "Mismatches",
    9: "Bulge Size",
    10: "Total",
    11: "PAM gen",
    12: "Variant Unique",
    13: "Samples",
    14: "Annotation Type",
    15: "Real Guide",
    16: "rsID",
    17: "AF",
    18: "Variant",
    19: "#Seq in cluster",
    20: "CFD",
    21: "CFD ref",
    22: "Highest CFD Risk Score",
}
# reference and non reference columns
COL_BOTH = [
    "Highest_CFD_Strand",
    "Chromosome",
    "Highest_CFD_start_coordinate",
    "Highest_CFD_aligned_spacer+PAM",
    "Highest_CFD_aligned_protospacer+PAM_REF",
    "Highest_CFD_aligned_protospacer+PAM_ALT",
    "Highest_CFD_mismatches",
    "Highest_CFD_bulges",
    "Highest_CFD_mismatches+bulges",
    "Highest_CFD_bulge_type",
    "Highest_CFD_PAM_gen",
    "Highest_CFD_score",
    "Highest_CFD_score_REF",
    "Highest_CFD_risk_score",
    "Not_found_in_REF",
    "Highest_CFD_variant_info_genome",
    "Highest_CFD_variant_MAF",
    "Highest_CFD_variant_rsID",
    "Highest_CFD_variant_samples",
    "Other_motifs",
    "Annotation_ENCODE",
]
# reference and non reference column types
COL_BOTH_TYPE = [
    "text",
    "text",
    "numeric",
    "text",
    "text",
    "text",
    "numeric",
    "numeric",
    "numeric",
    "text",
    "text",
    "numeric",
    "numeric",
    "numeric",
    "text",
    "text",
    "numeric",
    "text",
    "text",
    "numeric",
    "text",
]
# reference and non reference column renaming
COL_BOTH_RENAME = {
    0: "Highest_CFD_Strand",
    1: "Chromosome",
    2: "Highest_CFD_start_coordinate",
    3: "Highest_CFD_aligned_spacer+PAM",
    4: "Highest_CFD_aligned_protospacer+PAM_REF",
    5: "Highest_CFD_aligned_protospacer+PAM_ALT",
    6: "Highest_CFD_mismatches",
    7: "Highest_CFD_bulges",
    8: "Highest_CFD_mismatches+bulges",
    9: "Highest_CFD_bulge_type",
    10: "Highest_CFD_PAM_gen",
    11: "Highest_CFD_score",
    12: "Highest_CFD_score_REF",
    13: "Highest_CFD_risk_score",
    14: "Not_found_in_REF",
    15: "Highest_CFD_variant_info_genome",
    16: "Highest_CFD_variant_MAF",
    17: "Highest_CFD_variant_rsID",
    18: "Highest_CFD_variant_samples",
    19: "Other_motifs",
    37: "Annotation_ENCODE",
}
# genome database fields
GENOME_DATABASE = ["Reference", "Enriched", "Samples", "Dictionary", "Annotation"]
# results directory
RESULTS_DIR = "Results"
# assets directory
ASSETS_DIR = "assets"
# annotations directory
ANNOTATIONS_DIR = "Annotations"
# PAMs directory
PAMS_DIR = "PAMs"
# VCFs directory
VCFS_DIR = "VCFs"
# genomes directory
GENOMES_DIR = "Genomes"
# Post-process directory
POSTPROCESS_DIR = "PostProcess"
# Run parameters file
PARAMS_FILE = ".Params.txt"
# Log file
LOG_FILE = "log.txt"
# CRISPR guides file
GUIDES_FILE = ".guides.txt"
# samples files (LIST OF samplesID files, comprehensive samplesID file)
SAMPLES_FILE_LIST = ".samplesID.txt"
SAMPLES_ID_FILE = ".sampleID.txt"
# PAMs file
PAMS_FILE = ".pam.txt"
# email file
EMAIL_FILE = "email.txt"
# queue file
QUEUE_FILE = "queue.txt"
# data directory
DATA_DIR = "data"
# report images directory
IMGS_DIR = "imgs"
# guide column name
GUIDE_COLUMN = "Spacer+PAM"
# chromosome column name
CHR_COLUMN = "Chromosome"
# position column name
POS_COLUMN = "Start_coordinate_(highest_CFD)"
# mismatches column name
MM_COLUMN = "Mismatches_(highest_CFD)"
# bulges column name
BLG_COLUMN = "Bulges_(highest_CFD)"
# total column name
TOTAL_COLUMN = "Mismatches+bulges_(highest_CFD)"
# total column name with fewest_mm+bul
TOTAL_FEWEST_COLUMN = "Mismatches+bulges_(fewest_mm+b)"
# bulge type column name
BLG_T_COLUMN = "Bulge_type_(highest_CFD)"
# CFD score column name
CFD_COLUMN = "CFD_score_(highest_CFD)"
# CRISTA score column name
CRISTA_COLUMN = "CRISTA_score_(highest_CRISTA)"
# CFD risk score column name
RISK_COLUMN = "CFD_risk_score_(highest_CFD)"
# variant samples column name
SAMPLES_COLUMN = "Variant_samples_(highest_CFD)"
# variant CRISTA samples column name
SAMPLES_CRISTA_COLUMN = "Variant_samples_(highest_CRISTA)"
# variant fewest mm+b samples column name
SAMPLES_FEWEST_COLUMN = "Variant_samples_(fewest_mm+b)"
# variant genome CRISTA column name
VARIANTS_CRISTA = "Variant_info_genome_(highest_CRISTA)"
# variant genome CFD columns name
VARIANTS_CFD = "Variant_info_genome_(highest_CFD)"
# variant genome mm+b column name
VARIANTS_FEWEST = "Variant_info_genome_(fewest_mm+b)"
# results filtering criteria
FILTERING_CRITERIA = ["fewest", "CFD", "CRISTA"]
# filter mms + bulges
MMBULGES_FILTER = "fewest_mm+b"
# filter CFD
CFD_FILTER = "highest_CFD"
# filter CRISTA
CRISTA_FILTER = "highest_CRISTA"
# CRISPRme mail subject
MAIL_SUBJECT = "CRISPRme - Job completed"
# CRISPRme mail sender
MAIL_SENDER = "<SENDER OF RESULT MAIL>"
# SSL port (gmail account)
SSL_PORT = 465
# SpCas9 nuclease
CAS9 = "SpCas9"
# pandas series operator methods names
PANDAS_OPERATORS = ("eq", "ne", "lt", "le", "gt", "ge")
# job ID maximum length
JOBID_MAXLEN = 20
# maximum number of iterations to generate job ID
JOBID_ITERATIONS_MAX = 10
# allowed variants datasets (1000 genomes, human diversity project, custom data)
VARIANTS_DATA = ["1000G", "HGDP", "PV"]
# CRISPRme paper link
PAPER_LINK = "https://rdcu.be/c1GYQ"
# CRISPRme github page link
GITHUB_LINK = "https://github.com/pinellolab/CRISPRme"
# manual page image directory
MANUAL_IMGS = "manual_page-images"


def drop_columns(table: pd.DataFrame, filter_criterion: str) -> List[str]:
    """Recover the columns to drop from the result table.
    Empty DataFrames are allowed.
    ...

    Parameters
    ----------
    table : pd.DataFrame
        Results table
    filter_criterion : str
        Table filtering criterion

    Returns
    -------
    List[str]
        Columns to drop
    """

    if not isinstance(table, pd.DataFrame):
        raise TypeError(f"Expected {pd.DataFrame.__name__}, got {type(table).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if filter_criterion not in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filtering criterion ({filter_criterion})")
    drops = []
    if filter_criterion == FILTERING_CRITERIA[0]:  # mms + bulges
        drops = [
            col
            for col in table.columns.tolist()
            if (CFD_FILTER in col or CRISTA_FILTER in col)
        ]
    elif filter_criterion == FILTERING_CRITERIA[1]:  # CFD
        drops = [
            col
            for col in table.columns.tolist()
            if (MMBULGES_FILTER in col or CRISTA_FILTER in col)
        ]
    elif filter_criterion == FILTERING_CRITERIA[2]:  # CRISTA
        drops = [
            col
            for col in table.columns.tolist()
            if (MMBULGES_FILTER in col or CFD_FILTER in col)
        ]
    else:  # we should never go here
        raise ValueError(f"Wrong filtering criterion {filter_criterion}")
    assert bool(drops)
    return drops


def write_json(dropdown_value: str, job_id: str) -> None:
    """Write auxiliary file to keep track of filetring criterion
    when displaying tables in Summary by Mismatches and Bulges and
    Summary by Sample tabs.

    ...

    Parameters
    ----------
    dropdown_value : str
        Table filtering criterion
    job_id : str
        Unique job identifier

    Returns
    -------
    None
    """
    if not isinstance(dropdown_value, str):
        raise TypeError(f"Expected {str.__name__}, got {type(dropdown_value).__name__}")
    dropdown_json_file = os.path.join(current_working_directory, RESULTS_DIR, job_id, ".dropdown.json")
    try:
        handle = open(dropdown_json_file, mode="w")
        handle.write(f"{dropdown_value}")
    except OSError as e:
        raise e
    finally:
        handle.close()


def read_json(job_id: str) -> str:
    """Read the auxiliary file to recover the filtering criterion
    selected by the user with the dropdown.

    ...

    Parameters
    ----------
    job_id : str
        Unique job identifier

    Returns
    -------
    str
        Table filtering criterion
    """
    dropdown_json_file = os.path.join(current_working_directory, RESULTS_DIR, job_id, ".dropdown.json")
    if not os.path.isfile(dropdown_json_file):
        raise FileNotFoundError(f"Unable to locate {dropdown_json_file}")
    try:
        handle = open(dropdown_json_file, mode="r")
        while True:
            line = handle.readline().strip()
            if not line:
                break
            filter_criterion = line
    except OSError as e:
        raise e
    finally:
        handle.close()
    assert filter_criterion in FILTERING_CRITERIA
    return filter_criterion


def get_query_column(filter_criterion: str) -> Dict[str, str]:
    """Recover the names of the columns to display after in
    Summary by Mismatches/Bulges and Summary by Sample tabs.

    ...

    Parameters
    ----------
    filter_criterion : str
        Table filtering criterion

    Returns
    -------
    Dict[str, str]
        Columns to keep in the summary table after filtering
    """

    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if filter_criterion not in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filtering criterion ({filter_criterion})")
    query_columns = {
        "start": "Start_coordinate",
        "mm": "Mismatches",
        "bul": "Bulges",
        "bul_type": "Bulge_type",
        "sort": "",
        "samples": "",
    }
    if filter_criterion == FILTERING_CRITERIA[0]:  # fewest mm+bulges
        for key in query_columns:
            query_columns[key] = "_".join([query_columns[key], f"({MMBULGES_FILTER})"])
            query_columns["sort"] = TOTAL_FEWEST_COLUMN
            query_columns["samples"] = SAMPLES_FEWEST_COLUMN
    elif filter_criterion == FILTERING_CRITERIA[1]:  # cfd
        for key in query_columns:
            query_columns[key] = "_".join([query_columns[key], f"({CFD_FILTER})"])
            query_columns["sort"] = CFD_COLUMN
            query_columns["samples"] = SAMPLES_COLUMN
    elif filter_criterion == FILTERING_CRITERIA[2]:  # crista
        for key in query_columns:
            query_columns[key] = "_".join([query_columns[key], f"({CRISTA_FILTER})"])
            query_columns["sort"] = CRISTA_COLUMN
            query_columns["samples"] = SAMPLES_CRISTA_COLUMN
    else:
        raise ValueError
    return query_columns


def split_filter_part(filter_part: str) -> Tuple[str, str, str]:
    """Split the data table filter in its parts.

    ...

    Parameters
    ----------
    filter_part : str
        Filter

    Returns
    -------
    Tuple[str, str, str]
        Filter fields
    """

    if not isinstance(filter_part, str):
        raise TypeError(f"Expected {str.__name__}, got {type(filter_part).__name__}")
    for operator_type in operators:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                name = name_part[(name_part.find("{") + 1) : name_part.rfind("}")]
                value_part = value_part.strip()
                v0 = value_part[0]
                if v0 == value_part[-1] and v0 in ("'", '"', "`"):
                    value = value_part[1:-1].replace(("\\" + v0), v0)
                else:
                    try:
                        value = float(value_part)
                    except ValueError:
                        value = value_part
                # word operators need spaces after them in the filter string,
                # but we don't want these later
                return name, operator_type[0].strip(), value
    return [None] * 3


def generate_table(
    dataframe: pd.DataFrame,
    id_table: str,
    genome_type: str,
    guide: Optional[str] = "",
    job_id: Optional[str] = "",
    max_rows: Optional[int] = 2600,
) -> html.Table:
    """Generate a html table from a given pandas DataFrame.

    ...

    Parameters
    ----------
    dataframe : pd.DataFrame
        Input dataframe
    id_table : str
        HTML table identifier
    genome_type: str
        Genome type
    guide : str
        Guide
    job_id : str
        Unique job identifier
    max_rows : int
        Maximum number of rows to display

    Returns
    -------
    html.Table
        HTML table
    """

    if not isinstance(dataframe, pd.DataFrame):
        raise TypeError(
            f"Expected {type(pd.DataFrame).__name__}, got {type(dataframe).__name__}"
        )
    if not isinstance(id_table, str):
        raise TypeError(f"Expected {str.__name__}, got {type(id_table).__name__}")
    if not isinstance(genome_type, str):
        raise TypeError(f"Expected {str.__name__}, got {type(genome_type).__name__}")
    if not isinstance(guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide).__name__}")
    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if not isinstance(max_rows, int):
        raise TypeError(f"Expected {int.__name__}, got {type(max_rows).__name__}")
    # build table header
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
                html.Th(x, style={"vertical-align": "middle", "text-align": "center"})
                for x in ["Reference", "Variant", "Combined"]
            ]
        )
    )
    # add body to html table
    table_html = html.Table(
        header +
        # append body
        [
            html.Tr(
                [
                    (
                        html.Td(
                            html.A(
                                dataframe.loc[i, col],
                                href="".join(
                                    [
                                        "result?job=",
                                        f"{job_id}#{guide}new",
                                        dataframe.loc[i, "Bulge Type"],
                                        str(dataframe.loc[i, "Bulge Size"]),
                                        str(dataframe.loc[i, "Mismatches"]),
                                    ]
                                ),
                                target="_blank",
                            ),
                            style={"vertical-align": "middle", "text-align": "center"},
                        )
                        if col == ""
                        else html.Td(
                            dataframe.iloc[i][col],
                            style={"vertical-align": "middle", "text-align": "center"},
                        )
                    )
                    for col in dataframe.columns
                ]
            )
            for i in range(min(dataframe.shape[0], max_rows))
        ],
        style={"display": "inline-block"},
        id=id_table,
    )
    return table_html


def generate_table_samples(
    dataframe: pd.DataFrame,
    id_table: str,
    page: int,
    guide: Optional[str] = "",
    job_id: Optional[str] = "",
    max_rows: Optional[int] = 10,
) -> html.Table:
    """Generate a html table from a given pandas DataFrame.

    The table will be displayed when selecting the targets for a specific
    sample.

    ...

    Parameters
    ----------
    dataframe : pd.DataFrame
        Input dataframe
    id_table : str
        HTML table identifier
    page : int
        Current webpage
    guide : str
        Guide
    job_id : str
        Unique job identifier
    max_rows : int
        Maximum number of rows to display

    Returns
    -------
    html.Table
        HTML table
    """

    if not isinstance(dataframe, pd.DataFrame):
        raise TypeError(
            f"Expected {type(pd.DataFrame).__name__}, got {type(dataframe).__name__}"
        )
    if not isinstance(id_table, str):
        raise TypeError(f"Expected {str.__name__}, got {type(id_table).__name__}")
    if not isinstance(page, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page).__name__}")
    if not isinstance(guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide).__name__}")
    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if not isinstance(max_rows, int):
        raise TypeError(f"Expected {int.__name__}, got {type(max_rows).__name__}")
    if max_rows < 1:
        raise ValueError(f"Forbidden number of rows to display ({max_rows})")
    # force dataframe fields to be of str type
    dataframe = dataframe.astype(str)
    rows_remaining = len(dataframe) - (page - 1) * max_rows
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] +
        # Body
        [
            html.Tr(
                [
                    (
                        html.Td(
                            html.A(
                                dataframe.iloc[i + (page - 1) * max_rows][col],
                                href="".join(
                                    [
                                        "result?job=",
                                        job_id,
                                        "#",
                                        guide,
                                        "-Sample-",
                                        dataframe.iloc[i + (page - 1) * max_rows][
                                            "Sample"
                                        ],
                                    ]
                                ),
                                target="_blank",
                            )
                        )
                        if col == ""
                        else html.Td(dataframe.iloc[i + (page - 1) * max_rows][col])
                    )
                    for col in dataframe.columns
                ]
            )
            for i in range(min(rows_remaining, max_rows))
        ],
        style={"display": "inline-block"},
        id=id_table,
    )


def generate_table_position(
    dataframe: pd.DataFrame,
    id_table: str,
    page: int,
    mms: int,
    bulges: int,
    guide: Optional[str] = "",
    job_id: Optional[str] = "",
    max_rows: Optional[int] = 10,
):
    """Generate a html table from a given pandas DataFrame.

    The table will be displayed when selecting the targets found in a
    specific genomic region.

    ...

    Parameters
    ----------
    dataframe : pd.DataFrame
        Input dataframe
    id_table : str
        HTML table identifier
    page : int
        Current page
    mms : int
        Mismatches
    bulges : int
        Bulges
    guide : str
        Guide
    job_id : str
        Unique job identifier
    max_rows : int
        Maximum number of rows to display

    Returns
    -------
    html.Table
        HTML table
    """

    if not isinstance(dataframe, pd.DataFrame):
        raise TypeError(
            f"Expected {type(pd.DataFrame).__name__}, got {type(dataframe).__name__}"
        )
    if not isinstance(id_table, str):
        raise TypeError(f"Expected {str.__name__}, got {type(id_table).__name__}")
    if not isinstance(page, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page).__name__}")
    if not isinstance(mms, int):
        raise TypeError(f"Expected {int.__name__}, got {type(mms).__name__}")
    if not isinstance(bulges, int):
        raise TypeError(f"Expected {int.__name__}, got {type(bulges).__name__}")
    if not isinstance(guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide).__name__}")
    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if not isinstance(max_rows, int):
        raise TypeError(f"Expected {int.__name__}, got {type(max_rows).__name__}")
    if max_rows < 1:
        raise ValueError(f"Forbidden number of rows to display ({max_rows})")
    rows_remaining = dataframe.shape[0] - (page - 1) * max_rows
    # build table header
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
    # add mismatches to header
    mms_header = []
    for mm in range(mms + 1):
        mms_header.append(
            html.Th(
                f"{mm} MM",
                style={"vertical-align": "middle", "text-align": "center"},
            )
        )
    header.append(html.Tr(mms_header))
    # build table body
    data = []
    for i in range(min(rows_remaining, max_rows)):
        first_cells = [
            html.Td(
                dataframe.loc[(i + (page - 1) * max_rows), "Chromosome"],
                rowSpan=str(bulges + 1),
                style={"vertical-align": "middle", "text-align": "center"},
            ),
            html.Td(
                dataframe.loc[(i + (page - 1) * max_rows), "Position"],
                rowSpan=str(bulges + 1),
                style={"vertical-align": "middle", "text-align": "center"},
            ),
            html.Td(
                dataframe.loc[(i + (page - 1) * max_rows), "Best Target"],
                rowSpan=str(bulges + 1),
                style={"vertical-align": "middle", "text-align": "center"},
            ),
            html.Td(
                dataframe.loc[(i + (page - 1) * max_rows), "Min Mismatch"],
                rowSpan=str(bulges + 1),
                style={"vertical-align": "middle", "text-align": "center"},
            ),
            html.Td(
                dataframe.loc[(i + (page - 1) * max_rows), "Min Bulge"],
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
                dataframe.loc[(i + (page - 1) * max_rows), col],
                style={"vertical-align": "middle", "text-align": "center"},
            )
            for col in dataframe.columns[5 : 5 + mms + 1]
        ]
        data.append(
            html.Tr(
                first_cells
                + mm_cells
                + [
                    html.Td(
                        html.A(
                            "Show Targets",
                            href="".join(
                                [
                                    "result?job=",
                                    f"{job_id}#{guide}-Pos-",
                                    str(
                                        dataframe.loc[
                                            (i + (page - 1) * max_rows), "Chromosome"
                                        ]
                                    ),
                                    "-",
                                    str(
                                        dataframe.loc[
                                            (i + (page - 1) * max_rows), "Position"
                                        ]
                                    ),
                                ]
                            ),
                            target="_blank",
                        ),
                        rowSpan=str(bulges + 1),
                        style={"vertical-align": "middle", "text-align": "center"},
                    )
                ]
            )
        )
        for b in range(bulges):
            data.append(
                html.Tr(
                    [
                        html.Th(
                            f"{b + 1} Bulge",
                            style={"vertical-align": "middle", "text-align": "center"},
                        )
                    ]
                    + [
                        html.Td(dataframe.loc[(i + (page - 1) * max_rows), col])
                        for col in dataframe.columns[
                            5 + (b + 1) * (mms + 1) : 5 + (b + 1) * (mms + 1) + mms + 1
                        ]
                    ]
                )
            )
    return html.Table(header + data, style={"display": "inline-block"}, id=id_table)


def parse_contents(contents: str) -> bytearray:
    """Read the content of uploaded files and encode it into bits.

    ...

    Parameters
    ---------
    contents : str
        Contents to encode

    Returns
    -------
    bytearray
        byte-like object
    """

    if not isinstance(contents, str):
        raise TypeError(f"Expected {str.__name__}, got {type(contents).__name__}")
    content_type, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)  # decode data
    return decoded


def select_same_len_guides(guides: str) -> str:
    """If the user provides guides of different lengths, compute the length of
    the first given guide and keep only those with the same length.

    ...

    Parameters
    ----------
    guides : str
        Guides

    Returns
    -------
    str
        Selected guides
    """

    if not isinstance(guides, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guides).__name__}")
    length = len(guides.split("\n")[0])
    same_len_guides = [guide for guide in guides.split("\n") if len(guide) == length]
    same_len_guides = "\n".join(same_len_guides).strip()
    return same_len_guides


def friendly_pam_label(raw: str) -> str:
    """Human-readable label for a PAM/nuclease token.

    PAM filenames look like ``<len>bp-<motif>-<enzyme>`` (e.g. ``20bp-NGG-SpCas9``,
    ``23bp-TTTV-Cas12a``, ``20bp-NNN-NO-PAM``). Since the Cas-protein selector was
    removed, the PAM dropdown must be self-describing: show the enzyme and the PAM
    motif together (e.g. "SpCas9 · NGG"), with the pamless case spelled out. The
    dropdown *value* stays the raw token, so search behaviour is unaffected. Bare
    tokens (e.g. a plain "NO"/"SpCas9") are still handled gracefully.
    """

    if raw in ("NO", "NO-PAM"):
        return "No PAM (search all sites)"
    parts = raw.split("-")
    if len(parts) >= 3 and parts[0].endswith("bp"):
        motif = parts[1]
        enzyme = "-".join(parts[2:])
        if enzyme in ("NO-PAM", "NO") or motif == "N" * len(motif):
            return f"No PAM — search all sites ({motif})"
        return f"{enzyme} · {motif}"
    if raw.endswith("-NO-PAM"):
        prefix = raw.split("-")[0]  # e.g. "20bp"
        return f"{prefix} - No PAM (search all sites)"
    return raw


def get_available_PAM() -> List:
    """Recover the PAMs currently available in the /PAMs directory.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    List
        Available PAM files
    """

    pams_files = [
        f
        for f in os.listdir(os.path.join(current_working_directory, PAMS_DIR))
        if (
            not f.startswith(".")  # ignore hidden files
            and os.path.isfile(os.path.join(current_working_directory, PAMS_DIR, f))
        )
    ]
    # remove '.txt' from filenames
    pams_files = [f.replace(".txt", "") for f in pams_files]
    # skip temporary PAMs (used during dictionary updating). The label is
    # user-friendly; the value stays the raw filename so searches are unaffected.
    pams = [
        {"label": friendly_pam_label(pam), "value": pam}
        for pam in pams_files
        if "tempPAM" not in pam
    ]
    return pams


def get_available_CAS() -> List:
    """Recover the Cas proteins currently available in the /PAMs directory.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    List
        Availbale Cas proteins
    """

    cas_files = [
        f
        for f in os.listdir(os.path.join(current_working_directory, PAMS_DIR))
        if (
            not f.startswith(".")  # ignore hidden files
            and os.path.isfile(os.path.join(current_working_directory, PAMS_DIR, f))
        )
    ]
    # removed .txt from filenames
    cas_files = [f.replace(".txt", "") for f in cas_files]
    # skip temporary PAMs (used during dictionary updating)
    casprots = sorted(
        casprot.split(".")[0].split("-")[2]
        for casprot in cas_files
        if "tempPAM" not in casprot
    )
    # Collapse case-variant duplicates: the distributed PAMs can ship the same
    # nuclease twice differing only in case (e.g. 20bp-NGG-SpCas9.txt and
    # 20bp-NGG-spCas9.txt), which otherwise shows the nuclease twice in the
    # dropdown. Keep one spelling per nuclease (uppercase sorts first, so the
    # canonical "SpCas9" wins over "spCas9").
    seen = {}
    for casprot in casprots:
        seen.setdefault(casprot.lower(), casprot)
    casprots_data = [
        {"label": friendly_pam_label(casprot), "value": casprot}
        for casprot in sorted(seen.values(), key=str.lower)
    ]
    return casprots_data


_BUILTIN_VCF_DATASETS = ("1000G", "HGDP", "hg38_1000G", "hg38_HGDP")


def vcf_reference_genome(dataset: str) -> Optional[str]:
    """Best-effort reference genome a VCF dataset belongs to.

    A VCF is in a specific reference genome's coordinates and must only be paired
    with that genome. Determined by, in order: (1) an explicit marker file
    ``VCFs/<dataset>/.reference_genome`` written when the dataset is added; (2) an
    enriched genome ``Genomes/<G>+<dataset>`` already built for some installed G;
    (3) the built-in convention (1000G/HGDP are hg38); (4) a ``<G>_...`` dataset
    name prefix matching an installed genome. Returns None if it cannot be
    determined.
    """
    cwd = current_working_directory
    # marker lives *beside* the dataset dir (not inside), so it also covers
    # server folders registered via a symlink
    marker = os.path.join(cwd, VCFS_DIR, f".{dataset}.refgenome")
    if os.path.isfile(marker):
        try:
            with open(marker) as fh:
                g = fh.read().strip()
            if g:
                return g
        except OSError:
            pass
    genomes_dir = os.path.join(cwd, GENOMES_DIR)
    if os.path.isdir(genomes_dir):
        suffix = f"+{dataset}"
        for d in os.listdir(genomes_dir):
            if d.endswith(suffix) and not d.endswith("_INDELS"):
                return d[: -len(suffix)]
    if dataset in _BUILTIN_VCF_DATASETS:
        return "hg38"
    installed = {g["value"].replace(" ", "_") for g in get_available_genomes()}
    for g in installed:
        if dataset.startswith(g + "_"):
            return g
    return None


def get_custom_VCF(genome_value: str) -> List:
    """Recover the user's custom VCF datasets **that match the selected genome**.

    A VCF is only valid with the reference genome it was called against, so this
    lists custom datasets whose reference genome (:func:`vcf_reference_genome`)
    equals ``genome_value``. Datasets whose reference cannot be determined are
    still shown but flagged, so the user can verify. The built-in 1000G/HGDP
    datasets are excluded here (offered via the variant checklist instead).
    """
    if genome_value is not None:
        if not isinstance(genome_value, str):
            raise TypeError(
                f"Expected {str.__name__}, got {type(genome_value).__name__}"
            )
    genome_value = (genome_value or "").replace(" ", "_")
    vcf_dirs = [
        d
        for d in os.listdir(os.path.join(current_working_directory, VCFS_DIR))
        if (
            not d.startswith(".")  # ignore hidden directories
            and os.path.isdir(os.path.join(current_working_directory, VCFS_DIR, d))
        )
    ]
    vcfs = []
    for d in vcf_dirs:
        if d in _BUILTIN_VCF_DATASETS or "None" in d:
            continue
        ref = vcf_reference_genome(d)
        if ref is None:  # unknown reference -> show but flag for verification
            vcfs.append({"label": f"{d} (reference genome unverified)", "value": d})
        elif ref == genome_value:  # matches the selected genome
            vcfs.append({"label": d, "value": d})
        # else: belongs to a different genome -> hide (prevents wrong pairing)
    return vcfs


def get_available_genomes() -> List:
    """Recover genomes available in the /Genomes directory.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    List
        Available genomes
    """

    genomes = [
        d
        for d in os.listdir(os.path.join(current_working_directory, GENOMES_DIR))
        if os.path.isdir(os.path.join(current_working_directory, GENOMES_DIR, d))
    ]
    genomes = [g.replace("_", " ") for g in genomes]
    genomes_dirs = [
        {"label": d, "value": d} for d in genomes if ("+" not in d and "None" not in d)
    ]
    return genomes_dirs


def get_available_indexes() -> List:
    """Recover the USER-FACING precomputed CRISPRitz indexes under genome_library/.

    Only the logical variant/reference indexes are returned — NOT:
      * ``*_INDELS`` companions (the indel-genome index is an internal part of a
        variant index, built + used automatically; never a standalone choice), and
      * symlink aliases (naming-reconciliation links point at an index already
        listed — dedupe by real path so each index appears once).
    Labels are friendly (parsed from ``<motif>_<bMax+1>_<genome>[+<dataset>]``) so
    the list reads in the same vocabulary as the search form.

    Returns
    -------
    List
        ``{"label", "value"}`` dicts, one per real index (value = directory name).
    """

    idx_root = os.path.join(current_working_directory, "genome_library")
    if not os.path.isdir(idx_root):
        return []
    seen_real = set()
    indexes = []
    for d in sorted(os.listdir(idx_root)):
        p = os.path.join(idx_root, d)
        if d.startswith(".") or not os.path.isdir(p):
            continue
        if d.endswith("_INDELS"):  # internal indel companion — auto-used, not a choice
            continue
        real = os.path.realpath(p)  # collapse symlink aliases onto their target
        if real in seen_real:
            continue
        seen_real.add(real)
        # custom name if the user set one (general: any index may carry a
        # .display_label sidecar), else the parsed convention label
        custom = ""
        lf = os.path.join(p, ".display_label")
        if os.path.isfile(lf):
            try:
                custom = open(lf).read().strip()
            except OSError:
                custom = ""
        indexes.append({"label": custom or _friendly_index_label(d), "value": d})
    return indexes


def _friendly_index_label(name: str) -> str:
    """Human label for a genome_library index dir ``<motif>_<bMax+1>_<genome>[+<dataset>]``.

    e.g. ``NGG_3_hg38+hg38_1000G_HGDP`` -> "NGG · hg38 + 1000G_HGDP (≤2 bulges)";
    ``NNN_3_hg38+...`` adds "pamless — serves any PAM". Falls back to the raw name.
    """
    parts = name.split("_", 2)
    if len(parts) < 3 or not parts[1].isdigit():
        return name
    motif, nstr, rest = parts[0], parts[1], parts[2]
    bulges = int(nstr) - 1
    genome, dataset = (rest.split("+", 1) + [""])[:2]
    if dataset.startswith(genome + "_"):  # strip redundant genome prefix on the dataset
        dataset = dataset[len(genome) + 1:]
    label = motif
    if motif and set(motif) == {"N"}:
        label += " (pamless — serves any PAM)"
    label += f" · {genome}"
    if dataset:
        label += f" + {dataset}"
    label += f" (≤{bulges} bulge{'s' if bulges != 1 else ''})"
    return label


def pam_motif(pam_value: str) -> Optional[str]:
    """The PAM motif (e.g. 'NGG') for a PAM dropdown value, matching how index
    folders are named. Reads PAMs/<value>.txt and extracts the PAM substring the
    same way ``build-index-only`` / ``index-genome`` do. Returns None on error.
    """
    try:
        with open(os.path.join(current_working_directory, PAMS_DIR, pam_value + ".txt")) as fh:
            line = fh.readline()
        seq = line.split()[0]
        pos = int(line.split()[1])
        n = abs(pos)
        return seq[:n] if pos < 0 else seq[-n:]
    except (OSError, ValueError, IndexError):
        return None


def index_max_bulges(genome: str, pam_value: str, vcf: Optional[str] = None) -> int:
    """Max bulges an existing TST index supports for this genome+PAM(+VCF).

    Index folders are named ``<motif>_<N>_<genome>[+<vcf>]`` where N = bMax+1, so
    the usable bulge count is N-1. Returns 0 if no matching index exists (only a
    0-bulge, index-free search is possible). A variant search also needs the
    variant index, so callers should take the min with the reference result.
    """
    motif = pam_motif(pam_value)
    if not motif or not genome:
        return 0
    genome = genome.replace(" ", "_")
    lib = os.path.join(current_working_directory, "genome_library")
    if not os.path.isdir(lib):
        return 0
    tail = f"_{genome}+{vcf}" if vcf else f"_{genome}"
    best = 0
    for d in os.listdir(lib):
        if not os.path.isdir(os.path.join(lib, d)) or d.endswith("_INDELS"):
            continue
        if vcf:
            if not d.endswith(tail):
                continue
        else:
            if "+" in d or not d.endswith(tail):
                continue
        head = d[: -len(tail)]  # '<motif>_<N>'
        parts = head.rsplit("_", 1)
        if len(parts) != 2 or parts[0] != motif or not parts[1].isdigit():
            continue
        best = max(best, int(parts[1]))
    return max(0, best - 1)


def has_variant_index(genome: str, dataset: str) -> bool:
    """True if a variant TST index exists for genome+dataset (bulge-search ready)."""
    genome = (genome or "").replace(" ", "_")
    lib = os.path.join(current_working_directory, "genome_library")
    if not os.path.isdir(lib):
        return False
    marker = f"_{genome}+{dataset}"
    return any(
        marker in d
        and not d.endswith("_INDELS")
        and os.path.isdir(os.path.join(lib, d))
        for d in os.listdir(lib)
    )


def variant_dataset_present(genome: str, dataset: str) -> bool:
    """True if a variant dataset's data exists for a genome: the VCF folder, an
    enriched genome, or a built variant index. Used to only offer variant options
    that are actually usable for the selected genome."""
    cwd = current_working_directory
    genome = (genome or "").replace(" ", "_")
    if os.path.isdir(os.path.join(cwd, VCFS_DIR, dataset)) or os.path.isdir(
        os.path.join(cwd, VCFS_DIR, f"{genome}_{dataset}")
    ):
        return True
    if os.path.isdir(os.path.join(cwd, GENOMES_DIR, f"{genome}+{dataset}")):
        return True
    return has_variant_index(genome, dataset)


def get_all_vcf_datasets() -> List:
    """List every VCF dataset directory under VCFs/ (no genome filtering).

    Unlike :func:`get_custom_VCF` (which hides the built-in datasets and those
    matching the selected genome), this lists *all* installed datasets for the
    Settings data-manager table.
    """

    vcf_root = os.path.join(current_working_directory, VCFS_DIR)
    if not os.path.isdir(vcf_root):
        return []
    dirs = [
        d
        for d in os.listdir(vcf_root)
        if (
            not d.startswith(".")
            and os.path.isdir(os.path.join(vcf_root, d))
            and "None" not in d
        )
    ]
    return [{"label": d, "value": d} for d in sorted(dirs)]


def get_variant_dataset_options(genome: str) -> List:
    """Variant-dataset dropdown options for a selected genome.

    A single self-describing list driven by what is installed for THAT genome:
    "Reference only" is always first; then each built-in dataset present for the
    genome (1000G, HGDP), and a combined entry when more than one is available
    (e.g. "1000G + HGDP", value "1000G+HGDP"). For a genome with no variant data
    (e.g. a freshly added non-human assembly) only "Reference only" is offered, so
    the control never shows datasets that cannot apply. The combined value is
    split on '+' by the search wiring back into the individual datasets.
    """
    genome_norm = (genome or "").replace(" ", "_")
    options = [{"label": "Reference only (no variants)", "value": "ref"}]
    # Each option is ONE dataset = ONE enriched-genome/index = ONE search. A combined
    # panel (e.g. 1000G+HGDP) must be a single *merged* VCF dataset, not two datasets
    # searched separately and merged after -- the merged VCF is both what the shipped
    # combined index is built from and the only way cross-dataset haplotypes (a target
    # created by a 1000G variant next to an HGDP variant) are found. So we simply list
    # the dataset folders that exist for the genome; a merged panel appears once it is
    # registered as its own dataset (no synthetic "A+B" two-run option).
    labels = {"1000G": "1000 Genomes Project (1000G)", "HGDP": "HGDP"}
    for ds in ("1000G", "HGDP"):
        if variant_dataset_present(genome_norm, ds):
            options.append({"label": labels.get(ds, ds), "value": ds})
    # any additional installed dataset folder for this genome (e.g. a merged panel or
    # a custom VCF registered via Settings), matched by the genome token in its name
    seen = {o["value"] for o in options}
    for ds in get_all_vcf_datasets():
        val = ds["value"] if isinstance(ds, dict) else ds
        core = val[len(genome_norm) + 1:] if val.startswith(genome_norm + "_") else val
        if genome_norm and genome_norm in val and core not in seen:
            options.append({"label": core, "value": core})
            seen.add(core)
    return options


def get_annotation_options(genome: str) -> List:
    """Annotation dropdown options for a selected genome.

    Genome-driven, mirroring the variant selector: "No annotation" is always first;
    the built-in ENCODE cCREs + GENCODE bundle is an hg38 resource so it is offered
    only for hg38; and any installed custom annotation whose filename carries the
    genome token (e.g. ``...susScr11.bed``) is listed for that genome. Annotations
    that cannot apply to the selected genome are never shown.
    """
    genome_norm = (genome or "").replace(" ", "_")
    options = [{"label": "No annotation", "value": "none"}]
    if genome_norm == "hg38" and os.path.isfile(
        os.path.join(current_working_directory, ANNOTATIONS_DIR, "dhs+encode+gencode.hg38.bed")
    ):
        options.append({"label": "ENCODE cCREs + GENCODE genes (hg38)", "value": "EN"})
    for ann in get_custom_annotations():
        val = ann["value"] if isinstance(ann, dict) else ann
        # match the genome token in the filename; skip the built-in already covered
        if genome_norm and genome_norm in val and "encode" not in val.lower():
            options.append({"label": val, "value": val})
    return options


def get_pam_options(genome: str, variant_choice: Optional[str]) -> List:
    """PAM dropdown options for the current genome + variant selection.

    Reference-only searches can use ANY installed PAM (a PAM lacking a per-PAM
    index just runs index-free via ``-r``), so all PAMs are offered. When a variant
    dataset is included, a variant search needs a variant index for that PAM, so the
    list is restricted to PAMs that actually have one for genome+dataset -- except a
    pamless ``NNN`` variant index, which serves every PAM (the requested PAM is
    enforced by the post-search filter), so its presence unlocks all PAMs again.
    """
    all_pams = get_available_PAM()
    if variant_choice in (None, "", "ref"):
        return all_pams
    genome_norm = (genome or "").replace(" ", "_")
    datasets = [d for d in str(variant_choice).split("+") if d]
    lib = os.path.join(current_working_directory, "genome_library")
    if not os.path.isdir(lib) or not datasets:
        return all_pams
    # collect variant-index motifs available for EVERY selected dataset
    per_dataset_motifs = []
    for ds in datasets:
        # An index is named <motif>_<N>_<genome>+<vcf_folder>. The dropdown value
        # is the genome-stripped folder (e.g. "1000G_HGDP" for VCFs/hg38_1000G_HGDP,
        # "1000G" for VCFs/hg38_1000G), so match BOTH the short form
        # (<genome>+<ds>, e.g. hg38+1000G) and the full-folder form
        # (<genome>+<genome>_<ds>, e.g. hg38+hg38_1000G_HGDP) — otherwise a combined
        # panel's index (built under the full name) is missed and the PAM list
        # wrongly falls back to "all".
        # The vcf_folder is the SUFFIX after '+', so match on endswith (not
        # substring) — else "1000G" spuriously matches "1000G_HGDP" indexes and a
        # 1000G-alone search wrongly inherits the combined panel's pamless (=all
        # PAMs). Accept both the genome-stripped dropdown value ("+1000G") and the
        # full-folder form ("+hg38_1000G_HGDP").
        suffixes = (f"+{ds}", f"+{genome_norm}_{ds}")
        motifs = set()
        for d in os.listdir(lib):
            if d.endswith("_INDELS") or not os.path.isdir(os.path.join(lib, d)):
                continue
            if any(d.endswith(sfx) for sfx in suffixes):
                motifs.add(d.split("_")[0])  # <motif>_<N>_<genome>+<vcf_folder>
        per_dataset_motifs.append(motifs)
    usable = set.intersection(*per_dataset_motifs) if per_dataset_motifs else set()
    if any(m == "N" * len(m) for m in usable):  # pamless NNN index -> any PAM
        return all_pams
    return [p for p in all_pams if pam_motif(p["value"]) in usable] or all_pams


def get_custom_annotations() -> List:
    """Recover user's annotation data.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    List
        User's annotation data
    """

    annotation_data = glob(os.path.join(current_working_directory, ANNOTATIONS_DIR, "*.bed"))
    annotations = [
        {"label": ann.strip().split("/")[-1], "value": ann.strip().split("/")[-1]}
        for ann in annotation_data
        if (
            "encode" not in ann
            and "None" not in ann
            and "dummy" not in ann
            and "tmp" not in ann
        )
    ]
    return annotations


def _sort_bed(fname: str, outfname: str) -> str:
    """Sorts a BED file and writes the sorted output to a new file.

    Uses sort-bed to sort the input BED file and saves the result to the specified 
    output file. Raises an error if sorting fails.

    Args:
        fname: Path to the input BED file.
        outfname: Path where the sorted BED file will be written.

    Returns:
        The path to the sorted BED file.

    Raises:
        SystemExit: If sorting fails.
    """
    code = subprocess.call(f"sort-bed {fname} > {outfname}", shell=True)
    if code != 0:
        raise subprocess.SubprocessError("Sorting BED file failed")
    assert os.path.isfile(outfname)
    return outfname

def compress_file(fname: str) -> str:
    """Compresses a file using bgzip and returns the path to the compressed file.

    Uses bgzip to compress the specified file. The operation is idempotent and
    non-destructive: if the plain input file is already missing but the compressed
    ``.gz`` counterpart exists (e.g. because a previous submission already
    compressed a shared annotation file), the existing ``.gz`` is reused instead
    of failing. The source file is kept (``bgzip -k``) so shared annotation files
    are not consumed by the first job that touches them, and leftover ``.tmp``/
    lock cruft from interrupted runs is cleaned up defensively.

    Args:
        fname: Path to the file to be compressed. May already carry a ``.gz``
            suffix, in which case that path is returned unchanged.

    Returns:
        The path to the compressed file with a .gz extension.

    Raises:
        subprocess.SubprocessError: If compression fails.
    """
    # normalize plain/compressed paths so the function can be called with either
    if fname.endswith(".gz"):
        gz_path, plain_path = fname, fname[:-3]
    else:
        gz_path, plain_path = f"{fname}.gz", fname
    # defensively clean up cruft left behind by interrupted runs so a stale
    # temporary or lock file cannot corrupt the shared annotation directory
    for cruft in (f"{plain_path}.tmp", f"{gz_path}.tmp", f"{gz_path}.tbi.lock"):
        if os.path.isfile(cruft):
            subprocess.call(f"rm -f {cruft}", shell=True)
    # idempotent path: the plain file was already compressed by a previous run
    # (bgzip -f had removed it) - reuse the existing .gz instead of failing
    if not os.path.isfile(plain_path):
        if os.path.isfile(gz_path):
            return gz_path
        raise subprocess.SubprocessError(
            f"Cannot compress file, neither {plain_path} nor {gz_path} found"
        )
    # non-destructive compression: keep the source file (-k) and overwrite any
    # stale .gz (-f) so re-submitting the same job is safe
    code = subprocess.call(f"bgzip -k -f {plain_path}", shell=True)
    if code != 0:
        raise subprocess.SubprocessError("Compressing and indexing file failed")
    assert os.path.isfile(gz_path)
    return gz_path


def _mv_file(fname: str, outfname: str) -> str:
    """Renames or moves a file to a new location.

    Uses the mv command to move or rename the specified file. Raises an error if 
    the operation fails.

    Args:
        fname: Path to the source file.
        outfname: Path to the destination file.

    Returns:
        The path to the moved or renamed file.

    Raises:
        SystemExit: If the move or rename operation fails.
    """
    code = subprocess.call(f"mv {fname} {outfname}", shell=True)
    if code != 0:
        raise subprocess.SubprocessError("Renaming file failed")
    assert os.path.isfile(outfname)
    return outfname


def sort_annotation(annotationfile: str) -> str:
    """Sorts, compresses, and replaces a BED annotation file for downstream 
    analysis.

    Decompresses the input annotation file, sorts it, compresses it with bgzip, 
    and replaces the original file. Raises an error if any step fails.

    Args:
        annotationfile: Path to the input annotation file.

    Returns:
        The path to the sorted and compressed annotation file.

    Raises:
        SystemExit: If decompression, sorting, compression, or renaming fails.
    """
    if annotationfile.endswith(".gz"):
        gz_path, plain_path = annotationfile, annotationfile[:-3]
    else:
        gz_path, plain_path = f"{annotationfile}.gz", annotationfile
    decompressed_here = False
    if not os.path.isfile(plain_path):
        if os.path.isfile(gz_path):
            _decompress_file(gz_path, plain_path)
            decompressed_here = True
        else:
            raise subprocess.SubprocessError(
                f"Annotation file not found: {annotationfile}"
            )
    annotationfile_sorted = _sort_bed(plain_path, f"{plain_path}.tmp.sorted.bed")
    annotationfile_sorted_bgzip = compress_file(annotationfile_sorted)
    result = _mv_file(annotationfile_sorted_bgzip, gz_path)
    # clean up the intermediate sorted BED (compress_file keeps its source with
    # -k, so remove it here) to avoid leaving cruft in the shared annotation dir
    if os.path.isfile(annotationfile_sorted):
        subprocess.call(f"rm -f {annotationfile_sorted}", shell=True)
    if decompressed_here and os.path.isfile(plain_path):
        subprocess.call(f"rm {plain_path}", shell=True)
    return result


def _decompress_file(fname: str, outfname: str) -> str:
    """Decompresses a bgzipped file to a specified output file."""
    code = subprocess.call(f"gunzip -k -c {fname} > {outfname}", shell=True)
    if code != 0:
        raise subprocess.SubprocessError("Decompressing annotation file failed")
    assert os.path.isfile(outfname)
    return outfname
