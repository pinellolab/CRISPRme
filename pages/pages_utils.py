"""Define static variables and utilities functions used throughout CRISPRme's
webpages.
"""

from app import operators, WORKINGDIR

from typing import Dict, List, Optional, Tuple
from glob import glob

import dash_html_components as html
import pandas as pd

import base64
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
    dropdown_json_file = os.path.join(WORKINGDIR, RESULTS_DIR, job_id, ".dropdown.json")
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
    dropdown_json_file = os.path.join(WORKINGDIR, RESULTS_DIR, job_id, ".dropdown.json")
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
    if filter_criterion == FILTERING_CRITERIA[0]:
        for key in query_columns.keys():
            query_columns[key] = "_".join([query_columns[key], f"({MMBULGES_FILTER})"])
            query_columns["sort"] = TOTAL_FEWEST_COLUMN
            query_columns["samples"] = SAMPLES_FEWEST_COLUMN
    elif filter_criterion == FILTERING_CRITERIA[1]:
        for key in query_columns.keys():
            query_columns[key] = "_".join([query_columns[key], f"({CFD_FILTER})"])
            query_columns["sort"] = CFD_COLUMN
            query_columns["samples"] = SAMPLES_COLUMN
    elif filter_criterion == FILTERING_CRITERIA[2]:
        for key in query_columns.keys():
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
        for f in os.listdir(os.path.join(WORKINGDIR, PAMS_DIR))
        if (
            not f.startswith(".")  # ignore hidden files
            and os.path.isfile(os.path.join(WORKINGDIR, PAMS_DIR, f))
        )
    ]
    # remove '.txt' from filenames
    pams_files = [f.replace(".txt", "") for f in pams_files]
    # skip temporary PAMs (used during dictionary updating)
    pams = [{"label": pam, "value": pam} for pam in pams_files if "tempPAM" not in pam]
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
        for f in os.listdir(os.path.join(WORKINGDIR, PAMS_DIR))
        if (
            not f.startswith(".")  # ignore hidden files
            and os.path.isfile(os.path.join(WORKINGDIR, PAMS_DIR, f))
        )
    ]
    # removed .txt from filenames
    cas_files = [f.replace(".txt", "") for f in cas_files]
    # skip temporary PAMs (used during dictionary updating)
    casprots = [
        casprot.split(".")[0].split("-")[2]
        for casprot in cas_files
        if "tempPAM" not in casprot
    ]
    # remove potential duplicates
    casprots = set(casprots)
    casprots_data = [
        {"label": casprot, "value": casprot} for casprot in sorted(casprots)
    ]
    return casprots_data


def get_custom_VCF(genome_value: str) -> List:
    """Recover user's VCFs.

    ...

    Paramters
    ---------
    genome_value : str
        Genome

    Returns
    -------
    List
        User's VCFs.
    """

    if genome_value is not None:
        if not isinstance(genome_value, str):
            raise TypeError(
                f"Expected {str.__name__}, got {type(genome_value).__name__}"
            )
    vcf_dirs = [
        d
        for d in os.listdir(os.path.join(WORKINGDIR, VCFS_DIR))
        if (
            not d.startswith(".")  # ignore hidden directories
            and os.path.isdir(os.path.join(WORKINGDIR, VCFS_DIR, d))
        )
    ]
    genome_value = genome_value.replace(" ", "_")
    vcfs = [
        {"label": d, "value": d}
        for d in vcf_dirs
        if (
            "hg38_HGDP" not in d
            and "hg38_1000G" not in d
            and "None" not in d
            and genome_value not in d
        )
    ]
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
        for d in os.listdir(os.path.join(WORKINGDIR, GENOMES_DIR))
        if os.path.isdir(os.path.join(WORKINGDIR, GENOMES_DIR, d))
    ]
    genomes = [g.replace("_", " ") for g in genomes]
    genomes_dirs = [
        {"label": d, "value": d} for d in genomes if ("+" not in d and "None" not in d)
    ]
    return genomes_dirs


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

    annotation_data = glob(os.path.join(WORKINGDIR, ANNOTATIONS_DIR, "*.bed"))
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
