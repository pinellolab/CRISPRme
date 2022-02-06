"""Python script containing static variables used in CRISPRme's result page. 
"""

from typing import Dict, List
from app import current_working_directory

import pandas as pd

import os

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
GENOME_DATABASE = ["Reference", "Enriched",
                   "Samples", "Dictionary", "Annotation"]
# results directory
RESULTS_DIR = "Results"
# data directory
DATA_DIR = "data"
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
# results filtering criteria
FILTERING_CRITERIA = ["fewest", "CFD", "CRISTA"]
# filter mms + bulges
MMBULGES_FILTER = "fewest_mm+b"
# filter CFD
CFD_FILTER = "highest_CFD"
# filter CRISTA
CRISTA_FILTER = "highest_CRISTA"


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
        raise TypeError(
            f"Expected {pd.DataFrame.__name__}, got {type(table).__name__}")
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
        raise TypeError(
            f"Expected {str.__name__}, got {type(dropdown_value).__name__}")
    dropdown_json_file = os.path.join(
        current_working_directory, RESULTS_DIR, job_id, ".dropdown.json"
    )
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
    dropdown_json_file = os.path.join(
        current_working_directory, RESULTS_DIR, job_id, ".dropdown.json"
    )
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
    filter_crietrion : str
        Table filtering criterion
    
    Returns
    -------
    Dict[str, str]
        Columns to keep in the summary table after filtering
    """
    query_columns = {
        "start": "Start_coordinate",
        "mm": "Mismatches",
        "bul": "Bulges",
        "bul_type": "Bulge_type",
        "sort": "",
        "samples": ""
    }
    if filter_criterion == FILTERING_CRITERIA[0]:
        for key in query_columns.keys():
            query_columns[key] = "_".join(
                [query_columns[key], f"({MMBULGES_FILTER})"]
            )
            query_columns['sort'] = TOTAL_FEWEST_COLUMN
            query_columns['samples'] = SAMPLES_FEWEST_COLUMN
    elif filter_criterion == FILTERING_CRITERIA[1]:
        for key in query_columns.keys():
            query_columns[key] = "_".join(
                [query_columns[key], f"({CFD_FILTER})"]
            )
            query_columns['sort'] = CFD_COLUMN
            query_columns['samples'] = SAMPLES_COLUMN
    elif filter_criterion == FILTERING_CRITERIA[2]:
        for key in query_columns.keys():
            query_columns[key] = "_".join(
                [query_columns[key], f"({CRISTA_FILTER})"]
            )
            query_columns['sort'] = CRISTA_COLUMN
            query_columns['samples'] = SAMPLES_CRISTA_COLUMN
    else:
        raise ValueError
    return query_columns
