"""Canonical column schema for CRISPRme post-analysis intermediate tables.

DRAFT — foundation for issue #143 item 3 (a shared, single-source-of-truth schema
so the reporting layer stops relying on hardcoded, mutually-inconsistent column
indices). Not yet wired into the call sites; opened for discussion (see #143).

Background — three *different* layouts currently exist for the per-target
intermediate tables, which is the root cause of the `resultIntegrator.py`
`float(rsID)` crash in #143:

    A) new_simple_analysis.py:706  (initial write header, 20 cols)
       ... 15:rsID 16:AF 17:SNP 18:#Seq_in_cluster 19:Reference   (+CFD appended)
    B) new_simple_analysis.py:844  (later sed rewrite, 22 cols)
       ... 15:rsID 16:AF 17:SNP 18:Reference 19:CFD_ref 20:CFD 21:#Seq_in_cluster
    C) observed *.bestCFD_INDEL.txt.tmp header (22 cols)
       ... 16:rsID 17:AF ...   (Reference at position 4, shifting everything after)

`resultIntegrator.py` hardcodes indices that match (C) — target[16]=rsID,
target[17]=MAF, target[20]=CFD, target[21]=CFD_ref — so any row that does not
match (C) exactly (e.g. a multi-rsID `rsA;rsB` dbSNP id shifted by a
`bcftools norm -m-` split) lands a non-numeric value under a float() and crashes.

Proposed resolution (for discussion):
  1. Pick ONE canonical layout and emit it consistently from new_simple_analysis.py
     (drop the divergent sed rewrite), OR
  2. Have every consumer resolve columns BY NAME from the file's own header via
     `columns_from_header()` below, instead of positional indices.

Option (2) is the more robust long-term fix and is what this module enables.
"""

# Canonical per-sub-record column names (one block per scoring system: the
# per-target row concatenates the highest-CFD, fewest-mm+bulges, and CRISTA
# blocks). Names mirror new_simple_analysis.py's header vocabulary.
CANONICAL_COLUMNS = [
    "Bulge_type", "crRNA", "DNA", "Chromosome", "Position", "Cluster_Position",
    "Direction", "Mismatches", "Bulge_Size", "Total", "PAM_gen", "Var_uniq",
    "Samples", "Annotation_Type", "Real_Guide", "rsID", "AF", "SNP",
    "Reference", "CFD_ref", "CFD", "Seq_in_cluster",
]

# Columns the reporting layer parses as floats (MAF lists + CFD scores). Kept
# here so a single edit updates every guard/consumer. These are NAMES; resolve
# to indices per file via columns_from_header().
FLOAT_COLUMN_NAMES = ("AF", "CFD", "CFD_ref")


def columns_from_header(header_line):
    """Return {column_name: index} parsed from an actual header line.

    Lets consumers read columns by name instead of hardcoded positions, which is
    robust to the layout differences above. `header_line` may start with '#'.
    """
    names = header_line.lstrip("#").rstrip("\n").split("\t")
    return {name.lstrip("#").strip(): i for i, name in enumerate(names)}


def float_indices(header_line):
    """Indices (in this file's layout) of the columns that must parse as floats."""
    idx = columns_from_header(header_line)
    return tuple(idx[name] for name in FLOAT_COLUMN_NAMES if name in idx)
