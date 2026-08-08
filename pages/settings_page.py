"""Settings / Data-Manager page.

A graphical, non-expert way to add reference genomes, precomputed indexes, VCF
datasets, annotations and PAMs to the local CRISPRme data folder — so the search
form (which auto-discovers whatever lives under Genomes/, genome_library/, VCFs/,
Annotations/, PAMs/) picks them up with no command line.

Large data (multi-GB genomes/VCFs) never goes through an in-memory browser
upload: it is fetched server-side (UCSC by assembly, HuggingFace by name, an
explicit URL, or an already-staged server path). Only small files (annotations)
use a direct upload. Long operations run as detached jobs on a dedicated
single-slot executor so they never starve the two search slots, and their
progress is polled from the browser exactly like a search.

All mutations are local-mode only (disabled when app.ONLINE is True). Publishing
an index to the shared HuggingFace repo is maintainer-only (app.MAINTAINER_MODE).
"""

from app import (
    app,
    current_working_directory,
    settings_executor,
    ONLINE,
    MAINTAINER_MODE,
)
from .pages_utils import (
    get_available_genomes,
    get_available_indexes,
    get_available_PAM,
    get_available_CAS,
    get_all_vcf_datasets,
    get_custom_annotations,
)

from dash import Input, Output, State, html, dcc, no_update
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from flask import request

from typing import List, Optional
import subprocess
import tarfile
import string
import random
import base64
import shutil
import gzip
import os

SETTINGS_DIR = "Settings"


# ---------------------------------------------------------------------------
# Chunked / resumable browser upload (large local files)
# ---------------------------------------------------------------------------
# Large genome/VCF files cannot go through dcc.Upload (base64, in-memory). A
# small client-side helper (assets/chunked_upload.js) slices the chosen file and
# POSTs it here chunk-by-chunk; we append each chunk to a staging .part file and,
# on the final chunk, finalize it into the right data folder. Because chunks are
# appended by (name, index) the upload is resumable: re-sending from the last
# received index continues the same .part file.
def _uploads_dir() -> str:
    d = os.path.join(current_working_directory, SETTINGS_DIR, "uploads")
    os.makedirs(d, exist_ok=True)
    return d


def _finalize_upload(
    target: str, part_path: str, name: str, dataset: str, genome: str = ""
) -> str:
    """Move a completed upload into the correct data folder (unpacking FASTA)."""
    if target == "annotation":
        if not name.endswith(".bed"):
            raise ValueError("annotation upload must be a .bed file")
        dest_dir = os.path.join(current_working_directory, "Annotations")
        os.makedirs(dest_dir, exist_ok=True)
        os.replace(part_path, os.path.join(dest_dir, name))
        return name
    if target == "genome":
        assembly = name
        for ext in (".chromFa.tar.gz", ".fa.gz", ".fasta.gz", ".tar.gz", ".fa", ".fasta"):
            if assembly.endswith(ext):
                assembly = assembly[: -len(ext)]
                break
        err = _validate_name(assembly)
        if err:
            raise ValueError(err)
        dest_dir = os.path.join(current_working_directory, "Genomes", assembly)
        os.makedirs(dest_dir, exist_ok=True)
        if name.endswith((".tar.gz", ".tgz")):
            with gzip.open(part_path, "rb") as fin, tarfile.open(fileobj=fin, mode="r") as tar:
                tar.extractall(dest_dir)
            os.remove(part_path)
        elif name.endswith(".gz"):
            with gzip.open(part_path, "rb") as fin, open(
                os.path.join(dest_dir, f"{assembly}.fa"), "wb"
            ) as fout:
                shutil.copyfileobj(fin, fout)
            os.remove(part_path)
        else:
            os.replace(part_path, os.path.join(dest_dir, f"{assembly}.fa"))
        return assembly
    if target == "vcf":
        err = _validate_name(dataset)
        if err:
            raise ValueError(err)
        dest_dir = os.path.join(current_working_directory, "VCFs", dataset)
        os.makedirs(dest_dir, exist_ok=True)
        os.replace(part_path, os.path.join(dest_dir, name))
        _write_vcf_marker(dataset, genome)  # record the reference genome
        return dataset
    raise ValueError(f"unknown upload target {target!r}")


@app.server.route("/settings/upload-chunk", methods=["POST"])
def _settings_upload_chunk():
    if ONLINE:
        return ("data management is disabled on the public server", 403)
    target = request.headers.get("X-Target", "")
    raw_name = request.headers.get("X-File-Name", "")
    dataset = os.path.basename(request.headers.get("X-Dataset", "") or "")
    genome = request.headers.get("X-Genome", "") or ""
    name = os.path.basename(raw_name)
    try:
        idx = int(request.headers.get("X-Chunk-Index", "0"))
        total = int(request.headers.get("X-Total-Chunks", "1"))
    except ValueError:
        return ("bad chunk headers", 400)
    if target not in ("genome", "vcf", "annotation") or not name or ".." in name:
        return ("bad target or file name", 400)
    part = os.path.join(_uploads_dir(), name + ".part")
    with open(part, "wb" if idx == 0 else "ab") as fout:
        fout.write(request.get_data())
    if idx + 1 >= total:  # last chunk -> finalize
        try:
            result = _finalize_upload(target, part, name, dataset, genome)
        except Exception as e:
            if os.path.exists(part):
                os.remove(part)
            return (f"finalize failed: {e}", 500)
        return (f"complete:{result}", 200)
    return ("chunk received", 200)


# ---------------------------------------------------------------------------
# Background-job plumbing
# ---------------------------------------------------------------------------
def _new_job_id(prefix: str = "set") -> str:
    """A short unique job id, checked against existing Settings/ job dirs."""
    settings_root = os.path.join(current_working_directory, SETTINGS_DIR)
    os.makedirs(settings_root, exist_ok=True)
    while True:
        jid = prefix + "_" + "".join(
            random.choice(string.ascii_uppercase + string.digits) for _ in range(10)
        )
        if not os.path.isdir(os.path.join(settings_root, jid)):
            return jid


def _run_settings_job(cmd: str, jobdir: str, stage: str) -> None:
    """Run a data-management shell command, recording stage markers.

    Module-level (picklable) so it can be submitted to the ProcessPoolExecutor.
    Success/failure is decided by the process exit code — NOT by whether stderr
    is non-empty (subprocesses legitimately write warnings to stderr) — and
    recorded as ``<stage>\\tEnd`` / ``<stage>\\tFAILED`` in log.txt for the UI.
    """
    logtxt = os.path.join(jobdir, "log.txt")
    logerr = os.path.join(jobdir, "log_error.txt")
    logverb = os.path.join(jobdir, "log_verbose.txt")
    try:
        with open(logverb, "w") as vout, open(logerr, "w") as eout:
            rc = subprocess.call(cmd, shell=True, stdout=vout, stderr=eout)
        with open(logtxt, "a") as lt:
            lt.write(f"{stage}\tEnd\n" if rc == 0 else f"{stage}\tFAILED\n")
    except Exception as e:  # pragma: no cover - defensive
        with open(logerr, "a") as eout:
            eout.write(f"\n{e}\n")
        with open(logtxt, "a") as lt:
            lt.write(f"{stage}\tFAILED\n")


def launch_settings_job(argv: List[str], stage: str) -> str:
    """Launch ``crisprme.py <argv...> --path <cwd>`` as a detached job.

    Returns the job id; progress is polled by :func:`refresh_settings_job`.
    """
    job_id = _new_job_id()
    jobdir = os.path.join(current_working_directory, SETTINGS_DIR, job_id)
    os.makedirs(jobdir, exist_ok=True)
    with open(os.path.join(jobdir, "log.txt"), "w") as lt:
        lt.write(f"{stage}\tStart\n")
    # crisprme.py is on PATH inside the activated env; append --path so data
    # lands under the app's working directory.
    quoted = " ".join(_shlex_quote(a) for a in argv)
    cmd = f"crisprme.py {quoted} --path {_shlex_quote(current_working_directory)}"
    settings_executor.submit(_run_settings_job, cmd, jobdir, stage)
    return job_id


def _shlex_quote(s: str) -> str:
    import shlex

    return shlex.quote(str(s))


def _validate_name(name: str) -> Optional[str]:
    """Return an error string if a user-supplied data name is unsafe/unusable."""
    if not name or not name.strip():
        return "Please provide a name."
    name = name.strip()
    if name == "None" or any(c in name for c in ("/", "\\", " ", "+", "..")):
        return (
            "Name must not contain spaces, '+', 'None', '..' or path separators "
            "(these break auto-discovery)."
        )
    return None


# Benign stderr lines that must not mask the real error in the failure banner.
_BENIGN_ERR_PATTERNS = (
    "BiopythonWarning",
    "warnings.warn",
    "DtypeWarning",
    "SyntaxWarning",
    "FutureWarning",
    "DeprecationWarning",
    "importing Biopython",
    "This is bad practice",
    "pyproject.toml",
)


def _error_tail(logerr_path: str, n: int = 3) -> str:
    """Most relevant tail of a job's log_error.txt for the failure banner.

    Subprocesses emit benign warnings (Biopython import notes, pandas dtype
    warnings, ...) to stderr; showing only the last few raw lines can bury the
    real error under them. Prefer the last ``n`` NON-benign lines; fall back to
    the raw tail if everything looks benign.
    """
    if not os.path.isfile(logerr_path):
        return ""
    try:
        with open(logerr_path) as eh:
            lines = [ln.strip() for ln in eh.read().splitlines() if ln.strip()]
    except OSError:
        return ""
    if not lines:
        return ""
    meaningful = [ln for ln in lines if not any(p in ln for p in _BENIGN_ERR_PATTERNS)]
    return "  ".join((meaningful or lines)[-n:])


def _write_vcf_marker(dataset: str, genome: str) -> None:
    """Record which reference genome a VCF dataset belongs to (a marker file
    beside the dataset dir), so the search form only pairs it with that genome."""
    if not dataset or not genome:
        return
    vdir = os.path.join(current_working_directory, "VCFs")
    os.makedirs(vdir, exist_ok=True)
    try:
        with open(os.path.join(vdir, f".{dataset}.refgenome"), "w") as fh:
            fh.write(genome.replace(" ", "_"))
    except OSError:
        pass


# ---------------------------------------------------------------------------
# Disk-usage helpers
# ---------------------------------------------------------------------------
def _dir_size(path: str) -> int:
    """Total bytes under a path (following into subdirs, ignoring symlinks)."""
    total = 0
    if not os.path.exists(path):
        return 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            fp = os.path.join(root, f)
            try:
                if not os.path.islink(fp):
                    total += os.path.getsize(fp)
            except OSError:
                pass
    return total


def _fmt_size(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{int(n)} B" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


# HuggingFace catalog (what is available to download), fetched once per app run.
_HF_CATALOG: dict = {}


def _hf_catalog(component: str) -> List[dict]:
    """Cached list of HuggingFace-available items for a component ([] if offline)."""
    if component not in _HF_CATALOG:
        try:
            import sys

            pp = os.path.join(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "PostProcess"
            )
            if pp not in sys.path:
                sys.path.insert(0, pp)
            from crisprme_hf import list_available_downloads

            _HF_CATALOG[component] = list_available_downloads(component)
        except Exception:
            _HF_CATALOG[component] = []
    return _HF_CATALOG[component]


def _hf_options(component: str, installed: List[str]) -> List[dict]:
    """Dropdown options for a HuggingFace component, size-labelled + installed-marked."""
    inst = set(installed)
    opts = []
    for item in _hf_catalog(component):
        name = item["name"]
        label = f"{name} ({_fmt_size(item.get('size', 0))})"
        if name in inst:
            label += " — installed"
        opts.append({"label": label, "value": name})
    return opts


# Rough "space required" hints for the download actions (uncompressed on disk).
_SIZE_HINTS = {
    "genome": "a reference genome is typically ~1-3 GB",
    "index": "a genome-wide bulge index is typically ~5-10 GB",
    "vcf": "a population VCF set (e.g. 1000G) is typically ~15-20 GB",
}


def _render_storage() -> html.Div:
    """Per-category disk used by already-downloaded data, total, and free disk."""
    cats = [
        ("Genomes", "Genomes"),
        ("VCF datasets", "VCFs"),
        ("Indexes", "genome_library"),
        ("Annotations", "Annotations"),
    ]
    rows = []
    total = 0
    for label, d in cats:
        sz = _dir_size(os.path.join(current_working_directory, d))
        total += sz
        rows.append(html.Li(f"{label}: {_fmt_size(sz)}"))
    try:
        free = shutil.disk_usage(current_working_directory).free
        free_str = _fmt_size(free)
    except OSError:
        free_str = "unknown"
    return html.Div(
        [
            html.Ul(rows, style={"margin": "4px 0"}),
            html.Div(f"Total data used: {_fmt_size(total)}", style={"fontWeight": "bold"}),
            html.Div(f"Free disk available: {free_str}", style={"color": "#367"}),
        ],
        style={"fontSize": "0.9em"},
    )


# ---------------------------------------------------------------------------
# Installed-data tables
# ---------------------------------------------------------------------------
def _table(rows: List[str], empty: str) -> html.Div:
    if not rows:
        return html.Div(html.I(empty), style={"color": "#666", "padding": "6px"})
    return html.Ul([html.Li(r) for r in rows], style={"margin": "6px 0"})


def _render_all_tables() -> html.Div:
    genomes = [g["value"] for g in get_available_genomes()]
    indexes = [i["value"] for i in get_available_indexes()]
    vcfs = [v["value"] for v in get_all_vcf_datasets()]
    anns = [a["value"] for a in get_custom_annotations()]
    cas = [c["label"] for c in get_available_CAS()]
    pams = [p["value"] for p in get_available_PAM()]
    return html.Div(
        [
            html.Details(
                [html.Summary(f"Installed genomes ({len(genomes)})"), _table(genomes, "none yet")],
                open=True,
            ),
            html.Details(
                [html.Summary(f"Installed indexes ({len(indexes)})"), _table(indexes, "none yet")]
            ),
            html.Details(
                [html.Summary(f"Installed VCF datasets ({len(vcfs)})"), _table(vcfs, "none yet")]
            ),
            html.Details(
                [html.Summary(f"Installed annotations ({len(anns)})"), _table(anns, "none yet")]
            ),
            html.Details(
                [
                    html.Summary(f"Installed nucleases / PAMs ({len(pams)})"),
                    _table([f"{c}" for c in cas], "none yet"),
                    _table(pams, "none yet"),
                ]
            ),
        ]
    )


def _deletable_options() -> List:
    """All removable data items, tagged by type: value='<type>:<name>'."""
    opts = []
    for g in get_available_genomes():
        opts.append({"label": f"Genome: {g['value']}", "value": f"genome:{g['value']}"})
    for i in get_available_indexes():
        opts.append({"label": f"Index: {i['value']}", "value": f"index:{i['value']}"})
    for v in get_all_vcf_datasets():
        opts.append({"label": f"VCF dataset: {v['value']}", "value": f"vcf:{v['value']}"})
    for a in get_custom_annotations():
        opts.append({"label": f"Annotation: {a['value']}", "value": f"annotation:{a['value']}"})
    return opts


def _delete_path(kind: str, name: str) -> str:
    """Resolve the on-disk path for a '<type>:<name>' deletable item."""
    name = os.path.basename(name)  # never allow traversal
    roots = {
        "genome": ("Genomes", name),
        "index": ("genome_library", name),
        "vcf": ("VCFs", name),
        "annotation": ("Annotations", name),
    }
    if kind not in roots:
        raise ValueError(f"unknown data type {kind!r}")
    sub, leaf = roots[kind]
    return os.path.join(current_working_directory, sub, leaf)


# ---------------------------------------------------------------------------
# Small reusable UI helpers
# ---------------------------------------------------------------------------
def _add_card(title: str, blurb: str, body: List) -> dbc.Card:
    return dbc.Card(
        dbc.CardBody(
            [html.H5(title, className="card-title"), html.P(blurb, style={"color": "#555"})]
            + body
        ),
        style={"margin-bottom": "1rem"},
    )


# ---------------------------------------------------------------------------
# Page layout
# ---------------------------------------------------------------------------
def settings_page() -> List:
    """Build the Settings / Data-Manager layout (or a disabled notice online)."""
    if ONLINE:
        return [
            dbc.Container(
                dbc.Alert(
                    [
                        html.H4("Data management is disabled on the public server"),
                        html.P(
                            "Adding genomes, indexes, VCFs, annotations or PAMs "
                            "mutates the shared data folder, so it is only "
                            "available when you run CRISPRme locally (Docker or a "
                            "local install). See the Docker quickstart."
                        ),
                    ],
                    color="warning",
                ),
                style={"margin-top": "2rem"},
            )
        ]

    genome_sources = [
        {"label": "UCSC (by assembly name)", "value": "ucsc"},
        {"label": "HuggingFace (by assembly name)", "value": "hf"},
        {"label": "Direct URL (.fa.gz / .tar.gz)", "value": "url"},
    ]
    installed_genomes = get_available_genomes()
    installed_pams = get_available_PAM()

    # ---- Genomes -----------------------------------------------------------
    genome_card = _add_card(
        "Add a reference genome",
        "Download a reference assembly straight into your data folder. Example: "
        "the pig genome susScr11 from UCSC. Space required: "
        f"{_SIZE_HINTS['genome']} (check 'Free disk available' on the right).",
        [
            dbc.Row(
                [
                    dbc.Col(
                        dcc.Input(
                            id="genome-assembly-input",
                            placeholder="assembly, e.g. susScr11",
                            type="text",
                            style={"width": "100%"},
                        ),
                        width=4,
                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            id="genome-source-dropdown",
                            options=genome_sources,
                            value="ucsc",
                            clearable=False,
                        ),
                        width=4,
                    ),
                    dbc.Col(html.Button("Download genome", id="genome-add-btn"), width=4),
                ]
            ),
            dcc.Input(
                id="genome-url-input",
                placeholder="download URL (only for source = Direct URL)",
                type="text",
                style={"width": "100%", "margin-top": "0.5rem"},
            ),
            html.Div(
                dcc.Dropdown(
                    id="genome-hf-name",
                    options=_hf_options(
                        "genome", [g["value"] for g in get_available_genomes()]
                    ),
                    placeholder="…or pick a genome available on HuggingFace",
                ),
                style={"margin-top": "0.4rem"},
            ),
            html.Small(
                "Or upload a local genome file (.fa / .fa.gz / .tar.gz) — chunked, "
                "no browser size limit. The assembly name is taken from the file name."
            ),
            html.Div(
                className="crisprme-chunk-upload",
                style={"margin-top": "0.3rem"},
                **{"data-target": "genome"},
            ),
            html.Div(id="genome-feedback", style={"color": "#b00", "margin-top": "0.4rem"}),
        ],
    )

    # ---- Indexes -----------------------------------------------------------
    index_card = _add_card(
        "Add a precomputed index / build one",
        "Bulge-enabled searches need a genome index. Download a ready-made one "
        "from HuggingFace, or build one locally from an installed genome + PAM. "
        f"Space required: {_SIZE_HINTS['index']}.",
        [
            html.B("Download a prebuilt index (HuggingFace)"),
            dbc.Row(
                [
                    dbc.Col(
                        dcc.Dropdown(
                            id="index-hf-name",
                            options=_hf_options(
                                "index", [i["value"] for i in get_available_indexes()]
                            ),
                            placeholder="pick an index available on HuggingFace",
                        ),
                        width=8,
                    ),
                    dbc.Col(html.Button("Download index", id="index-hf-btn"), width=4),
                ]
            ),
            html.Hr(),
            html.B("Build an index locally"),
            dbc.Row(
                dbc.Col(
                    dcc.Dropdown(
                        id="index-build-vcf",
                        options=[{"label": "(reference only — no variants)", "value": ""}]
                        + get_all_vcf_datasets(),
                        value="",
                        clearable=False,
                    ),
                    width=6,
                ),
                style={"margin-bottom": "0.4rem"},
            ),
            html.Small(
                "Pick a VCF dataset to pre-build a variant-aware index (genome "
                "enrichment + SNP/indels indexing) so the first variant search is "
                "fast — this build is slower and uses more disk."
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.Small("Genome"),
                            dcc.Dropdown(
                                id="index-build-genome",
                                options=installed_genomes,
                                placeholder="genome",
                            ),
                        ],
                        width=3,
                    ),
                    dbc.Col(
                        [
                            html.Small("PAM / nuclease"),
                            dcc.Dropdown(
                                id="index-build-pam",
                                options=installed_pams,
                                placeholder="PAM",
                            ),
                        ],
                        width=3,
                    ),
                    dbc.Col(
                        [
                            html.Small("DNA bulges"),
                            dcc.Input(
                                id="index-build-bdna",
                                type="number",
                                value=1,
                                min=0,
                                max=2,
                                style={"width": "100%"},
                            ),
                        ],
                        width=2,
                    ),
                    dbc.Col(
                        [
                            html.Small("RNA bulges"),
                            dcc.Input(
                                id="index-build-brna",
                                type="number",
                                value=1,
                                min=0,
                                max=2,
                                style={"width": "100%"},
                            ),
                        ],
                        width=2,
                    ),
                    dbc.Col(
                        [html.Small(" "), html.Button("Build", id="index-build-btn")],
                        width=2,
                    ),
                ],
                align="start",
            ),
            html.Small(
                "The index only needs bulge counts (max DNA/RNA bulges it will "
                "support); mismatches are set later, at search time."
            ),
            html.Div(id="index-feedback", style={"color": "#b00", "margin-top": "0.4rem"}),
        ]
        + _publish_index_section(),
    )

    # ---- VCF datasets ------------------------------------------------------
    vcf_card = _add_card(
        "Add a VCF dataset",
        "Variant datasets are large, so they are fetched server-side (from "
        "HuggingFace, a URL, or an existing folder on this machine) or chunk-"
        "uploaded. CRISPRme works with the files kept COMPRESSED — provide "
        "bgzip-compressed, per-chromosome files (the chromosome in the file name, "
        f"e.g. mydata.chr1.vcf.gz); no flat .vcf. Space: {_SIZE_HINTS['vcf']}.",
        [
            dbc.Row(
                [
                    dbc.Col(
                        dcc.Input(
                            id="vcf-name",
                            placeholder="dataset name, e.g. my_pig_vcf",
                            type="text",
                            style={"width": "100%"},
                        ),
                        width=4,
                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            id="vcf-source-dropdown",
                            options=[
                                {"label": "HuggingFace dataset", "value": "hf"},
                                {"label": "Register a server folder", "value": "path"},
                            ],
                            value="hf",
                            clearable=False,
                        ),
                        width=4,
                    ),
                    dbc.Col(html.Button("Add VCF dataset", id="vcf-add-btn"), width=4),
                ]
            ),
            html.Div(
                [
                    html.Small(
                        "Reference genome this VCF is called against (required — a "
                        "VCF only works with its matching reference):"
                    ),
                    dcc.Dropdown(
                        id="vcf-ref-genome",
                        options=installed_genomes,
                        placeholder="reference genome (e.g. hg38, susScr11)",
                    ),
                ],
                style={"margin-top": "0.4rem"},
            ),
            html.Div(
                dcc.Dropdown(
                    id="vcf-hf-name",
                    options=_hf_options(
                        "vcf", [v["value"] for v in get_all_vcf_datasets()]
                    ),
                    placeholder="…or pick a VCF dataset available on HuggingFace (source = HuggingFace)",
                ),
                style={"margin-top": "0.4rem"},
            ),
            dcc.Input(
                id="vcf-path-input",
                placeholder="absolute path to a folder of .vcf.gz (for 'Register a server folder')",
                type="text",
                style={"width": "100%", "margin-top": "0.5rem"},
            ),
            html.Small(
                "Or upload a local .vcf.gz into the dataset named above — chunked, "
                "no browser size limit."
            ),
            html.Div(
                className="crisprme-chunk-upload",
                style={"margin-top": "0.3rem"},
                **{
                    "data-target": "vcf",
                    "data-dataset-input": "vcf-name",
                    "data-genome-input": "vcf-ref-genome",
                },
            ),
            html.Div(id="vcf-feedback", style={"color": "#b00", "margin-top": "0.4rem"}),
        ],
    )

    # ---- Annotations -------------------------------------------------------
    annotation_card = _add_card(
        "Add an annotation (BED)",
        "Annotation BED files are small, so you can upload one directly.",
        [
            dcc.Upload(
                id="annotation-upload",
                children=html.Div(["Drag and drop or ", html.A("select a .bed file")]),
                style={
                    "width": "100%",
                    "height": "60px",
                    "lineHeight": "60px",
                    "borderWidth": "1px",
                    "borderStyle": "dashed",
                    "borderRadius": "5px",
                    "textAlign": "center",
                },
                multiple=False,
            ),
            html.Div(id="annotation-feedback", style={"color": "#b00", "margin-top": "0.4rem"}),
        ],
    )

    # ---- PAMs --------------------------------------------------------------
    pam_card = _add_card(
        "Define a nuclease / PAM",
        "PAM files are tiny, so just describe the nuclease and it is written "
        "directly (e.g. name SpCas9, length 20, motif NGG, 3').",
        [
            dbc.Row(
                [
                    dbc.Col(
                        dcc.Input(
                            id="pam-name",
                            placeholder="nuclease name, e.g. SpCas9",
                            type="text",
                            style={"width": "100%"},
                        ),
                        width=3,
                    ),
                    dbc.Col(
                        dcc.Input(
                            id="pam-length",
                            type="number",
                            value=20,
                            min=1,
                            max=30,
                            style={"width": "100%"},
                        ),
                        width=2,
                    ),
                    dbc.Col(
                        dcc.Input(
                            id="pam-motif",
                            placeholder="motif, e.g. NGG",
                            type="text",
                            style={"width": "100%"},
                        ),
                        width=3,
                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            id="pam-position",
                            options=[
                                {"label": "3' (Cas9-like)", "value": "3"},
                                {"label": "5' (Cas12a-like)", "value": "5"},
                            ],
                            value="3",
                            clearable=False,
                        ),
                        width=2,
                    ),
                    dbc.Col(html.Button("Add PAM", id="pam-add-btn"), width=2),
                ]
            ),
            html.Div(id="pam-feedback", style={"margin-top": "0.4rem"}),
        ],
    )

    # ---- Remove data -------------------------------------------------------
    delete_card = _add_card(
        "Remove installed data",
        "Delete a downloaded reference to free disk space. Anything you remove "
        "here can be downloaded again later.",
        [
            dbc.Row(
                [
                    dbc.Col(
                        dcc.Dropdown(
                            id="delete-item",
                            options=_deletable_options(),
                            placeholder="select data to delete",
                        ),
                        width=8,
                    ),
                    dbc.Col(
                        html.Button("Delete", id="delete-btn", style={"color": "#b00"}),
                        width=4,
                    ),
                ]
            ),
            dcc.ConfirmDialog(
                id="delete-confirm",
                message="Delete this data? You can download it again later.",
            ),
            html.Div(id="delete-feedback", style={"margin-top": "0.4rem"}),
        ],
    )

    return [
        dbc.Container(
            [
                html.H2("Settings / Data Manager", style={"margin-top": "1rem"}),
                html.P(
                    "Add data to your local CRISPRme folder. New data appears in "
                    "the search form automatically."
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                genome_card,
                                index_card,
                                vcf_card,
                                annotation_card,
                                pam_card,
                                delete_card,
                            ],
                            width=8,
                        ),
                        dbc.Col(
                            [
                                html.H5("Storage"),
                                html.Div(id="settings-storage", children=_render_storage()),
                                html.Hr(),
                                html.H5("Installed data"),
                                html.Div(id="settings-tables-container", children=_render_all_tables()),
                                html.Hr(),
                                html.H5("Activity"),
                                html.Div(
                                    id="settings-progress",
                                    children=html.I("No operation running."),
                                ),
                            ],
                            width=4,
                        ),
                    ]
                ),
                dcc.Store(id="settings-active-job"),
                dcc.Interval(id="settings-check", interval=3000, disabled=True),
            ],
            fluid=True,
            style={"margin-top": "1rem"},
        )
    ]


def _publish_index_section() -> List:
    """Maintainer-only 'publish index to HuggingFace' controls (else nothing)."""
    if not MAINTAINER_MODE:
        return []
    return [
        html.Hr(),
        html.B("Publish an index to HuggingFace (maintainers)"),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Dropdown(
                        id="index-publish-name",
                        options=get_available_indexes(),
                        placeholder="index to publish",
                    ),
                    width=8,
                ),
                dbc.Col(html.Button("Publish", id="index-publish-btn"), width=4),
            ]
        ),
    ]


# ---------------------------------------------------------------------------
# Launch callbacks (all guarded against ONLINE mode)
# ---------------------------------------------------------------------------
def _start(job_id: str):
    """Standard 3-tuple for a launched job: (store, interval-enabled, feedback)."""
    return job_id, False, ""


@app.callback(
    [
        Output("settings-active-job", "data", allow_duplicate=True),
        Output("settings-check", "disabled", allow_duplicate=True),
        Output("genome-feedback", "children"),
    ],
    [Input("genome-add-btn", "n_clicks")],
    [
        State("genome-assembly-input", "value"),
        State("genome-source-dropdown", "value"),
        State("genome-url-input", "value"),
        State("genome-hf-name", "value"),
    ],
    prevent_initial_call=True,
)
def add_genome(n, assembly, source, url, hf_name):
    if n is None or ONLINE:
        raise PreventUpdate
    if hf_name:  # a genome picked from the HuggingFace catalog takes precedence
        assembly, source = hf_name, "hf"
    err = _validate_name(assembly or "")
    if err:
        return no_update, no_update, err
    argv = ["download", "--what", "genome", "--ref", assembly.strip(), "--source", source]
    if source == "url":
        if not url:
            return no_update, no_update, "Direct URL source requires a URL."
        argv += ["--url", url.strip()]
    return _start(launch_settings_job(argv, "Download genome"))


@app.callback(
    [
        Output("settings-active-job", "data", allow_duplicate=True),
        Output("settings-check", "disabled", allow_duplicate=True),
        Output("index-feedback", "children"),
    ],
    [Input("index-hf-btn", "n_clicks")],
    [State("index-hf-name", "value")],
    prevent_initial_call=True,
)
def add_index_hf(n, name):
    if n is None or ONLINE:
        raise PreventUpdate
    # name comes from the HF catalog dropdown; index names may contain '+'
    # (variant indexes), so only guard against path-traversal / spaces here.
    if not name:
        return no_update, no_update, "Pick an index from the HuggingFace list."
    if any(c in name for c in ("/", "\\", " ", "..")):
        return no_update, no_update, "Invalid index name."
    argv = ["download", "--what", "index", "--index-name", name.strip()]
    return _start(launch_settings_job(argv, "Download index"))


@app.callback(
    [
        Output("settings-active-job", "data", allow_duplicate=True),
        Output("settings-check", "disabled", allow_duplicate=True),
        Output("index-feedback", "children", allow_duplicate=True),
    ],
    [Input("index-build-btn", "n_clicks")],
    [
        State("index-build-genome", "value"),
        State("index-build-pam", "value"),
        State("index-build-bdna", "value"),
        State("index-build-brna", "value"),
        State("index-build-vcf", "value"),
    ],
    prevent_initial_call=True,
)
def build_index(n, genome, pam, bdna, brna, vcf):
    if n is None or ONLINE:
        raise PreventUpdate
    if not genome or not pam:
        return no_update, no_update, "Select an installed genome and PAM."
    genome_dir = os.path.join(current_working_directory, "Genomes", genome.replace(" ", "_"))
    pam_file = os.path.join(current_working_directory, "PAMs", f"{pam}.txt")
    argv = [
        "build-index-only",
        "--genome",
        genome_dir,
        "--pam",
        pam_file,
        "--bDNA",
        str(int(bdna or 0)),
        "--bRNA",
        str(int(brna or 0)),
    ]
    stage = "Build index"
    if vcf:  # pre-build a variant-aware index (enrichment + SNP/indels indexing)
        argv += ["--vcf", os.path.join(current_working_directory, "VCFs", vcf)]
        stage = "Build variant index"
    return _start(launch_settings_job(argv, stage))


@app.callback(
    [
        Output("settings-active-job", "data", allow_duplicate=True),
        Output("settings-check", "disabled", allow_duplicate=True),
        Output("vcf-feedback", "children"),
    ],
    [Input("vcf-add-btn", "n_clicks")],
    [
        State("vcf-name", "value"),
        State("vcf-source-dropdown", "value"),
        State("vcf-path-input", "value"),
        State("vcf-hf-name", "value"),
        State("vcf-ref-genome", "value"),
    ],
    prevent_initial_call=True,
)
def add_vcf(n, name, source, path, hf_name, ref_genome):
    if n is None or ONLINE:
        raise PreventUpdate
    if source == "hf" and hf_name:  # a dataset picked from the HuggingFace catalog
        name = hf_name
    err = _validate_name(name or "")
    if err:
        return no_update, no_update, err
    if not ref_genome:
        return (
            no_update,
            no_update,
            "Select the reference genome this VCF is called against.",
        )
    name = name.strip()
    # record the reference genome so the search form only pairs this VCF with it
    _write_vcf_marker(name, ref_genome)
    if source == "hf":
        argv = ["download", "--what", "vcf", "--dataset", name]
        return _start(launch_settings_job(argv, "Download VCF"))
    # register an existing server folder by symlinking it under VCFs/<name>
    if not path or not os.path.isdir(path):
        return no_update, no_update, "Provide an absolute path to an existing folder."
    try:
        dest = os.path.join(current_working_directory, "VCFs", name)
        if not os.path.exists(dest):
            os.symlink(os.path.abspath(path), dest)
    except OSError as e:
        return no_update, no_update, f"Could not register folder: {e}"
    return (
        no_update,
        no_update,
        f"Registered VCF dataset '{name}' for {ref_genome}.",
    )


@app.callback(
    Output("annotation-feedback", "children"),
    [Input("annotation-upload", "contents")],
    [State("annotation-upload", "filename")],
    prevent_initial_call=True,
)
def add_annotation(contents, filename):
    if contents is None or ONLINE:
        raise PreventUpdate
    if not filename or not filename.endswith(".bed"):
        return "Please upload a .bed file."
    err = _validate_name(filename[:-4])
    if err:
        return err
    try:
        _content_type, content_string = contents.split(",", 1)
        data = base64.b64decode(content_string)
        os.makedirs(os.path.join(current_working_directory, "Annotations"), exist_ok=True)
        dest = os.path.join(current_working_directory, "Annotations", filename)
        with open(dest, "wb") as fout:
            fout.write(data)
    except Exception as e:
        return f"Upload failed: {e}"
    return f"Added annotation '{filename}'. Reload the page to see it in the list."


@app.callback(
    Output("pam-feedback", "children"),
    [Input("pam-add-btn", "n_clicks")],
    [
        State("pam-name", "value"),
        State("pam-length", "value"),
        State("pam-motif", "value"),
        State("pam-position", "value"),
    ],
    prevent_initial_call=True,
)
def add_pam(n, name, length, motif, position):
    if n is None or ONLINE:
        raise PreventUpdate
    err = _validate_name(name or "")
    if err:
        return err
    name = name.strip()
    if "-" in name:
        return "Nuclease name must not contain a hyphen."
    if not motif or any(c not in "ACGTRYSWKMBDHVN" for c in motif.upper()):
        return "Motif must be IUPAC nucleotide letters (e.g. NGG, TTTV)."
    motif = motif.upper()
    try:
        glen = int(length)
    except (TypeError, ValueError):
        return "Guide length must be a number."
    # compose the CRISPRme PAM token exactly like the built-in PAM files
    if position == "3":  # PAM at 3' end (Cas9-like): guide then motif
        token = "N" * glen + motif
        pos = len(motif)
    else:  # PAM at 5' end (Cas12a-like): motif then guide
        token = motif + "N" * glen
        pos = -len(motif)
    filename = f"{glen}bp-{motif}-{name}.txt"
    try:
        os.makedirs(os.path.join(current_working_directory, "PAMs"), exist_ok=True)
        with open(os.path.join(current_working_directory, "PAMs", filename), "w") as fout:
            fout.write(f"{token} {pos}\n")
    except OSError as e:
        return f"Could not write PAM: {e}"
    return html.Span(
        f"Added PAM '{filename}'. Reload the page to use it in a search.",
        style={"color": "green"},
    )


if MAINTAINER_MODE:

    @app.callback(
        [
            Output("settings-active-job", "data", allow_duplicate=True),
            Output("settings-check", "disabled", allow_duplicate=True),
            Output("index-feedback", "children", allow_duplicate=True),
        ],
        [Input("index-publish-btn", "n_clicks")],
        [State("index-publish-name", "value")],
        prevent_initial_call=True,
    )
    def publish_index_cb(n, name):
        if n is None or ONLINE:
            raise PreventUpdate
        if not name:
            return no_update, no_update, "Select an index to publish."
        index_dir = os.path.join(current_working_directory, "genome_library", name)
        argv = ["publish-index", "--index", index_dir]
        return _start(launch_settings_job(argv, "Publish index"))


# ---------------------------------------------------------------------------
# Progress polling
# ---------------------------------------------------------------------------
@app.callback(
    [
        Output("settings-progress", "children"),
        Output("settings-check", "disabled"),
        Output("settings-tables-container", "children"),
        Output("settings-storage", "children"),
    ],
    [Input("settings-check", "n_intervals")],
    [State("settings-active-job", "data")],
    prevent_initial_call=True,
)
def refresh_settings_job(n, job_id):
    if not job_id:
        raise PreventUpdate
    logtxt = os.path.join(current_working_directory, SETTINGS_DIR, job_id, "log.txt")
    if not os.path.isfile(logtxt):
        raise PreventUpdate
    with open(logtxt) as handle:
        log = handle.read()
    # find the stage name from the Start marker
    stage = "Operation"
    for line in log.splitlines():
        if "\tStart" in line:
            stage = line.split("\t")[0]
            break
    if f"{stage}\tEnd" in log:
        # success: stop polling and refresh the installed-data tables + storage
        return (
            html.Div(f"{stage}: done.", style={"color": "green"}),
            True,
            _render_all_tables(),
            _render_storage(),
        )
    if f"{stage}\tFAILED" in log:
        logerr = os.path.join(current_working_directory, SETTINGS_DIR, job_id, "log_error.txt")
        tail = _error_tail(logerr)
        return (
            html.Div(
                [html.Div(f"{stage}: failed.", style={"color": "#b00"}), html.Small(tail)]
            ),
            True,
            no_update,
            no_update,
        )
    # still running
    return (
        html.Div(
            [
                dbc.Spinner(size="sm"),
                html.Span(f"  {stage}: running... (this can take a while)"),
            ]
        ),
        False,
        no_update,
        no_update,
    )


# ---------------------------------------------------------------------------
# Delete installed data
# ---------------------------------------------------------------------------
@app.callback(
    Output("delete-confirm", "displayed"),
    [Input("delete-btn", "n_clicks")],
    [State("delete-item", "value")],
    prevent_initial_call=True,
)
def ask_delete(n, item):
    if n is None or ONLINE or not item:
        raise PreventUpdate
    return True  # open the confirmation dialog


@app.callback(
    [
        Output("delete-feedback", "children"),
        Output("settings-tables-container", "children", allow_duplicate=True),
        Output("settings-storage", "children", allow_duplicate=True),
    ],
    [Input("delete-confirm", "submit_n_clicks")],
    [State("delete-item", "value")],
    prevent_initial_call=True,
)
def do_delete(n, item):
    if not n or ONLINE or not item:
        raise PreventUpdate
    kind, _sep, name = item.partition(":")
    try:
        path = _delete_path(kind, name)
    except ValueError as e:
        return html.Span(str(e), style={"color": "#b00"}), no_update, no_update
    if not os.path.exists(path) and not os.path.islink(path):
        return (
            html.Span(f"{item} not found (already removed?).", style={"color": "#b00"}),
            _render_all_tables(),
            _render_storage(),
        )
    freed = 0
    try:
        if os.path.islink(path):
            os.unlink(path)  # a registered server folder: only drop the symlink
        elif os.path.isdir(path):
            freed = _dir_size(path)
            shutil.rmtree(path)
        else:
            freed = os.path.getsize(path)
            os.remove(path)
    except OSError as e:
        return html.Span(f"Delete failed: {e}", style={"color": "#b00"}), no_update, no_update
    return (
        html.Span(
            f"Deleted {item} (freed {_fmt_size(freed)}). You can download it again "
            "later.",
            style={"color": "green"},
        ),
        _render_all_tables(),
        _render_storage(),
    )
