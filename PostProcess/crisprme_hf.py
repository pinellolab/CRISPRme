"""HuggingFace-backed fast download/upload of CRISPRme reference data and
precomputed indexes.

CRISPRme's reference bundle (genome, annotations, PAMs, sample-ID files, variant
VCFs) and its precomputed CRISPRitz indexes can be hosted on a HuggingFace
dataset repository and fetched over HF's CDN, which is typically far faster and
more reliable than the original UCSC/FTP sources. This module is the thin client
for that: it resolves the target repo, downloads requested components into the
canonical CRISPRme directory layout (decompressing where needed), and — on the
publish side — uploads a locally built index back to the repo.

Design notes
------------
* Everything network-facing goes through ``huggingface_hub``. The dependency is
  imported lazily so that importing this module (e.g. for the pure-logic unit
  tests) never requires it; only the functions that actually touch HF do.
* Public downloads need no token. Uploads (``publish_index``) need a write token,
  read from ``--token`` or the ``HF_TOKEN`` / ``HUGGING_FACE_HUB_TOKEN`` env vars.
* The remote layout mirrors the local one so mapping is trivial:

      genomes/<ref>/<chrom>.fa[.gz]      -> Genomes/<ref>/<chrom>.fa
      annotations/<file>                 -> Annotations/<file>
      PAMs/<file>                        -> PAMs/<file>
      samplesIDs/<file>                  -> samplesIDs/<file>
      VCFs/<dataset>/<file>.vcf.gz       -> VCFs/<dataset>/<file>.vcf.gz
      indexes/<index_name>.tar.gz        -> genome_library/<index_name>/

  Large FASTA is stored gzipped and decompressed after download; VCFs are kept
  bgzipped (CRISPRme consumes them that way); indexes travel as a single
  ``.tar.gz`` and are unpacked into ``genome_library/``.
"""

from typing import Dict, List, Optional
import os
import sys
import io
import json
import gzip
import shutil
import tarfile
from datetime import datetime, timezone

# Default HuggingFace dataset repo. Overridable per-invocation with --hf-repo or
# the CRISPRME_HF_REPO environment variable (e.g. switch to a pinellolab org
# repo once it exists).
DEFAULT_HF_REPO = "lucapinello/crisprme-data"

# Components understood by the download command and their remote sub-paths.
# This layout is the single source of truth for the HF dataset tree; the
# transparent setup/complete_test fast-path (utils.hf_fetch) targets the same
# paths (genomes/<ref>/, vcfs/<dataset>/, samplesIDs/, annotations/, indexes/).
_COMPONENT_PREFIXES = {
    "genome": "genomes",
    "annotations": "annotations",
    "pams": "pams",
    "samples": "samplesIDs",
    "vcf": "vcfs",
    "index": "indexes",
}

# Local destination directory for each component (relative to the working dir).
_COMPONENT_LOCALDIR = {
    "genome": "Genomes",
    "annotations": "Annotations",
    "pams": "PAMs",
    "samples": "samplesIDs",
    "vcf": "VCFs",
    "index": "genome_library",
}


def resolve_repo(cli_repo: Optional[str] = None) -> str:
    """Resolve the HuggingFace repo id: CLI value > env var > built-in default."""
    if cli_repo:
        return cli_repo
    return os.environ.get("CRISPRME_HF_REPO", DEFAULT_HF_REPO)


def resolve_token(cli_token: Optional[str] = None) -> Optional[str]:
    """Resolve an HF token for uploads.

    Order: explicit ``--token`` > ``HF_TOKEN`` > ``HUGGING_FACE_HUB_TOKEN`` >
    a cached ``huggingface-cli login`` token. The cached-login fallback means
    that on a machine already logged into HuggingFace, ``publish-index`` works
    with no extra flags or env vars (matching how other HF tools behave).
    """
    tok = (
        cli_token
        or os.environ.get("HF_TOKEN")
        or os.environ.get("HUGGING_FACE_HUB_TOKEN")
    )
    if tok:
        return tok
    try:  # cached `huggingface-cli login` token, if any (no network)
        from huggingface_hub import get_token

        return get_token()
    except Exception:  # huggingface_hub absent or no cached token
        return None


def _require_hf():
    """Import huggingface_hub lazily, with an actionable error if it is missing."""
    try:
        import huggingface_hub  # noqa: F401

        return huggingface_hub
    except ImportError as e:  # pragma: no cover - exercised only without the dep
        raise ImportError(
            "The 'huggingface_hub' package is required for HuggingFace "
            "downloads/uploads but is not installed. Install it with "
            "'conda install -c conda-forge huggingface_hub' or "
            "'pip install huggingface_hub'."
        ) from e


def list_available_downloads(
    component: str, repo: Optional[str] = None
) -> List[Dict]:
    """List items available to download from the HuggingFace dataset repo.

    Returns ``[{"name": str, "size": int}, ...]`` for the given component:
    ``genome``/``vcf`` -> the assembly/dataset folders under ``genomes/`` or
    ``vcfs/`` (size = sum of that folder's files); ``index`` -> the
    ``indexes/<name>.tar.gz`` archives. Returns ``[]`` on any error (offline,
    missing repo, dependency absent) so callers can fall back gracefully.
    """
    if component not in ("genome", "vcf", "index"):
        return []
    try:
        hf = _require_hf()
        repo = resolve_repo(repo)
        api = hf.HfApi()
        prefix = _COMPONENT_PREFIXES[component]
        if component == "index":
            items = []
            for entry in api.list_repo_tree(
                repo, repo_type="dataset", path_in_repo=prefix
            ):
                size = getattr(entry, "size", None)
                path = getattr(entry, "path", "")
                if size is not None and path.endswith(".tar.gz"):
                    name = os.path.basename(path)[: -len(".tar.gz")]
                    items.append({"name": name, "size": size})
            return sorted(items, key=lambda d: d["name"])
        # genome / vcf: folders under prefix/, sum their files' sizes
        sizes: Dict[str, int] = {}
        for entry in api.list_repo_tree(
            repo, repo_type="dataset", path_in_repo=prefix, recursive=True
        ):
            size = getattr(entry, "size", None)
            path = getattr(entry, "path", "")
            if size is None:  # a folder entry, not a file
                continue
            rel = path[len(prefix) + 1 :]
            top = rel.split("/")[0]
            if top:
                sizes[top] = sizes.get(top, 0) + size
        return [{"name": k, "size": v} for k, v in sorted(sizes.items())]
    except Exception:
        return []


def decompress_gz(path: str, remove_source: bool = True) -> str:
    """Gunzip ``path`` (``*.gz`` -> ``*``) and return the decompressed path.

    Used for FASTA that is stored gzipped on HF to save space/bandwidth but must
    be plain text for CRISPRme/CRISPRitz. VCFs are intentionally NOT run through
    this (they must stay bgzipped).
    """
    if not path.endswith(".gz"):
        return path
    out = path[:-3]
    with gzip.open(path, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    if remove_source:
        os.remove(path)
    return out


def _hf_snapshot(repo: str, allow_patterns: List[str], local_dir: str,
                 token: Optional[str] = None) -> str:
    """Download the files matching ``allow_patterns`` from ``repo`` into
    ``local_dir`` (flat, following the repo tree). Returns ``local_dir``."""
    hf = _require_hf()
    os.makedirs(local_dir, exist_ok=True)
    hf.snapshot_download(
        repo_id=repo,
        repo_type="dataset",
        allow_patterns=allow_patterns,
        local_dir=local_dir,
        token=token,
    )
    return local_dir


def download_component(
    component: str,
    workdir: str,
    repo: Optional[str] = None,
    ref: str = "hg38",
    dataset: Optional[str] = None,
    index_name: Optional[str] = None,
    token: Optional[str] = None,
) -> str:
    """Download one CRISPRme reference component from HF into the canonical
    local layout under ``workdir``.

    Args:
        component: one of genome|annotations|pams|samples|vcf|index.
        workdir: working directory the CRISPRme dir-tree lives under.
        repo: HF repo id (defaults via :func:`resolve_repo`).
        ref: reference-genome name (for ``genome``/``index``), e.g. "hg38".
        dataset: variant dataset name for ``vcf`` (e.g. "1000G", "HGDP").
        index_name: precomputed index directory name for ``index``
            (e.g. "NGG_2_hg38").

    Returns:
        The local destination directory populated by the download.
    """
    if component not in _COMPONENT_PREFIXES:
        raise ValueError(
            f"Unknown component '{component}'. Expected one of: "
            f"{', '.join(sorted(_COMPONENT_PREFIXES))}."
        )
    repo = resolve_repo(repo)
    remote_prefix = _COMPONENT_PREFIXES[component]
    local_dir = os.path.join(workdir, _COMPONENT_LOCALDIR[component])

    if component == "genome":
        patterns = [f"{remote_prefix}/{ref}/*"]
        staging = os.path.join(workdir, ".hf_stage_genome")
        _hf_snapshot(repo, patterns, staging, token)
        src = os.path.join(staging, remote_prefix, ref)
        if not os.path.isdir(src):
            shutil.rmtree(staging, ignore_errors=True)
            raise ValueError(
                f"No genome '{ref}' found in HuggingFace repo '{repo}' "
                f"(looked for '{remote_prefix}/{ref}/'). Check --ref, or the repo "
                f"may not have this genome uploaded yet."
            )
        dest = os.path.join(local_dir, ref)
        os.makedirs(dest, exist_ok=True)
        for fn in sorted(os.listdir(src)):
            moved = shutil.move(os.path.join(src, fn), os.path.join(dest, fn))
            if moved.endswith(".gz"):
                decompress_gz(moved)
        shutil.rmtree(staging, ignore_errors=True)
        return dest

    if component == "vcf":
        if not dataset:
            raise ValueError("--dataset is required to download a 'vcf' component")
        patterns = [f"{remote_prefix}/{dataset}/*"]
        staging = os.path.join(workdir, ".hf_stage_vcf")
        _hf_snapshot(repo, patterns, staging, token)
        src = os.path.join(staging, remote_prefix, dataset)
        if not os.path.isdir(src):
            shutil.rmtree(staging, ignore_errors=True)
            raise ValueError(
                f"No VCF dataset '{dataset}' found in HuggingFace repo '{repo}' "
                f"(looked for '{remote_prefix}/{dataset}/'). Check --dataset, or "
                f"the repo may not have this dataset uploaded yet."
            )
        dest = os.path.join(local_dir, dataset)
        os.makedirs(dest, exist_ok=True)
        for fn in sorted(os.listdir(src)):
            shutil.move(os.path.join(src, fn), os.path.join(dest, fn))  # keep bgzip
        shutil.rmtree(staging, ignore_errors=True)
        return dest

    if component == "index":
        if not index_name:
            raise ValueError("--index-name is required to download an 'index' component")
        patterns = [f"{remote_prefix}/{index_name}.tar.gz"]
        staging = os.path.join(workdir, ".hf_stage_index")
        _hf_snapshot(repo, patterns, staging, token)
        tarball = os.path.join(staging, remote_prefix, f"{index_name}.tar.gz")
        if not os.path.isfile(tarball):
            shutil.rmtree(staging, ignore_errors=True)
            raise ValueError(
                f"No index '{index_name}' found in HuggingFace repo '{repo}' "
                f"(looked for '{remote_prefix}/{index_name}.tar.gz'). Check "
                f"--index-name, or build it locally with 'crisprme.py "
                f"build-index-only' — the repo may not have it uploaded yet."
            )
        os.makedirs(local_dir, exist_ok=True)
        with tarfile.open(tarball) as tf:
            # surface the provenance manifest (if present) without dropping it
            # into genome_library/; extract only the index folder members
            members = [m for m in tf.getmembers() if m.name != "manifest.json"]
            tf.extractall(local_dir, members=members)  # -> genome_library/<index_name>/
            try:
                mf = tf.extractfile("manifest.json")
                if mf is not None:
                    meta = json.load(mf)
                    sys.stderr.write(
                        f"Downloaded index {index_name} "
                        f"(built {meta.get('created_at', 'unknown')})\n"
                    )
            except (KeyError, json.JSONDecodeError):
                pass  # no/!invalid manifest — proceed, the index itself is what matters
        shutil.rmtree(staging, ignore_errors=True)
        return os.path.join(local_dir, index_name)

    # flat components: annotations, pams, samples
    patterns = [f"{remote_prefix}/*"]
    staging = os.path.join(workdir, f".hf_stage_{component}")
    _hf_snapshot(repo, patterns, staging, token)
    src = os.path.join(staging, remote_prefix)
    if not os.path.isdir(src):
        shutil.rmtree(staging, ignore_errors=True)
        raise ValueError(
            f"No '{component}' files found in HuggingFace repo '{repo}' "
            f"(looked for '{remote_prefix}/'). The repo may not have this "
            f"component uploaded yet."
        )
    os.makedirs(local_dir, exist_ok=True)
    for fn in sorted(os.listdir(src)):
        moved = shutil.move(os.path.join(src, fn), os.path.join(local_dir, fn))
        if component == "annotations" and moved.endswith(".gz") and not moved.endswith(".bed.gz"):
            decompress_gz(moved)
    shutil.rmtree(staging, ignore_errors=True)
    return local_dir


def publish_index(
    index_dir: str,
    repo: Optional[str] = None,
    token: Optional[str] = None,
) -> str:
    """Tar a locally built genome_library index and upload it to HF.

    Args:
        index_dir: path to the ``genome_library/<index_name>`` directory.
        repo: HF repo id (defaults via :func:`resolve_repo`).
        token: HF write token (defaults via :func:`resolve_token`).

    Returns:
        The remote path (``indexes/<index_name>.tar.gz``) the index was uploaded to.
    """
    index_dir = os.path.abspath(index_dir.rstrip("/"))
    if not os.path.isdir(index_dir):
        raise ValueError(f"Index directory {index_dir} does not exist")
    index_name = os.path.basename(index_dir)
    repo = resolve_repo(repo)
    token = resolve_token(token)
    if not token:
        raise ValueError(
            "An HF write token is required to publish an index. Provide --token "
            "or set HF_TOKEN (never commit it)."
        )
    hf = _require_hf()
    # self-describing provenance travels inside the tarball (adopted from the
    # #141 design): PAM / bulge / genome parsed from the folder name
    # <PAM>_<bMax+1>_<ref>, plus a UTC timestamp.
    manifest = {"name": index_name, "created_at": datetime.now(timezone.utc).isoformat()}
    parts = index_name.rsplit("_", 2)
    if len(parts) == 3:
        manifest["pam"], manifest["index_bmax"], manifest["genome"] = (
            parts[0],
            parts[1],
            parts[2],
        )
    tarball = f"{index_dir}.tar.gz"
    with tarfile.open(tarball, "w:gz") as tf:
        tf.add(index_dir, arcname=index_name)  # -> genome_library/<name>/ on unpack
        manifest_bytes = json.dumps(manifest, indent=2).encode()
        info = tarfile.TarInfo(name="manifest.json")
        info.size = len(manifest_bytes)
        tf.addfile(info, io.BytesIO(manifest_bytes))
    remote_path = f"indexes/{index_name}.tar.gz"
    api = hf.HfApi()
    api.upload_file(
        path_or_fileobj=tarball,
        path_in_repo=remote_path,
        repo_id=repo,
        repo_type="dataset",
        token=token,
    )
    os.remove(tarball)
    return remote_path
